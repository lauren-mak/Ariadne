//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_data.hpp"
#include "io/reads_io/read_processor.hpp"
#include "valid_kmer_generator.hpp"

#include "io/reads_io/ireadstream.hpp"
#include "config_struct_hammer.hpp"

#include "data_structures/mph_index/kmer_index_builder.hpp"

#include "io/kmers_io/kmer_iterator.hpp"
#include "utils/adt/bf.hpp"

using namespace hammer;

class BufferFiller;

struct KMerComparator {
    bool operator()(const KMer &l, const KMer &r) const {
      for (size_t i = 0; i < KMer::DataSize ; ++i) {
        if (l.data()[i] != r.data()[i]) {
          return (l.data()[i] < r.data()[i]);
        }
      }

      return false;
    }
};


class HammerFilteringKMerSplitter : public KMerSortingSplitter<hammer::KMer> {
 public:
  typedef std::function<bool(const KMer&)> KMerFilter;

  HammerFilteringKMerSplitter(std::string &work_dir,
                              KMerFilter filter = [](const KMer&) { return true; })
      : KMerSortingSplitter<hammer::KMer>(work_dir, hammer::K),
      filter_(std::move(filter)) {}

  path::files_t Split(size_t num_files) override;

 private:
  KMerFilter filter_;

  friend class BufferFiller;
};

class BufferFiller {
  size_t processed_;
  HammerFilteringKMerSplitter &splitter_;

 public:
  BufferFiller(HammerFilteringKMerSplitter &splitter)
      : processed_(0), splitter_(splitter) {}

  size_t processed() const { return processed_; }

  bool operator()(const Read &r) {
    int trim_quality = cfg::get().input_trim_quality;

    // FIXME: Get rid of this
    Read cr = r;
    size_t sz = cr.trimNsAndBadQuality(trim_quality);

    #pragma omp atomic
    processed_ += 1;

    if (sz < hammer::K)
      return false;

    unsigned thread_id = omp_get_thread_num();
    ValidKMerGenerator<hammer::K> gen(cr);
    bool stop = false;
    for (; gen.HasMore(); gen.Next()) {
      KMer seq = gen.kmer();
      if (!splitter_.filter_(seq))
        continue;

      stop |= splitter_.push_back_internal( seq, thread_id);
      stop |= splitter_.push_back_internal(!seq, thread_id);
    }

    return stop;
  }
};

path::files_t HammerFilteringKMerSplitter::Split(size_t num_files) {
  unsigned nthreads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
  size_t reads_buffer_size = cfg::get().count_split_buffer;

  INFO("Splitting kmer instances into " << num_files << " buckets. This might take a while.");

  path::files_t out = PrepareBuffers(num_files, nthreads, reads_buffer_size);

  size_t n = 15;
  BufferFiller filler(*this);
  const auto& dataset = cfg::get().dataset;
  for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
    INFO("Processing " << *I);
    ireadstream irs(*I, cfg::get().input_qvoffset);
    while (!irs.eof()) {
      hammer::ReadProcessor rp(nthreads);
      rp.Run(irs, filler);
      DumpBuffers(out);
      VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

      if (filler.processed() >> n) {
        INFO("Processed " << filler.processed() << " reads");
        n += 1;
      }
    }
  }
  INFO("Processed " << filler.processed() << " reads");

  return out;
}

static inline void Merge(KMerStat &lhs, const KMerStat &rhs) {
  lhs.set_count(lhs.count() + rhs.count());
  lhs.total_qual *= rhs.total_qual;
  lhs.qual += rhs.qual;
}

static void PushKMer(KMerData &data,
                     KMer kmer, const unsigned char *q, double prob) {
  size_t idx = data.checking_seq_idx(kmer);
  if (idx == -1ULL)
      return;
  KMerStat &kmc = data[idx];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, (float)prob, q));
  kmc.unlock();
}

static void PushKMerRC(KMerData &data,
                       KMer kmer, const unsigned char *q, double prob) {
  unsigned char rcq[K];

  // Prepare RC kmer with quality.
  kmer = !kmer;
  for (unsigned i = 0; i < K; ++i)
    rcq[K - i - 1] = q[i];

  size_t idx = data.checking_seq_idx(kmer);
  if (idx == -1ULL)
      return;
  KMerStat &kmc = data[idx];
  kmc.lock();
  Merge(kmc,
        KMerStat(1, (float)prob, rcq));
  kmc.unlock();
}

class KMerDataFiller {
  KMerData &data_;

 public:
  KMerDataFiller(KMerData &data)
      : data_(data) {}

  bool operator()(const Read &r) {
    int trim_quality = cfg::get().input_trim_quality;

    // FIXME: Get rid of this
    Read cr = r;
    size_t sz = cr.trimNsAndBadQuality(trim_quality);

    if (sz < hammer::K)
      return false;

    ValidKMerGenerator<hammer::K> gen(cr);
    const char *q = cr.getQualityString().data();
    while (gen.HasMore()) {
      KMer kmer = gen.kmer();
      const unsigned char *kq = (const unsigned char*)(q + gen.pos() - 1);

      PushKMer(data_, kmer, kq, 1 - gen.correct_probability());
      PushKMerRC(data_, kmer, kq, 1 - gen.correct_probability());

      gen.Next();
    }

    return false;
  }
};

class KMerMultiplicityCounter {
  bf::bitcounting_bloom_filter<KMer, 2> bf_;
  size_t processed_;

  public:
  KMerMultiplicityCounter(size_t size)
      : bf_([](const KMer &k, uint64_t seed) { return k.GetHash((uint32_t)seed); },
            4 * size),
        processed_(0) {}

  ~KMerMultiplicityCounter() {}

  bool operator()(const Read &r) {
      int trim_quality = cfg::get().input_trim_quality;

      // FIXME: Get rid of this
      Read cr = r;
      size_t sz = cr.trimNsAndBadQuality(trim_quality);

#   pragma omp atomic
      processed_ += 1;

      if (sz < hammer::K)
          return false;

      ValidKMerGenerator<hammer::K> gen(cr);
      for (; gen.HasMore(); gen.Next()) {
          KMer kmer = gen.kmer();

          bf_.add(kmer);
          bf_.add(!kmer);
      }

      return (processed_ % (256*1024) == 0 && processed_ > 1024);
  }

  size_t count(const KMer &k) const {
      return bf_.lookup(k);
  }

  size_t processed() const {
      return processed_;
  }
};

void KMerDataCounter::BuildKMerIndex(KMerData &data) {
  // Build the index
  std::string workdir = cfg::get().input_working_dir;

  // Optionally perform a filtering step
  size_t kmers = 0;
  std::string final_kmers;
  if (cfg::get().count_filter_singletons) {
      INFO("Filtering singleton k-mers");

      size_t buffer_size = cfg::get().count_split_buffer;
      if (buffer_size == 0)
          buffer_size = 512 * 1024 * 1024;
      KMerMultiplicityCounter mcounter(buffer_size);

      size_t n = 15;
      for (const auto &reads : cfg::get().dataset.reads()) {
          INFO("Processing " << reads);
          ireadstream irs(reads, cfg::get().input_qvoffset);
          while (!irs.eof()) {
              hammer::ReadProcessor rp(omp_get_max_threads());
              rp.Run(irs, mcounter);
              VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");

              if (mcounter.processed() >> n) {
                  INFO("Processed " << mcounter.processed() << " reads");
                  n += 1;
              }
          }
      }
      INFO("Total " << mcounter.processed() << " reads processed");

      // FIXME: Reduce code duplication
      HammerFilteringKMerSplitter splitter(workdir,
                                           [&] (const KMer &k) { return mcounter.count(k) > 1; });
      KMerDiskCounter<hammer::KMer> counter(workdir, splitter);

      kmers = KMerIndexBuilder<HammerKMerIndex>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter, /* save final */ true);
      final_kmers = counter.GetFinalKMersFname();
  } else {
      HammerFilteringKMerSplitter splitter(workdir);
      KMerDiskCounter<hammer::KMer> counter(workdir, splitter);

      kmers = KMerIndexBuilder<HammerKMerIndex>(workdir, num_files_, omp_get_max_threads()).BuildIndex(data.index_, counter, /* save final */ true);
      final_kmers = counter.GetFinalKMersFname();
  }


  // Check, whether we'll ever have enough memory for running BH and bail out earlier
  double needed = 1.25 * (double)kmers * (sizeof(KMerStat) + sizeof(hammer::KMer));
  if (needed > (double) get_memory_limit())
      FATAL_ERROR("The reads contain too many k-mers to fit into available memory. You need approx. "
                  << needed / 1024.0 / 1024.0 / 1024.0
                  << "GB of free RAM to assemble your dataset");

  {
    INFO("Arranging kmers in hash map order");
    data.kmers_.set_size(kmers);
    data.kmers_.set_data(new hammer::KMer::DataType[kmers * hammer::KMer::GetDataSize(hammer::K)]);

    unsigned nthreads = std::min(cfg::get().count_merge_nthreads, cfg::get().general_max_nthreads);
    auto kmers_its = io::make_kmer_iterator<hammer::KMer>(final_kmers, hammer::K, 16*nthreads);

#   pragma omp parallel for num_threads(nthreads) schedule(guided)
    for (size_t i = 0; i < kmers_its.size(); ++i) {
        auto &kmer_it = kmers_its[i];
        for (; kmer_it.good(); ++kmer_it) {
            size_t kidx = data.index_.seq_idx(hammer::KMer(hammer::K, *kmer_it));
            memcpy(data.kmers_[kidx].data(), *kmer_it, hammer::KMer::TotalBytes);
        }
    }

    unlink(final_kmers.c_str());
  }
}

void KMerDataCounter::FillKMerData(KMerData &data) {
  // Now use the index to fill the kmer quality information.
  INFO("Collecting K-mer information, this takes a while.");
  data.data_.resize(data.kmers_.size());

  KMerDataFiller filler(data);
  const auto& dataset = cfg::get().dataset;
  for (auto I = dataset.reads_begin(), E = dataset.reads_end(); I != E; ++I) {
    INFO("Processing " << *I);
    ireadstream irs(*I, cfg::get().input_qvoffset);
    hammer::ReadProcessor rp(omp_get_max_threads());
    rp.Run(irs, filler);
    VERIFY_MSG(rp.read() == rp.processed(), "Queue unbalanced");
  }

  INFO("Collection done, postprocessing.");

  size_t singletons = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    VERIFY(data[i].count());

    // Make sure all the kmers are marked as 'Bad' in the beginning
    data[i].mark_bad();

    if (data[i].count() == 1)
      singletons += 1;
  }

  INFO("There are " << data.size() << " kmers in total. "
       "Among them " << singletons << " (" <<  100.0 * (double)singletons / (double)data.size() << "%) are singletons.");
}