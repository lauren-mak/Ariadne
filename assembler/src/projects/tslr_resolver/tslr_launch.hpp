#pragma once

#include <barcode_map_construction.hpp>
#include <tslr_resolver.hpp>
#include <projects/spades/pair_info_count.hpp>
#include <modules/stages/simplification.hpp>
#include <modules/pipeline/genomic_info_filler.hpp>
#include <projects/spades/gap_closer.hpp>
#include "projects/spades/distance_estimation.hpp"

namespace spades {

    void run_tslr_resolver(const std::string& path_to_tslr_dataset, const std::string& path_to_reference) {
        INFO("Starting from stage " << cfg::get().entry_point.c_str());

        debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                                cfg::get().tmp_dir,
                                                cfg::get().ds.reads.lib_count(),
                                                cfg::get().ds.reference_genome,
                                                cfg::get().flanking_range,
                                                cfg::get().pos.max_mapping_gap,
                                                cfg::get().pos.max_gap_diff);
        StageManager manager({cfg::get().developer_mode,
                              cfg::get().load_from,
                              cfg::get().output_saves});
        conj_gp.kmer_mapper.Attach();
        manager.add(new debruijn_graph::Construction())
                .add(new debruijn_graph::GenomicInfoFiller());
        if (!cfg::get().ts_res.ideal_reads) {
            manager.add(new debruijn_graph::GapClosing("early_gapcloser"))
                   .add(new debruijn_graph::Simplification)
                   .add(new debruijn_graph::GapClosing("late_gapcloser"))
                   .add(new debruijn_graph::SimplificationCleanup());
        }

        manager.add(new BarcodeMapConstructionStage(cfg::get().K, path_to_tslr_dataset))
                .add(new debruijn_graph::PairInfoCount())
                .add(new debruijn_graph::DistanceEstimation())
                .add(new TslrResolverStage(cfg::get().K, cfg::get().output_dir + "resolver_output.fasta", path_to_reference));
        INFO("Output directory: " << cfg::get().output_dir);

        manager.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("TSLR resolver finished.");
    }
} //spades
