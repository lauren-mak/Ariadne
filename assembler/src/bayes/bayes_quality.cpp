/*
 * bayes_quality.cpp
 *
 *  Created on: 23.03.2011
 *      Author: snikolenko
 */

#include "bayes_quality.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include "quality.hpp"
#include "read.hpp"
#include "ireadstream.hpp"

using namespace std;

LOGGER("b");

namespace bayes_quality {

bool BayesQualityGenome::isAvailable(size_t readno, size_t j, const PSeq & curpseq) {
	if (readno < 0) return true;
	size_t mid = readMids[readno];
	if (!(availableReadsHead_[curpseq][readno])) return false;
	PSeq ps = PSeq(genome_.Subseq(j + mid - INS, j + mid - INS +PREPROCESS_SEQ_LENGTH));
	for (size_t testj = j + mid - INS; testj < j + mid + DEL + 1; ++testj) {
		if (availableReadsMid_[ps][readno]) return true;
		ps = ps << genome_[testj+1];
	}
	return false;
}

int BayesQualityGenome::simpleMatch(const LASeq & seq, const LASeq & gen) {
	int res = 0;
	for (size_t i=0; i<LA_SIZE; ++i) {
		if (seq[i] != gen[i]) ++res;
		if (res > LA_ERRORS) break;
	}
	return res;
}

MatchResults BayesQualityGenome::ProcessOneReadBQ(const Sequence & seq, const QVector & qvec, size_t readno, const vector<size_t> & indices) {
	MatchResults mr_;
	// compute total quality of the read
	mr_.totalq_ = accumulate(qvec.begin(), qvec.end(), 0);
	
	size_t seqsize = seq.size();
	
	// initialize temporary vectors of possible start positions
	vector<QVector*> qv(INS+DEL+1);
	size_t qvsize_ = gensize_ - seqsize + 1 + INS;
	if (indices.size() > 0) qvsize_ = indices.size();
	for (size_t i=0; i < INS+DEL+1; ++i) {
		qv[i] = new QVector(qvsize_);
		fill(qv[i]->begin(), qv[i]->end(), VERYLARGEQ);
	}
	
	INFO("thread " << omp_get_thread_num() << ": " << seq.str());
	// { ostringstream os; for (size_t i=0; i<indices.size(); ++i) os << indices[i] << " "; INFO(os.str()); }

	#ifdef USE_PREPROCESSING
	PSeq curpseq = PSeq(genome_.Subseq(0,PREPROCESS_SEQ_LENGTH));
	#endif
	
	#ifdef USE_LOOKAHEAD
	LASeq laseq  = LASeq(genome_.Subseq(0, LA_SIZE));
	LASeq laread = LASeq(seq.Subseq(0, LA_SIZE));
	size_t mid = seqsize / 2;
	if (seqsize - mid < LA_SIZE) mid = seqsize - LA_SIZE;
	LASeq larmid = LASeq(seq.Subseq(mid, mid+LA_SIZE));
	vector<LASeq> lavec;
	for (size_t k=mid-DEL; k<=mid+INS; ++k) {
		lavec.push_back(LASeq(genome_.Subseq(k, k+LA_SIZE)));
	}
	#endif
	
	char curgc, curgs;
	bool badPreprocPosition = false; bool goodLAPosition = true;
	
	// dynamic programming
	for (size_t jind = 0; jind < qvsize_; ++jind) {
		size_t j = jind; if (indices.size() > 0) j = indices[jind];
		if (j > gensize_ - seqsize - INS) continue;
		// INFO("thread " << omp_get_thread_num() << ": " << j << " size=" << genome_.size() << " indices " << j+INS+1);
		++totalPos_;
		
		#ifdef USE_PREPROCESSING
		// is the read already ruled out during preprocessing?
		badPreprocPosition = !isAvailable(readno, j, curpseq);
		curpseq = curpseq << genome_[j+PREPROCESS_SEQ_LENGTH];
		#endif
		
		#ifdef USE_LOOKAHEAD
		goodLAPosition = (simpleMatch(laread, laseq) <= LA_ERRORS);
		laseq = laseq << genome_[j+LA_SIZE];
		for (size_t k=0; k<=DEL+INS+1; ++k) {
			goodLAPosition = goodLAPosition || (simpleMatch(larmid, lavec[k]) <= LA_ERRORS);
			if (j+mid-DEL+k+LA_SIZE < gensize_) lavec[k] = lavec[k] << genome_[j+mid-DEL+k+LA_SIZE];
			else goodLAPosition = true; // if we are at the very end already, we should skip this check
		}
		#endif

		if (badPreprocPosition || !goodLAPosition) {
			continue;
		}
		
		++totalGood_;

		for (size_t i = 0; i < seqsize+1; ++i) {			
			// genome is over
			if ( i + j > gensize_ ) {
				for (size_t s=0; s < seqsize-i+1; ++s) qv[s]->at(jind) = VERYLARGEQ;
				continue;
			}
						
			if ( i == 0 ) {  // first column
				qv[0]->at(jind) = 0;
				for (size_t s=1; s <= INS; ++s)	qv[s]->at(jind) = INSERT_Q*s;
				for (size_t s=INS+1; s < INS+DEL+1; ++s) qv[s]->at(jind) = DELETE_Q*(s-INS);
				continue;
			}

			// now we are inside the matrix
			curgs = seq[i-1];
			curgc = genome_[i+j-1]; // just a hack to avoid extra Sequence::[]
			if ( curgs != curgc ) qv[0]->at(jind) += qvec[i-1];
			
			// processing inserts
			for (size_t s=1; s <= INS; ++s) {
				if (i+s > seqsize) { continue; } // the read is over, too many inserts
				
				// case 1: quality value accumulated so far with s inserts + current mismatch value (down-right)
				int curq = qv[s]->at(jind) + (seq[i-1+s] != curgc) * qvec[i-1+s];
				
				// case 2: quality value accumulated for s-1 inserts + insert value (down)
				int prevq = qv[s-1]->at(jind) + INSERT_Q;
				
				// to get the result, we add up the corresponding values and take the logarithm
				qv[s]->at(jind) = AddTwoQualityValues(curq, prevq);
			}
			
			// processing deletes
			for (size_t s=INS+1; s < INS+DEL+1; ++s) {
				if (i+j+(s-INS) > gensize_) { qv[s]->at(jind) = VERYLARGEQ; continue; } // the genome is over, too many deletes

				// case 1: quality value accumulated so far with s-INS deletes + current mismatch value (down-right)
				int curq = qv[s]->at(jind) + (curgs != genome_[i-1+j+s-INS]) * qvec[i-1];

				// case 2: quality value accumulated for s-INS-1 deletes + delete value (right)
				int prevq = (s==INS+1) ? (qv[0]->at(jind) +   DELETE_Q) 
									   : (qv[s-1]->at(jind) + DELETE_Q);
									   
				// to get the result, we add up the corresponding values and take the logarithm
				qv[s]->at(jind) = AddTwoQualityValues(curq, prevq);
			}
		}
	}
	
	mr_.bestq_ = mr_.totalq_;
	mr_.inserts_ = 0; mr_.deletes_ = 0;
	mr_.index_ = 0;

	// now add them all up
	double res = 0;
	for (size_t i=0; i < INS+DEL+1; ++i) {
		for (size_t jind = 0; jind < qvsize_; ++jind) {
			res += pow(10, -(qv[i]->at(jind))/10.0);
			if (qv[i]->at(jind) < mr_.bestq_) {
				mr_.bestq_ = (int)(qv[i]->at(jind));
				mr_.index_ = jind; if (indices.size() > 0) mr_.index_ = indices[jind];
				if (i <= INS) { mr_.inserts_ = i; mr_.deletes_ = 0; }
				else { mr_.inserts_ = 0; mr_.deletes_ = i-INS; }
			}
		}
	}
	
	// find and pretty print the best match
	// to do so, we run dynamic programming once again (since we haven't been storing the entire matrix)
	vector<QVector *> m(seqsize+1);
	for (size_t i = 0; i < seqsize+1; ++i) {
		m[i] = new QVector(seqsize + mr_.deletes_ + 1);
		fill(m[i]->begin(), m[i]->end(), 0);
	}
	// first row
	for (size_t j = 0; j < seqsize + mr_.deletes_ + 1; ++j) {
		if (j > mr_.deletes_) { m[0]->at(j) = VERYLARGEQ; }
		else {
			m[0]->at(j) = DELETE_Q*j;
		}
	}
	// first column
	for (size_t i = 1; i < seqsize + 1; ++i) {
		if (i > mr_.inserts_) { m[i]->at(0) = VERYLARGEQ; }
		else {
			m[i]->at(0) = INSERT_Q*i;
		}
	}
	// other elements
	for (size_t i = 1; i < seqsize + 1; ++i) {
		for (size_t j = 1; j < seqsize + mr_.deletes_ + 1; ++j) {
			//cout << i << " " << j << " " << index_+j;
			if (i > mr_.inserts_ + j) { m[i]->at(j) = VERYLARGEQ; continue; }
			if (j > mr_.deletes_ + i) { m[i]->at(j) = VERYLARGEQ; continue; }
			// if (j > gensize_ - mr_.index_ + 1) { m[i]->at(j) = VERYLARGEQ; continue; }
			m[i]->at(j) = AddThreeQualityValues(
				(mr_.index_+j <= gensize_) ? (m[i-1]->at(j-1) + (genome_[mr_.index_+j-1] != seq[i-1])*qvec[i-1]) : VERYLARGEQ,
				              (mr_.deletes_ > 0) ? (m[i]->at(j-1)   + DELETE_Q)                                        : VERYLARGEQ,
				              (mr_.inserts_ > 0) ? (m[i-1]->at(j)   + INSERT_Q)                                        : VERYLARGEQ );
			//cout << " mat\n";
		}
	}
	mr_.match_.clear();
	ostringstream osm, osmr, osmp;
	size_t ifirst_ = seqsize;
	size_t j = min_element(m[ifirst_]->begin(), m[ifirst_]->end()) - m[ifirst_]->begin();
	for (size_t i = ifirst_; i >= 1; --i) {
		// go up
		if ((mr_.inserts_ > 0) && ( j == 0 || ((m[i-1]->at(j) < m[i-1]->at(j-1)) && (m[i-1]->at(j) < m[i]->at(j-1))) )) {
			mr_.match_.push_back( 2 ); osm << '-'; osmp << ' '; osmr << nucl(seq[i-1]); continue;
		}
		// or left
		if ((mr_.deletes_ > 0) && ((m[i]->at(j-1) < m[i-1]->at(j) && m[i]->at(j-1) < m[i-1]->at(j-1)) )) {
			mr_.match_.push_back( 3 );
			osm << nucl(genome_[j-1+mr_.index_]);
			osmp << ' '; osmr << '-'; ++i; --j;
			continue;
		}
		// or up-left
		mr_.match_.push_back( (int)(genome_[j-1+mr_.index_] == seq[i-1]) );
		osm << nucl(genome_[j-1+mr_.index_]); osmr << nucl(seq[i-1]);
		if (genome_[j-1+mr_.index_] == seq[i-1]) osmp << "|"; else osmp << " "; --j;
	}
	// the first character is always here
	
	reverse(mr_.match_.begin(), mr_.match_.end());
	mr_.matchreadstr_ = osmr.str(); reverse(mr_.matchreadstr_.begin(), mr_.matchreadstr_.end());
	mr_.matchstr_ = osm.str(); reverse(mr_.matchstr_.begin(), mr_.matchstr_.end());
	mr_.matchprettystr_ = osmp.str(); reverse(mr_.matchprettystr_.begin(), mr_.matchprettystr_.end());
	
	for (size_t i = 0; i < seqsize+1; ++i) {
		m[i]->clear(); delete m[i];
	}
	m.clear();
	for (size_t i=0; i < INS+DEL+1; ++i) {
		qv[i]->clear();  delete qv[i];
	}
	qv.clear();
	mr_.prob_ = res;
	return mr_;
}

MatchResults BayesQualityGenome::ReadBQPreprocessed(const Read & r, size_t readno, size_t readssize, const BowtieResults & bwres) {
	QVector q(r.getQualityString().size());
	copy(r.getQualityString().begin(), r.getQualityString().end(), q.begin());
	Sequence s = r.getSequence();
	Sequence cs = !s;
	
	vector<size_t> bwmatches, cbwmatches;
	
	#ifdef USE_BOWTIE
	size_t rdind = bwres.size();
	for (size_t i=0; i<bwres.size(); ++i) {
		if (bwres[i].size() > 0 && bwres[i][0][BWMAP_IND_RDNO] == readno) {
			rdind = i; break;
		}
	}
	if (rdind < bwres.size()) {
		for (size_t j=0; j<bwres[rdind].size(); ++j) {
			// INFO(bwres[rdind][j][BWMAP_IND_REVC] << " " << bwres[rdind][j][BWMAP_IND_GIND] << " " << bwres[rdind][j][BWMAP_IND_RIND]);
			if (bwres[rdind][j][BWMAP_IND_REVC] == 0) {
				for (size_t match = bwres[rdind][j][BWMAP_IND_GIND]-bwres[rdind][j][BWMAP_IND_RIND] - INS;
						match < bwres[rdind][j][BWMAP_IND_GIND]-bwres[rdind][j][BWMAP_IND_RIND] + DEL + 1;
						++match)
				bwmatches.push_back( match );

				bwmatches.push_back(bwres[rdind][j][BWMAP_IND_GIND] - bwres[rdind][j][BWMAP_IND_RIND] );
				// INFO(r.getSequence());
				// INFO(genome_.Subseq(bwres[rdind][j][BWMAP_IND_GIND] - bwres[rdind][j][BWMAP_IND_RIND],
				//				 	bwres[rdind][j][BWMAP_IND_GIND] - bwres[rdind][j][BWMAP_IND_RIND] + r.getSequence().size()));
			} else {
				for (size_t match = bwres[rdind][j][BWMAP_IND_GIND]-(s.size() - BOWTIE_SEED_LENGTH - bwres[rdind][j][BWMAP_IND_RIND]) - INS;
							match < bwres[rdind][j][BWMAP_IND_GIND]-(s.size() - BOWTIE_SEED_LENGTH - bwres[rdind][j][BWMAP_IND_RIND]) + DEL + 1;
							++match)
					cbwmatches.push_back( match );
				// INFO(!(r.getSequence()));
				// INFO(genome_.Subseq( (bwres[rdind][j][BWMAP_IND_GIND] - (s.size() - BOWTIE_SEED_LENGTH - bwres[rdind][j][BWMAP_IND_RIND]) ),
				//					  (bwres[rdind][j][BWMAP_IND_GIND] - (s.size() - BOWTIE_SEED_LENGTH - bwres[rdind][j][BWMAP_IND_RIND]) + r.getSequence().size())));
			}
		}
	}
	sort( bwmatches.begin(),  bwmatches.end());  bwmatches.erase(unique(  bwmatches.begin(),   bwmatches.end()),  bwmatches.end());
	sort(cbwmatches.begin(), cbwmatches.end()); cbwmatches.erase(unique( cbwmatches.begin(),  cbwmatches.end()), cbwmatches.end());
	#endif

	MatchResults mr1, mr2;

	#ifdef USE_BOWTIE
	if (bwmatches.size() > 0 || cbwmatches.size() == 0) {
		mr1 = ProcessOneReadBQ(s, q, readno, bwmatches);
		//INFO(s);
		INFO("Read as it is gives result " << mr1.prob_);
		if (mr1.bestq_ <= GOOD_ENOUGH_Q) {
			return mr1;
		}
	} else INFO("No matches for read as it is.");
	#else
	mr1 = ProcessOneReadBQ(s, q, readno, bwmatches);
	//INFO(s);
	INFO("Read as it is gives result " << mr1.prob_);
	if (mr1.bestq_ <= GOOD_ENOUGH_Q) {
		return mr1;
	}
	#endif

	
	#ifdef USE_BOWTIE
	if (cbwmatches.size() > 0 || bwmatches.size() == 0) {
		mr2 = ProcessOneReadBQ(!s, q, readno + readssize, cbwmatches);
		//INFO(!s);
		INFO("Complementary read gives result " << mr2.prob_);
		if (mr2.prob_ > mr1.prob_) {
			return mr2;
		}
		else return mr1;
	} else {
		INFO("No matches for complementary read.");
		return mr1;
	}
	#else 
	mr2 = ProcessOneReadBQ(!s, q, readno + readssize, cbwmatches);
	//INFO(!s);
	INFO("Complementary read gives result " << mr2.prob_);
	if (mr2.prob_ > mr1.prob_) {
		return mr2;
	}
	else return mr1;
	#endif
}

void BayesQualityGenome::PreprocessOneMapOneReadWithShift(const PSeq &ps, PreprocReadsMap & rm, int mapno, PreprocMap & map, size_t readno, size_t reads_size, size_t shift) {	
	
	// see if we can find the same
	PreprocReadsMap::iterator itrm = rm.find(ps);
	if (itrm != rm.end()) {  // found it!
		if (itrm->second.second == 0) { // found from heads
			for (PreprocMap::iterator it = map.begin(); it != map.end(); ++it) {
				it->second[readno] = availableReadsHead_[ps][itrm->second.first];
			}
		}
		else if (itrm->second.second == 1) { // found from mids
			for (PreprocMap::iterator it = map.begin(); it != map.end(); ++it) {
				it->second[readno] = availableReadsMid_[ps][itrm->second.first];
			}
		}
	}

	// if we cannot, we add ps to rm
	rm.insert(make_pair(ps, make_pair(readno, mapno)));

	for (PreprocMap::iterator it = map.begin(); it != map.end(); ++it) {
		int inters = 0;
		for (size_t k = 0; k < PREPROCESS_SEQ_LENGTH; ++k) {
			if (it->first[k] == ps[k]) {
				++inters;
				if (inters >= NEED_INTERSECTION) { it->second[readno] = true; break; }
			}
		}
	}
}

void BayesQualityGenome::PreprocessReads(const vector<Read> &reads) {
	availableReadsHead_.clear();
	availableReadsMid_.clear();
	fillAvailableReads(2*reads.size(), "");
	readMids.clear(); for (size_t i = 0; i < 2*reads.size(); ++i) readMids.push_back(0);
	
	// this map serves as an example to make preprocessing more efficient
	PreprocReadsMap rm;

	for (size_t i = 0; i < reads.size(); ++i) {
		Sequence s = reads[i].getSequence();
		size_t mid = s.size() / 2;
		if (s.size() - mid < PREPROCESS_SEQ_LENGTH) mid = s.size() - PREPROCESS_SEQ_LENGTH;
		readMids[i] = mid; readMids[i+reads.size()] = mid;
		
		PSeq   ps(s.Subseq(0, PREPROCESS_SEQ_LENGTH));
		PSeq  cps((!s).Subseq(0, PREPROCESS_SEQ_LENGTH));
		PSeq  psm(s.Subseq(mid, mid+PREPROCESS_SEQ_LENGTH));
		PSeq cpsm((!s).Subseq(mid, mid+PREPROCESS_SEQ_LENGTH));

		PreprocessOneMapOneReadWithShift( ps, rm, 0, availableReadsHead_, i, 			  reads.size(), 0);
		PreprocessOneMapOneReadWithShift(cps, rm, 0, availableReadsHead_, i+reads.size(), reads.size(), 0);
				
		PreprocessOneMapOneReadWithShift( psm, rm, 1, availableReadsMid_, i, 			  reads.size(), mid);
		PreprocessOneMapOneReadWithShift(cpsm, rm, 1, availableReadsMid_, i+reads.size(), reads.size(), mid);
	}
}

void BayesQualityGenome::ProcessReads(const vector<Read> &reads) {
	#ifdef USE_PREPROCESSING
	INFO("Preprocessing...");
	PreprocessReads(reads);
	INFO("  ...done.");
	#endif

	#ifdef WRITE_STATSFILE
	ofstream os(STATSFILENAME, ios::out);
	if (os) {
		os << "# The order of reads matches the input file" << endl;
		os << "#          sumq  bestloc  bestq  replacements  inserts  dels  matchstring" << endl;
	} else {
		cerr << "Cannot write to statistics file\n";
		abort();
	}
	#endif

	for (size_t i=0; i<reads.size(); ++i) {
	
		BowtieResults bwres;
		MatchResults mr = ReadBQPreprocessed(reads[i], i, BIGREADNO, bwres);
		# pragma omp critical
		{
			INFO(mr.matchreadstr_);
			INFO(mr.matchprettystr_);
			INFO(mr.matchstr_);
			ostringstream m; for (size_t j=0; j< mr.match_.size(); ++j) m << mr.match_[j]; INFO(m.str());
			INFO(mr.prob_ << ", best: " << mr.bestq_ << "/" << mr.totalq_ << " at " << mr.index_ << " with " << mr.inserts_ << " inserts and " << mr.deletes_ << " deletes");
			#ifdef WRITE_STATSFILE
			size_t replacements = 0; for (size_t j=0; j< mr.match_.size(); ++j) if (mr.match_[j] == 0) ++replacements;
			os << setw(7) << i << ":" << setw(15) << mr.prob_ << "\t" << mr.index_ << "\t" << mr.bestq_ << "\t" << replacements << "\t" << mr.inserts_ << "\t" << mr.deletes_ << "\t" << m.str().data() << endl;
			#endif
		}

		/*double res = ReadBQPreprocessed(reads[i], i, reads.size());

		INFO(LastMatchReadString());
		INFO(LastMatchPrettyString());
		INFO(LastMatchString());
		ostringstream m; for (size_t i=0; i< LastMatch().size(); ++i) m << LastMatch()[i]; INFO(m.str());
		INFO(res << ", best: " << LastMatchQ() << "/" << LastTotalQ() << " at " << LastMatchIndex() << " with " << LastMatchInserts() << " inserts and " << LastMatchDeletes() << " deletes");
		#ifdef WRITE_STATSFILE
		size_t replacements = 0; for (size_t i=0; i< LastMatch().size(); ++i) if (LastMatch()[i] == 0) ++replacements;
		os << setw(15) << res << "\t" << LastMatchIndex() << "\t" << LastMatchQ() << "\t" << replacements << "\t" << LastMatchInserts() << "\t" << LastMatchDeletes() << "\t" << m.str().data() << endl;
		#endif*/
	}
	
	#ifdef WRITE_STATSFILE
	os.close();
	#endif
}

void BayesQualityGenome::ProcessReads(const char *filename) {
	#ifdef USE_PREPROCESSING
	INFO("Preprocessing...");
	PreprocessReads(reads);
	INFO("  ...done.");
	#endif

	#ifdef WRITE_STATSFILE
	ofstream os(STATSFILENAME, ios::out);
	if (os) {
		os << "# The order of reads matches the input file" << endl;
		os << "#          sumq  bestloc  bestq  replacements  inserts  dels  matchstring" << endl;
	} else {
		cerr << "Cannot write to statistics file\n";
		abort();
	}
	#endif

	ireadstream ifs(filename);
	assert(ifs.is_open());
	Read r;
	int batch_no = 1;
	size_t readno = 0;
	// skip several reads -- needed to continue the runs
	while (!ifs.eof()) {
		if (readno > SKIP_READS) break;
		ifs >> r;
		r.trimNs();
		++readno;
	}
	
	while (!ifs.eof()) {
		// fill up an array of reads
		vector<Read> v;
		for (size_t i=0; i<READ_BATCH; ++i) {
			ifs >> r;
			r.trimNs();
			v.push_back(r);
			if (ifs.eof()) break;
		}
		INFO("Batch " << batch_no << " of size " << v.size() << " is ready.");
		++batch_no;
		
		#ifdef USE_BOWTIE
		INFO("writing read parts");
		string tmpname = writeReadPartsForBowtie(v, readno);
		INFO(tmpname);
		char *tmpname2 = strdup("/tmp/tmpfileXXXXXX");
		int fd2 = mkstemp(tmpname2);
		assert(fd2>=0);
		string cmd = bowtieCmd_ + " " + tmpname + " > " + tmpname2;
		INFO(cmd.data());
		system(cmd.data());
		BowtieResults bwres = readBowtieResults(fd2, v.size(), readno);
		bwres_ = bwres;
		cmd = "rm -rf "; cmd.append(tmpname);
		system(cmd.data());
		cmd = "rm -rf "; cmd.append(tmpname2);
		system(cmd.data());
		#endif
		
		vector<MatchResults> mrv(v.size());
		
		INFO("Entering parallel section");
		omp_set_num_threads(THREADS_NUM);
		#pragma omp parallel for shared(os, readno, mrv, bwres) private(r)
		for (int i=0; i<v.size(); ++i) {
			r = v[i];
			if (r.size() < MIN_READ_SIZE) {
				mrv[i].prob_ = -1;
				continue;
			}
		
			INFO("Hello from thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " read no " << readno+i);
			MatchResults mr = ReadBQPreprocessed(r, readno+i, BIGREADNO, bwres);
			mrv[i] = mr;

			INFO(mr.matchreadstr_);
			INFO(mr.matchprettystr_);
			INFO(mr.matchstr_);
			ostringstream m; for (size_t j=0; j< mr.match_.size(); ++j) m << mr.match_[j]; INFO(m.str());
			INFO(mr.prob_ << ", best: " << mr.bestq_ << "/" << mr.totalq_ << " at " << mr.index_ << " with " << mr.inserts_ << " inserts and " << mr.deletes_ << " deletes");
		}
		
		#ifdef WRITE_STATSFILE
		for (size_t i=0; i<v.size(); ++i) {
			if (mrv[i].prob_ < 0) {
				os << setw(7) << readno+i << " skipped: only " << setw(2) << r.size() << " known letters at the beginning" << endl;
				continue;
			}
			ostringstream m; for (size_t j=0; j< mrv[i].match_.size(); ++j) m << mrv[i].match_[j];
			size_t replacements = 0; for (size_t j=0; j< mrv[i].match_.size(); ++j) if (mrv[i].match_[i] == 0) ++replacements;
			os << setw(7) << readno+i << ":" << setw(15) << mrv[i].prob_ << "\t" << mrv[i].index_ << "\t" << mrv[i].bestq_ << "\t" << replacements << "\t" << mrv[i].inserts_ << "\t" << mrv[i].deletes_ << "\t" << m.str().data() << endl;	
		}
		#endif
		
		readno += v.size();
	}
	#ifdef WRITE_STATSFILE
	os.close();
	#endif
}

float BayesQualityGenome::AddTwoQualityValues(float curq, float prevq) {
	// there are a few shortcuts here
	if ( curq < prevq - 1 ) return curq;
	else if ( curq > prevq + 1 ) return prevq;
	else if ( curq == prevq + 1 ) return prevq - TENLOGELEVENTENTH;
	else if ( curq == prevq - 1 ) return curq - TENLOGELEVENTENTH;
	else if ( curq == prevq ) return prevq - TENLOGTWO;
	return min(curq, prevq);
}

int BayesQualityGenome::AddTwoQualityValues(int curq, int prevq) {
	// there are a few shortcuts here
	if ( curq < prevq - 1 ) return curq;
	else if ( curq > prevq + 1 ) return prevq;
	else if ( curq == prevq + 1 ) return prevq - TENLOGELEVENTENTH;
	else if ( curq == prevq - 1 ) return curq - TENLOGELEVENTENTH;
	else if ( curq == prevq ) return prevq - TENLOGTWO;
	return min(curq, prevq);
}

float BayesQualityGenome::AddThreeQualityValues(float diagq, float leftq, float rightq) {
	return AddTwoQualityValues(diagq, AddTwoQualityValues(leftq, rightq));
}

int BayesQualityGenome::AddThreeQualityValues(int diagq, int leftq, int rightq) {
	return AddTwoQualityValues(diagq, AddTwoQualityValues(leftq, rightq));
}



void BayesQualityGenome::fillAvailableReads(size_t vecSize, const string & s) {
	if (s.size() == PREPROCESS_SEQ_LENGTH) {
		vector<bool> newVec(vecSize);
		fill(newVec.begin(), newVec.end(), false);
		availableReadsHead_.insert(make_pair(PSeq(s), newVec));
		availableReadsMid_.insert(make_pair(PSeq(s), newVec));
		return;
	}
	fillAvailableReads(vecSize, (s + 'A'));
	fillAvailableReads(vecSize, (s + 'C'));
	fillAvailableReads(vecSize, (s + 'G'));
	fillAvailableReads(vecSize, (s + 'T'));	
}

void BayesQualityGenome::writeReadPartsForBowtie(const vector<Read> & v, int fd, size_t readno) {
	// ofstream os; os.attach(fd);
	FILE *f = fdopen(fd, "w");
	for (size_t i = 0; i < v.size(); ++i) {
		if (v[i].size() < MIN_READ_SIZE) continue;
		for (size_t j = 0; j < BOWTIE_SEED_NUM; ++j) {
			size_t ind = (j * (v[i].size() - BOWTIE_SEED_LENGTH) ) / BOWTIE_SEED_NUM;
			if (ind + BOWTIE_SEED_LENGTH > v[i].size()) ind = v[i].size() - BOWTIE_SEED_LENGTH;
			string qual = v[i].getQualityString().substr(ind, BOWTIE_SEED_LENGTH);
			for (size_t k=0; k<qual.size(); ++k) {
				qual[k] = (char)(qual[k]+Read::PHRED_OFFSET);
			}
			fprintf(f, "@%u %u\n%s\n+\n%s\n", i+readno, ind, v[i].getSequenceString().substr(ind, BOWTIE_SEED_LENGTH).data(), qual.data());
			// os << "@" << i+readno << " " << ind << endl <<  << endl << "+" << endl;
			// os << endl;
		}
	}
	fclose(f);
}

BowtieResults BayesQualityGenome::readBowtieResults(int fd, size_t noofreads, size_t readno) {
	//ifstream is; is.attach(fd);
	FILE *f = fdopen(fd, "r");
	BowtieResults res(noofreads);
	char buf[1024];
	char rev;
	unsigned int i, j, ind;
	while (!feof(f)) {
		if ( fgets(buf, 1024, f) == NULL && feof(f) ) break; 
		if (strlen(buf) == 0) break;
		sscanf(buf, "%u %u\t%c\t%*s\t%u\t%*s\n", &i, &j, &rev, &ind);
		vector<size_t> curv(BWMAP_IND_GIND+1);
		if (rev == '+') curv[BWMAP_IND_REVC] = 0; else curv[BWMAP_IND_REVC] = 1;
		curv[BWMAP_IND_RDNO] = i;
		curv[BWMAP_IND_RIND] = j;
		curv[BWMAP_IND_GIND] = ind;
		res[i-readno].push_back(curv);
	}
	fclose(f);
	
	return res;
}


}

