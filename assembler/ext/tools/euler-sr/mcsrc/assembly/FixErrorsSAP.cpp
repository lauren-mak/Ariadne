/***************************************************************************
 * Title:          FixErrorsSAP.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "hash/HashUtils.h"
#include "ListSpectrum.h"
#include "BitSpectrum.h"
#include "MultTuple.h"
#include "NumericTuple.h"
#include "IntegralTupleStatic.h"


#include <vector>
#include <iostream>
//#include <hash_map>
#include <ext/hash_map>
#include <map>


class Stats {
public:
  ssize_t numEdge;
  ssize_t numIns;
  ssize_t numDel;
  ssize_t numMut;
  ssize_t numNoSolid;
  ssize_t numNoPathFound;
  ssize_t numMultiplePaths;
  ssize_t numErrorAtEnd;
  Stats() {
    Reset();
  }
  void Reset() {
    numEdge = numIns = numDel = numMut = 0;
    numNoSolid = 0;
    numNoPathFound = 0;
    numMultiplePaths = 0;
		numErrorAtEnd = 0;
  }
	Stats &Append(const Stats &s) {
    numEdge += s.numEdge;
    numIns += s.numIns;
    numDel += s.numDel;
    numMut += s.numMut;
    numNoSolid += s.numNoSolid;
    numNoPathFound  += s.numNoPathFound;
    numMultiplePaths += s.numMultiplePaths;
		numErrorAtEnd += s.numErrorAtEnd;
    return *this;
	}
  Stats &operator+=(const Stats &s) {
		Append(s);
		return *this;
  }
	friend std::ostream &operator<<(std::ostream &strm, const Stats &s) {
		strm << s.numEdge <<" " << s.numIns 
				 <<" " << s.numDel << " " 
				 << s.numMut <<" " 
				 << s.numNoSolid << " " << s.numNoPathFound 
				 << " " << s.numMultiplePaths << std::endl;
		return strm;
	}
};

class SEdge {
public:
  ssize_t prevNuc;
  char prevLevel;
  ssize_t prevPosition;
	MultTuple tuple;
  ssize_t score;
  static int tupleSize;
  SEdge() {
    prevNuc = 0;
    prevLevel = 0;
    prevPosition = 0;
    tuple = MultTuple("");
    score = 0;
  }
  SEdge& operator=(const SEdge& e) {
		if (this != &e) {
			prevNuc = e.prevNuc;
			prevLevel = e.prevLevel;
			prevPosition = e.prevPosition;
			CopyTuple(e.tuple);
			score = e.score;
		}
    return *this;
  }
  void CopyTuple(const std::string &t) {
    // in case the type of t changes, encapsulate the copy function
    tuple = t;
  }
};

int SEdge::tupleSize = 0;
class FixParams {
public:
  ssize_t maxGap;
  ssize_t gapOpen;
  ssize_t gapExtend;
  ssize_t misMatch;
  ssize_t scoreThreshold;
  ssize_t span;
  ssize_t maxTrim;
  ssize_t edgeLimit;
};

typedef std::map<MultTuple, SEdge> Cell;
typedef Cell::iterator EdgeIterator;


typedef Cell* Column;
typedef Column* Grid;
typedef std::vector<Grid> Cube;
void ReverseSeq(DNASequence &seq, DNASequence &rev);
void PatchSeq(DNASequence &seq, ssize_t pos, DNASequence &patch, ssize_t replaceLength);
void CreateGrid(ssize_t dim, Grid& grid);
void DeleteGrid(ssize_t dim, Grid &grid);

template <typename T_Spectrum>
ssize_t FindSolidPosition(DNASequence &seq, T_Spectrum &spectrum, ssize_t span, 
											ssize_t &lastSolidPos, ssize_t& pos);


void StoreEdge(ssize_t prevPos, ssize_t prevLevel, ssize_t prevNuc, MultTuple prevTuple, 
							Cell& cell, MultTuple tuple, ssize_t score);

ssize_t FindMinimumScoreEdge(Cube &matrix, FixParams &p, EdgeIterator &edgeIt, 
												 ssize_t &minN, ssize_t &minK, ssize_t& numMin, ssize_t &minScore, 
												 MultTuple &minTuple, ssize_t &minTrim, ssize_t &trim);

void Backtrack(Cube &matrix, ssize_t pos, ssize_t level, ssize_t nuc, MultTuple tuple, 
							DNASequence &seq, Stats &fixStats, FixParams &params,
							ssize_t &firstEdit, ssize_t &end);

template <typename T_Spectrum>
ssize_t SolidifyRead(DNASequence &seq, T_Spectrum &spectrum, 
								 FixParams &params, Stats &fixStats);

template<typename T_Spectrum>
ssize_t SolidifyUntilFixed(DNASequence &seq, T_Spectrum &spectrum, 
											 int tupleSize, FixParams &params,
											 DNASequence &fixedSeq, ssize_t &replaceLength, 
											 Stats &fixStats,
											 ssize_t &lastEdit, ssize_t &end);

ssize_t IsValid(MultTuple &tuple);

IntMatrix ScoreMat;

void PrintUsage() {
	std::cout << "usage: fixErrorsSAP readsFile spectrumFile tupleSize outputFile [options] " << std::endl;
	std::cout << "options: " << std::endl;
	std::cout << " -minMult   The minimum multiplicity k-mer to consider solid " << std::endl
						<< " -maxGap    The maximum gap to consider (insertions or deletions). " << std::endl
						<< "           Increasing maxGap may substantially increase the run time. " << std::endl
						<< " -gapOpen   Not implemented. " << std::endl
						<< " -gapExtend (1) Cost to extend a gap " << std::endl
						<< " -misMatch  (1) Cost to mutate a nucleotide " << std::endl
						<< " -span n    minimum span needed for solid region" << std::endl
						<< " -discardFile \"\" Name of file to output unfixable reads " << std::endl
						<< "             null value means keep discards in original file " << std::endl
						<< " -maxScore  (infinity) Maximum score to permit extension of graph " << std::endl
						<< " -maxTrim   The maximum amount to clip off the ends of reads when the fix is " 
						<< std::endl
						<< " -edgeLimit Only find fixes if they are more than edgeLimit away from the edge" 
						<< std::endl
						<< " -readFixFile filename  Print the numbers of fixes (indel, mutation) to filename " 
						<< std::endl
						<< "            ambiguous or impossible at the end, but complete in the middle. " 
						<< std::endl
						<< " -startScore score (3) Start fixing reads with 'score' max score.  This "<< std::endl
						<< "              helps fixing reads that only have 1 or 2 errors." << std::endl
						<< " -stepScore  step (2) Fix reads starting at 'startScore', and increase " << std::endl
						<< "              the search space by 'step' mutations until 'maxScore' is reached."
						<< std::endl;
	std::cout << " -spectrumType [full|concise|numeric] The type of file that the spectrum is."<<std::endl
						<< "                Currently there is no auto-detect, and so if a concise"<<std::endl
						<< "                spectrum is used it must be specified here" << std::endl
						<< " -lockFile file Wait on lock 'file' to read the spectrum." << std::endl
						<< "                This is helpful when many processes are reading the."<<std::endl
						<< "                same file at the same time off an NFS." << std::endl
						<< "                Not implemented." << std::endl; // TODO: Why is it here but not implemented?
	


}
int main(int argc, char* argv[]) {
  if (argc < 4) {
    PrintUsage();
    exit(0);
  }
  int tupleSize;
  int argi = 1;
  std::string readsFile        = argv[argi++];
  std::string spectrumFileName = argv[argi++];
  tupleSize = atoi(argv[argi++]);
  std::string outputFileName   = argv[argi++];
  std::string discardFileName  = "";
  std::string readFixFile      = "";
	std::string spectrumType = "full";
  FixParams params;
  params.gapOpen = 4;
  params.gapExtend = 1;
  params.maxGap = 7;
  params.span = 4;
  params.maxTrim = 7;
  params.edgeLimit = 3;
  params.scoreThreshold = 100000;

  ssize_t minMult;
  minMult = 20;
  params.misMatch = 1;
  ssize_t startScore = 3;
  ssize_t stepScore  = 2;
  while (argi < argc) {
    if (strcmp(argv[argi], "-startScore") == 0) {
      startScore = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-stepScore") == 0) {
      stepScore = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-maxTrim") == 0 ) {
      ++argi;
      params.maxTrim = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-minMult") == 0 ) {
      ++argi;
      minMult = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-span") == 0 ) {
      ++argi;
      params.span = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-maxGap") == 0) {
      ++argi;
      params.maxGap = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-gapOpen") == 0 ) {
      ++argi;
      params.gapOpen = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-gapExtend") == 0 ) {
      ++argi;
      params.gapExtend = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-misMatch") == 0 ) {
      ++argi;
      params.misMatch = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-maxScore") == 0 ) {
      ++argi;
      params.scoreThreshold = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-discardFile") == 0 ) {
      ++argi;
      discardFileName = argv[argi];
    }
    else if (strcmp(argv[argi], "-edgeLimit") == 0 ) {
      ++argi;
      params.edgeLimit = atoi(argv[argi]);
    }
    else if (strcmp(argv[argi], "-readFixFile") == 0 ) {
      ++argi;
      readFixFile = argv[argi];
    }
		else if (strcmp(argv[argi], "-spectrumType") == 0) {
			spectrumType = argv[++argi];
		}
    else {
      PrintUsage();
      std::cout << "Bad option: " << argv[argi] << std::endl;
			exit(1);
    }
    ++argi;
  }

	if (spectrumType != "full" and 
			spectrumType != "concise" and
			spectrumType != "numeric") {
		PrintUsage();
		std::cout << "Invalid spectrum type: " << spectrumType << std::endl;
		exit(1);
	}


	std::string reportFileName = FormReportName(readsFile);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);



  ssize_t i, j;
  
  CreateMatrix(ScoreMat, 5, 5);
  for (i = 0; i < 4; i++ ) {
    for (j = 0; j < 4; j++ ) {
      if (i != j) {
				ScoreMat[i][j] = params.misMatch;
      }
      else {
				ScoreMat[i][j] = 0;
      }
    }
  }
  for (i = 0; i < 4; i++) {
    ScoreMat[i][4] = params.misMatch;
    ScoreMat[4][i] = params.misMatch;
  }
  

  std::vector<DNASequence*> reads;
  std::ifstream seqIn;
  openck(readsFile, seqIn, std::ios::in, report);

  std::ofstream seqOut, discardOut, readFixOut;
  openck(outputFileName, seqOut, std::ios::out, report);

  if (discardFileName != "") 
    openck(discardFileName, discardOut, std::ios::out, report);
  if (readFixFile != "") 
    openck(readFixFile, readFixOut, std::ios::out, report);

  std::cout << "getting reads "; std::cout.flush();
  DNASequence *read;
  while (SeqReader::GetSeq(seqIn, read, SeqReader::noConvert)) {
    reads.push_back(read);
  }
  std::cout << "done." << reads.size() << std::endl;
  std::cout << "reading tuples >= " << minMult << " ... "; std::cout.flush();
  ListSpectrum<MultTuple> spectrum;

	ListSpectrum<NumericTuple> numericSpectrum;
	// Configure the bit spectrum
	
	BitSpectrum<MultTuple> bitSpectrum(tupleSize);
	bitSpectrum.FindOnlySolid();

	//	Spectrum<MultTuple>* spectrumPtr;


	if (spectrumType == "full") {
		spectrum.tupleSize = tupleSize;
		//		spectrumPtr = &spectrum;
		spectrum.Read(spectrumFileName, minMult);
	}
	else if (spectrumType == "concise") {
		//		spectrumPtr = &bitSpectrum;
		bitSpectrum.Read(spectrumFileName, minMult);
	}
	else if (spectrumType == "numeric") {
		//		spectrumPtr = &numericSpectrum;
		numericSpectrum.Read(spectrumFileName, minMult);
	}

	/*	spectrumPtr->Read(spectrumFileName, minMult);
	if (spectrumPtr->size() == 0) {
		std::cout << "WARNING: no tuples above " << minMult<< " were found." << std::endl;
		std::cout << "    Error correction is not possible." << std::endl;
		exit(1);
	}
  std::cout << "done. " << spectrumPtr->size() << std::endl;
	*/
	SEdge::tupleSize = tupleSize;
	//	Tuple::tupleSize = tupleSize;
	Tuple::SetTupleSize(tupleSize);

  //UNUSED// ssize_t s;
  //UNUSED// HashValue hashValue;

  Stats stats, totalStats;
  // Now try to fix all reads.
  ssize_t r;
  ssize_t numFixed = 0;
  std::vector<char> fixed;
  fixed.resize(reads.size());
  for (r = 0; r < reads.size(); r++) {
    fixed[r] = 0;
  }
  ssize_t maxScore = params.scoreThreshold;
  ssize_t score;
  //UNUSED// ssize_t step;
  ssize_t nFixed; 
	ssize_t readFixed = 0;
  for (score = startScore; score <= maxScore; score+= stepScore ) {
    params.scoreThreshold = score;
    nFixed = 0; 
    for (r = 0; r < reads.size(); r++) {
      stats.Reset();
			
      if (r % 1000 == 999) {
				std::cout << ".";
				std::cout.flush();
      }
      if (r % 50000 == 49999 )
				std::cout << " " << totalStats ;
						 
      if (fixed[r])
				continue;
			readFixed = 0;
			if (spectrumType == "full" and SolidifyRead(*reads[r], spectrum, params, stats)) 
				readFixed = 1;
			if (spectrumType == "concise" and SolidifyRead(*reads[r], bitSpectrum, params, stats)) 
				readFixed = 1;
			if (spectrumType == "numeric" and SolidifyRead(*reads[r], numericSpectrum, params, stats)) 
				readFixed = 1;
			
			if (readFixed) {
				fixed[r] = 1;
        nFixed++;
			}
			totalStats.Append(stats);
    }
    std::cout << "threshold: " << score << " fixed: " << nFixed << std::endl;
  }
  
  for (r = 0; r < reads.size(); r++ ) {
    if (fixed[r] == 0) {
      if (discardFileName != "") {
				reads[r]->_ascii = 1;
				reads[r]->PrintSeq(discardOut);
				discardOut << std::endl;
      }
    }
    else {
      reads[r]->PrintSeq(seqOut);
      seqOut << std::endl;
			if (readFixFile != "" ) {
				readFixOut << reads[r]->namestr << " ";
				readFixOut << stats.numEdge << " edges " << stats.numMut << " mutations " << stats.numIns 
									 << " insertions and " << stats.numDel << " deletions." << std::endl;
			}
      totalStats += stats;
      numFixed++;
    }
  }
  std::cout << "------------------------------------------------------------" << std::endl;
  std::cout << "fixed: " << numFixed << std::endl;
  std::cout << "# mutations:  " << totalStats.numMut << std::endl;
  std::cout << "# insertions: " << totalStats.numIns << std::endl;
  std::cout << "# deletions:  " << totalStats.numDel << std::endl;
  std::cout << "Sequences were not fixed because of: " << std::endl;
  std::cout << "  " << totalStats.numNoSolid << " sequences did not have any solid tuples." << std::endl;
  std::cout << "  " << totalStats.numMultiplePaths << " sequences had multiple paths." << std::endl;
  std::cout << "  " << totalStats.numNoPathFound << " sequences had no valid support paths "<< std::endl;
	std::cout << "  " << totalStats.numErrorAtEnd << " could not be fixed due to an error within " << params.edgeLimit << " of an end "
						<< std::endl;
  for (r = 0; r < reads.size(); r++ ) {
    reads[r]->Reset();
    delete reads[r];
  }

	EndReport(report);
	report.close();
	return 0;
}


void CreateGrid(ssize_t dim, Grid& grid) {
  grid = new Column[dim+1];
  ssize_t d;
  for (d = 0; d <= dim; d++ ) {
    grid[d] = new Cell[4];
  }
}

void DeleteGrid(ssize_t dim, Grid &grid) {
  ssize_t d;
  for (d = 0; d <= dim; d++ ) {
    delete[] grid[d];
  }
  delete[] grid;
}

ssize_t IsValid(Tuple &tuple) {
	ssize_t i;
	for (i = 0; i < tuple.size(); i++) {
		if (numeric_nuc_index[(unsigned char) tuple[i]] >= 4) 
			return 0;
	}
	return 1;
}


template <typename T_Spectrum>
ssize_t FindSolidPosition(DNASequence &seq, T_Spectrum &spectrum, ssize_t span,
											ssize_t &lastSolidPos, ssize_t& pos) {
  MultTuple tuple;
  ssize_t solidSpanFound = 0;
  lastSolidPos = -1;
#ifdef VERBOSE
  std::cout << "Finding a solid position " << std::endl;
#endif
  for (pos = 0; pos < seq.length - spectrum.tupleSize + 1; pos++) {
		tuple.assign((char*) &(seq.seq[pos]));
		//    if (GetHashValue(seq, pos, spectrum.tupleSize, tuple)) {
		if (tuple.Valid()) {
      if (spectrum.FindTuple(tuple) != -1) {
				// found a solid tuple
				if (lastSolidPos == -1) {
#ifdef VERBOSE
					std::cout << "found solid at : " << pos << " " << tuple << std::endl;
#endif
					lastSolidPos = pos;
				}
				else {
#ifdef VERBOSE
					std::cout << "still solid at: " << pos << " " << tuple << std::endl;
#endif
					// if the span is long enough, return this position
					if (pos - lastSolidPos >= span) {
						solidSpanFound = 1;
					}
					// oherwise, check more tuples to see if they are solid
				}
      }
      else {
				// found an erroneous tuple
//				std::cout << "found an error at: " << pos << std::endl;
				if (lastSolidPos > -1 and pos - lastSolidPos - 1 >= span) {
					// If this tuple is erroneous, but we are currently in a valid span
					// of solid tuples, return the previous position, as this does not have an
					// error.
					pos -= 1;
					return 1;
				}
				else {
					//	  std::cout << "found an error at " << pos << " and last solid at: " << lastSolidPos << std::endl;
					lastSolidPos = -1;
				}
      }
    } else {
      // found an invalid (masked, etc.) tuple, keep looking
      if (lastSolidPos > -1 and pos - lastSolidPos -1 >= span) {
				pos -=1;
				return 1;
      }
      else 
				lastSolidPos = -1;
    }
  }
	pos--;
  return solidSpanFound;
}

template <typename T_Spectrum>
ssize_t SolidifyUntilFixed(DNASequence &seq, T_Spectrum &spectrum, 
											 FixParams &params, DNASequence &fixedSeq, 
											 ssize_t &replaceLength, Stats &stats, ssize_t &lastEdit, ssize_t &end) {

  stats.Reset();
  // the reads starts having errorst after solidpos. fix them.
  ssize_t startPos;
  ssize_t fixed;
  ssize_t pos; // for iterating overpositions
  ssize_t k; // for iterating over gaps
  ssize_t n; // for iterating over nucleotides
  fixed = 0;
  ssize_t i;
  replaceLength = 0;
  Grid grid;
	#ifdef VERBOSE
	std::cout << "fixing: " << std::endl;
	seq._ascii = 1;
	seq.PrintSeq(std::cout);
	std::cout << std::endl;
	#endif
  // start at 'pos' and fix until a tuple is found that corresponds
  // to the read
  ssize_t init = 0;
  SEdge edge;
  MultTuple tuple;
	tuple.assign((char*)seq.seq);
  startPos = pos = 0 ;
  Cube cube;
  ssize_t newEdgeCreated = 1;
  MultTuple newTuple;
  newTuple.reserve(spectrum.tupleSize);
  while (pos < seq.length - spectrum.tupleSize + 1 and 
				 fixed == 0 and
				 newEdgeCreated) {
		//    std::cout << ".";
		//   std::cout.flush();
    newEdgeCreated = 0;
    CreateGrid(params.maxGap, grid);
    cube.push_back(grid);

		// for each gap
		for (k = 0; k < params.maxGap; k++) {
			// foreach nucleotide
			if (init == 0) {
				// we do not want to compute values for all 4 nucleotides
				// on the first iteration.

				// this is the first iteration, so we need to initialize
				// the only edge.  We assume that this points to nothing.
				newEdgeCreated = 1;
				edge.prevNuc = -1;
				edge.prevLevel = -1;
				edge.prevPosition = -1;
				edge.score = 0;
				//if (!GetHashValue(seq, pos, spectrum.tupleSize, tuple)) {
				if (!tuple.Valid()) {
					std::cout << "Error, can only solidify starting on valid tuples (internal error)"
										<< std::endl;
					exit(1);
				}
				edge.CopyTuple(tuple);
				init = 1;

				// store the edge here.
				grid[0][numeric_nuc_index[(unsigned char) edge.tuple[spectrum.tupleSize-1]]][edge.tuple] = edge;
#ifdef VERBOSE
				std::cout << "created a cell for column 0 nucleotide: " 
									<< (edge.tuple[spectrum.tupleSize-1]) << " " << edge.tuple <<std::endl;
#endif
			}
			else {
				// assign cell values for each nucleotide
				for ( n = 0; n < 4; n++ ) {
					// try to find incoming edges to each nucleotide.

					// Try to fix using mutations and deletions.
					// want to start at the position before pos, and 
					// work back up until either maxgap
					EdgeIterator edgeIt, edgeEnd;
					ssize_t m;
					unsigned char prevNuc;
					ssize_t score;
	   
					for (m = 1; m <= std::min(params.maxGap, pos - startPos); m++ ) {
	     
						// reference cell at position [pos -solidPos]
						//                   level k
						//                   nucleotide n
						// and find all edges going into it

						// Consider edges from all previous nucleotides into this one.
						for (prevNuc = 0; prevNuc < 4 ; prevNuc++) {
							// Look at all edges stored at the previous nucleotide 
							// to see if the k-1 tuple overlaps with a 
							for (edgeIt = cube[pos - startPos - m][k][prevNuc].begin();
									 edgeIt != cube[pos - startPos - m][k][prevNuc].end();
									 edgeIt++) {
								Concatenate( edgeIt->first, n, newTuple);
#ifdef VERBOSE
								std::cout << "from mp " << pos - startPos - m << " del: " << m-1 << ", level " << k << ", n " << nuc_char[prevNuc] << "  to: p " 
													<< pos - startPos << ", n " << nuc_char[n] << " size: " 
													<< cube[pos - startPos - m][k][prevNuc].size() 			   
													<< " " << edgeIt->first << " " << newTuple << std::endl;
#endif
								/*
									std::cout << newTuple << std::endl;
								*/
								if (spectrum.FindTuple(newTuple) != -1) {
									// create a new edge here.
									// compare the nucleotide at 'pos' with what we are 
									// changing to at 'n'
									if (unmasked_nuc_index[seq.seq[pos + spectrum.tupleSize-1]] < 4) 
										score = ScoreMat[n][unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]]] + 
											edgeIt->second.score + params.gapExtend * (m-1);
									else {
										/*
											std::cout << "using 'n' score because of" <<unmasked_nuc_index[seq.seq[pos + spectrum.tupleSize-1]] << std::endl;
										*/
										score = 1 + edgeIt->second.score + params.gapExtend * (m-1);
									}
#ifdef VERBOSE
									std::cout << "allowed extension: " << pos << " d: " << m << " " 
														<< (ssize_t) seq.seq[pos+ spectrum.tupleSize-1] << "(" 
														<< (ssize_t) unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]] << ") -> " << n << " score: "<< score 
														<< " mutation " <<  (ssize_t) ScoreMat[n][unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]]] 
														<< " gap: " << params.gapExtend * (m-1) << " m - 1: " << m-1 << " " << score << std::endl;
									std::cout << "prev score: " << edgeIt->second.score << std::endl;
#endif
									// Store a new edge if this score is better than
									// all other edges
									// and edge is defined as a 'prev tuple' 'cur tuple' edge
									// for each nucleotide.
		   
									// so we need to pass in the set at the current position
									// the current level 'm', the current nucleotide 'n'
									// and try to connect it to the tuple 'newTuple'
									if (score < params.scoreThreshold) {
#ifdef VERBOSE
										std::cout << "storing edge into cube: p: " << pos - startPos << " level: " << k << " nuc: " << n << std::endl;
#endif
										newEdgeCreated = 1;
										StoreEdge(pos - startPos - m, k, prevNuc, edgeIt->first,
															cube[pos - startPos][k][n], newTuple, score);
										stats.numEdge++;
									}
								}
								/*
									else {
		     
									std::cout << "did not find " << newTuple << std::endl;
									}
								*/
							}
						}
					} // done fixing deletions

					if (k > 0) {
						for (prevNuc = 0; prevNuc < 4; prevNuc++ ) {
							for (edgeIt = cube[pos - startPos][k-1][prevNuc].begin();
									 edgeIt != cube[pos - startPos][k-1][prevNuc].end();
									 ++edgeIt) {
								Concatenate(edgeIt->first, n, newTuple);
#ifdef VERBOSE
								std::cout << "from kp " << pos - startPos  << ", k " << k-1 << ", n " << nuc_char[prevNuc] << "  to: p " 
													<< pos - startPos << ", k " << k << ", n " << nuc_char[n] << " size: " 
													<< cube[pos - startPos][k][prevNuc].size() << " " << edgeIt->first << " " << newTuple << std::endl;
#endif
								if (spectrum.FindTuple(newTuple) != -1) {
									if (unmasked_nuc_index[seq.seq[pos+ spectrum.tupleSize-1]] < 4) 
										score = ScoreMat[n][unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]]] + 
											edgeIt->second.score + params.gapExtend;
									else
										score = 10000 + edgeIt->second.score + params.gapExtend;
#ifdef VERBOSE
									std::cout << "allowed insert-extension: " << pos << " d: " << m << " " 
														<< (ssize_t) seq.seq[pos+ spectrum.tupleSize-1] << "(" 
														<< (ssize_t) unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]] << ") -> " << n << " score: "<< score 
														<< " mutation " <<  (ssize_t) ScoreMat[n][unmasked_nuc_index[seq.seq[pos+spectrum.tupleSize-1]]] 
														<< " gap: " << params.gapExtend * (m-1) << " m - 1: " << m-1 << " " << score << std::endl;
									std::cout << "prev score: " << edgeIt->second.score << std::endl;
#endif

									if (score < params.scoreThreshold) {
#ifdef VERBOSE
										std::cout << "storing insert edge into cube: p: " << pos - startPos << " level: " << k << " nuc: " << n << std::endl;
#endif
										newEdgeCreated = 1;
										StoreEdge(pos - startPos, k-1, prevNuc, edgeIt->first,
															cube[pos - startPos][k][n], newTuple, score);
										stats.numEdge++;
									}
								}
								else {
									/*
										std::cout << "did not ifind "<< newTuple << std::endl;
									*/
								}
							}
						}
					}
				}
      }
    }
    pos++;
  }

  //std::cout << std::endl;
  EdgeIterator minEdgeIt;
  ssize_t minN, minK;
  ssize_t numMin, minScore;
  minScore = 0;
  numMin   = 0;
  ssize_t success = 0;
  ssize_t trim;
  MultTuple minTuple;
  ssize_t minTrim;
#ifdef VERBOSE
	std::cout << "nec: " << newEdgeCreated << " cs:"  << cube.size() << " last: " << cube.size() + spectrum.tupleSize - 1 
						<< " " << seq.length - params.maxTrim << std::endl;
#endif
  if (newEdgeCreated == 1 or 
      (newEdgeCreated == 0 and cube.size() + spectrum.tupleSize - 1 > seq.length - params.maxTrim)) {
    // newEdgeCreated == 1 means that this ended with an edge that 'explains' the last nucleotide in 'seq'.
    // newEdgeCreated == 0 means that the method had to bail out since no possible explanations were 
    //                     available for seq at (pos + tupleSize - 1).  However it is possible that 
    //                     some amount of the sequence was fixed, so we want to keep that, and discard 
    //                     the remainder of the read.  If the last position fixed (startPos + cube.size()) 
    //                     is within params.maxTrim of the end then we can just use that.
		#ifdef VERBOSE
    std::cout << "finding min scoring edge, with cube of size : " << cube.size() << " and edge: "<< newEdgeCreated << std::endl;
		#endif
    if (newEdgeCreated == 0) {
      DeleteGrid(params.maxGap, cube[cube.size()-1]);
      cube.pop_back();
    }
    ssize_t numMin;
    if ( (numMin = FindMinimumScoreEdge(cube, params, minEdgeIt, 
																				minN, minK, numMin, minScore, minTuple, 
																				minTrim, trim)) == 1 ) {
      // make sure we finished by creating an edge, otherwise 
      // no good edges were found.
      
			//          std::cout << "found a minimum scoring edge at: " << minN << " k: " << minK << " score: " << minEdgeIt->second.score << std::endl;
			//           std::cout << "backtracking of length: " << cube.size()-1 - trim-minK << std::endl;
      Backtrack(cube, cube.size()-1-trim-minK, minK, minN, minTuple, fixedSeq, stats, params, lastEdit, end);
			#ifdef VERBOSE
      std::cout << "got fix: ";
      fixedSeq._ascii = 1;
      fixedSeq.PrintSeq(std::cout);
      std::cout << std::endl;
			#endif
      success = 1;
      replaceLength = seq.length - spectrum.tupleSize;

    }
    else {
			//    std::cout << "There were " << numMin << " possible fixes for this sequence " << std::endl;
      stats.numMultiplePaths++;
    }
  }
  else {
    if (newEdgeCreated == 0) {
      stats.numNoPathFound++;
    }
		else {
			stats.numErrorAtEnd++;
		}
  }
  if (success == 1) {
    // std::cout << "fixed: " << seq.namestr << std::endl;
  }
  else {
    // std::cout << "no fix for " << seq.namestr << std::endl;
  }
  for (i = 0; i < cube.size(); i++ ) {
    DeleteGrid(params.maxGap, cube[i]);
  }
  return success;
}

ssize_t FindMinimumScoreEdge(Cube &matrix, FixParams &params, 
												 EdgeIterator &minEdge, 
												 ssize_t &minN, ssize_t &minK, ssize_t& numMin, ssize_t &minScore,
												 MultTuple &minTuple, ssize_t &minTrim,
												 ssize_t &trim) {
  numMin = 0;
  minScore = 999999999;
  ssize_t last = matrix.size() - 1;
  // look over all nucleotides
  ssize_t n, k;
  EdgeIterator edgeIt, edgeEnd;
  minN = -1; minK = -1;
  numMin = 0;
  trim = 0;
  while (numMin != 1 and 
				 trim < params.maxTrim) {
    for (k = 0; k < params.maxGap; k++) {
      for (n = 0; n < 4; n++ ) {
				if (last - trim - k >= 0) {
					for (edgeIt = matrix[last-trim-k][k][n].begin();
							 edgeIt != matrix[last-trim-k][k][n].end();
							 ++edgeIt) {
						if (edgeIt->second.score < minScore) {
							minScore = edgeIt->second.score;
							minN     = n;
							minK     = k;
							minEdge  = edgeIt;
							numMin   = 1;
							minTuple = edgeIt->first;
							minTrim  = trim;
						}
						else if (edgeIt->second.score == minScore) {
							++numMin;
						}
					}
				}
      }
    }
    if (numMin != 1)
      ++trim;
  }
  return numMin;
}

void Backtrack(Cube &matrix, ssize_t pos, ssize_t level, ssize_t nuc, MultTuple tuple, 
							DNASequence &seq, Stats &stats, FixParams &params, 
							ssize_t &firstEdit, ssize_t &end) {

  // Trace a path in 'matrix' starting at the tuple referenced
  // at the cell matrix[pos][level][nuc].

  if (pos < 0) 
    return; // return 0;
  std::vector<unsigned char> newSeq;
  assert(matrix[pos][level][nuc].find(tuple) != matrix[pos][level][nuc].end());
  SEdge edge;
  ssize_t length = 0;
  ssize_t trim = 0;
	//UNUSED// ssize_t firstOk= -1;
	firstEdit = -1;
	end = pos;
  while (matrix[pos][level][nuc][tuple].prevNuc != -1) {
    length++;
    newSeq.push_back(nuc_char[nuc]);
    edge = matrix[pos][level][nuc][tuple];

    // there must be a link to the previous cell
    assert(matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc].find(edge.tuple) !=
					 matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc].end());

    // store some statistics about the fix
    ssize_t score = matrix[pos][level][nuc][tuple].score;
    ssize_t prevScore = matrix[edge.prevPosition][edge.prevLevel][(unsigned char) edge.prevNuc][edge.tuple].score;
    
    if (score > prevScore and length < params.maxTrim) {
      trim = length;
    }
		if (score > prevScore and firstEdit == -1 ){
			firstEdit = pos;
			//			std::cout << "first edit: " << firstEdit << " end: " << end << " span: " << end - firstEdit << std::endl;
		}
    
    if (edge.prevPosition == pos - 1 and
				prevScore < score)
      ++stats.numMut;
    else if (edge.prevLevel < level) {
      ++stats.numIns;
    }
    else if ( edge.prevPosition < pos - 1 )
      ++stats.numDel;
    
    pos   = edge.prevPosition;;
    level = edge.prevLevel;
    nuc   = edge.prevNuc;
    tuple = edge.tuple;

  }
  if (trim > 0) {
    //    std::cout << "trimmed " << trim << std::endl;
  }

  seq.Reset(newSeq.size() - trim);
  ssize_t i;
  for (i = newSeq.size()-1 ; i >= trim; i--) {
    seq.seq[newSeq.size() - i - 1] = newSeq[i];
  }
  
}

void StoreEdge(ssize_t prevPos, ssize_t prevLevel, ssize_t prevNuc, MultTuple prevTuple, 
							 Cell& cell, MultTuple tuple, ssize_t score) {
  // What the heck is going on here?  Why didn't I comment????
  // Ok, a cell is a collection of tuples that end in the same nucleotide.
  // I do this so that I can have a fixed number of cells at each iteration of error correction.
  // There may be many different tuples that end in the same iteration
  // First try and locate the tuple in the cell.  If it is not 
  // in the cell, append it.
  // If it is in the cell, if the score to reach this tuple is better than
  // the previous score, use the current path.  Otherwise, use the previous path.

  Cell::iterator cellIt;
  cellIt = cell.find(tuple);
  if (cellIt == cell.end() ) {
    cell[tuple].prevNuc = prevNuc;
    cell[tuple].prevLevel = prevLevel;
    cell[tuple].prevPosition   = prevPos;
    cell[tuple].CopyTuple(prevTuple);
		assert(cell[tuple].tuple.size() == cell[tuple].tuple.tupleSize);
    cell[tuple].score     = score;
  }
  else {
    if (cellIt->second.score > score) {
      cell[tuple].prevNuc      = prevNuc;
      cell[tuple].prevLevel    = prevLevel;
      cell[tuple].prevPosition = prevPos;
      cell[tuple].CopyTuple(prevTuple);
      cell[tuple].score        = score;
    }
  }
}


template <typename T_Spectrum>
ssize_t SolidifyRead(DNASequence &seq, T_Spectrum &spectrum, 
								 FixParams &params, Stats &stats) {

  DNASequence toFix, solidSeq;
  ssize_t replacedLength;
  ssize_t solidPos;
  ssize_t lastSolidPos;
  Stats readStats;
  /*
    std::cout << "fixing full read " << std::endl;
    seq.PrintSeq(std::cout);
    std::cout << std::endl;
  */
  std::string namestr= seq.namestr;
  if (FindSolidPosition(seq, spectrum, params.span, lastSolidPos, solidPos) ) {
    //    std::cout << "first solid: " << lastSolidPos << " sp: " << solidPos << std::endl;
    // make sure there are errors to correct
//		std::cout << lastSolidPos << " " << solidPos << " " << seq.length - spectrum.tupleSize << std::endl;
    if (lastSolidPos == 0 and solidPos + spectrum.tupleSize == seq.length) {
      // this sequence is ok, probably just return
#ifdef VERBOSE
      std::cout << "The full sequence is solid." << std::endl;
#endif
      return 1;
    }
    if (solidPos < seq.length - spectrum.tupleSize) {
      toFix.namestr = seq.namestr;
      toFix.seq = &seq.seq[solidPos];
      toFix.length = seq.length - solidPos ;
      toFix._ascii = seq._ascii;
#ifdef VERBOSE
      std::cout << "first: " << lastSolidPos << " solid: " << solidPos << std::endl;
			//      std::cout << "fixing: " << std::endl;
			//      toFix.PrintSeq(std::cout);
			//      std::cout << std::endl;
#endif
			ssize_t lastEdit, end;
      if (SolidifyUntilFixed(toFix, spectrum, params, solidSeq, replacedLength, readStats, lastEdit, end) == 0) {
				//				std::cout << " not fixable in forward direction " << std::endl;
				stats.Append(readStats);
				return 0;
      }
			if (end - lastEdit < params.edgeLimit) {
				stats.numErrorAtEnd++;
				return 0;
			}

#ifdef VERBOSE
      std::cout << "solid pos: " << solidPos + spectrum.tupleSize << " " << solidSeq.length << " " << replacedLength << std::endl;
      std::cout << "solidified sequence "<< std::endl;
      solidSeq._ascii = 1;
      solidSeq.PrintSeq(std::cout);
      std::cout << std::endl;
      std::cout << "patching " << solidPos + spectrum.tupleSize << " solidseqlen: "<<solidSeq.length << " replacing len: " << replacedLength << std::endl;
#endif

      PatchSeq(seq, solidPos + spectrum.tupleSize, solidSeq, replacedLength);
      solidSeq.Reset();
      seq._ascii = 1;
#ifdef VERBOSE
			//////////////////////////////////////////////////
      std::cout << "newly solid seq: " << std::endl;
      seq.PrintSeq(std::cout);
      std::cout << std::endl;
			//////////////////////////////////////////////////
#endif



    }
    if (lastSolidPos > 0) {
      DNASequence reverse;
      MakeRC(seq, reverse);

			//////////////////////////////////////////////////
#ifdef VERBOSE
			std::cout << "last solid pos: " << lastSolidPos << std::endl;
			std::cout << "fixing seq at pos: " << reverse.length - lastSolidPos - 1 - spectrum.tupleSize << std::endl;
			std::cout << "of length: " << lastSolidPos + spectrum.tupleSize << " seq: " << reverse.length << std::endl;
#endif
			//////////////////////////////////////////////////


      toFix.seq = &reverse.seq[reverse.length - lastSolidPos - spectrum.tupleSize];
      toFix.length = lastSolidPos + spectrum.tupleSize;
      toFix._ascii = 1;

      
			//////////////////////////////////////////////////
#ifdef VERBOSE
      std::cout << "fixing: " << std::endl;
      toFix.PrintSeq(std::cout);
      std::cout << std::endl;
      std::cout << "reverse to fix: " << toFix.length << std::endl;
#endif
			//////////////////////////////////////////////////

			ssize_t lastEdit, end;
      if (SolidifyUntilFixed(toFix, spectrum, params, 
														 solidSeq, replacedLength, readStats, lastEdit, end) == 0) {
				//        std::cout << "not fixable in reverse direction " << std::endl;
				stats += readStats;
				return 0;
      }
			if (end - lastEdit < params.edgeLimit) {
				stats.numErrorAtEnd++;
				return 0;
			}
      // patch in the error
      solidSeq._ascii = 1;

			//////////////////////////////////////////////////
#ifdef VERBOSE
      std::cout << "patching: " << reverse.length - lastSolidPos << " " << solidSeq.length << " " << replacedLength << std::endl;
      std::cout << "fixed seq (rc): " << std::endl;
      solidSeq.PrintSeq(std::cout);
      std::cout << std::endl;
#endif
			//////////////////////////////////////////////////

      PatchSeq(reverse, reverse.length - lastSolidPos, solidSeq, replacedLength);
      //      std::cout << "rc fixed seq: " << std::endl;
      solidSeq.Reset();
      reverse._ascii = 1;
      /*
				reverse.PrintSeq(std::cout);
				std::cout << std::endl;
      */
      MakeRC(reverse, seq);
      seq.namestr = reverse.namestr;
    }
    stats += readStats;
    seq._ascii = 1;
    seq.namestr = namestr;
    return 1;
  }
  else {
    // No solid position found
    stats.numNoSolid++;
    return 0;
  }
  /*
    std::cout << "the fixed sequence is; " << std::endl;
    seq.PrintSeq(std::cout);
    std::cout << std::endl;
  */
}

void PatchSeq(DNASequence &seq, ssize_t pos, DNASequence &patch, ssize_t replaceLength) {
  DNASequence patchedSeq;
  patchedSeq.Reset(seq.length - replaceLength + patch.length);

  patchedSeq.namestr = seq.namestr;
  // copy the unchanged segments
  patchedSeq._ascii = 0;
  ssize_t i;

  // Write nonsense into the sequence to make sure I patch every position.
  for (i = 0; i < patchedSeq.length; i++ )
    patchedSeq.seq[i] = 255;

#ifdef VERBOSE
  std::cout << "ps init: " << std::endl;
  patchedSeq.PrintSeq(std::cout);
  std::cout << std::endl;
#endif

  // The sequence is ok up to 'pos'
  for (i = 0; i < pos; i++ ){
    patchedSeq.seq[i] = seq.seq[i];
  }

#ifdef VERBOSE
  std::cout << "ps now: " << std::endl;
  patchedSeq._ascii = 1;
  patchedSeq.PrintSeq(std::cout);
  std::cout << std::endl;
  std::cout << "sl: " << seq.length << " pos: " << pos << " rl: " << replaceLength << " pl " << patch.length << std::endl;
#endif

  for (i = 0; i < seq.length - pos - replaceLength; i++)
    patchedSeq.seq[pos + patch.length + i] = seq.seq[pos + replaceLength + i];

#ifdef VERBOSE
  std::cout << "ps now (with end) " << std::endl;
  patchedSeq.PrintSeq(std::cout);
  std::cout << std::endl;
#endif

  // copy the new segment
  for (i = 0; i < patch.length; i++) 
    patchedSeq.seq[pos + i] = patch.seq[i];

#ifdef VERBOSE
  std::cout << "copied patch: " << std::endl;
  patch.PrintSeq(std::cout);
  std::cout << std::endl;
  std::cout << "result: "<< std::endl;
  patchedSeq.PrintSeq(std::cout);
  std::cout << std::endl;
#endif
  seq = patchedSeq;
  patchedSeq.Reset();
}



void ReverseSeq(DNASequence &seq, DNASequence &rev) {
  // not reverse complement, just reverse
  ssize_t i;
  rev.Reset(seq.length);
  for (i = 0; i < seq.length; i++ ) {
    rev.seq[rev.length - i - 1] = seq.seq[i];
  }
}
