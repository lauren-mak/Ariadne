/***************************************************************************
 * Title:          enum2aln.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "EnumUtils.h"
#include "StripGen.h"
#include "TupleLib.h"
#include "SeqReader.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>

void PrintUsage();
void initenv(int argc, char *argv[], 
	     std::string &inpfile, 
	     std::string &refSeqFile,
	     std::string &qrySeqFile,
	     std::string &outfile,
	     ssize_t &lineLength);


int main(int argc, char* argv[]) {
  std::string inputFileName, outputFileName, refSeqFileName, qrySeqFileName;
  inputFileName = outputFileName = refSeqFileName = qrySeqFileName = "";
  std::ofstream outfile;
  std::ostream *output;
  ssize_t lineLength;

  lineLength = 60;
  initenv(argc, argv, inputFileName, refSeqFileName, qrySeqFileName, outputFileName, lineLength);
  if (inputFileName == "" ||
      refSeqFileName == "" ||
      qrySeqFileName == "") {
    std::cout << "You must specify a value for every argument." << std::endl;
    PrintUsage();
  }

  std::ifstream inFile;
  inFile.open(inputFileName.c_str());
  if ( ! inFile.good() ) {
    std::cout << "Could not open infile " << std::endl;
    exit(0);
  }

  // Get the sequences
  std::ifstream refFile;
  refFile.open(refSeqFileName.c_str());
  if ( ! refFile.good() ) {
    std::cout << "Could not open reffile " << std::endl;
    exit(0);
  }

  std::ifstream qryFile;
  qryFile.open(qrySeqFileName.c_str());
  if ( ! qryFile.good() ) {
    std::cout << "Could not open qryfile " << std::endl;
    exit(0);
  }

  if (outputFileName == "")
    output = &std::cout;
  else {
    outfile.open(outputFileName.c_str());
    output = &outfile;
  }

  DNASequence refSeq, qrySeq;

  SeqReader seqReader(&refFile);

  // Read the ref sequence (first contig)
  seqReader.GetSeq(refSeq, SeqReader::noConvert);

  // read the qry sequence (first contig)
  seqReader.SetSource(&qryFile);
  seqReader.GetSeq(qrySeq, SeqReader::noConvert);
  
  // For holding descripion of mapping
  ssize_t *enumerations, *refLocations, *qryLocations, *size, *stripsMerged;
  ssize_t length;
  // Read the mapping (called read strips for historic reasons).
  ReadStrips(inFile,
	     enumerations, refLocations, refLocations,
	     qryLocations, qryLocations, size,
	     stripsMerged, length);

  // No strips input
  if (length == 0) {
    std::cout << "No sequences were read." << std::endl;
    return 0;
  }

  PrintAlignment(output,
		 refSeq,
		 qrySeq,
		 enumerations, refLocations, qryLocations, length);
}




void initenv(int argc, char *argv[], 
	     std::string &inpfile, 
	     std::string &refSeqFile,
	     std::string &qrySeqFile,
	     std::string &outfile,
	     ssize_t &lineLength) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "e:a:r:q:l:")) != EOF) {
    switch(copt) {
    case 'e':
      inpfile = optarg;
      continue;
    case 'a':
      outfile = optarg;
      continue;
    case 'r':
      refSeqFile = optarg;
      continue;
    case 'q':
      qrySeqFile = optarg;
      continue;
    case 'l':
      lineLength = atoi(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  if (inpfile == "") {
    PrintUsage();
    exit(1);
  }
}

void PrintUsage() {
  std::cout << "Create an alignment view from the enumeration of two sequences " <<std::endl;
  std::cout << " generated by tupal." << std::endl;
  std::cout << "usage:  enum2aln  -e enumeration  -r refSeq -q qrySeq -a alignment(output)"<< std::endl;
  std::cout << "   -e enumeration of qrySeq versus refSeq " << std::endl;
  std::cout << "   -r reference dna sequence. " << std::endl;
  std::cout << "   -q qrySeqFileName " << std::endl;
  std::cout << "   -a alignment-file (stdout) " << std::endl;
}

