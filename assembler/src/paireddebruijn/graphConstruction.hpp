#ifndef GRAPHCONSTRUCTION_H_
#define GRAPHCONSTRUCTION_H_
#include "simpleGraph.hpp"
#include "pairedGraph.hpp"
#include "common.hpp"
#include "graphVisualizer.hpp"

using namespace paired_assembler;

void constructGraph();
edgesMap sequencesToMap(string parsed_k_sequence, bool usePaired = true);
void createVertices(gvis::GraphScheme<int> &g, edgesMap &edges, verticesMap &verts, longEdgesMap &longEdges);
int expandRight(edgesMap &edges, verticesMap &verts, ll &finishKmer, Sequence* &finishSeq);
int expandLeft(edgesMap &edges, verticesMap &verts, ll &startKmer, Sequence* &startSeq);
int checkUniqueWayLeft(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayLeft(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int checkUniqueWayRight(edgesMap &edges, ll finishKmer, Sequence *finishSeq);
int goUniqueWayRight(edgesMap &edges, ll &finishKmer, Sequence* &finishSeq);
int storeVertex(gvis::GraphScheme<int> &g, verticesMap &verts, ll newKmer, Sequence* newSeq);
void resetVertexCount();
void expandDefinite(verticesMap &verts, longEdgesMap &longEdges);
void outputLongEdges(longEdgesMap &longEdges);
void traceReads(verticesMap &verts, longEdgesMap &longEdges);


#endif /*GRAPHCONSTRUCTION_H_*/
