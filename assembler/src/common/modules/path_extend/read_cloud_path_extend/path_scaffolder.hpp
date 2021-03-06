//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/path_extender.hpp"
#include "common/pipeline/graph_pack.hpp"

namespace path_extend {

class StartFinder {
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef unordered_map<ScaffoldVertex, ScaffoldVertex> TransitionMap;

    const Graph &g_;
 public:
    StartFinder(const Graph &g);

    std::unordered_set<ScaffoldVertex> GetStarts(const TransitionMap &transition_map) const;
};

class PathScaffolder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

    const debruijn_graph::conj_graph_pack& gp_;
    const ScaffoldingUniqueEdgeStorage& unique_storage_;
    size_t small_path_length_threshold_;
    size_t large_path_length_threshold_;

 public:
    PathScaffolder(const conj_graph_pack &gp_,
                   const ScaffoldingUniqueEdgeStorage &unique_storage_,
                   size_t small_path_length_threshold_,
                   size_t large_path_length_threshold_);

    void MergePaths(const PathContainer &paths) const;

 private:

    void MergeUnivocalEdges(const vector<ScaffoldEdge> &scaffold_edges) const;

    void ExtendPathAlongConnections(const ScaffoldVertex& start,
                                    const std::unordered_map<ScaffoldVertex, ScaffoldVertex> &merge_connections,
                                    const std::unordered_map<ScaffoldVertex, size_t> &start_to_length) const;

    DECL_LOGGER("PathScaffolder");
};

}