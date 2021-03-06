//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "reference_path_index.hpp"
#include "common/barcode_index/cluster_storage/barcode_cluster.hpp"

namespace path_extend {
namespace validation {

class PathClusterValidator {
    ReferencePathIndex ref_path_index_;

 public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::vector<ScaffoldVertex> SimplePath;

 public:
    PathClusterValidator(const ReferencePathIndex &ref_path_index);

    bool IsCorrect(const cluster_storage::Cluster &cluster) const;

    bool IsCorrect(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;

    bool IsCovered(const cluster_storage::Cluster &cluster) const;

    bool IsCovered(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;

    bool IsCovered(const scaffold_graph::ScaffoldVertex &vertex) const;

    void PrintRefIndexInfo(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;

    boost::optional<SimplePath> GetReferencePath(const set<ScaffoldVertex> &vertices) const;
};

}
}