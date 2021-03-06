//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "transition_extractor.hpp"

namespace path_extend {
namespace validation {
class SimpleTransitionGraphValidator {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    const ContigTransitionStorage reference_transition_storage_;

 public:
    SimpleTransitionGraphValidator(const ContigTransitionStorage &reference_transition_storage);

    boost::optional<vector<ScaffoldVertex>> GetCorrectPath(const SimpleTransitionGraph& graph,
                                                           ScaffoldVertex source, ScaffoldVertex sink) const;


    DECL_LOGGER("SimpleTransitionGraphValidator");
};

class SimpleTransitionGraphValidatorConstructor {
    const conj_graph_pack& gp_;

 public:
    SimpleTransitionGraphValidatorConstructor(const conj_graph_pack &gp);

    SimpleTransitionGraphValidator GetValidator(const string &path_to_reference) const;
};
}
}