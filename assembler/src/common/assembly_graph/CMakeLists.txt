############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(assembly_graph CXX)

add_library(assembly_graph STATIC
            components/connected_component.cpp paths/bidirectional_path.cpp paths/bidirectional_path_io/io_support.cpp paths/bidirectional_path_io/bidirectional_path_output.cpp graph_support/scaff_supplementary.cpp ../modules/alignment/edge_index_refiller.cpp graph_support/coverage_uniformity_analyzer.cpp contracted_graph/contracted_graph.cpp contracted_graph/contracted_graph_builder.cpp contracted_graph/contracted_graph_helper.cpp contracted_graph/graph_condensation.cpp contracted_graph/contracted_statistics.cpp contracted_graph/contracted_statistics.hpp)
target_link_libraries(assembly_graph hattrie)
