############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(path_extend CXX)

add_library(path_extend STATIC pe_config_struct.cpp
                    pipeline/launch_support.cpp
                    pipeline/launcher.cpp
                    pipeline/extenders_logic.cpp
                    scaffolder2015/extension_chooser2015.cpp
                    scaffolder2015/scaffold_graph.cpp
                    scaffolder2015/scaffold_graph_constructor.cpp
                    scaffolder2015/scaffold_graph_visualizer.cpp
                    scaffolder2015/connection_condition2015.cpp 
                    scaffolder2015/path_polisher.cpp
                    scaffolder2015/scaffold_vertex.cpp
                    read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.cpp
                    read_cloud_path_extend/conjugate_score_extractor.hpp
                    read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_construction_pipeline.cpp
                    read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage.cpp
                    read_cloud_path_extend/validation/scaffold_graph_validation.cpp
                    read_cloud_path_extend/validation/transition_extractor.cpp
                    read_cloud_path_extend/scaffold_graph_extractor.cpp
                    read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.cpp
                    read_cloud_path_extend/scaffold_graph_dijkstra.hpp
                    read_cloud_path_extend/transitions/transitions.cpp
                    read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.cpp
                    read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.cpp
                    read_cloud_path_extend/intermediate_scaffolding/pe_extraction.cpp
                    read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp
                    read_cloud_path_extend/intermediate_scaffolding/predicate_builders.cpp
                    read_cloud_path_extend/path_extend_dijkstras.hpp
                    read_cloud_path_extend/path_scaffolder.cpp
                    read_cloud_path_extend/scaffold_graph_construction/containment_index_threshold_finder.cpp
                    read_cloud_path_extend/read_cloud_polisher_support.cpp
                    read_cloud_path_extend/scaffold_graph_construction/construction_callers.cpp
                    read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.cpp
                    read_cloud_path_extend/scaffold_graph_gap_closer/cloud_scaffold_graph_gap_closer.cpp
                    read_cloud_path_extend/scaffold_graph_gap_closer/scaffold_graph_gap_closer.cpp
                    read_cloud_path_extend/fragment_model/secondary_stats_estimators.cpp
                    read_cloud_path_extend/extender_support/candidate_selectors.cpp
                    read_cloud_path_extend/extender_support/entry_collectors.cpp
                    read_cloud_path_extend/statistics/path_cluster_statistics.cpp
                    read_cloud_path_extend/validation/transition_subgraph_validation.cpp
        read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_path_cleaner.cpp
                    read_cloud_path_extend/validation/reference_path_index.cpp
                    read_cloud_path_extend/validation/path_cluster_validation.cpp
                    read_cloud_path_extend/statistics/cloud_check_statistics.cpp
                    read_cloud_path_extend/statistics/component_validation.cpp
                    read_cloud_path_extend/statistics/long_edge_dataset.cpp
                    read_cloud_path_extend/statistics/perfect_transitions.cpp
                    read_cloud_path_extend/statistics/perfect_clouds.cpp
                    read_cloud_path_extend/statistics/coverage_break_analysis.cpp
                    read_cloud_path_extend/statistics/path_scaffolder_analyzer.cpp
                    read_cloud_path_extend/contracted_graph_scaffolding/contracted_graph_scaffolder.cpp
                    read_cloud_path_extend/contracted_graph_scaffolding/contracted_graph_scaffolder.hpp
                    read_cloud_path_extend/contracted_graph_scaffolding/contracted_gfa_writer.cpp
                    read_cloud_path_extend/contracted_graph_scaffolding/contracted_gfa_writer.hpp
                    read_cloud_path_extend/statistics/short_edge_dataset.cpp
                    read_cloud_path_extend/statistics/short_edge_dataset.hpp read_cloud_path_extend/statistics/split_index_statistics.cpp read_cloud_path_extend/statistics/split_index_statistics.hpp)

target_link_libraries(path_extend assembly_graph ssw)


