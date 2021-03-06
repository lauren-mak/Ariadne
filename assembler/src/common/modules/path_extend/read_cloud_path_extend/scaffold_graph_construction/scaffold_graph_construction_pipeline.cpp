//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_construction_pipeline.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "containment_index_threshold_finder.hpp"
#include "read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {

ScaffolderParams ScaffolderParamsConstructor::ConstructScaffolderParams(
        const Graph& g, size_t min_length, cluster_model::ClusterStatisticsExtractor primary_extractor) const {
    size_t length_threshold = min_length;
    size_t tail_threshold = min_length;
    size_t count_threshold = cfg::get().ts_res.scaff_con.count_threshold;
    double split_procedure_strictness = cfg::get().ts_res.scaff_con.split_procedure_strictness;
    size_t transitive_distance_threshold = cfg::get().ts_res.scaff_con.transitive_distance_threshold;

    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t min_length_for_barcode_collection = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;

    const double initial_distance_percentile = cfg::get().ts_res.scaff_con.cluster_length_percentile;
    const double score_percentile = cfg::get().ts_res.scaff_con.score_percentile;

    size_t initial_distance = primary_extractor.GetLengthPercentile(initial_distance_percentile);
    auto score_estimation_params = GetScoreEstimationParams(g, primary_extractor, score_percentile,
                                                            initial_distance_percentile, min_length);

    ScaffolderParams result(length_threshold, tail_threshold, count_threshold,
                            connection_length_threshold, connection_count_threshold,
                            initial_distance, split_procedure_strictness,
                            transitive_distance_threshold,
                            min_length_for_barcode_collection, score_estimation_params);
    return result;
}

LongEdgePairGapCloserParams ScaffolderParamsConstructor::ConstructGapCloserParamsFromMainParams(
        const ScaffolderParams &params,
        const Graph &g,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t unique_edge_length,
        size_t max_threads) const {
    const double SHORT_EDGE_THRESHOLD = 0.01;
//    auto threshold_estimator_params = params.score_estimation_params_;
//    ShortEdgeScoreThresholdEstimatorFactory threshold_estimator_factory(g, barcode_extractor,
//                                                                        threshold_estimator_params.training_edge_length_threshold_,
//                                                                        unique_edge_length,
//                                                                        threshold_estimator_params.max_cluster_gap_,
//                                                                        threshold_estimator_params.score_percentile_,
//                                                                        max_threads);
//    auto score_threshold_estimator = threshold_estimator_factory.GetThresholdEstimator();
//    //fixme magic
//    auto connection_score_threshold = score_threshold_estimator->GetThreshold() * 0.5;
//    if (connection_score_threshold < MIN_SHORT_EDGE_THRESHOLD) {
//        INFO("Estimated threshold is too small, using " << MIN_SHORT_EDGE_THRESHOLD << " as default threshold");
//        connection_score_threshold = MIN_SHORT_EDGE_THRESHOLD;
//    }
    double connection_score_threshold = SHORT_EDGE_THRESHOLD;
    
    double relative_coverage_threshold = cfg::get().ts_res.scaff_con.relative_coverage_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t tail_threshold = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
    bool normalize_using_cov = true;
    return LongEdgePairGapCloserParams(connection_count_threshold, tail_threshold, connection_score_threshold,
                                       relative_coverage_threshold, connection_length_threshold,
                                       normalize_using_cov);
}
ScaffolderParams::ScoreEstimationParams ScaffolderParamsConstructor::GetScoreEstimationParams(
        const Graph &g, cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor,
        double score_percentile, double cluster_length_percentile, size_t block_length) const {
    size_t max_training_gap = cluster_statistics_extractor.GetLengthPercentile(cluster_length_percentile);
//    double length_normalizer = static_cast<double>(cfg::get().ts_res.long_edge_length_max_upper_bound);
//    const double block_normalizer = 2.0 / 3.0;
//    double distance_normalizer = (static_cast<double>(block_length) + length_normalizer * block_normalizer) / length_normalizer;
//    max_training_gap = static_cast<size_t>(distance_normalizer * static_cast<double>(max_training_gap));
    size_t min_training_length = 2 * block_length + max_training_gap;
    INFO("Max training gap: " << max_training_gap);
    INFO("Min training length: " << min_training_length);
    ScaffolderParams::ScoreEstimationParams score_estimation_params(score_percentile, max_training_gap,
                                                                    min_training_length);
    return score_estimation_params;
}

ScaffoldGraphConstructionPipeline::ScaffoldGraphConstructionPipeline(
    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> initial_constructor_, const Graph &g,
    const path_extend::ScaffolderParams &params)
    : initial_constructor_(initial_constructor_), construction_stages_(),
      intermediate_results_(), g_(g), params_(params) {}

void ScaffoldGraphConstructionPipeline::Run() {
    auto initial_graph_ptr = make_shared<ScaffoldGraph>(g_);
    INFO("Constructing initial scaffold graph");
    initial_graph_ptr = initial_constructor_->Construct();
    INFO("Constructed initial graph");
    INFO(initial_graph_ptr->VertexCount() << " vertices and " << initial_graph_ptr->EdgeCount()
                                          << " edges in initial graph");
    intermediate_results_.push_back({initial_graph_ptr, "Initial graph"});
    for (const auto &stage: construction_stages_) {
        INFO("Starting " << stage->getName());
        auto constructor = stage->GetScaffoldGraphConstuctor(params_, *(intermediate_results_.back().first));
        auto next_graph = constructor->Construct();
        intermediate_results_.push_back({next_graph, stage->getName()});
        INFO(next_graph->VertexCount() << " vertices and " << next_graph->EdgeCount() << " edges in current graph");
        path_extend::ScaffoldGraphExtractor extractor;
        auto univocal_edges = extractor.ExtractUnivocalEdges(*next_graph);
        INFO(univocal_edges.size() << " univocal edges in current graph");
    }
}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraph> ScaffoldGraphConstructionPipeline::GetResult() const {
    return intermediate_results_.back().first;
}
void ScaffoldGraphConstructionPipeline::AddStage(shared_ptr<IterativeScaffoldGraphConstructorCaller> stage) {
    construction_stages_.push_back(stage);
}
vector<std::pair<shared_ptr<path_extend::scaffold_graph::ScaffoldGraph>, string>> ScaffoldGraphConstructionPipeline::GetIntermediateResults() const {
    return intermediate_results_;
}

ScaffoldGraphPipelineConstructor::ScaffoldGraphPipelineConstructor(const conj_graph_pack &gp) : gp_(gp) {}
shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> ScaffoldGraphPipelineConstructor::ConstructSimpleEdgeIndex(
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices,
        ScaffoldGraphPipelineConstructor::BarcodeIndexPtr barcode_extractor,
        const ScaffolderParams &params, size_t max_threads) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = params.tail_threshold_;
    const size_t length_threshold = params.min_length_for_barcode_collection_;
    const size_t count_threshold = params.count_threshold_;
    auto tail_threshold_getter = make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads, scaffold_vertices);
    auto scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    return scaffold_index_extractor;
}

ScaffoldGraphConstructionPipeline BasicScaffoldGraphPipelineConstructor::ConstructPipeline(
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    INFO("Constructing scaffold graph with length threshold " << min_length_);
    path_extend::ScaffolderParamsConstructor params_constructor;

    path_extend::cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor(gp_.read_cloud_distribution_pack);
    auto params = params_constructor.ConstructScaffolderParams(gp_.g, min_length_, cluster_statistics_extractor);
    auto initial_constructor =
        make_shared<path_extend::scaffold_graph::UniqueScaffoldGraphConstructor>(gp_.g,
                                                                                 unique_storage_,
                                                                                 scaffold_vertices,
                                                                                 params.initial_distance_,
                                                                                 max_threads_);
    auto iterative_constructor_callers = ConstructStages(params, scaffold_vertices);
    INFO("Created " << iterative_constructor_callers.size() << " stages");
    ScaffoldGraphConstructionPipeline pipeline(initial_constructor, gp_.g, params);
    for (const auto stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    return pipeline;
}
BasicScaffoldGraphPipelineConstructor::BasicScaffoldGraphPipelineConstructor(
        const conj_graph_pack &gp,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        size_t max_threads_,
        size_t min_length_)
    : ScaffoldGraphPipelineConstructor(gp),
      unique_storage_(unique_storage),
      barcode_extractor_(barcode_extractor_),
      max_threads_(max_threads_),
      min_length_(min_length_) {}

ScaffoldGraphConstructionPipeline GapScaffoldGraphPipelineConstructor::ConstructPipeline(
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    INFO("Constructing scaffold graph with length threshold " << min_length_);
    path_extend::ScaffolderParamsConstructor params_constructor;
    path_extend::cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor(gp_.read_cloud_distribution_pack);
    auto params = params_constructor.ConstructScaffolderParams(gp_.g, min_length_, cluster_statistics_extractor);
    auto initial_constructor = GetInitialConstructor(params, scaffold_vertices);
    auto iterative_constructor_callers = ConstructStages(params, scaffold_vertices);
    INFO("Created " << iterative_constructor_callers.size() << " stages");
    ScaffoldGraphConstructionPipeline pipeline(initial_constructor, gp_.g, params);
    for (const auto& stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    return pipeline;
}
GapScaffoldGraphPipelineConstructor::GapScaffoldGraphPipelineConstructor(
        const conj_graph_pack &gp,
        const ScaffoldingUniqueEdgeStorage& unique_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        size_t max_threads_,
        size_t min_length_)
    : ScaffoldGraphPipelineConstructor(gp),
      unique_storage_(unique_storage),
      barcode_extractor_(barcode_extractor_),
      max_threads_(max_threads_),
      min_length_(min_length_) {}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GapScaffoldGraphPipelineConstructor::GetInitialConstructor(
        path_extend::ScaffolderParams params,
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    auto scaffold_index_extractor = ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);
    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(gp_.g, scaffold_index_extractor);
    vector<ScaffoldVertex> scaff_vertex_vector;
    std::copy(scaffold_vertices.begin(), scaffold_vertices.end(), back_inserter(scaff_vertex_vector));
    //fixme move somewhere
    const double score_threshold = 0.095;
    INFO("Setting containment index threshold to " << score_threshold);

    auto initial_constructor =
        make_shared<path_extend::scaffold_graph::ScoreFunctionScaffoldGraphConstructor>(gp_.g,
                                                                                        scaffold_vertices,
                                                                                        score_function,
                                                                                        score_threshold,
                                                                                        max_threads_);
    return initial_constructor;
}
FullScaffoldGraphPipelineConstructor::FullScaffoldGraphPipelineConstructor(
        const conj_graph_pack &gp,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        size_t max_threads_,
        size_t min_length_)
    : BasicScaffoldGraphPipelineConstructor(gp, unique_storage, barcode_extractor_, max_threads_, min_length_) {}

vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> FullScaffoldGraphPipelineConstructor::ConstructStages(
        path_extend::ScaffolderParams params,
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    auto scaffold_index_extractor = ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);

    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    iterative_constructor_callers.push_back(make_shared<BarcodeScoreConstructorCaller>(gp_.g, barcode_extractor_,
                                                                                       scaffold_index_extractor,
                                                                                       max_threads_));
    iterative_constructor_callers.push_back(
        make_shared<BarcodeConnectionConstructorCaller>(gp_.g, barcode_extractor_, scaffold_index_extractor,
                                                        unique_storage_, max_threads_));
//    bool scaffolding_mode = false;
//    iterative_constructor_callers.push_back(
//        make_shared<CompositeConnectionConstructorCaller>(gp_, barcode_extractor_, scaffold_index_extractor,
//                                                          unique_storage_, max_threads_, scaffolding_mode));

    //fixme works only with 2 graphs pipeline
    const size_t min_pipeline_length = cfg::get().ts_res.long_edge_length_lower_bound;
    bool launch_full_pipeline = min_length_ > min_pipeline_length;

    if (launch_full_pipeline) {
        const double EDGE_LENGTH_FRACTION = 0.5;
        auto fraction_tail_threshold_getter =
            make_shared<barcode_index::FractionTailThresholdGetter>(gp_.g, EDGE_LENGTH_FRACTION);
        auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               count_threshold, length_threshold,
                                                                               max_threads_, scaffold_vertices);
        auto split_scaffold_index_extractor =
            make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
        iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g,
                                                                                        split_scaffold_index_extractor,
                                                                                        max_threads_));
        iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, max_threads_));
    }
    return iterative_constructor_callers;
}
BinningScaffoldGraphPipelineConstructor::BinningScaffoldGraphPipelineConstructor(
        const conj_graph_pack &gp,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        size_t max_threads_,
        size_t min_length_)
    : GapScaffoldGraphPipelineConstructor(gp, unique_storage, barcode_extractor_, max_threads_, min_length_) {}

vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> BinningScaffoldGraphPipelineConstructor::ConstructStages(
        path_extend::ScaffolderParams params,
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    auto scaffold_index_extractor = ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    return iterative_constructor_callers;
}
MergingScaffoldGraphPipelineConstructor::MergingScaffoldGraphPipelineConstructor(
        const conj_graph_pack &gp,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        size_t max_threads_,
        size_t min_length_)
    : GapScaffoldGraphPipelineConstructor(gp, unique_storage, barcode_extractor_, max_threads_, min_length_) {}

vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> MergingScaffoldGraphPipelineConstructor::ConstructStages(
        path_extend::ScaffolderParams params,
        const set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t length_threshold = params.min_length_for_barcode_collection_;
    const size_t count_threshold = params.count_threshold_;
    auto scaffold_index_extractor = ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);

    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    //fixme works only with 2 graphs pipeline
//    const size_t min_pipeline_length = cfg::get().ts_res.long_edge_length_lower_bound;
//    bool launch_full_pipeline = min_length_ > min_pipeline_length;

    INFO("Min scaffold vertex length: " << min_length_);
    const double EDGE_LENGTH_FRACTION = 0.5;
    auto fraction_tail_threshold_getter =
        make_shared<barcode_index::FractionTailThresholdGetter>(gp_.g, EDGE_LENGTH_FRACTION);
    auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               count_threshold, length_threshold,
                                                                               max_threads_, scaffold_vertices);
    auto split_scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
    iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g,
                                                                                        split_scaffold_index_extractor,
                                                                                        max_threads_));
    iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, max_threads_));
    return iterative_constructor_callers;
}
} //path_extend

