//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "distribution_extractor.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "common/barcode_index/cluster_storage/initial_cluster_storage_builder.hpp"

namespace path_extend {
namespace cluster_model {

class MinTrainingLengthEstimator {
    const Graph& g_;
    size_t min_total_length_;
    size_t optimal_total_length_;
    size_t min_edges_;

 public:
    MinTrainingLengthEstimator(const Graph &g, size_t min_total_length, size_t optimal_total_length, size_t min_edges)
        : g_(g), min_total_length_(min_total_length), optimal_total_length_(optimal_total_length), min_edges_(min_edges) {}

    boost::optional<size_t> EstimateTrainingLength() const {
        size_t min_length = 5000;
        vector<size_t> edge_length_initial_list;
        boost::optional<size_t> result;
        omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(g_);
        for (const auto& edge: edge_it_helper) {
            if (g_.length(edge) >= min_length) {
                edge_length_initial_list.push_back(g_.length(edge));
            }
        }
        if (edge_length_initial_list.size() < min_edges_) {
            return result;
        }
        vector<std::pair<size_t, size_t>> length_rev_cumulative_list;
        std::sort(edge_length_initial_list.begin(), edge_length_initial_list.end(), std::greater<size_t>());
        size_t current_sum = 0;
        for (const auto& length: edge_length_initial_list) {
            current_sum += length;
            length_rev_cumulative_list.emplace_back(length, current_sum);
        }
        size_t total_long_length = current_sum;
        if (total_long_length < min_total_length_) {
            return result;
        }
        result = length_rev_cumulative_list.back().first;
        auto it = length_rev_cumulative_list.begin();
        current_sum = 0;
        size_t optimal_total_length = std::min(optimal_total_length_, total_long_length / 2);
        while (current_sum <= optimal_total_length and it != length_rev_cumulative_list.end()) {
            result = (*it).first;
            current_sum = (*it).second;
            ++it;
        }
        INFO("Estimated training length: " << result.get());
        return result;
    }
};

class MinTrainingLengthEstimatorHelper {
 public:

    boost::optional<size_t> EstimateTrainingLength(const Graph &g) const {
        const size_t min_training_edges = cfg::get().ts_res.min_training_edges;
        const size_t min_training_total_length = cfg::get().ts_res.min_training_total_length;
        const size_t optimal_training_total_length = cfg::get().ts_res.optimal_training_total_length;
        MinTrainingLengthEstimator training_length_estimator(g, min_training_total_length,
                                                             optimal_training_total_length, min_training_edges);
        return training_length_estimator.EstimateTrainingLength();
    }
};

template<class ClusterStatType>
struct DistributionStatistics {
  ClusterStatType min_;
  ClusterStatType max_;
  ClusterStatType mean_;
  ClusterStatType median_;

  DistributionStatistics(ClusterStatType min_, ClusterStatType max_, ClusterStatType mean_, ClusterStatType median_) :
      min_(min_), max_(max_), mean_(mean_), median_(median_) {}

  friend std::ostream& operator<< (std::ostream& stream, const DistributionStatistics<ClusterStatType> &statistics)
  {
      stream << "Min: " << statistics.min_ << std::endl
             << "Max: " << statistics.max_ << std::endl
             << "Mean: " << statistics.mean_ << std::endl
             << "Median: " << statistics.median_ << std::endl;
      return stream;
  }
};

struct StatisticsPack {
  DistributionStatistics<ClusterLength> length_statistics_;
  DistributionStatistics<ClusterCoverage> coverage_statistics_;

  StatisticsPack(const DistributionStatistics<ClusterLength> &length_statistics_,
                 const DistributionStatistics<ClusterCoverage> &coverage_statistics_) : length_statistics_(
      length_statistics_), coverage_statistics_(coverage_statistics_) {}
};

typedef std::map<size_t, StatisticsPack> StatisticsContainer;

class StatisticsSerializer {

 public:
    void PrintDistributions(const string& path, const DistributionPack& distributions) {
        auto length_stream = std::ofstream(fs::append_path(path, "length"));
        auto coverage_stream = std::ofstream(fs::append_path(path, "coverage"));
        length_stream << distributions.length_distribution_;
        coverage_stream << distributions.coverage_distribution_;
    }


    void PrintStatistics(const StatisticsContainer& statistics, const string& output_path) {
        auto stats_stream = ofstream(fs::append_path(output_path, "statistics"));
        for (const auto& entry: statistics) {
            stats_stream << "Distance " << entry.first << ':' << std::endl;
            stats_stream << "Length" << std::endl;
            stats_stream << entry.second.length_statistics_ << std::endl;
            stats_stream << "Coverage" << std::endl;
            stats_stream << entry.second.coverage_statistics_;
            stats_stream << std::endl << "------------" << std::endl;
        }
    }
};

class ClusterDistributionExtractor {
 private:
    const conj_graph_pack& gp_;
    size_t min_read_threshold_;
    size_t min_edge_length_;
    size_t min_cluster_offset_;
    size_t max_threads_;
 public:
    ClusterDistributionExtractor(const conj_graph_pack &gp_,
                                 size_t min_read_threshold_,
                                 size_t min_edge_length_,
                                 size_t min_cluster_offset_,
                                 size_t max_threads_)
        : gp_(gp_),
          min_read_threshold_(min_read_threshold_),
          min_edge_length_(min_edge_length_),
          min_cluster_offset_(min_cluster_offset_),
          max_threads_(max_threads_) {}

    DistributionPack GetDistributionsForDistance(size_t distance_threshold) {
        auto cluster_storage = GetInitialClusterStorage(distance_threshold);

//fixme remove code duplication by using boost::variant or smh else
//fixme lambdas everywhere
        auto cluster_predicate = [this](const cluster_storage::Cluster& cluster) {
          VERIFY_MSG(cluster.Size() == 1, "Cluster covers multiple edges!");
          path_extend::scaffold_graph::EdgeGetter edge_getter;
          auto map_info = cluster.GetMappings()[0];
          EdgeId cluster_edge = edge_getter.GetEdgeFromScaffoldVertex(map_info.GetEdge());
          size_t edge_length = this->gp_.g.length(cluster_edge);
          size_t left_pos = map_info.GetLeft();
          size_t right_pos = map_info.GetRight();
          return left_pos >= min_cluster_offset_ and right_pos + min_cluster_offset_ <= edge_length;
        };

        std::function<boost::optional<ClusterLength>(const cluster_storage::Cluster&)> length_extractor =
            [cluster_predicate](const cluster_storage::Cluster &cluster) {
                boost::optional<ClusterLength> result;
                if (cluster_predicate(cluster)) {
                    result = cluster.GetSpan();
                }
                return result;
            };
        std::function<boost::optional<ClusterCoverage>(const cluster_storage::Cluster&)> coverage_extractor =
            [cluster_predicate](const cluster_storage::Cluster &cluster) {
                boost::optional<ClusterCoverage> result;
                if (cluster_predicate(cluster)) {
                    result = cluster.GetCoverage();
                }
                return result;
            };
        SimpleDistributionExtractor distribution_extractor;

        DEBUG("Constructed extractors");

        auto length_distribution = distribution_extractor.ExtractDistribution(cluster_storage, length_extractor);
        auto coverage_distribution = distribution_extractor.ExtractDistribution(cluster_storage, coverage_extractor);

        DEBUG("Extracted distributions");

        DistributionPack result(length_distribution, coverage_distribution);
        cluster_storage.Clear();
        return result;
    }

    DistributionPack GetClusterDistributions() {
        INFO("Extracting read cloud cluster statistics");
        const size_t min_distance = 5000;
        const size_t max_distance = 41000;
        const size_t distance_step = 5000;
        vector<size_t> distances;
        for (size_t distance = min_distance; distance <= max_distance; distance += distance_step) {
            distances.push_back(distance);
        }
        StatisticsContainer statistics_container;
        for (size_t distance: distances) {
            DEBUG("Getting distributions for distance " << distance);
            auto distributions = GetDistributionsForDistance(distance);
            auto length_statistics = GetDistributionStatistics(distributions.length_distribution_);
            auto coverage_statistics = GetDistributionStatistics(distributions.coverage_distribution_);
            StatisticsPack statistics(length_statistics, coverage_statistics);
            statistics_container.insert({distance, statistics});
            DEBUG("Length mean: " << length_statistics.mean_);
            DEBUG("Length median: " << length_statistics.median_);
            DEBUG("Coverage mean: " << coverage_statistics.mean_);
            DEBUG("Coverage median: " << coverage_statistics.median_);
        }

        size_t optimal_distance = EstimateDistance(statistics_container);
        INFO("Estimated gap within cluster: " << optimal_distance);
        auto statistics = statistics_container.at(optimal_distance);
        INFO("Estimated mean cluster length: " << statistics.length_statistics_.mean_);
        INFO("Estimated median cluster length: " << statistics.length_statistics_.median_);
        INFO("Estimated median cluster coverage: " << statistics.coverage_statistics_.median_)
        return GetDistributionsForDistance(optimal_distance);
    }

 private:
    template<class T>
    DistributionStatistics<T> GetDistributionStatistics(SimpleDistribution<T>& distribution) const {
        T min_element = std::numeric_limits<T>::max();
        T max_element = 0;
        T element_sum = 0;
        if (not distribution.is_sorted()) {
            distribution.sort();
        }
        for (const T& element: distribution) {
            min_element = std::min(min_element, element);
            max_element = std::max(max_element, element);
            element_sum += element;
        }
        T mean = element_sum / static_cast<T>(distribution.size());
        T median = distribution.at(distribution.size() / 2);

        DistributionStatistics<T> result(min_element, max_element, mean, median);
        return result;
    }

    cluster_storage::ClusterStorage GetInitialClusterStorage(size_t distance_threshold) {
        auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
        omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(gp_.g);
        std::set<scaffold_graph::ScaffoldVertex> long_edges;
        for (const auto& edge: edge_it_helper) {
            if (gp_.g.length(edge) >= min_edge_length_) {
                long_edges.insert(edge);
            }
        }
        DEBUG("Found " << long_edges.size() << " long edges in the graph");
        auto edge_cluster_extractor =
            make_shared<cluster_storage::AccurateEdgeClusterExtractor>(gp_.g, barcode_extractor,
                                                                       distance_threshold, min_read_threshold_);
        cluster_storage::EdgeInitialClusterStorageBuilder initial_builder(gp_.g, edge_cluster_extractor, long_edges,
                                                                          distance_threshold, min_read_threshold_,
                                                                          max_threads_);
        auto initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();
        return initial_cluster_storage.get_cluster_storage();
    }

    size_t EstimateDistance(const StatisticsContainer& statistics) const {
        const double mean_growth_relative_threshold = 1.05;
        size_t prev_mean = 0;
        size_t current_distance = 0;
        VERIFY(statistics.size() > 0);
        for (const auto& entry: statistics) {
            current_distance = entry.first;
            size_t current_length_mean = entry.second.length_statistics_.mean_;
            size_t length_mean_threshold = static_cast<size_t>(static_cast<double>(prev_mean) * mean_growth_relative_threshold);
            if (prev_mean != 0 and current_length_mean <= length_mean_threshold) {
                INFO("Estimated fragment length median: " << current_length_mean);
                return current_distance;
            }
            prev_mean = current_length_mean;
        }
        return current_distance;
    }

    DECL_LOGGER("ClusterDistributionExtractor");
};

class ClusterStatisticsExtractor {
    DistributionPack cluster_distributions_;

 public:
    explicit ClusterStatisticsExtractor(const DistributionPack &cluster_distributions) :
        cluster_distributions_(cluster_distributions) {}

    DistributionPack GetDistributionPack() {
        return cluster_distributions_;
    }

    size_t GetLengthPercentile(double percent) {
        return GetPercentile(cluster_distributions_.length_distribution_, percent);
    }

    double GetCoveragePercentile(double percent) {
        return GetPercentile(cluster_distributions_.coverage_distribution_, percent);
    }

 private:
    template<class T>
    T GetPercentile(SimpleDistribution<T>& distribution, double percent) {
        if (not distribution.is_sorted()) {
            distribution.sort();
        }
        size_t index = static_cast<size_t>(static_cast<double>(distribution.size()) * percent);
        return distribution.at(index);
    }
};

class ClusterStatisticsExtractorHelper {
    const conj_graph_pack& gp_;
    const size_t max_threads_;

 public:
    ClusterStatisticsExtractorHelper(const conj_graph_pack &gp, size_t max_threads)
        : gp_(gp), max_threads_(max_threads) {}

    ClusterStatisticsExtractor GetStatisticsExtractor() const {
        const size_t DEFAULT_TRAINING_LENGTH = 10000;
        const size_t MIN_CLUSTER_OFFSET = 10000;
        const size_t LENGTH_TO_OFFSET = 5;
        const size_t MIN_READ_THRESHOLD = 5;

        MinTrainingLengthEstimatorHelper training_length_estimator_helper;
        const auto min_training_length_result = training_length_estimator_helper.EstimateTrainingLength(gp_.g);
        size_t min_training_length = DEFAULT_TRAINING_LENGTH;
        if (min_training_length_result.is_initialized()) {
            min_training_length = min_training_length_result.get();
        }

        const size_t min_cluster_offset = std::min(MIN_CLUSTER_OFFSET, min_training_length / LENGTH_TO_OFFSET);

        path_extend::cluster_model::ClusterDistributionExtractor distribution_analyzer(gp_,
                                                                                       MIN_READ_THRESHOLD,
                                                                                       min_training_length,
                                                                                       min_cluster_offset,
                                                                                       cfg::get().max_threads);
        auto cluster_distribution_pack = distribution_analyzer.GetClusterDistributions();
        cluster_model::ClusterStatisticsExtractor primary_parameters_extractor(cluster_distribution_pack);
        return primary_parameters_extractor;
    }
};
}
}