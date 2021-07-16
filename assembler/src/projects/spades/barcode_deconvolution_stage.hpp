#ifndef BARCODE_DECONVOLUTION_STAGE_HPP
#define BARCODE_DECONVOLUTION_STAGE_HPP

#include "assembly_graph/paths/path_processor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"
#include "pipeline/stage.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/memory_limit.hpp"
#include <algorithm>
#include <cstdlib>
#include <map>
#include <numeric>
#include <omp.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>

namespace debruijn_graph {

    /*
     Minimal required definitions for the read cloud deconvolution module. 
    */

    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;

    class BarcodeDeconvolutionStage : public spades::AssemblyStage {
    public:
        BarcodeDeconvolutionStage() : AssemblyStage("Ariadne", "barcode_deconvolution") {}
        void run(conj_graph_pack &gp, const char*);
    };

    /*
     Tools to process barcode strings and single reads.
    */

    inline int CalcCutoff(vector<int> cloud_sizes, float size_cutoff)
    {
        size_t size = cloud_sizes.size();
        if (size == 0){
            return 0;  // Undefined, really.
        }
        else{
            sort(cloud_sizes.begin(), cloud_sizes.end());
            int cutoff_index = (int)(size * size_cutoff);
            return cloud_sizes[cutoff_index];
        }
    }

    inline int SmallestVectorIndex(std::vector<std::vector<std::tuple<std::string, std::string, std::string, int>>>& v) {
        int smallest = INT_MAX;
        int index = -1;
        for ( size_t i = 0; i < v.size(); ++i ) {
            if ( v[i].size() < smallest ) {
                smallest = v[i].size();
                index = i;
            }
        }
        return index;
    }

    inline bool check_forward(int i) {
        return i == 0 || i % 2 == 0;
    }

    std::ofstream MakeOutputStream(std::string suffix) {
        std::ofstream stream;
        const size_t bufsize = 1024*1024;
        char buf[bufsize];
        stream.rdbuf()->pubsetbuf(buf, bufsize);
        std::string prefix = cfg::get().output_dir + std::to_string(cfg::get().search_distance);
        stream.open(prefix + suffix, std::ofstream::out);
        return stream;
    }

    inline std::string GetTenXBarcodeFromRead(const io::PairedRead &read) {
        std::string delimiter = "BX:Z:";
        std::string delimiter2 = "-1";
        size_t start_pos = read.first().name().find(delimiter);
        size_t delimiter_size = delimiter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.first().name().substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(delimiter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    inline std::string GetTenXBarcodeFromRead(const std::string &read) {
        std::string delimiter = "BX:Z:";
        std::string delimiter2 = "-1";
        size_t start_pos = read.find(delimiter);
        size_t delimiter_size = delimiter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(delimiter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    void MakeRead(io::SingleRead& read,
                  debruijn_graph::conj_graph_pack& gp,
                  std::vector<std::vector<std::vector<int>>>& connected_reads,
                  std::vector<std::vector<std::pair<int, MappingPath<EdgeId>>>>& mapping_record,
                  std::vector<std::vector<std::tuple<std::string, std::string, std::string, int>>>& read_record,
                  int read_index)
    // Store the reads in the analysis structures (tmp_connected and tmp_mapping) and the reporting structure (tmp_reads).
    {
        auto mapper = MapperInstance(gp);
        auto path = mapper->MapRead(read);
        connected_reads.back().emplace_back( std::vector<int>() );
        mapping_record.back().emplace_back( read_index, path );
        if ( check_forward(read_index) ) {
            connected_reads.back().back().push_back( read_index + 1 ); // Connect first read to second read.
            read_record.back().emplace_back( std::make_tuple( read.name(), read.GetSequenceString(),
                                                              read.GetPhredQualityString(), 1) );
        } else {
            connected_reads.back().back().push_back( read_index - 1 ); // Connect second read to first read.
            read_record.back().emplace_back( std::make_tuple( read.name(), Reverse(Complement(read.GetSequenceString())),
                                                              Reverse(read.GetPhredQualityString()), 2) );
        }
    }

    std::vector<VertexId> VerticesReachedFrom(VertexId& start_vertex,
                                              debruijn_graph::conj_graph_pack &gp, int edge_size)
    // Find all vertices reachable within the search distance starting from the 3'-most vertex of the read's mapping path.
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);

        return bounded_dijkstra.ReachedVertices();
    }
    std::vector<VertexId> ConjugateVerticesReachedFrom(VertexId& start_vertex,
                                                       debruijn_graph::conj_graph_pack &gp, int edge_size)
    // Find all vertices reachable within the search distance with the 3'-most vertex of the read's mapping path as the destination.
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);

        return bounded_dijkstra.ReachedVertices();
    }

    inline std::string GetUpdatedReadName(std::string& read, size_t cloud_number)
    // Updates original read cloud number to enhanced read cloud number.
    {
        size_t start_pos = read.find("-1"); // The original cloud number that comes with the barcodes.
        if (start_pos != string::npos) {
            return read.substr(0, start_pos) + "-" + std::to_string(cloud_number);
        } else {
            // If no barcode, get rid of the '_RC' or '_SUBSTR' that cloudSPAdes tags on and report plain read-name.
            start_pos = read.find("_RC");
            std::string tmp = read.substr(0, start_pos);
            start_pos = tmp.find("_SUBSTR");
            return tmp.substr(0, start_pos);
        }
    }

    inline std::string MakeFastQString(std::string& name, std::string& seq, std::string& qual)
    // Concatenates read information into FastQ format.
    {
        std::string read_string = "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
        return read_string;
    }

    void ReformatYAMLs()
    // Move old input_dataset.yaml to a different name, and write a new version of it so that
    // barcode_index_construction.cpp will pick up on the enhanced FastQs instead of the original FastQs.
    {
        std::string out_dir = cfg::get().output_base;
        std::string old_name_str = out_dir + "input_dataset.yaml";
        char old_name[old_name_str.length() + 1];
        strcpy(old_name, old_name_str.c_str());
        std::string new_name_str = out_dir + "input_dataset_original.yaml";
        char new_name[new_name_str.length() + 1];
        strcpy(new_name, new_name_str.c_str());
        std::rename(old_name, new_name);

        std::ofstream stream;
        stream.open(old_name_str, std::ofstream::out);
        std::string prefix = cfg::get().output_dir + std::to_string(cfg::get().search_distance);
        stream << "- \"left reads\":\n" << "  - \"" << prefix << ".R1.fastq\"\n  \"orientation\": \"fr\"\n"
               << "  \"right reads\":\n" << "  - \"" << prefix << ".R2.fastq\"\n  \"type\": \"clouds10x\"";
        stream.close();
    }

    /*
     Workhorse read information processing functions.
    */

    int CountReads(const lib_t& lib_10x)
    // Count number of reads in the dataset.
    {
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;
        int num_reads_total = 0;

        while ( !stream->eof() ) {
            *stream >> read;
            num_reads_total += 2;
        }
        stream->close();
        INFO(num_reads_total << " reads to process");
        return num_reads_total;
    }

    void OutputReads(std::vector<std::tuple<std::string, std::string, std::string, int>>& master_record,
                     std::ofstream& fastq_stream_forward,
                     std::ofstream& fastq_stream_reverse)
    // Report read clouds.
    {
        std::map<std::string,std::tuple<std::string, std::string, std::string>> holding_space;
        std::string fastq_string_forward = "";
        std::string fastq_string_reverse = "";
        while ( !master_record.empty() ){
            auto curr_read = master_record.front();
            auto read_name = std::get<0>(curr_read);
            if ( holding_space.find(read_name) == holding_space.end() ) { // Current read is first read in the set.
                holding_space[read_name] = std::make_tuple( std::get<0>(curr_read), std::get<1>(curr_read), std::get<2>(curr_read) );
            } else { // Current read is second read in the set. Output both.
                auto prev_read = holding_space[read_name];
                if (std::get<3>(curr_read) == 2) {
                    fastq_string_forward += MakeFastQString(std::get<0>(prev_read), std::get<1>(prev_read), std::get<2>(prev_read));
                    fastq_string_reverse += MakeFastQString(std::get<0>(curr_read), std::get<1>(curr_read), std::get<2>(curr_read));
                } else {
                    fastq_string_forward += MakeFastQString(std::get<0>(curr_read), std::get<1>(curr_read), std::get<2>(curr_read));
                    fastq_string_reverse += MakeFastQString(std::get<0>(prev_read), std::get<1>(prev_read), std::get<2>(prev_read));
                }
                holding_space.erase(read_name);
            }
            master_record.erase(master_record.begin());
        }
        fastq_stream_forward << fastq_string_forward << flush;
        fastq_stream_reverse << fastq_string_reverse << flush;
    }

    int LoadReads(std::vector<std::vector<std::vector<int>>>& connected_reads,
                  std::vector<std::vector<std::pair<int, MappingPath<EdgeId>>>>& mapping_record,
                  std::vector<std::vector<std::tuple<std::string, std::string, std::string, int>>>& read_record,
                  debruijn_graph::conj_graph_pack& graph_pack,
                  const lib_t& lib_10x, int num_reads_start, int num_reads_goal,
                  std::ofstream& fastq_stream_forward,
                  std::ofstream& fastq_stream_reverse)
    // Divides reads from library in PairedEasyStream into the original read clouds.
    {
        // Tools to locate reads along the compacted de Bruijn graph.
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;

        // Initialize barcode string and loop-tracker values.
        std::string current_barcode;
        int num_cloud_reads = 0;
        int num_reads_total = 0;
        // std::vector<std::tuple<std::string, std::string, std::string, int>> unbarcoded_record; 

        while ( !stream->eof() ) {
            *stream >> read;
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if ( barcode_string.empty() ) { // If these reads do not have barcodes, skip.
                num_reads_total += 2;
                continue;
            }
            if ( num_reads_total + 2 <= num_reads_start ) { // If these reads have already been loaded, skip.
                num_reads_total += 2;
                continue;
            }
            if ( num_reads_total + 2 > num_reads_goal ) { // If loading two more reads puts the number over the limit, stop.
                break;
            }
            num_reads_total += 2;
            // if ( barcode_string.empty() ) { // This read is not barcoded. Output straight away.
            //     std::string first_in_pair = read.first().name(); 
            //     std::string second_in_pair = read.second().name(); 
            //     unbarcoded_record.emplace_back( std::make_tuple( GetUpdatedReadName(first_in_pair, 0), 
            //                                                      read.first().GetSequenceString(),
            //                                                      read.first().GetPhredQualityString(), 1) );
            //     unbarcoded_record.emplace_back( std::make_tuple( GetUpdatedReadName(second_in_pair, 0), 
            //                                                      Reverse(Complement(read.second().GetSequenceString())),
            //                                                      Reverse(read.second().GetPhredQualityString()), 2) );
            // } else {
            if ( barcode_string != current_barcode ){ // This barcoded read belongs to the next read-cloud. Start new storage objects.
                connected_reads.emplace_back( std::vector<std::vector<int>>() );
                mapping_record.emplace_back( std::vector<std::pair<int, MappingPath<EdgeId>>>() );
                read_record.emplace_back( std::vector<std::tuple<std::string, std::string, std::string, int>>() );
                num_cloud_reads = 0;
            }
            // For each pair of reads, map them to the assembly graph and extract their FastQ information.
            MakeRead(read.first(), graph_pack, connected_reads, mapping_record, read_record, num_cloud_reads);
            MakeRead(read.second(), graph_pack, connected_reads, mapping_record, read_record, num_cloud_reads + 1);

            num_cloud_reads += 2;
            current_barcode = barcode_string;
            // }
        }
        stream->close();
        // INFO(unbarcoded_record.size() << " unbarcoded reads, now outputting");
        // OutputReads(unbarcoded_record, fastq_stream_forward, fastq_stream_reverse);
        int num_reads_loaded = num_reads_total - num_reads_start;
        INFO(num_reads_loaded << " reads loaded");
        return num_reads_total;
    }

    void ClusterReads(std::vector<std::vector<int>>& tmp_connected,
                      std::vector<std::pair<int, MappingPath<EdgeId>>>& tmp_mapping,
                      debruijn_graph::conj_graph_pack &gp)
    // Generates all connections between reads in the same read cloud by constructing a Djikstra graph from accessible vertices and traversing it
    // to find other reads with mapping paths containing those vertices.
    {
        int search_dist = cfg::get().search_distance;
        for (size_t i = 0; i < tmp_mapping.size(); ++i) { // For each read i...
            auto read_1 = tmp_mapping[i].second;
            if (read_1.size()){  // Check to see if the first read has a mapping path.
                EdgeId read_1_edge_end = read_1.edge_at(read_1.size() - 1);
                EdgeId read_1_edge_start = read_1.edge_at(0);
                for(size_t j = i + 1; j < tmp_mapping.size(); ++j) { // ...check to see if all other reads j can be connected to it.
                    auto read_2 = tmp_mapping[j].second;
                    if ( read_2.size() && !(check_forward(i) && j == i + 1) ){  // Check to see if read j has a mapping path and isn't the paired read.
                        EdgeId read_2_edge_end = read_2.edge_at(read_2.size() - 1);
                        EdgeId read_2_edge_start = read_2.edge_at(0);

                        if (read_1_edge_end == read_2_edge_start ||
                            read_1_edge_end == gp.g.conjugate(read_2_edge_end)) { // Does read j start on the same edge that read i ends?
                            long int read_distance = read_2.start_pos() - read_1.end_pos();
                            if (std::abs(read_distance) < search_dist && read_distance >= 0){
                                tmp_connected[ tmp_mapping[i].first ].push_back( tmp_mapping[j].first );
                                tmp_connected[ tmp_mapping[j].first ].push_back( tmp_mapping[i].first );
                            }
                        }

                        else if (read_1_edge_start == read_2_edge_end ||
                                 read_1_edge_start == gp.g.conjugate(read_2_edge_start)) { // Does read j end on the same edge that read i starts?
                            int path2_conj_end = gp.g.length(read_2_edge_end) - read_2.end_pos();
                            long int read_distance = path2_conj_end - read_1.end_pos();
                            if (std::abs(read_distance) < search_dist && read_distance >= 0) {
                                tmp_connected[ tmp_mapping[i].first ].push_back( tmp_mapping[j].first );
                                tmp_connected[ tmp_mapping[j].first ].push_back( tmp_mapping[i].first );
                            }
                        }

                        else { // Otherwise, is the read between read j traversable within the Djikstra graph?

                            // Look for the set of all vertices reachable from read 1's 3'-most vertex.
                            VertexId startVertex = gp.g.EdgeEnd(read_1.back().first); // 3'-most vertex of the read's mapping path.
                            std::vector<VertexId> reached_vertices;
                            int endDist = gp.g.length(read_1_edge_end) - read_1.end_pos(); // Distance between end of read and 3'-most vertex.
                            int reducedEndDist = search_dist - endDist; // The amount of distance to search starting from the 3'-most vertex.
                            if (reducedEndDist > 0){
                                reached_vertices = VerticesReachedFrom(startVertex, gp, reducedEndDist); // Find the list of vertices that can be reached with a read of MappingPath<EdgeId>*.;
                            }
                            std::sort(reached_vertices.begin(), reached_vertices.end()); // A sorted list of VertexIDs, which are just integers.

                            // Look for the set of all vertices that can reach read 1's 5'-most vertex.
                            VertexId conjStartVertex = gp.g.EdgeStart(read_1.front().first); // 5'-most vertex of the read's mapping path.
                            std::vector<VertexId> conjugate_reached_vertices;
                            int startDist = read_1.mapping_at(0).mapped_range.start_pos;
                            int reducedStartDist = search_dist - startDist;
                            if (reducedStartDist > 0){
                                conjugate_reached_vertices = ConjugateVerticesReachedFrom(conjStartVertex, gp, reducedStartDist);
                            }
                            std::sort(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end());

                            // Check each end-vertex of read 2 and its reverse complement to see if it connects to read 1.
                            for (size_t k = 0; k < read_2.size(); ++k) {
                                auto read_2_vertex = gp.g.EdgeEnd(read_2[k].first);
                                // This is the end vertex of an edge that the second read sits on. Assuming the read is going in the 'forward' direction.
                                if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), read_2_vertex) ||
                                    std::binary_search(conjugate_reached_vertices.begin(),
                                                       conjugate_reached_vertices.end(), read_2_vertex)) {
                                    tmp_connected[tmp_mapping[i].first].push_back(tmp_mapping[j].first);
                                    tmp_connected[tmp_mapping[j].first].push_back(tmp_mapping[i].first);
                                    break;
                                } else {
                                    read_2_vertex = gp.g.conjugate(gp.g.EdgeEnd(read_2[k].first));
                                    // This is the end vertex of an edge that the reverse conjugate of the second read sits on.
                                    if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), read_2_vertex) ||
                                        std::binary_search(conjugate_reached_vertices.begin(),
                                                           conjugate_reached_vertices.end(), read_2_vertex)) {
                                        tmp_connected[tmp_mapping[i].first].push_back(tmp_mapping[j].first);
                                        tmp_connected[tmp_mapping[j].first].push_back(tmp_mapping[i].first);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<std::tuple<std::string, std::string, std::string, int>> MakeEnhanced(std::vector<std::vector<int>>& tmp_connected,
                     std::vector<std::tuple<std::string, std::string, std::string, int>>& tmp_reads, int cloud_size_filter,
                     std::ofstream& enh_cloud_stats)
    // Cluster reads that are in the same enhanced read by depth-first traversing through sets of connected reads.
    {
        std::vector<std::tuple<std::string, std::string, std::string, int>> master_record;
        std::string barcode = GetTenXBarcodeFromRead(std::get<0>(tmp_reads[0]));

        if ( tmp_reads.size() >= cloud_size_filter ) { // If the original read cloud is larger than the mean...
            std::vector<int> cloud_queue;
            std::vector<int> visited_reads;
            std::vector<std::vector<std::tuple<std::string, std::string, std::string, int>>> enhanced_record;
            std::vector<std::tuple<std::string, std::string, std::string, int>> unenhanced_record;

            for (size_t i = 0; i < tmp_connected.size(); ++i) {
                auto read_find = std::find(visited_reads.begin(), visited_reads.end(), i);
                if (read_find == visited_reads.end()) { // If this read hasn't already been grouped, process it.
                    size_t num_reads = 0; // The number of reads in the enhanced cloud.
                    visited_reads.push_back(i); // Add read index to vector of already-grouped reads.
                    cloud_queue.push_back(i);
                    std::vector<std::tuple<std::string, std::string, std::string, int>> temporary_enhanced;
                    do {
                        int &cloud_index = cloud_queue.front(); // Get the unique index of read j.
                        temporary_enhanced.emplace_back(tmp_reads[cloud_index]);
                        ++num_reads;
                        for (int j = 0; j <
                                        tmp_connected[cloud_index].size(); ++j) // Do this for every read k that read j is connected to.
                        {
                            auto current_index = std::find(visited_reads.begin(), visited_reads.end(),
                                                           tmp_connected[cloud_index][j]);
                            if (current_index == visited_reads.end()) { // If the current read has not been grouped...
                                auto additional_read = tmp_connected[cloud_index][j];
                                cloud_queue.emplace_back(additional_read); // ...add it to the cloud queue
                                visited_reads.emplace_back(
                                        additional_read); // ...add it to the vector of grouped reads.
                            }
                        }
                        cloud_queue.erase(cloud_queue.begin());
                    } while (!cloud_queue.empty());
                    if (num_reads == 2) { // Just a pair of reads. Add to the un-enhanced read cloud.
                        unenhanced_record.insert(std::end(unenhanced_record), std::begin(temporary_enhanced),
                                                 std::end(temporary_enhanced));
                    } else {
                        enhanced_record.emplace_back(temporary_enhanced);
                    }
                }
            }

            if (unenhanced_record.size() > 0) {
                if (unenhanced_record.size() == 2) { // If the un-enhanced cloud is just a pair...
                    int smallest_index = SmallestVectorIndex(enhanced_record);
                    for (auto &read : enhanced_record[smallest_index]) { // ...add unenhanced reads to the smallest enhanced cloud
                        unenhanced_record.push_back(read);
                    }
                    enhanced_record.erase(enhanced_record.begin() + smallest_index);
                }
                enh_cloud_stats << barcode << "," << 0 << "," << unenhanced_record.size()
                                << std::endl; // Reporting un-enhanced read cloud information.
                for (auto &read : unenhanced_record) {
                    master_record.emplace_back(std::make_tuple(GetUpdatedReadName(std::get<0>(read), 0),
                                                               std::get<1>(read), std::get<2>(read), std::get<3>(read)));
                }
            }
            for (size_t i = 0; i < enhanced_record.size(); ++i) {
                int enh_index = i + 1;
                enh_cloud_stats << barcode << "," << enh_index << "," << enhanced_record[i].size()
                                << std::endl; // Reporting un-enhanced read cloud information.
                for (auto &read : enhanced_record[i]) {
                    master_record.emplace_back(std::make_tuple(GetUpdatedReadName(std::get<0>(read), enh_index),
                                                               std::get<1>(read), std::get<2>(read), std::get<3>(read)));
                }
            }
        } else { // Otherwise, add the whole original read cloud to the output stream.
            for ( auto& read : tmp_reads ) {
                master_record.emplace_back(std::make_tuple( GetUpdatedReadName(std::get<0>(read), 0),
                                                                std::get<1>(read), std::get<2>(read), std::get<3>(read)) );
            }
            // if ( barcode.empty() ) barcode = "NA"; // If this is an unbarcoded read-pair, set it to NA.
            enh_cloud_stats << barcode << "," << 0 << "," << master_record.size() << std::endl; // Reporting un-enhanced read cloud information.
        }
        return master_record;
    }

    /*
     main() function.
     */

    void BarcodeDeconvolutionStage::run(conj_graph_pack &gp, const char*)
    // Run the read cloud deconvolution module processReads() for every read library in the dataset.
    {
        gp.EnsureIndex();
        if (!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
        INFO("Read cloud deconvolution starting");
        config::dataset& dataset_info = cfg::get_writable().ds;
        lib_t& lib_10x = dataset_info.reads[0];

        // Output stream objects.
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        std::ofstream fastq_stream_forward = MakeOutputStream(".R1.fastq");
        std::ofstream fastq_stream_reverse = MakeOutputStream(".R2.fastq");
        std::ofstream stat_stream = MakeOutputStream(".summary.csv");

        // Collection of read information from the original read clouds
        std::vector<std::vector<std::vector<int>>> connected_reads;
        std::vector<std::vector<std::pair<int, MappingPath<EdgeId>>>> mapping_record;
        std::vector<std::vector<std::tuple<std::string, std::string, std::string, int>>> read_record;

        // Setting memory-based limits on read-chunks loaded at once.
        int num_loadable_reads = (int)((utils::get_free_memory() * 1.0) / 850 ); // Available memory / estimated memory per read in bytes
        int num_reads_section_end = 0;
        int num_reads_section_goal = num_loadable_reads;
        int num_reads_total = CountReads(lib_10x);
        int cloud_size_filter = cfg::get().size_cutoff;
        INFO("Cutoff for original cloud size is " << cloud_size_filter << " reads");

        while (num_reads_section_end < num_reads_total ) {
            INFO("Loading original read clouds from library");
            num_reads_section_end = LoadReads(connected_reads, mapping_record, read_record, gp, lib_10x, num_reads_section_end, num_reads_section_goal, fastq_stream_forward, fastq_stream_reverse);
            size_t num_original_clouds = mapping_record.size();
            INFO(num_original_clouds << " original read clouds loaded");
#pragma omp parallel for shared(connected_reads, mapping_record, read_record) schedule(dynamic, 1) num_threads(cfg::get().max_threads)
            for (size_t c = 0; c < num_original_clouds; ++c) {
                if (c % 100000 == 0) {
                    INFO(c << ": Deconvolving and reporting " << mapping_record[c].size() << " reads on thread " << omp_get_thread_num());
                }
                std::string barcode_string = GetTenXBarcodeFromRead(std::get<0>(read_record[c][0]));
                if ( barcode_string != "NA" && mapping_record[c].size() >= cloud_size_filter ) {
                    ClusterReads(connected_reads[c], mapping_record[c], gp);
                }
#pragma omp critical
                {
                    std::vector<std::tuple<std::string, std::string, std::string, int>> master_record = MakeEnhanced(connected_reads[c], read_record[c], cloud_size_filter, stat_stream);
                    OutputReads(master_record, fastq_stream_forward, fastq_stream_reverse);
                }
            }
            connected_reads.clear();
            mapping_record.clear();
            read_record.clear();
            num_reads_section_goal = num_reads_section_end + num_loadable_reads;
        }
        fastq_stream_forward.close();
        fastq_stream_reverse.close();
        stat_stream.close();
        // ReformatYAMLs();
        INFO("Read cloud deconvolution finished");
    }

}

#endif