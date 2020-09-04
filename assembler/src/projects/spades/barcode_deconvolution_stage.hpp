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
#include <numeric>
#include <omp.h>
#include <string>
#include <sstream>
#include <unordered_map>
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

    inline int SmallestVectorIndex(std::vector<std::vector<std::tuple<std::string, std::string, std::string>>>& v) {
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
                  std::vector<std::vector<std::tuple<std::string, std::string, std::string>>>& read_record,
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
                                                              read.GetPhredQualityString() ) );
        } else {
            connected_reads.back().back().push_back( read_index - 1 ); // Connect second read to first read.
            read_record.back().emplace_back( std::make_tuple( read.name(), Complement(read.GetSequenceString()),
                                                              Reverse(read.GetPhredQualityString()) ) );
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
            // If no barcode, get rid of the '_RC' that metaSPAdes tags on and report plain read-name.
            start_pos = read.find("_RC");
            return read.substr(0, start_pos);
        }
    }

    inline std::string MakeFastQString(std::string& name, std::string& seq, std::string& qual)
    // Concatenates read information into FastQ format.
    {
        std::string read_string = "@" + name + '\n' + seq + "\n+\n" + qual;
        return read_string;
    }

    /*
     Workhorse read information processing functions.
    */

    int SetCloudFilter(debruijn_graph::conj_graph_pack& graph_pack,
                        const lib_t& lib_10x, int& num_reads_total)
    // Calculate the mean size of the original read clouds to serve as a filter for processing read clouds later.
    {
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;

        // Initialize barcode string and loop-tracker values.
        std::string current_barcode;
        int num_cloud_reads = 0;
        std::vector<int> cloud_sizes;

        while ( !stream->eof() ) {
            *stream >> read;
            num_reads_total += 2;
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if ( !barcode_string.empty() ){
                if ( barcode_string != current_barcode && !current_barcode.empty() ){
                    cloud_sizes.push_back(num_cloud_reads);
                    num_cloud_reads = 0;
                }
                num_cloud_reads += 2;
                current_barcode = barcode_string;
            }
        }
        cloud_sizes.push_back(num_cloud_reads);
        stream->close();
        double cloud_mean_size = std::accumulate(cloud_sizes.begin(), cloud_sizes.end(), 0) / cloud_sizes.size();
        INFO(num_reads_total << " reads to process");
        return (int)cloud_mean_size;
    }

    int LoadReads(std::vector<std::vector<std::vector<int>>>& connected_reads,
                  std::vector<std::vector<std::pair<int, MappingPath<EdgeId>>>>& mapping_record,
                  std::vector<std::vector<std::tuple<std::string, std::string, std::string>>>& read_record, std::vector<std::tuple<std::string, std::string, std::string>>& unbarcoded_record,
                  debruijn_graph::conj_graph_pack& graph_pack,
                  const lib_t& lib_10x, int num_reads_start, int num_reads_goal)
    // Divides reads from library in PairedEasyStream into the original read clouds.
    {
        // Tools to locate reads along the compacted de Bruijn graph.
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;

        // Initialize barcode string and loop-tracker values.
        std::string current_barcode;
        int num_cloud_reads = 0;
        int num_reads_total = 0;

        while ( !stream->eof() ) {
            *stream >> read;
            num_reads_total += 2;
            if (num_reads_total <= num_reads_start) {
                continue;
            }
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if ( !barcode_string.empty() ){
                if ( barcode_string != current_barcode ){ // This read belongs to the next read-cloud. Start new storage objects.
                    if ( num_reads_total >= num_reads_goal ) {
                        break;
                    }
                    connected_reads.emplace_back( std::vector<std::vector<int>>() );
                    mapping_record.emplace_back( std::vector<std::pair<int, MappingPath<EdgeId>>>() );
                    read_record.emplace_back( std::vector<std::tuple<std::string, std::string, std::string>>() );
                    num_cloud_reads = 0;
                }
                // For each pair of reads, map them to the assembly graph and extract their FastQ information.
                MakeRead(read.first(), graph_pack, connected_reads, mapping_record, read_record, num_cloud_reads);
                MakeRead(read.second(), graph_pack, connected_reads, mapping_record, read_record, num_cloud_reads + 1);

                num_cloud_reads += 2;
                current_barcode = barcode_string;
            }
            else { // This read does not have a barcode, and will be reported at the end.
                unbarcoded_record.emplace_back( std::make_tuple( read.first().name(), read.first().GetSequenceString(),
                                                                 read.first().GetPhredQualityString() ) );
                std::string reverse_name = read.second().name();
                size_t start_pos = reverse_name.find("_RC");
                unbarcoded_record.emplace_back( std::make_tuple( reverse_name.substr(0, start_pos), Complement(read.second().GetSequenceString()),
                                                                 Reverse(read.second().GetPhredQualityString()) ) );
            }
        }
        stream->close();
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
                            read_1_edge_end == gp.g.conjugate(read_2_edge_end) || // Does read j start on the same edge that read i ends?
                            read_1_edge_start == read_2_edge_end ||
                            read_1_edge_start == gp.g.conjugate(read_2_edge_start) // Does read j end on the same edge that read i starts?
                                ) { // During previous tests, too few reads are being connected. Test to see if this fixes it.
                            tmp_connected[ tmp_mapping[i].first ].push_back( tmp_mapping[j].first );
                            tmp_connected[ tmp_mapping[j].first ].push_back( tmp_mapping[i].first );
                        } else { // Otherwise, is the read between read j traversable within the Djikstra graph?

                            // Look for the set of all vertices reachable from read 1's 3'-most vertex.
                            VertexId startVertex = gp.g.EdgeEnd(read_1.back().first); // 3'-most vertex of the read's mapping path.
                            std::vector<VertexId> reached_vertices;
                            int endDist = gp.g.length(read_1_edge_end) - read_1.end_pos(); // Distance between end of read and 3'-most vertex.
                            int reducedEndDist = cfg::get().barcode_distance - endDist; // The amount of distance to search starting from the 3'-most vertex.
                            if (reducedEndDist > 0){
                                reached_vertices = VerticesReachedFrom(startVertex, gp, reducedEndDist); // Find the list of vertices that can be reached with a read of MappingPath<EdgeId>*.;
                            }
                            std::sort(reached_vertices.begin(), reached_vertices.end()); // A sorted list of VertexIDs, which are just integers.

                            // Look for the set of all vertices that can reach read 1's 5'-most vertex.
                            VertexId conjStartVertex = gp.g.EdgeStart(read_1.front().first); // 5'-most vertex of the read's mapping path.
                            std::vector<VertexId> conjugate_reached_vertices;
                            int startDist = read_1.mapping_at(0).mapped_range.start_pos;
                            int reducedStartDist = cfg::get().barcode_distance - startDist;
                            if (reducedStartDist > 0){
                                conjugate_reached_vertices = ConjugateVerticesReachedFrom(conjStartVertex, gp, reducedStartDist);
                            }
                            std::sort(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end());

                            // Evaluate the 5' and 3'-most vertices of read 2.
                            VertexId forward_end_three = gp.g.EdgeEnd(read_2_edge_end);
                            // Assuming read is same-strand. The 3'-most vertex of the last edge that read 2 sits on.
                            VertexId reverse_end_three = gp.g.conjugate(gp.g.EdgeStart(read_2_edge_start));
                            // Assuming read is opposite-strand. The 5'-most vertex of the first edge that read 2 sits on.
                            VertexId forward_end_five = gp.g.EdgeStart(read_2_edge_start);
                            // Assuming read is same-strand. The 5'-most vertex of the last edge that read 2 sits on.
                            VertexId reverse_end_five = gp.g.conjugate(gp.g.EdgeEnd(read_2_edge_end));
                            // Assuming read is opposite-strand. The 3'-most vertex of the last edge that read 2 sits on.

                            if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), forward_end_five) ||
                                std::binary_search(reached_vertices.begin(), reached_vertices.end(), reverse_end_five) ||
                                // If the 5'-most vertex can be found in the Djikstra graph, then add an edge between the two reads.
                                std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), forward_end_three) ||
                                std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), forward_end_three)) {
                                // If the 3'-most vertex can be found in the conjugate Djikstra graph, then add an edge between the two reads.
                                tmp_connected[ tmp_mapping[i].first ].push_back( tmp_mapping[j].first );
                                tmp_connected[ tmp_mapping[j].first ].push_back( tmp_mapping[i].first );
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<std::tuple<std::string, std::string, std::string>> MakeEnhanced(std::vector<std::vector<int>>& tmp_connected,
                     std::vector<std::tuple<std::string, std::string, std::string>>& tmp_reads, int cloud_size_filter,
                     std::ofstream& enh_cloud_stats)
    // Cluster reads that are in the same enhanced read by depth-first traversing through sets of connected reads.
    {
        std::vector<std::tuple<std::string, std::string, std::string>> master_record;
        std::string barcode = GetTenXBarcodeFromRead(std::get<0>(tmp_reads[0]));

        if ( tmp_reads.size() >= cloud_size_filter ) { // If the original read cloud is larger than the mean...
            std::vector<int> cloud_queue;
            std::vector<int> visited_reads;
            std::vector<std::vector<std::tuple<std::string, std::string, std::string>>> enhanced_record;
            std::vector<std::tuple<std::string, std::string, std::string>> unenhanced_record;

            for (size_t i = 0; i < tmp_connected.size(); ++i) {
                auto read_find = std::find(visited_reads.begin(), visited_reads.end(), i);
                if (read_find != visited_reads.end()) { // If this read has already been grouped, skip.
                    continue;
                } else {
                    size_t num_reads = 0; // The number of reads in the enhanced cloud.
                    visited_reads.push_back(i); // Add read index to vector of already-grouped reads.
                    cloud_queue.push_back(i);
                    std::vector<std::tuple<std::string, std::string, std::string>> temporary_enhanced;
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

            if (unenhanced_record.size() == 2) { // If the un-enhanced cloud is just a pair...
                if (enhanced_record.size() == 1) { // ...and there is only one enhanced cloud, no enhancement done.
                    unenhanced_record.insert(std::end(unenhanced_record), std::begin(enhanced_record[0]),
                                             std::end(enhanced_record[0]));
                    enhanced_record.clear();
                } else if (enhanced_record.size() >
                           1) { // ...and there are multiple clouds, add unenhanced reads to the smallest enhanced cloud.
                    for (auto &read : unenhanced_record) {
                        enhanced_record[SmallestVectorIndex(enhanced_record)].push_back(read);
                    }
                    unenhanced_record.clear();
                }
            } // If there are no pairs that were unenhanced, nothing happens.

            enh_cloud_stats << barcode << "," << 0 << "," << unenhanced_record.size()
                            << std::endl; // Reporting un-enhanced read cloud information.
            for (auto &read : unenhanced_record) {
                master_record.emplace_back(std::make_tuple(GetUpdatedReadName(std::get<0>(read), 0),
                                                           std::get<1>(read), std::get<2>(read)));
            }
            for (size_t i = 0; i < enhanced_record.size(); ++i) {
                int enh_index = i + 1;
                enh_cloud_stats << barcode << "," << enh_index << "," << enhanced_record[i].size()
                                << std::endl; // Reporting un-enhanced read cloud information.
                for (auto &read : enhanced_record[i]) {
                    master_record.emplace_back(std::make_tuple(GetUpdatedReadName(std::get<0>(read), enh_index),
                                                               std::get<1>(read), std::get<2>(read)));
                }
            }
        } else { // Otherwise, add the whole original read cloud to the output stream.
            for ( auto& read : tmp_reads ) {
                master_record.emplace_back(std::make_tuple( GetUpdatedReadName(std::get<0>(read), 0),
                                                                std::get<1>(read), std::get<2>(read)) );
            }
            enh_cloud_stats << barcode << "," << 0 << "," << master_record.size() << std::endl; // Reporting un-enhanced read cloud information.
        }
//        enh_cloud_stats.flush();
        return master_record;
    }

    void OutputReads(std::vector<std::tuple<std::string, std::string, std::string>>& master_record,
                     std::ofstream& fastq_stream_forward,
                     std::ofstream& fastq_stream_reverse)
    // Report read clouds.
    {
        int direction = 0;
        while ( !master_record.empty() ){
            auto read = master_record.front();
            if ( check_forward(direction) ) { // If first read in pair, put in forward file.
                fastq_stream_forward << MakeFastQString(std::get<0>(read), std::get<1>(read), std::get<2>(read)) << std::endl;
            } else { // If second read in pair, put in reverse complement file.
                fastq_stream_reverse << MakeFastQString(std::get<0>(read), std::get<1>(read), std::get<2>(read)) << std::endl;
            }
            ++direction;
            master_record.erase(master_record.begin());
        }
//        fastq_stream_forward.flush();
//        fastq_stream_reverse.flush();
    }

    /*
     main() function.
     */

    void BarcodeDeconvolutionStage::run(conj_graph_pack &gp, const char*)
    // Run the read cloud deconvolution module processReads() for every read library in the dataset.
    {
        gp.EnsureIndex();
        if (!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
        INFO("There are " << gp.g.size() << " vertices in the assembly graph");
        INFO("Read cloud deconvolution starting");
        config::dataset& dataset_info = cfg::get_writable().ds;
        lib_t& lib_10x = dataset_info.reads[0];

        // Output stream objects.
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        std::ofstream fastq_stream_forward, fastq_stream_reverse, stat_stream;
        const size_t bufsize = 1024*1024;
        char buf[bufsize];
        fastq_stream_forward.rdbuf()->pubsetbuf(buf, bufsize);
        fastq_stream_reverse.rdbuf()->pubsetbuf(buf, bufsize);
        stat_stream.rdbuf()->pubsetbuf(buf, bufsize);
        std::string file_name = cfg::get().output_dir + std::to_string(cfg::get().barcode_distance);
        fastq_stream_forward.open (file_name + ".R1.fastq", std::ofstream::out);
        fastq_stream_reverse.open (file_name + ".R2.fastq", std::ofstream::out);
        stat_stream.open (file_name + ".summary.csv", std::ofstream::out);

        // Collection of read information from the original read clouds
        std::vector<std::vector<std::vector<int>>> connected_reads;
        std::vector<std::vector<std::pair<int, MappingPath<EdgeId>>>> mapping_record;
        std::vector<std::vector<std::tuple<std::string, std::string, std::string>>> read_record;
        std::vector<std::tuple<std::string, std::string, std::string>> unbarcoded_record;

        // Setting memory-based limits on read-chunks loaded at once.
        int num_loadable_reads = (int)((utils::get_free_memory() * 1.0) / 1300 ); // Available memory / estimated memory per read in bytes
        int free_mem = utils::get_free_memory() / ( 1024*1024*1024 );
        INFO( free_mem << "GB available, chunks of approximately " << num_loadable_reads << " reads processed at once");
        int num_reads_section_end = 0;
        int num_reads_section_goal = num_loadable_reads;
        int num_reads_total;
        int cloud_size_filter = SetCloudFilter(gp, lib_10x, num_reads_total);
        INFO("Mean original cloud size is " << cloud_size_filter << " reads");

        while (num_reads_section_end < num_reads_total ) {
            INFO("Loading original read clouds from library");
            num_reads_section_end = LoadReads(connected_reads, mapping_record, read_record, unbarcoded_record, gp, lib_10x, num_reads_section_end, num_reads_section_goal);
            size_t num_original_clouds = mapping_record.size();
            INFO(num_original_clouds << " original read clouds loaded");
#pragma omp parallel for shared(connected_reads, mapping_record, read_record) schedule(dynamic, 1) num_threads(cfg::get().max_threads)
            for (size_t c = 0; c < num_original_clouds; ++c) {
                if (c % 100000 == 0) {
                    INFO(c << ": Deconvolving and reporting " << mapping_record[c].size() << " reads on thread " << omp_get_thread_num());
                }
                if ( mapping_record[c].size() >= cloud_size_filter ) {
                    ClusterReads(connected_reads[c], mapping_record[c], gp);
                }
#pragma omp critical
                {
                    std::vector<std::tuple<std::string, std::string, std::string>> master_record = MakeEnhanced(connected_reads[c], read_record[c], cloud_size_filter, stat_stream);
                    OutputReads(master_record, fastq_stream_forward, fastq_stream_reverse);
                }
            }
            connected_reads.clear();
            mapping_record.clear();
            read_record.clear();
            num_reads_section_goal = num_reads_section_end + num_loadable_reads;
        }

        OutputReads(unbarcoded_record, fastq_stream_forward, fastq_stream_reverse);
        stat_stream << "NA," << 0 << "," << unbarcoded_record.size() << std::endl; // Reporting un-enhanced read cloud information.
        fastq_stream_forward.close();
        fastq_stream_reverse.close();
        stat_stream.close();
        INFO("Read cloud deconvolution finished");
    }

}

#endif