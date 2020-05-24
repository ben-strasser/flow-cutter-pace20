#include "greedy_order.h"
#include "list_graph.h"
#include "node_flow_cutter.h"
#include "separator.h"
#include "tree_depth_decomposition.h"

#include "bfs_split_separator.h"

#include <limits>
#include <signal.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#ifdef PARALLELIZE
#include <atomic>
#include <omp.h>
#endif

#include <sys/time.h>
#include <unistd.h>
using namespace std;

bool print_status = false;
bool print_verbose_status = false;

ArrayIDIDFunc tail, head;
const char* volatile best_decomposition = 0;
int best_tree_depth = numeric_limits<int>::max();

void ignore_return_value(int) {}

unsigned long long program_start_milli_time;

unsigned long long get_milli_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (unsigned long long)(tv.tv_sec) * 1000 + (unsigned long long)(tv.tv_usec) / 1000;
}

void test_new_elimination_order(std::string name, const ArrayIDIDFunc& order)
{
    if(order.preimage_count() == 0){
        if (print_verbose_status) {
            string msg = name+" was aborted after " + to_string(get_milli_time() - program_start_milli_time) + "\n";
            ignore_return_value(write(STDERR_FILENO, msg.data(), msg.length()));
        }
        return;
    }
    ArrayIDFunc<int> parent = compute_parent_array_from_elimination_order(tail, head, order);
    int depth = compute_tree_depth_of_parent_array(parent);
    if (depth < best_tree_depth) {
        char* new_decomposition;
        {
            std::string s = format_parent_array(parent, depth);
            new_decomposition = new char[s.length() + 1];
            memcpy(new_decomposition, s.c_str(), s.length() + 1);
        }
#ifdef PARALLELIZE
#pragma omp critical
#endif
        {
            if (depth < best_tree_depth) {
                delete[] best_decomposition;
                best_decomposition = new_decomposition;
                new_decomposition = nullptr;
                best_tree_depth = depth;
                if (print_status) {
                    string msg = "depth " + to_string(best_tree_depth) + " found after " + to_string(get_milli_time() - program_start_milli_time) + " ms by " + name + "\n";
                    ignore_return_value(write(STDERR_FILENO, msg.data(), msg.length()));
                }
            }
        }
        delete[] new_decomposition;
    }else{
        if (print_verbose_status) {
            string msg = "not better depth " + to_string(depth) + " found after " + to_string(get_milli_time() - program_start_milli_time) + " ms by " + name + "\n";
            ignore_return_value(write(STDERR_FILENO, msg.data(), msg.length()));
        }
    }
}

char no_decomposition_message[] = "programm was aborted before any decomposition was computed\n";

#ifdef PARALLELIZE
volatile atomic_flag only_one_thread_in_signal_handler = ATOMIC_FLAG_INIT;
#endif

void signal_handler(int)
{
#ifdef PARALLELIZE
    while (only_one_thread_in_signal_handler.test_and_set()) {
    }
#endif

    const char* x = best_decomposition;
    if (x != 0)
        ignore_return_value(write(STDOUT_FILENO, x, strlen(x)));
    else if (print_status) {
        ignore_return_value(write(STDOUT_FILENO, no_decomposition_message,
            sizeof(no_decomposition_message)));
    }

    _Exit(EXIT_SUCCESS);
}

int main(int argc, char* argv[])
{
    signal(SIGTERM, signal_handler);
    signal(SIGINT, signal_handler);
    signal(SIGSEGV, signal_handler);

    int random_seed = 0;

    try {
        {
            string input_file_name = "-";
            for (int i = 1; i < argc; ++i) {
                if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")){
                    char msg[] = "Computes a tree depth decomposition given a graph. The graph is read by default from stdin in the PACE 2020 graph format. The output is written to stdout. Status messages can be written to stderr, if requested. By default, nothing is written to stderr. The program supports the following options:\n"
                    "  -h,--help Print this message and do nothing else\n"
                    "  --status  Print a message to stderr each time a\n"
                    "            better tree depth decomposition is found\n"
                    "  --verbose Print a message to stderr each time a tree\n"
                    "            depth decomposition was computed independent\n"
                    "            of whether it is better than the best one \n"
                    "            found so far\n"
                    "  -i <file> Instead of reading the graph from stdin,\n"
                    "            read it from <file>\n"
                    "  -s <seed> Use <seed> as seed for the random number\n"
                    "            generator. This must be an integer. The\n"
                    "            default seed is 0.\n";
                    ignore_return_value(write(STDERR_FILENO, msg, sizeof(msg)-1));
                    return 1;
                } else if (!strcmp(argv[i], "--verbose")) {
                    print_verbose_status = true;
                    print_status = true;
                } else if (!strcmp(argv[i], "--status")) {
                    print_status = true;
                } else if (!strcmp(argv[i], "-i") && i != argc - 1) {
                    ++i;
                    input_file_name = argv[i];
                } else if (!strcmp(argv[i], "-s") && i != argc - 1) {
                    ++i;
                    random_seed = atoi(argv[i]);
                }
            }

            auto g = uncached_load_pace_graph(input_file_name);
            tail = std::move(g.tail);
            head = std::move(g.head);
        }

        if (print_status) {
            string msg = "node_count = " + to_string(tail.image_count()) + " arc_count = " + to_string(tail.preimage_count()) + "\n";
            ignore_return_value(write(STDERR_FILENO, msg.data(), msg.length()));
        }

        if (print_status)
            program_start_milli_time = get_milli_time();

        const int node_count = tail.image_count();
                        
        #ifdef PARALLELIZE
        #pragma omp parallel
        #endif
        {
            try {
                std::minstd_rand rand_gen;
                rand_gen.seed(random_seed
                    #ifdef PARALLELIZE
                    + omp_get_thread_num()
                    #endif
                );

                #ifdef PARALLELIZE
                #pragma omp sections nowait
                #endif
                {
                    #ifdef PARALLELIZE
                    #pragma omp section
                    #endif
                    {
                        test_new_elimination_order("greedy order", compute_greedy_order(tail, head));

                        if(20*best_tree_depth > node_count)
                            test_new_elimination_order("refined bfs split in nested dissection", compute_tree_depth_order(
                                tail, head, 
                                [&](const ArrayIDIDFunc& tail, const ArrayIDIDFunc& head, int max_size) {
                                    return compute_separator_by_running_bfs(tail, head, max_size, rand_gen);
                                },
                                best_tree_depth - 1)
                            );
                    }

                    #ifdef PARALLELIZE
                    #pragma omp section
                    #endif
                    {
                        flow_cutter::Config config;
                        config.random_seed = rand_gen();
                        config.cutter_count = 0;
                        config.pierce_rating = flow_cutter::Config::PierceRating::max_target_minus_source_hop_dist;
                        config.max_cut_size = node_count;
                        test_new_elimination_order(
                            "edge flowcutter cutter_count=1 distant-source-target pierce_rating=max_target_minus_source_hop_dist random_seed=" + config.get("random_seed"),
                            compute_tree_depth_order(
                                tail, head,
                                flow_cutter::FastComputeSeparator(config),
                                best_tree_depth - 1));
                    }
                    #ifdef PARALLELIZE
                    #pragma omp section
                    #endif
                    { 
                        flow_cutter::Config config;
                        config.random_seed = rand_gen();
                        config.cutter_count = 0;
                        config.pierce_rating = flow_cutter::Config::PierceRating::max_target_minus_source_hop_dist;
                        config.max_cut_size = node_count;
                        test_new_elimination_order(
                            "flowcutter cutter_count=1 distant-source-target pierce_rating=max_target_minus_source_hop_dist random_seed=" + config.get("random_seed"),
                            compute_tree_depth_order(
                                tail, head,
                                flow_cutter::ComputeSeparator(config),
                                best_tree_depth - 1)
                        );
                    }
                }

                {
                    flow_cutter::Config config;
                    config.max_cut_size = node_count;

                    for (int i = 0;; ++i) {
                        
                        if (i % 2 == 0) {
                            config.pierce_rating = flow_cutter::Config::PierceRating::
                                max_target_minus_source_hop_dist;
                        } else {
                            config.pierce_rating = flow_cutter::Config::PierceRating::random;
                        }
                        
                        config.cutter_count = 1;
                        config.random_seed = rand_gen();
                        test_new_elimination_order(
                            "flowcutter"
                            " cutter_count="
                                + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                            compute_tree_depth_order(
                                tail, head,
                                flow_cutter::ComputeSeparator(config),
                                best_tree_depth - 1));


                        config.cutter_count = 2;
                        config.random_seed = rand_gen();
                        test_new_elimination_order(
                            "flowcutter"
                            " cutter_count="
                                + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                            compute_tree_depth_order(
                                tail, head,
                                flow_cutter::ComputeSeparator(config),
                                best_tree_depth - 1));

                        if((i%3)>0){
                                config.cutter_count = 3;
                                config.random_seed = rand_gen();
                                test_new_elimination_order(
                                    "flowcutter"
                                    " cutter_count="
                                        + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                                    compute_tree_depth_order(
                                        tail, head,
                                        flow_cutter::ComputeSeparator(config),
                                        best_tree_depth - 1));
                        }

                        if((i%20)>15){
                                config.cutter_count = 20;
                                config.random_seed = rand_gen();
                                test_new_elimination_order(
                                    "flowcutter"
                                    " cutter_count="
                                        + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                                    compute_tree_depth_order(
                                        tail, head,
                                        flow_cutter::ComputeSeparator(config),
                                        best_tree_depth - 1));
                        }

                        if((i%50)>30){
                                config.cutter_count = 40;
                                config.random_seed = rand_gen();
                                test_new_elimination_order(
                                    "flowcutter"
                                    " cutter_count="
                                        + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                                    compute_tree_depth_order(
                                        tail, head,
                                        flow_cutter::ComputeSeparator(config),
                                        best_tree_depth - 1));
                        }

                        if((i%100)>98){
                                config.cutter_count = 80;
                                config.random_seed = rand_gen();
                                test_new_elimination_order(
                                    "flowcutter"
                                    " cutter_count="
                                        + config.get("cutter_count") + " pierce_rating=" + config.get("pierce_rating") + " random_seed=" + config.get("random_seed"),
                                    compute_tree_depth_order(
                                        tail, head,
                                        flow_cutter::ComputeSeparator(config),
                                        best_tree_depth - 1));
                        }
                    }
                }
            }
            catch (...)
            {
                signal_handler(0);
            }
        }
    }
    catch (...)
    {
        signal_handler(0);
    }
}
