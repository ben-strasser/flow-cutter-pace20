#ifndef SEPARATOR_H
#define SEPARATOR_H

#include "optimize_separator.h"
#include "back_arc.h"
#include "distant_node.h"
#include "flow_cutter.h"
#include "flow_cutter_config.h"
#include "min_max.h"
#include "multi_arc.h"
#include "node_flow_cutter.h"
#include "tiny_id_func.h"

namespace flow_cutter {

class FastComputeSeparator {
public:
    FastComputeSeparator(Config config)
        : config(config)
    {
    }

    template <class Tail, class Head>
    std::vector<int> operator()(const Tail& tail, const Head& head, int max_separator_size) const
    {
        const int node_count = tail.image_count();
        const int arc_count = tail.preimage_count();
        (void)arc_count;

        auto out_arc = invert_sorted_id_id_func(tail);
        auto back_arc = compute_back_arc_permutation(tail, head);

        auto graph = flow_cutter::make_graph(
            make_const_ref_id_id_func(tail),
            make_const_ref_id_id_func(head),
            make_const_ref_id_id_func(back_arc),
            ConstIntIDFunc<1>(arc_count), // capacity
            make_const_ref_id_func(out_arc));

        Config my_config = config;
        my_config.max_cut_size = max_separator_size;

        auto cutter = make_simple_cutter(graph, my_config);
        std::vector<SourceTargetPair> pairs;
        if (config.cutter_count > 0)
            pairs = select_random_source_target_pairs(node_count, config.cutter_count, config.random_seed);
        else
            pairs = { compute_distant_node_pair(tail, head) };

        cutter.init(pairs, config.random_seed);

        std::vector<int> separator;

        for (;;) {
            int cut_size = cutter.get_current_cut().size();

            if (cut_size > max_separator_size)
                break;

            int small_side_size = cutter.get_current_smaller_cut_side_size();

            if (3 * small_side_size > (node_count - cut_size)) {
                separator = cutter.get_current_cut();
                for (int& s : separator)
                    s = head(s);
                std::sort(separator.begin(), separator.end());
                separator.erase(std::unique(separator.begin(), separator.end()), separator.end());
                separator = remove_nodes_from_separator_as_long_as_result_is_balanced(tail, head, std::move(separator));
                break;
            }

            if (!cutter.advance())
                break;
        }

        return separator;
    }

private:
    Config config;
};

class ComputeSeparator {
public:
    explicit ComputeSeparator(Config config)
        : config(config)
    {
    }

    template <class Tail, class Head>
    std::vector<int> operator()(const Tail& tail, const Head& head, int max_separator_size) const
    {

        const int node_count = tail.image_count();
        const int arc_count = tail.preimage_count();
        (void)arc_count;

        auto out_arc = invert_sorted_id_id_func(tail);
        auto back_arc = compute_back_arc_permutation(tail, head);

        auto expanded_graph = expanded_graph::make_graph(
            make_const_ref_id_id_func(tail), make_const_ref_id_id_func(head),
            make_const_ref_id_id_func(back_arc), make_const_ref_id_func(out_arc));

        Config my_config = config;
        my_config.max_cut_size = max_separator_size;

        auto cutter = make_simple_cutter(expanded_graph, my_config);
        std::vector<SourceTargetPair> pairs;
        if (config.cutter_count > 0)
            pairs = select_random_source_target_pairs(node_count, config.cutter_count, config.random_seed);
        else
            pairs = { compute_distant_node_pair(tail, head) };

        cutter.init(expanded_graph::expand_source_target_pair_list(pairs),
            config.random_seed);

        int balance_num, balance_div;
        switch (((unsigned)config.random_seed * (unsigned)node_count) % 3) {
        case 0:
            balance_num = 1;
            balance_div = 3;
        case 1:
            balance_num = 2;
            balance_div = 5;
        case 2:
            balance_num = 1;
            balance_div = 4;
        }

        double best_score = std::numeric_limits<double>::max();
        std::vector<int> separator;

        for (;;) {
            int cut_size = cutter.get_current_cut().size();
            int small_side_size = cutter.get_current_smaller_cut_side_size();

            if (cut_size > max_separator_size)
                break;

            

            if (balance_div * small_side_size > balance_num * (node_count - cut_size)) {

                double score = (double)cut_size / (double)small_side_size;
                if (score < best_score) {
                    separator = remove_nodes_from_separator_as_long_as_result_is_balanced(tail, head, expanded_graph::extract_original_separator(tail, head, cutter).sep);
                    score = (double)separator.size() / (double)small_side_size;
                    best_score = score;
                }

                double potential_best_next_score = (double)(cut_size + 1) / (double)(expanded_graph::expanded_node_count(node_count) / 2);
                if (potential_best_next_score >= best_score)
                    break;
            }

            if (!cutter.advance())
                break;
        }

        return separator;
    }

private:
    Config config;
};
} // namespace flow_cutter

#endif
