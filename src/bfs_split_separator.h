#ifndef BFS_SPLIT_SEPARATOR_H
#define BFS_SPLIT_SEPARATOR_H

#include "array_id_func.h"
#include "greedy_order.h"
#include "heap.h"
#include "id_func.h"
#include "id_multi_func.h"
#include "permutation.h"
#include "tiny_id_func.h"
#include "tree_root.h"
#include <vector>

class ActiveNodeSet{
public:
    ActiveNodeSet(const ArrayIDIDMultiFunc& successor, const BitIDFunc& side):
        successor(successor), side(side),
        active_node_list(successor.preimage_count()),
        active_node_count(0),
        is_active_node(successor.preimage_count()){
        is_active_node.fill(false);

        for(int x=0; x<successor.preimage_count(); ++x){
            for(int y:successor(x)){
                if(side(x) != side(y)){
                    active_node_list[active_node_count++]=x;
                    is_active_node.set(x, true);
                    break;
                }
            }
        }
    }

    void activate(int x){
        if(!is_active_node(x)){
            active_node_list[active_node_count++] = x;
            is_active_node.set(x, true);
        }
    }

    void notify_node_side_change(int x){
        activate(x);
        for(int y:successor(x))
            activate(y);
    }

    template<class Callback>
    void forall_active_nodes(const Callback&callback){
        int i = 0;
        while(i < active_node_count){
            int x = active_node_list[i];
            bool is_still_active = false;

            for(int y:successor(x)){
                if(side(x) != side(y)){
                    is_still_active = true;
                    break;
                }
            }
            if(!is_still_active){
                active_node_list[i] = active_node_list[--active_node_count];
                is_active_node.set(x, false);
            }else{
                callback(x);
                ++i;
            }
        }
    }

    template<class RandGen>
    void shuffle_active_nodes(RandGen&rand_gen){
        std::shuffle(active_node_list.begin(), active_node_list.begin()+active_node_count, rand_gen);
    }

    bool is_active(int x)const{
        return is_active_node(x);
    }

private:
    const ArrayIDIDMultiFunc& successor;
    const BitIDFunc& side;

    ArrayIDFunc<int>active_node_list;
    int active_node_count;
    BitIDFunc is_active_node;

};

template <class RandGen, class ShouldMoveSide>
int move_nodes(ActiveNodeSet&active_node_set, const ArrayIDIDMultiFunc& successor, BitIDFunc& side, int* side_size, RandGen& rand_gen, const ShouldMoveSide& should_move_side)
{
    int move_count = 0;

    active_node_set.shuffle_active_nodes(rand_gen);
    active_node_set.forall_active_nodes(
        [&](int x){
            unsigned my_side = side(x);
            unsigned other_side = 1 ^ my_side;

            int my_side_size = 0;
            int other_side_size = 0;
            for (int y : successor(x)) {
                if (my_side == side(y))
                    ++my_side_size;
                else
                    ++other_side_size;
            }

            if (should_move_side(my_side_size, other_side_size, side_size[my_side], side_size[other_side], 1, other_side)) {
                --side_size[my_side];
                side.set(x, other_side);
                active_node_set.notify_node_side_change(x);
                ++side_size[other_side];
                ++move_count;
            }
        }
    );

    return move_count;
}

template <class RandGen, class ShouldMoveSide>
int move_edges(ActiveNodeSet&active_node_set, const ArrayIDIDMultiFunc& successor, BitIDFunc& side, int* side_size, RandGen& rand_gen, const ShouldMoveSide& should_move_side)
{
    const int node_count = successor.preimage_count();

    auto compute_cut_size = [&] {
        int s = 0;
        for (int x = 0; x < node_count; ++x)
            for (int y : successor(x))
                if (side(x) != side(y))
                    ++s;

        return s / 2;
    };
    (void)compute_cut_size;

    int move_count = 0;

    BitIDFunc is_neighbor_of_x(node_count);
    is_neighbor_of_x.fill(false);

    active_node_set.shuffle_active_nodes(rand_gen);
    active_node_set.forall_active_nodes(
        [&](int x){

            for (int y : successor(x))
                is_neighbor_of_x.set(y, true);

            for (int y : successor(x)) {
                unsigned my_side = side(x);
                unsigned other_side = 1 ^ my_side;
                if ((x < y || !active_node_set.is_active(y)) && side(y) == my_side) {

                    int my_side_size = 0;
                    int other_side_size = 0;
                    for (int z : successor(x)) {

                        if (my_side == side(z))
                            ++my_side_size;
                        else
                            ++other_side_size;
                    }
                    for (int z : successor(y)) {

                        if (!is_neighbor_of_x(z)) {
                            if (my_side == side(z))
                                ++my_side_size;
                            else
                                ++other_side_size;
                        }
                    }

                    my_side_size -= 2;

                    if (should_move_side(my_side_size, other_side_size, side_size[my_side], side_size[other_side], 2, other_side)) {
                        assert(side_size[0] + side_size[1] == node_count);
                        side_size[my_side] -= 2;
                        side.set(x, other_side);
                        side.set(y, other_side);
                        active_node_set.notify_node_side_change(x);
                        active_node_set.notify_node_side_change(y);
                        side_size[other_side] += 2;
                        assert(side_size[0] + side_size[1] == node_count);
                        ++move_count;
                    }
                }
            }

            for (int y : successor(x))
                is_neighbor_of_x.set(y, false);
        }
    );
    return move_count;
}

template <class RandGen, class ShouldMoveSide>
int move_nodes_and_edges(ActiveNodeSet&active_node_set, const ArrayIDIDMultiFunc& successor, BitIDFunc& side, int* side_size, RandGen& rand_gen, const ShouldMoveSide& should_move_side)
{
    return move_nodes(active_node_set, successor, side, side_size, rand_gen, should_move_side) + move_edges(active_node_set, successor, side, side_size, rand_gen, should_move_side);
}

template <class RandGen>
void optimize_cut(const ArrayIDIDMultiFunc& successor, BitIDFunc& side, RandGen& rand_gen)
{
    const int node_count = successor.preimage_count();

    ActiveNodeSet active_node_set(successor, side);

    auto compute_cut_size = [&] {
        int s = 0;
        for (int x = 0; x < node_count; ++x)
            for (int y : successor(x))
                if (side(x) != side(y))
                    ++s;

        return s / 2;
    };
    (void)compute_cut_size;

    int side_size[2] = { 0, 0 };

    for (int x = 0; x < node_count; ++x)
        ++side_size[side(x)];

    auto move_if_cut_decrease = [&](
                                    int my_side_neighbor_count, int other_side_neighbor_count,
                                    int my_side_node_count, int other_side_node_count,
                                    int obj_size, int other_side) {
        return other_side_neighbor_count > my_side_neighbor_count && 3 * other_side_node_count + obj_size < 2 * node_count;
    };

    auto move_if_cut_decrease_or_balance_improvement = [&](
                                                           int my_side_neighbor_count, int other_side_neighbor_count,
                                                           int my_side_node_count, int other_side_node_count,
                                                           int obj_size, int other_side) {
        return move_if_cut_decrease(my_side_neighbor_count, other_side_neighbor_count, my_side_node_count, other_side_node_count, obj_size, other_side)
            || (my_side_neighbor_count == other_side_neighbor_count && other_side_node_count + obj_size < my_side_node_count);
    };

    auto move_to_0_if_same_cut_and_balance_fulfilled = [&](
                                                           int my_side_neighbor_count, int other_side_neighbor_count,
                                                           int my_side_node_count, int other_side_node_count,
                                                           int obj_size, int other_side) {
        return my_side_neighbor_count == other_side_neighbor_count && other_side == 0 && 3 * other_side_node_count + obj_size < 2 * node_count;
    };

    const int round_count = 8;

    auto decrease_cut_size = [&] {
        for (int round = 0; round < round_count; ++round)
            if (move_nodes_and_edges(active_node_set, successor, side, side_size, rand_gen, move_if_cut_decrease) < 10)
                break;
    };

    auto balance_cut = [&] {
        for (int round = 0; round < round_count; ++round)
            if (move_nodes_and_edges(active_node_set, successor, side, side_size, rand_gen, move_if_cut_decrease_or_balance_improvement) < 10)
                break;
    };

    auto shrink_side1 = [&] {
        for (int round = 0; round < round_count; ++round)
            if (move_nodes_and_edges(active_node_set, successor, side, side_size, rand_gen, move_to_0_if_same_cut_and_balance_fulfilled) < 10)
                break;
    };

    decrease_cut_size();
    balance_cut();
    for (int i = 0; i < 20; ++i) {
        shrink_side1();
        decrease_cut_size();
        balance_cut();
    }
}

std::vector<int> convert_cut_to_balanced_separator_or_no_separator(const ArrayIDIDFunc& tail, const ArrayIDIDFunc& head, const BitIDFunc& side)
{
    const int arc_count = tail.preimage_count();
    const int node_count = tail.image_count();

    std::vector<int> separator;

    BitIDFunc in_separator(node_count);
    in_separator.fill(false);

    int side_size[2] = { 0, 0 };
    for (int x = 0; x < node_count; ++x)
        ++side_size[side(x)];

    for (int xy = 0; xy < arc_count; ++xy) {
        int x = tail(xy);
        int y = head(xy);
        if (side(x) == 0 && side(y) == 1) {
            if (!in_separator(x) && !in_separator(y)) {
                if (side_size[0] < side_size[1]) {
                    separator.push_back(y);
                    in_separator.set(y, true);
                    --side_size[1];
                } else {
                    separator.push_back(x);
                    in_separator.set(x, true);
                    --side_size[0];
                }
            }
        }
    }

    if (3 * std::min(side_size[0], side_size[1]) < node_count - (int)separator.size())
        separator.clear();

    return separator;
}

template <class RandGen>
std::vector<int> compute_separator_by_running_bfs(const ArrayIDIDFunc& tail, const ArrayIDIDFunc& head, int max_size, RandGen& rand_gen)
{
    const int node_count = tail.image_count();

    std::vector<int> queue(node_count);
    int queue_begin, queue_end;

    BitIDFunc was_pushed(node_count);

    auto empty = [&] { return queue_begin == queue_end; };
    auto push = [&](int x) {queue[queue_end++]=x; was_pushed.set(x, true); };
    auto pop = [&] { return queue[queue_begin++]; };
    auto reset = [&] {queue_begin = 0; queue_end = 0; was_pushed.fill(false); };

    BitIDFunc side(node_count);

    ArrayIDIDMultiFunc successor = compute_successor_function(tail, head);

    std::vector<int> best_separator;

    const int round_count = 15;
    for (int round = 0; round < round_count; ++round) {
        reset();
        int s, t;
        do {
            s = rand_gen() % node_count;
            t = rand_gen() % node_count;
        } while (s == t);

        side.set(s, 0);
        side.set(t, 1);

        push(s);
        push(t);

        while (!empty()) {
            int x = pop();
            for (int y : successor(x)) {
                if (!was_pushed(y)) {
                    side.set(y, side(x));
                    push(y);
                }
            }
        }

        optimize_cut(successor, side, rand_gen);

        std::vector<int> separator = convert_cut_to_balanced_separator_or_no_separator(tail, head, side);
        if (best_separator.empty() || separator.size() < best_separator.size()) {
            best_separator = std::move(separator);
        }
    }

    if ((int)best_separator.size() > max_size)
        best_separator.clear();

    return best_separator;
}

#endif
