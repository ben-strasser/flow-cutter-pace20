#include "greedy_order.h"
#include "array_id_func.h"
#include "heap.h"
#include "id_func.h"
#include "id_multi_func.h"
#include "permutation.h"
#include "tiny_id_func.h"
#include "tree_root.h"
#include <vector>

namespace {

ArrayIDFunc<std::vector<int>> build_dyn_array(const ArrayIDIDFunc& tail,
    const ArrayIDIDFunc& head)
{
    const int node_count = tail.image_count();
    const int arc_count = tail.preimage_count();

    ArrayIDFunc<std::vector<int>> neighbors(node_count);

    for (int i = 0; i < arc_count; ++i)
        neighbors[tail(i)].push_back(head(i));

    for (int i = 0; i < node_count; ++i)
        std::sort(neighbors[i].begin(), neighbors[i].end());

    return neighbors; // NVRO
}

template <class Iter1, class Iter2, class Iter3, class T>
Iter3 set_union_and_remove_element(Iter1 a, Iter1 a_end, Iter2 b, Iter2 b_end,
    Iter3 out, const T& remove_element1,
    const T& remove_element2)
{
    while (a != a_end && b != b_end) {
        if (*a < *b) {
            if (*a != remove_element1 && *a != remove_element2)
                *out++ = *a;
            ++a;
        } else if (*a > *b) {
            if (*b != remove_element1 && *b != remove_element2)
                *out++ = *b;
            ++b;
        } else if (*a == *b) {
            if (*a != remove_element1 && *a != remove_element2)
                *out++ = *a;
            ++b;
            ++a;
        }
    }

    while (a != a_end) {
        if (*a != remove_element1 && *a != remove_element2)
            *out++ = *a;
        ++a;
    }

    while (b != b_end) {
        if (*b != remove_element1 && *b != remove_element2)
            *out++ = *b;
        ++b;
    }

    return out;
}

std::vector<int> contract_node(ArrayIDFunc<std::vector<int>>& graph, int node)
{

    for (int x : graph(node)) {
        std::vector<int> tmp;
        set_union_and_remove_element(graph(node).begin(), graph(node).end(),
            graph(x).begin(), graph(x).end(),
            std::back_inserter(tmp), node, x);
        graph[x].swap(tmp);
    }

    return std::move(graph[node]);
}

}

ArrayIDIDFunc compute_greedy_order(const ArrayIDIDFunc& tail,
    const ArrayIDIDFunc& head)
{
    const int node_count = tail.image_count();

    auto g = build_dyn_array(tail, head);

    min_id_heap<int> q(node_count);

    for (int x = 0; x < node_count; ++x)
        q.push(x, g(x).size());

    ArrayIDIDFunc order(node_count, node_count);
    int next_pos = 0;

    ArrayIDFunc<int> level(node_count);
    level.fill(0);

    while (!q.empty()) {
        auto x = q.pop();

        order[next_pos++] = x;
        if (150 < g(x).size()) {
            while (!q.empty()) {
                order[next_pos++] = q.pop();
            }
            break;
        }

        for (auto y : contract_node(g, x)) {
            if (level[y] < level[x] + 1)
                level[y] = level[x] + 1;
            q.push_or_set_key(y, 8 * level[y] + g(y).size());
        }
    }

    return order; // NVRO
}
