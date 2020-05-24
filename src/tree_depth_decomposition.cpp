#include "tree_depth_decomposition.h"
#include <sstream>

int compute_tree_depth_of_parent_array(const ArrayIDFunc<int>& parent)
{
    int node_count = parent.preimage_count();

    ArrayIDFunc<int> depth_of_node(node_count);
    depth_of_node.fill(-1);

    int tree_depth = 0;

    for (int x = 0; x < node_count; ++x) {
        if (depth_of_node[x] == -1) {
            int depth_of_x = 1;

            {
                int y = x;
                for (;;) {
                    int z = parent[y];
                    if (z == tree_root)
                        break;
                    if (depth_of_node[z] != -1) {
                        depth_of_x += depth_of_node[z];
                        break;
                    }
                    y = z;
                    ++depth_of_x;
                }
            }
            {
                int y = x;
                int depth_of_y = depth_of_x;
                for (;;) {
                    depth_of_node[y] = depth_of_y;
                    int z = parent[y];
                    if (z == tree_root)
                        break;
                    if (depth_of_node[z] != -1)
                        break;
                    y = z;
                    --depth_of_y;
                }
            }
            max_to(tree_depth, depth_of_x);
        }
    }
    return tree_depth;
}

std::string format_parent_array(const ArrayIDFunc<int>& parent, int depth)
{
    int node_count = parent.preimage_count();

    std::ostringstream out;
    out << depth << '\n';
    for (int i = 0; i < node_count; ++i) {
        static_assert(tree_root + 1 == 0, "");
        out << parent[i] + 1 << '\n';
    }
    return out.str();
}

ArrayIDFunc<int> compute_parent_array_from_elimination_order(
    const ArrayIDIDFunc& tail, const ArrayIDIDFunc& head,
    const ArrayIDIDFunc& order)
{
    const int node_count = tail.image_count();
    const int arc_count = tail.preimage_count();

    ArrayIDIDFunc rank = inverse_permutation(order);
    auto order_by_rank = [&](int l, int r) { return rank[l] < rank[r]; };

    ArrayIDFunc<int> parent(node_count);

    std::vector<std::vector<int>> must_be_ancestors(node_count);

    for (int xy = 0; xy < arc_count; ++xy) {
        int x = tail(xy);
        int y = head(xy);
        if (rank[x] < rank[y]) {
            must_be_ancestors[x].push_back(y);
        }
    }

    for (int x = 0; x < node_count; ++x) {
        std::sort(must_be_ancestors[x].begin(), must_be_ancestors[x].end(),
            order_by_rank);
    }

    for (int i = 0; i < node_count; ++i) {
        int x = order(i);
        assert(std::is_sorted(must_be_ancestors[x].begin(),
            must_be_ancestors[x].end(), order_by_rank));
        if (must_be_ancestors[x].empty()) {
            parent[x] = tree_root;
        } else {
            int p = must_be_ancestors[x].front();
            assert(rank[x] < rank[p]);
            parent[x] = p;
            std::vector<int> m(must_be_ancestors[x].size() - 1 + must_be_ancestors[p].size());
            std::merge(must_be_ancestors[x].begin() + 1, must_be_ancestors[x].end(),
                must_be_ancestors[p].begin(), must_be_ancestors[p].end(),
                m.begin(), order_by_rank);
            m.erase(std::unique(m.begin(), m.end()), m.end());
            m.swap(must_be_ancestors[p]);
            if (!must_be_ancestors[p].empty())
                assert(rank[p] < rank[must_be_ancestors[p].front()]);
        }
        must_be_ancestors[x] = std::vector<int>();
    }

    return parent;
}

int compute_tree_depth_of_order(const ArrayIDIDFunc& tail,
    const ArrayIDIDFunc& head,
    const ArrayIDIDFunc& elimination_order)
{
    return compute_tree_depth_of_parent_array(compute_parent_array_from_elimination_order(
        tail, head, elimination_order));
};
