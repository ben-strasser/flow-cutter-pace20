#ifndef OPTIMIZE_SEPARATOR_H
#define OPTIMIZE_SEPARATOR_H

#include "union_find.h"
#include "id_multi_func.h"
#include <vector>

template<class Tail, class Head>
std::vector<int>remove_nodes_from_separator_as_long_as_result_is_balanced(const Tail&tail, const Head&head, std::vector<int> separator){
        const int node_count = tail.image_count();        
        const int arc_count = tail.preimage_count();

        BitIDFunc in_separator(node_count);
        in_separator.fill(false);
        for(int x:separator)
                in_separator.set(x, true);

        UnionFind uf(node_count);
        for(int xy=0; xy<arc_count; ++xy){
                int x=tail(xy), y=head(xy);
                if(!in_separator(x) && !in_separator(y))
                        uf.unite(x, y);
        }

        BitIDFunc was_representative_counted(node_count);
        was_representative_counted.fill(false);

        ArrayIDIDMultiFunc successor = compute_successor_function(tail, head);

        auto compute_component_size_if_node_removed_from_separator = [&](int x){
                int comp_size = 1;
                for(int y:successor(x)){
                        if(!in_separator(y)){
                                int r = uf(y);
                                if(!was_representative_counted(r)){
                                        comp_size += uf.component_size(r);
                                        was_representative_counted.set(r, true);
                                }
                        }
                }
                for(int y:successor(x)){
                        if(!in_separator(y)){
                                was_representative_counted.set(uf(y), false);
                        }
                }
                return comp_size;
        };

        auto remove_node_from_separator = [&](int x){
                for(int y:successor(x)){
                        uf.unite(y, x);
                }
                in_separator.set(x, false);
        };

        int out = 0;
        const int separator_size = separator.size();
        for(int in=0; in<separator_size; ++in){

                if(3*compute_component_size_if_node_removed_from_separator(separator[in]) <= 2*node_count){
                        remove_node_from_separator(separator[in]);
                }else{
                        separator[out] = separator[in];
                        ++out;
                }
        }
        separator.erase(separator.begin()+out, separator.end());

        return separator;
}

#endif
