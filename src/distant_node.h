#ifndef DISTANT_NODE_H
#define DISTANT_NODE_H

#include "array_id_func.h"
#include "flow_cutter.h"
#include "id_multi_func.h"
#include "tiny_id_func.h"

template <class Tail, class Head>
flow_cutter::SourceTargetPair compute_distant_node_pair(const Tail& tail, const Head& head)
{
    const int node_count = tail.image_count();
    auto successor = compute_successor_function(tail, head);

    BitIDFunc was_pushed(node_count);
    was_pushed.fill(false);
    was_pushed.set(0, true);

    ArrayIDFunc<int> queue(node_count);
    int queue_begin, queue_end;

    auto reset_queue = [&](int x) {
        was_pushed.fill(false);
        was_pushed.set(x, true);
        queue_end = 1;
        queue_begin = 0;
        queue[0] = x;
    };

    auto push = [&](int x) {
        assert(queue_end != node_count);
        was_pushed.set(x, true);
        queue[queue_end++] = x;
    };

    auto pop = [&] {
        assert(queue_begin != queue_end);
        return queue[queue_begin++];
    };

    auto bfs = [&] {
        int last;
        do {
            int x = pop();
            last = x;
            for (int y : successor(x)) {
                if (!was_pushed(y)) {
                    push(y);
                }
            }
        } while (queue_begin != queue_end);
        return last;
    };

    reset_queue(0);
    reset_queue(bfs());

    flow_cutter::SourceTargetPair pair;
    pair.source = bfs();
    reset_queue(pair.source);
    pair.target = bfs();

    return pair;
}

#endif
