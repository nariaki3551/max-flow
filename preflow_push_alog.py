from sys import argv
from random import choice
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

PLOT = False
DFS = False
GAP = False
FREEZE = False
imgix = 0
relabel_operation = 0
saturating_flow = 0
nonsaturation_flow = 0

def usage():
    usage=f"""python3 {__file__} graph_data_file s t [-p | --plot] [-d | --dfs] [-g | --gap]
    s: source node
    t: sink   node
    [option]
    -d --dfs : set initial label by DFS
    -g --gap : operate gap-relabeling during algorithm
    -f --freeze: utilize the freezing algorithm
    -p --plot: draw graph during algorithm
              (You need to create "figs" directory for saveing figures)"""
    print(usage)

def main(graph_data_file, s, t):
    oriG = read_graph_file(graph_data_file)
    if PLOT: draw_graph(oriG, title='original network', node_label='node_name')

    print('GENERIC PREFLOW PUSH ALGORITHM')
    oriG = generic_preflow_push_algo(oriG)
    if PLOT: draw_graph(oriG, title='maximum flow', node_label='node_name', last=True)

    print()
    print('saturating_flow:', saturating_flow)
    print('nonsaturation_flow:', nonsaturation_flow)
    print('relabel_operation:', relabel_operation)


def generic_preflow_push_algo(oriG):
    global relabel_operation, saturating_flow, nonsaturation_flow

    def flow(u, v, f):
        print('push {} -> {}: {}'.format(u, v, f))
        oriG[u][v]['preflow']  += f
        resG[u][v]['preflow']  -= f
        resG[v][u]['preflow']  += f
        oriG.node[u]['excess'] -= f
        oriG.node[v]['excess'] += f

    excess = lambda i: oriG.node[i]['excess']
    d      = lambda i: oriG.node[i]['d']

    # generate residul network
    resG = generate_residul(oriG)

    # set distance lable
    set_ini_d(oriG, resG)

    # set excess
    set_ini_excess(oriG)

    # set preflow
    set_ini_preflow(oriG, resG)

    print('preprocess')
    # preprocess distance label
    N = len(oriG.nodes())
    oriG.node[s]['d'] = N
    resG.node[s]['d'] = N

    # preprocess preflow
    active_nodes = set()
    frozen_nodes = set()
    for v in oriG[s]:
        f = oriG[s][v]['capacity']
        flow(s, v, f)
        if v != t:
            active_nodes.add(v)
    if PLOT: draw_graph(oriG, title='preprocess', node_label='d')

    if DFS:
        # DFS labeling
        print('DFS levling')
        dfs_set_label(oriG, resG)
        if PLOT: draw_graph(oriG, title='DFS labeling', node_label='d')

    while active_nodes:
        i = active_nodes.pop() # choice an active node
        admissible_arcs = get_admissible_arcs(resG, i)
        print(f'active node {i} admissible arcs {admissible_arcs}')

        if admissible_arcs:
            v = choice(admissible_arcs)
            # push operation
            delta = push_units(oriG, resG, i, v)
            if (i, v) in oriG.edges():
                flow(i, v,  delta)
            else:
                flow(v, i, -delta)
            if excess(i):
                active_nodes.add(i)
            if v not in active_nodes and v not in [s, t]:
                active_nodes.add(v)
            if PLOT:
                tmp = active_nodes | {i} if excess(i) else active_nodes
                title='push {} → {}: {} units'.format(i, v, delta)
                draw_graph(oriG, active_nodes=tmp,
                    frozen_nodes=frozen_nodes, title=title, node_label='d')
        else:
            # gap-relabeling
            # Whether or not a node of label g dose not exists by labeling i
            g = d(i)
            if GAP and all(d(j) != g for j in oriG.node() if j != i):
                # operate gep-relabeling
                print(f'gap-relabeling operation g = {g}')
                for j in oriG.node():
                    if g < d(j) < N or j == i:
                        print(f'gap-relabeled {j}: {N}')
                        oriG.node[j]['d'] = resG.node[j]['d'] = N
                        if FREEZE and excess(j):
                            print('freeze {}: {}'.format(j, d(j)))
                            frozen_nodes.add(j)
                            if j in active_nodes:
                                active_nodes.remove(j)
                        elif excess(j):
                            active_nodes.add(j)
                        relabel_operation += 1
                        if PLOT:
                            title = f'gap-relabeled {j}: {N}'
                            draw_graph(oriG, active_nodes=active_nodes,
                                frozen_nodes=frozen_nodes, title=title, node_label='d')
            else:
                # normal relabel operation
                oriG.node[i]['d'] = resG.node[i]['d'] = relabel(resG, i)
                if FREEZE and d(i) >= N:
                    frozen_nodes.add(i)
                    print('freeze {}: {}'.format(i, oriG.node[i]['d']))
                    if PLOT:
                        title = 'freeze {}: {}'.format(i, oriG.node[i]['d'])
                        draw_graph(oriG, active_nodes=active_nodes,
                            frozen_nodes=frozen_nodes, title=title, node_label='d')
                else:
                    active_nodes.add(i)
                    print('relabeld {}: {}'.format(i, oriG.node[i]['d']))
                    if PLOT:
                        title = 'relabeld {}: {}'.format(i, oriG.node[i]['d'])
                        draw_graph(oriG, active_nodes=active_nodes,
                            frozen_nodes=frozen_nodes, title=title, node_label='d')

    if FREEZE and frozen_nodes:
        print('FREEZE operation')
        print('frozen_nodes', frozen_nodes)
        from collections import OrderedDict

        # label the node to the distance from s
        dfs_lable_from_s(oriG, resG)
        sdepth = lambda node: oriG.node[node]['sdepth']

        def iter_can_back(u, visited):
            children = [v for v in set(resG[u]) - set(visited)
                if u in oriG[v] and oriG[v][u]['preflow'] > 0]
            for i in sorted(children, key=lambda v: (v in {s} | frozen_nodes, oriG[v][u]['preflow'], -sdepth(v)), reverse=True):
                yield i


        if PLOT:
            title='label on depth from s'
            draw_graph(oriG, active_nodes=active_nodes,
                frozen_nodes=frozen_nodes, title=title, node_label='sdepth')
        
        # reverse flow from the frozen point faster than s
        for frozen_node in sorted(frozen_nodes, key=sdepth, reverse=True):
            frozen_nodes.remove(frozen_node)
            excess_value = excess(frozen_node)
            print(f'frozen_node {frozen_node} excess_value {excess_value}')
            # find the path from fronzen node to s or other fronzen node by DFS and flow
            # continue this process until there is no excess
            
            visited = OrderedDict.fromkeys([frozen_node])
            stack = [iter_can_back(frozen_node, visited)]

            while excess_value:
                children = stack[-1]
                child = next(children, None)

                if child is None:
                    stack.pop()
                    visited.popitem()
                elif oriG[child][list(visited)[-1]]['preflow'] == 0:
                    continue
                elif len(visited) < N-1:
                    if child in ({s} | frozen_nodes):
                        path = list(visited) + [child]
                        print('found path', path)
                        min_cap = min(oriG[v][u]['preflow'] for u, v in zip(path, path[1:]))
                        delta = min(excess_value, min_cap) # max flow value 
                        excess_value -= delta
                        for u, v in zip(path, path[1:]):
                            if (u, v) in oriG.edges():
                                flow(u, v,  delta)
                            else:
                                flow(v, u, -delta)
                            saturating_flow += 1
                            if PLOT:
                                tmp = frozen_nodes | {frozen_node} if excess_value else frozen_nodes
                                title='push {} → {}: {} units'.format(u, v, delta)
                                draw_graph(oriG, frozen_nodes=tmp, title=title, node_label='sdepth')

                        if excess_value > 0:
                            # back to the edge whose weight is delta
                            pivot_ix = 0
                            for u, v in zip(path, path[1:]):
                                if oriG[v][u]['preflow'] == 0:
                                    break
                                pivot_ix += 1
                            # len(path) == pivot_ixになるまでpopする
                            for _ in range(len(path)-pivot_ix-2):
                                stack.pop()
                                visited.popitem()

                    elif child not in visited:
                        visited[child] = None
                        stack.append(iter_can_back(child, visited))
    
    return oriG


def dfs_lable_from_s(oriG, resG):
    for u in oriG:
        if u == s:
            oriG.node[s]['sdepth'] = 0
        else:
            oriG.node[u]['sdepth'] = float('inf')
    labeld_nodes = {s}
    depth = 0
    stack = [s]
    next_stack = []
    while stack:
        depth += 1
        for u in stack:
            for v in oriG[u]:
                if oriG[u][v]['preflow'] > 0 and v not in labeld_nodes:
                    oriG.node[v]['sdepth'] = resG.node[v]['sdepth'] = depth
                    labeld_nodes.add(v)
                    next_stack.append(v)
        stack = next_stack
        next_stack = []


def dfs_set_label(oriG, resG):
    labeld_nodes = {s, t}
    depth = 0
    stack = [t]
    next_stack = []
    while stack:
        depth += 1
        for u in stack:
            for v in resG[u]:
                if v not in labeld_nodes:
                    oriG.node[v]['d'] = resG.node[v]['d'] = depth
                    labeld_nodes.add(v)
                    next_stack.append(v)
        stack = next_stack
        next_stack = []


def push_units(oriG, resG, i, v):
    global saturating_flow, nonsaturation_flow
    i_excess = oriG.node[i]['excess']
    res_flow = resG[i][v]['preflow']
    if res_flow < i_excess:
        saturating_flow += 1
    else:
        nonsaturation_flow += 1
    return min(i_excess, res_flow)


def relabel(resG, i):
    global relabel_operation
    relabel_operation += 1
    adj_i = [v for v in resG[i] if resG[i][v]['preflow'] > 0]
    return min([resG.node[v]['d'] for v in adj_i]) + 1


def get_admissible_arcs(resG, i):
    # device
    if s in resG[i]:
        if resG[i][s]['preflow'] > 0:
            if resG.node[i]['d'] == resG.node[s]['d'] + 1:
                return [s]
    admissible_arcs = []
    for v in set(resG[i]) - {s}:
        if resG[i][v]['preflow'] > 0:
            if resG.node[i]['d'] == resG.node[v]['d'] + 1:
                admissible_arcs.append(v)
    return admissible_arcs


def generate_residul(oriG):
    resG = nx.DiGraph()
    for n in oriG:
        pos = oriG.node[n]['pos']
        resG.add_node(n, pos=pos)
    for u, v in oriG.edges():
        capacity = oriG[u][v]['capacity']
        if (u, v) not in resG:
            resG.add_edge(u, v, capacity=capacity)
            resG.add_edge(v, u, capacity=capacity)
        else:
            resG[u][v]['capacity'] += capacity
            resG[v][u]['capacity'] += capacity
    return resG


def set_ini_d(oriG, resG):
    for n in oriG:
        oriG.node[n]['d'] = 0
        resG.node[n]['d'] = 0


def set_ini_preflow(oriG, resG):
    for u, v in oriG.edges():
        oriG[u][v]['preflow'] = 0
        resG[u][v]['preflow'] = resG[u][v]['capacity']
        resG[v][u]['preflow'] = 0


def set_ini_excess(oriG):
    for n in oriG:
        oriG.node[n]['excess'] = 0


def draw_graph(oriG, active_nodes=set(), frozen_nodes=set(), title=None, node_label='node_name', last=False):
    global imgix

    draw_nodes = nx.draw_networkx_nodes
    draw_edges = nx.draw_networkx_edges
    draw_labels = nx.draw_networkx_labels
    draw_edge_labels = nx.draw_networkx_edge_labels

    pos = {i: oriG.node[i]['pos'] for i in oriG}

    fig, ax = plt.subplots(figsize=(12, 8))
    # fig, ax = plt.subplots(figsize=(6, 4))
    if title is not None:
        ax.set_title(title)

    # draw nodes
    non_active_nodes = set(oriG) - active_nodes - frozen_nodes
    draw_nodes(oriG, pos=pos, nodelist=non_active_nodes, node_color='skyblue', ax=ax)
    draw_nodes(oriG, pos=pos, nodelist=frozen_nodes, node_color='grey', ax=ax)
    draw_nodes(oriG, pos=pos, nodelist=[s], node_color='blue', ax=ax)
    draw_nodes(oriG, pos=pos, nodelist=[t], node_color='darkgreen', ax=ax)
    for active_node in active_nodes:
        excess_ratio = calc_excess_ratio(oriG, active_node)
        draw_nodes(oriG, pos=pos, nodelist=[active_node], node_color='red', alpha=excess_ratio, ax=ax)
    nls = {i: oriG.node[i][node_label] for i in oriG}
    draw_labels(oriG, pos=pos, labels=nls, ax=ax)

    # draw edges
    draw_edges(oriG, pos=pos, edgelist=oriG.edges(), alpha=0.2)
    for u, v in oriG.edges():
        if 'preflow' in oriG[u][v]:
            width = oriG[u][v]['preflow']
        else:
            width = 1
        draw_edges(oriG, pos=pos, edgelist=[(u, v)], width=width, ax=ax)
    els = dict()
    for u, v in oriG.edges():
        if 'preflow' in oriG[u][v]:
            preflow  = oriG[u][v]['preflow']
        else:
            preflow = 0
        capacity = oriG[u][v]['capacity']
        els[u, v] = f'{int(preflow)}/{int(capacity)}'
    draw_edge_labels(oriG, pos=pos, edge_labels=els, ax=ax)

    # erase locator
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())

    if last:
        plt.savefig('figs/Last.png')
    else:
        plt.savefig('figs/Image{0:04d}.png'.format(imgix)); imgix = imgix + 1

    plt.close()


def calc_excess_ratio(oriG, i):
    to_i_nodes = [u for u in oriG if i in oriG[u]]
    max_excess = sum([oriG[u][i]['capacity'] for u in to_i_nodes])
    return oriG.node[i]['excess'] / max_excess


def read_graph_file(graph_data_file):
    G = nx.DiGraph()

    f = open(graph_data_file, 'r')
    for line in f.readlines():
        line = line.strip()
        if not line:
            continue
        if line == 'nodes':
            sigh = 'node'
            continue
        if line == 'edges':
            sigh = 'edge'
            continue
        if sigh == 'node':
            node, x, y = line.split()
            node, x, y = int(node), float(x), float(y)
            G.add_node(node, pos=(x, y), node_name=node)
        if sigh == 'edge':
            u, v, capacity = line.split()
            u, v, capacity = int(u), int(v), float(capacity)
            G.add_edge(u, v, capacity=capacity)
    f.close()

    return G

if __name__ == '__main__':
    if len(argv) >= 4:
        graph_data_file = argv[1]
        s = int(argv[2])
        t = int(argv[3])
        if '-d' in argv or '--dfs' in argv:
            DFS = True
            print('option DFS from OFF to ON')
        if '-g' in argv or '--gap' in argv:
            GAP = True
            print('option GAP from OFF to ON')
        if '-f' in argv or '--freeze' in argv:
            FREEZE = True
            print('option FREEZE from OFF to ON')
        if '-p' in argv or '--plot' in argv:
            PLOT = True
            print('option PLOT from OFF to ON')
        main(graph_data_file, s, t)
    else:
        usage()
