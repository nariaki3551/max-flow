from sys import argv
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

PLOT = False
imgix = 0
relabel_operation = 0
saturating_flow = 0
nonsaturation_flow = 0

def usage():
    usage=f"""python3 {__file__} graph_data_file s t [-p | --plot]
    s: source node
    t: sink   node
    [option]
    -p --plot: draw graph during algorithm
              (You need to create "figs" directory for saveing figures)"""
    print(usage)

def main(graph_data_file, s, t):
    oriG = read_graph_file(graph_data_file)
    if PLOT: draw_graph(oriG, title='original network')

    oriG = highest_label_preflow_push_algo(oriG)
    if PLOT: draw_graph(oriG, title='maximum flow', last=True)

    print()
    print('saturating_flow:', saturating_flow)
    print('nonsaturation_flow:', nonsaturation_flow)
    print('relabel_operation:', relabel_operation)


def highest_label_preflow_push_algo(oriG):

    def flow(u, v, f):
        print('flow {} -> {}: {}'.format(u, v, f))
        oriG[u][v]['preflow']  += f
        resG[u][v]['preflow']  -= f
        resG[v][u]['preflow']  += f
        oriG.node[u]['excess'] -= f
        oriG.node[v]['excess'] += f

    excess = lambda i: oriG.node[i]['excess']

    # generate residul network
    resG = generate_residul(oriG)

    # set distance lable
    oriG, resG = set_ini_d(oriG, resG)

    # set excess
    oriG = set_ini_excess(oriG)

    # set preflow
    oriG, resG = set_ini_preflow(oriG, resG)

    print('preprocess')
    # preprocess distance label
    num_nodes = len(oriG.nodes())
    oriG.node[s]['d'] = num_nodes
    resG.node[s]['d'] = num_nodes

    # preprocess preflow
    active_LIST = {i: set() for i in range(2*len(oriG.nodes()))}
    active_nodes = set()
    for v in oriG[s]:
        f = oriG[s][v]['capacity']
        flow(s, v, f)
        if v != t:
            active_LIST[0].add(v)
            active_nodes.add(v)
    if PLOT: draw_graph(oriG, title='preprocess')


    print('HIGHEST LABEL PREFLOW PUSH ALGORITHM')
    highest_level = 0
    while active_nodes:
        while not active_LIST[highest_level]:
            highest_level = highest_level - 1
        i = active_LIST[highest_level].pop()
        active_nodes.remove(i)
        admissible_arcs = get_admissible_arcs(resG, i)

        while excess(i) and admissible_arcs:
            v = admissible_arcs.pop()
            # push operation
            delta = push_units(oriG, resG, i, v)
            if (i, v) in oriG.edges():
                flow(i, v,  delta)
            else:
                flow(v, i, -delta)
            # add active_LIST newly active node
            if v not in active_nodes and v not in [s, t]:
                active_nodes.add(v)
                active_LIST[oriG.node[v]['d']].add(v)
                highest_level = max(highest_level, oriG.node[v]['d'])
            if PLOT:
                tmp = active_nodes | {i} if excess(i) else active_nodes
                title='flow {} → {}: {} units'.format(i, v, delta)
                draw_graph(oriG, active_nodes=tmp, title=title)

        if excess(i):
            # relabel operation
            oriG.node[i]['d'] = resG.node[i]['d'] = relabel(resG, i)
            print('relabeld {}: {}'.format(i, oriG.node[i]['d']))
            # add active_LIST relabeled node
            active_nodes.add(i)
            active_LIST[oriG.node[i]['d']].add(i)
            highest_level = max(highest_level, oriG.node[i]['d'])
            if PLOT:
                draw_graph(oriG, active_nodes=active_nodes, title='relabeld {}: {}'.format(i, oriG.node[i]['d']))

    return oriG


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
    admissible_arcs = []
    for v in resG[i]:
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
    return oriG, resG


def set_ini_preflow(oriG, resG):
    for u, v in oriG.edges():
        oriG[u][v]['preflow'] = 0
        resG[u][v]['preflow'] = resG[u][v]['capacity']
        resG[v][u]['preflow'] = 0
    return oriG, resG


def set_ini_excess(oriG):
    for n in oriG:
        oriG.node[n]['excess'] = 0
    return oriG


def draw_graph(oriG, active_nodes=set(), title=None, last=False):
    global imgix

    draw_nodes = nx.draw_networkx_nodes
    draw_edges = nx.draw_networkx_edges
    draw_labels = nx.draw_networkx_labels
    draw_edge_labels = nx.draw_networkx_edge_labels

    pos = {i: oriG.node[i]['pos'] for i in oriG}

    fig, ax = plt.subplots(figsize=(12, 8))
    if title is not None:
        ax.set_title(title)

    # draw nodes
    nonactive_nodes = set(oriG.nodes()) - set(active_nodes)
    draw_nodes(oriG, pos=pos, nodelist=nonactive_nodes, node_color='skyblue', ax=ax)
    draw_nodes(oriG, pos=pos, nodelist=[s], node_color='blue', ax=ax)
    draw_nodes(oriG, pos=pos, nodelist=[t], node_color='darkgreen', ax=ax)
    for active_node in active_nodes:
        excess_ratio = calc_excess_ratio(oriG, active_node)
        draw_nodes(oriG, pos=pos, nodelist=[active_node], node_color='red', alpha=excess_ratio, ax=ax)
    nls = {i: oriG.node[i]['d'] if 'd' in oriG.node[i] else i for i in oriG}
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
        els[u, v] = '{}/{}'.format(preflow, capacity)
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
            G.add_node(node, pos=(x, y))
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
        if '-p' in argv or '--plot' in argv:
            PLOT = True
        main(graph_data_file, s, t)
    else:
        usage()
