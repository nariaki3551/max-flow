from random import randint
import networkx as nx
import matplotlib.pyplot as plt

N = 8
M = 40


G = nx.DiGraph()
G.add_nodes_from(range(N))
while len(G.edges()) < M:
    u = randint(0, N-1)
    v = randint(0, N-1)
    if u == v or (u, v) in G.edges():
        continue
    G.add_edge(u, v, weight=randint(1, 10))
print('nodes')
pos = nx.spring_layout(G)
for i in range(N):
    print(i, *pos[i])
print('edges')
for u, v in G.edges():
    cap = G[u][v]['weight']
    print(u, v, cap)
# nx.draw(G)
# plt.show()
