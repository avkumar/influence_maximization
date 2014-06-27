import networkx as nx
import matplotlib.pyplot as plt
red=nx.random_lobster(100,0.9,0.9)
G=nx.sedgewick_maze_graph()
#nx.draw(red)
#nx.draw_random(red)
nx.draw_circular(G)
#nx.draw_spectral(G)


plt.show(red)
