import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


# we have number, connections, and then remodel
G = nx.watts_strogatz_graph(n=2000,k=20,p=0.00)
nx.draw_circular(G)
plt.show()
print(nx.average_clustering(G))
