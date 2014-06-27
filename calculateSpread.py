import random
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from itertools import combinations 

from itertools import permutations 


def power_set(List):
    subs = [list(j) for i in range(len(List)) for j in combinations(List, i+1)]
    return subs

def init(G):
    newG = G
    for i in range(1, len(G.edge)):
        for j in range(1, len(G.edge[i])):
            newG[i][j] = random.uniform(0, 1)
            print newG[i][j]
    return newG

def getGraph(G):
    newG = G 
    for i in newG.edges_iter():
        print i
    return newG
#	for i in range(1, len(G.edge)):
#		for j in range(1, len(G.edge[i])):
#		    if (G[i][j]]+ random.uniform(0, 1) > 1):
#                newG[i][j] = 1


def computeSpread(num_iterations):
    for m in range(1, num_iterations):
        G_m = getGraph(G)
        calculateValueOfCoalitions(G_m)
        


def calculateValueOfCoalitions(G_m):
    subs = power_set(G.nodes())
    print subs
#    for i in G_m.node():





if __name__ == "__main__":
    G = nx.read_gml('karate.gml')
    G_Prob = init(G)
    print G_Prob
    num_iterations=1
#    computeSpread(num_iterations)
#    calculateValueOfCoalitions(G) 



