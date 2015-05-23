'''
This file runs Dijkstra's algorithm on a given graph for all nodes
'''

import numpy as np

def dijkstra(numNodes, adj, source, distances, maxDistance, D):
    def greater(a, b):
        if distances[source][a] < distances[source][b]:
            return -1 
        if distances[source][a] > distances[source][b]:
            return 1
        return a < b
    distances [source][source] = 0.0
    heap = [i for i in range(numNodes)]
    heap = sorted(heap, cmp=greater)
    while len(heap) :
        u = heap[0]
        heap.remove(u)
        if maxDistance !=-1 and distances[source][u]  > maxDistance :   #If the exploration distance crosses max distance then we stop
            break
        if distances[source][u] != 'inf':    
            if adj.__contains__(u):
                for v in adj[u]:
                        if  (distances[source][u] + adj[u][v] < distances[source][v]):
                            heap.remove(v)
                            distances[source][v] = distances[source][u] + adj[u][v]
                            heap.append(v)
                            heap = sorted(heap, cmp=greater)

        D[source].append(u)
    return distances,D 


def dijkstraDistances(numNodes, adj, maxDistance):
    print "maxdistance is ", maxDistance
    distances = np.empty((numNodes, numNodes,))
    distances[:] = 'inf'
    D = {}
    for vert in range(numNodes):
        D[vert] = []
        distances,D = dijkstra(numNodes, adj, vert, distances, maxDistance, D)
#        print distances
    return distances,D
