from heapq import heappush, heappop
import numpy as np

def dijkstra(numNodes, adj, source, distances, maxDistance):
#    print adj
    distances [source][source] = 0.0
#    maxDistance = 1
    heap = []
    visited = [False]*numNodes
    heappush(heap, source)
    while len(heap) :
        u = heap[0]
        heap.remove(u)
        if distances[source][u] != 'inf':    
            if adj.__contains__(u):
                for v in adj[u]:
                    if not visited[v]:
                        if  (distances[source][u] + adj[u][v] < distances[source][v]):
                            distances[source][v] = distances[source][u] + adj[u][v]
                            heappush(heap, v)
        visited[u] = True
    return distances 





def dijkstraDistances(numNodes, adj, maxDistance):
    print "maxdistance is ", maxDistance
    distances = np.full((numNodes, numNodes),'inf')
    D = {}
    for vert in range(numNodes):
        D[vert] = []
        distances = dijkstra(numNodes, adj, vert, distances, maxDistance)
        for i in range(numNodes):
            if distances[vert][i] <= maxDistance:
                D[vert].append(i)
    return distances, D
