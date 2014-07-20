from heapq import heappush, heappop
import numpy as np

def dijkstra(numNodes, adj, source, distances, maxDistance, D):
#    print adj
    distances [source][source] = 0.0
#    maxDistance = 1
    heap = []
    visited = [False]*numNodes
    heappush(heap, source)
    while len(heap) :
        u = heap[0]
        heap.remove(u)
#        print "removed", u
        if maxDistance !=-1 and distances[source][u] > maxDistance:
            break
        if adj.__contains__(u):
            for v in adj[u]:
#                print "adj of ",u,"and ",v," is",adj[u][v]
                if not visited[v]:
 #                   print "inside adj of ",u,"and ",v," is",adj[u][v]
  #                  print distances[source][u], adj[u][v], distances[source][v]
                    if (distances[source][u] + adj[u][v] < distances[source][v]):
                        distances[source][v] = distances[source][u] + adj[u][v]
                        heappush(heap, v)
        visited[u] = True
        D[source].append(u)
    return distances,D 





def dijkstraDistances(numNodes, adj, maxDistance):
    distances = np.full((numNodes, numNodes),'inf')
    D = {}
    for vert in range(numNodes):
        D[vert] = []
        distances,D = dijkstra(numNodes, adj, vert, distances, maxDistance, D)
#        print distances
    return distances,D
