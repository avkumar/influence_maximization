

def computeSpread(numNodes,  adj, selectedNodes, Wcutoff, numIterations):
    active = [False]* numNodes
    Weights = [0]* numNodes
    for node in selectedNodes:
#       print node
       active [node] = True

    newActiveNodes = list(selectedNodes)
    numActiveNodesAtIter = [0]*numIterations

    for iter in range(numIterations):
        lastIterActivatedNodes = list(newActiveNodes)
        newActiveNodes = []
        for node in  lastIterActivatedNodes:
            if adj.__contains__(node):
                for u in adj[node]:
                    if not active[u]:
                        Weights[u] += adj[u][node] 
                        if Weights[u] > Wcutoff[u]:
                            active[u] = True
                            newActiveNodes.append(u)          
        numActiveNodesAtIter[iter] = sum(active)                    
    return numActiveNodesAtIter
