'''
This file compares the influence spread among different centrality measures while reading from json files

'''



import numpy as np
import operator 
import json
from pprint import pprint
import yaml



def readUndirectedGraph(filename, isWeighted) :
    weight = 1.0
    adj= {}
    with open(filename) as f:
       numNodes, numEdges = f.readline().split(',')
       numNodes , numEdges = int (numNodes), int (numEdges)

       for i  in range(0,int (numEdges)):
            if (isWeighted):
                node1, node2, weight = f.readline().split(',')
            else :
                node1, node2 = f.readline().split(',')
            node1, node2, weight = int (node1), int (node2), float (weight)    

            if not adj.__contains__(node1):
                adj[node1] = {}
            if node2 in adj[node1]:     
                adj[node1][node2] += float(weight)
            else:
                adj[node1][node2] = float(weight)
                    
            if not adj.__contains__(node2):
                adj[node2] = {}
            if node1 in adj[node2]:
                adj[node2][node1] += float(weight)
            else:    
                adj[node2][node1] = float(weight)

    return (numNodes, numEdges, adj)




def readRandomGraph(size):
    adj = {}
    for i in range(size):
        adj[i] = {}
        for j in range(i + 1, size):
            if not adj.__contains__(j):
                adj[j] = {}
            weight = random.randint(0,5)
            adj[i][j] = weight 
            adj[j][i] = weight
            
    return adj       
            
 
def getWcutoff(numNodes, adj, const, frac):
    alpha = [0]*numNodes
    for vert in range(numNodes):
        if adj.__contains__(vert):
            for u in adj[vert]:
                alpha[vert] += adj[vert][u]
    print alpha            
    return (np.array(const)+ np.array(alpha)* frac)           


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

if __name__ == "__main__":
  measures=[
  'eigen',
  'degree',
  'indegree',
  'betweenness',
  'closeness',
  'outdegree']
  mapping ={'power.json':'power.txt', 'hep_data.json':'hep_data', 'astro-ph.json':'astro-ph.txt'}
  files = ['power.json', 'hep_data.json', 'astro-ph.json']
  for file in files:
    numNodes, numEdges, adj = readUndirectedGraph(mapping[file], True)
    for measure in measures:
      json_data = open(file) 
      data = yaml.load(json_data)
      dict = {}
      dict = data[measure]
      print dict
      sorted_dict = sorted(dict.iteritems(), key=operator.itemgetter(1), reverse = True)
      print "sorted_dict", sorted_dict
      selectedNodes = []
      print sorted_dict 
  
      selectedNodes = [(i[0],) for i in sorted_dict] 
      weights =[0.1, 0.25, 0.5, 0.75]
      for w in weights:
        filename = 'spread'+str(file)+str(measure)+str(w)
        Wcutoff = getWcutoff(numNodes, adj, 0 , w)
        f = open(filename,'a')
        for i in range(len(selectedNodes)):     
            numSteps = 0.2* numNodes 
            print >> f, computeSpread(numNodes, adj, selectedNodes[0:i],  Wcutoff, numSteps)
      json_data.close()




