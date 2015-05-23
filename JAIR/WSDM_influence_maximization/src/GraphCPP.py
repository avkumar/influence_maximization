

'''
This example uses our methods of computation of Shapley value and Banzhaf index for all the five games
'''



#! /usr/bin/env python
import time
import exactMethods
import mcMethods
from exactMethods import computeShapExactGame1
from mcMethods import computeMCGame1
from exactMethods import computeBanzExactGame1
from mcMethods import computeMCGame5
from exactMethods import computeBanzBruteForceGame5
from exactMethods import computeShapBruteForceGame5
from computeInfluence import computeSpread
from exactMethods import computeBanzExactGame2
from mcMethods import computeMCGame2
import numpy as np
from dijkstra import dijkstraDistances
from exactMethods import computeBanzExactGame3
from mcMethods import computeMCGame3
from exactMethods import computeBanzExactGame4
from mcMethods import computeMCGame4


def readUndirectedGraph(filename, isWeighted) :
    weight = 1.0
    adj= {}
    with open(filename) as f:
        numNodes, numEdges = f.readline().split(',')
        numNodes , numEdges = int (numNodes), int (numEdges)
        for i in range(0,int (numEdges)):
            if (isWeighted):
                print i 
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

def readDirectedGraph(filename, isWeighted) :
    weight = 1.0
    adj= {}
    with open(filename) as f:
        numNodes, numEdges = f.readline().split(',')
        numNodes , numEdges = int (numNodes), int (numEdges)
        for i in range(0,int (numEdges)):
            if (isWeighted):
                node1, node2, weight = f.readline().split(',')
            else :
                node1, node2 = f.readline().split(',')
            node1, node2, weight = int (node1), int (node2), float (weight)    

            if not adj.__contains__(node1):
                adj[node1] = {}
            adj[node1][node2] = float(weight)
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
            
                        

def computeAvgDist(numNodes, distances):
    avgDis = 0 
    count = 0
    for i in range(numNodes):
        for j in range(numNodes):
            if distances[i][j] != 'inf' :
                avgDis += distances[i][j]
                count += 1
    avgDis /= float (count)
    return avgDis


def getWcutoff (numNodes):
    Wcutoff = [0]*numNodes
    return Wcutoff

def rankNodes(powerIndices):
    idx = sorted(xrange(len(powerIndices)), key=powerIndices.__getitem__)
    idx.reverse()
    return idx

def getThreshold(numNodes, adj, const, frac):
    threshold = [0]*numNodes
    for vert in range(numNodes):
        if adj.__contains__(vert):
            threshold[vert] = (const + frac * (len(adj[vert])) ) 
            if threshold[vert] - int (const + frac * (len(adj[vert])) ) > 0:
                threshold[vert] = int (const + frac * (len(adj[vert])) ) + 1 
            else:
                threshold[vert] = const + frac * (len(adj[vert]))   
                    
    return threshold

def getWcutoff(numNodes, adj, const, frac):
    alpha = [0]*numNodes
    for vert in range(numNodes):
        if adj.__contains__(vert):
            for u in adj[vert]:
                alpha[vert] += adj[vert][u]
    print alpha            
    return (np.array(const)+ np.array(alpha)* frac)           


def getDcutoff(numNodes, adj, d_cutoff):
    print "d_cutoff is", d_cutoff
    distances, D = dijkstraDistances(numNodes, adj, d_cutoff)
    print distances
    return distances, D


def f1(x):
    if x is 'inf':
        return 0
    return float (1.0/(1 + x))

def f2(x):
    if x is 'inf':
        return 0
    return 1.0/(1 + pow(2.0, x))
        

def compareBanzGame1(numNodes, adj, maxIterations):
    start_time = time.time()
    banzExactGame1 = computeBanzExactGame1(numNodes, adj) 
    print "banzExactGame1 is ", banzExactGame1 
    start_banzMCGame1 = time.time() 
    print "computation time for banzExactGame1 is ", start_banzMCGame1 - start_time 
    banzMCGame1, banzMCForIter, timeMCForIter =  computeMCGame1(numNodes, adj, maxIterations, 'banz')
    print "banzMCGame1 is ", banzMCGame1
    print "computation time for banzMCGame1 is", time.time() - start_banzMCGame1
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame1) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame1) + np.array([0.00001]))))
    print marginalError     
    print timeMCForIter
    return marginalError, timeMCForIter  



def compareShapGame1(numNodes, adj, maxIterations):
    start_time = time.time()
    shapExactGame1 = computeShapExactGame1(numNodes, adj) 
    print "shapExactGame1 is ", shapExactGame1 
    start_shapMCGame1 = time.time() 
    print "computation time for shapExactGame1 is ", start_shapMCGame1 - start_time 
    shapMCGame1, shapMCForIter, timeMCForIter =  computeMCGame1(numNodes, adj, maxIterations, 'shap')
    print "shapMCGame1 is ", shapMCGame1
    print "computation time for shapMCGame1 is", time.time() - start_shapMCGame1
    marginalError = []
    for iter in range(len(shapMCForIter)):
        marginalError.append(max(np.absolute(np.array(shapMCGame1) -np.array(shapMCForIter[iter]))/ (np.array(shapMCGame1) + np.array([0.00001]))))
    print marginalError     
    print timeMCForIter
    return marginalError, timeMCForIter  



def compareBanzGame2(numNodes,adj,numThreshold, maxIterations):
    start_time = time.time()
    banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshold) 
    print "banzExactGame2 is ", banzExactGame2 
    start_banzMCGame2 = time.time()
    print "computation time for banzExactGame2 is ", start_banzMCGame2 - start_time 
    banzMCGame2, banzMCForIter, timeMCForIter =  computeMCGame2(numNodes, adj, maxIterations, numThreshold, 'banz')
    print "banzMCGame2 is", banzMCGame2
    print "computation time for banzMCGame2 is", time.time() - start_banzMCGame2
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame2) - np.array(banzMCForIter[iter]))/ (np.array(banzMCGame2) + np.array([0.00001]))))
    return marginalError, timeMCForIter  


def compareBanzGame3(numNodes, adj, D, maxIterations):
    start_time = time.time()
    banzExactGame3 = computeBanzExactGame3(numNodes, adj, D) 
    print "banzExactGame3 is ", banzExactGame3 
    start_banzMCGame3 = time.time()
    print "computation time for banzExactGame3 is ", start_banzMCGame3 - start_time 
    banzMCGame3, banzMCForIter, timeMCForIter =  computeMCGame3(numNodes, adj, maxIterations, D, 'banz')
    print "banzMCGame3 is", banzMCGame3
    print "computation time for banzMCGame3 is", time.time() - start_banzMCGame3
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame3) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame3) + np.array([0.00001]))))
    print marginalError     
    return marginalError, timeMCForIter  


def compareBanzGame4(numNodes, adj, payoffFunction, distances, D, maxIterations):
    start_time = time.time()
    filename = '../results/compareBanzGame4_Err'+str(numNodes)+str(maxIterations)+str(payoffFunction)    
    f = open(filename, 'a')
    banzExactGame4 = computeBanzExactGame4(numNodes, adj, payoffFunction ,distances, D) 
    print >> f, "banzExactGame4 is,"
    print >> f,  banzExactGame4 
    start_banzMCGame4 = time.time()
    print >> f, "computation time of banzExactGame4 is,"
    print >> f,  start_banzMCGame4 - start_time 
    banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, payoffFunction, distances, 'D',  'banz')
    print >> f, "banzMCGame4 is"
    print >> f,  banzMCGame4
    print >> f,  time.time() - start_banzMCGame4
    f.close
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame4) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame4) + np.array([0.00001]))))
    print marginalError     
    filename = '../results/compareBanzGame4_marginalErr'+ str(numNodes)+str(payoffFunction)
    f = open(filename, 'a')
    print >> f, marginalError
    print >> f, timeMCForIter 
    f.close


def compareShapAndBanz(numNodes, adj, maxIterations, Wcutoff, numSteps, selectedNodes, cutoffFraction):
    shapMCGame5, shapMCForIter, timeMCForIter = computeMCGame5(numNodes, adj, maxIterations, Wcutoff, 'shap')
    banzMCGame5, banzMCForIter, timeMCForIter =  computeMCGame5(numNodes, adj, maxIterations,  Wcutoff , 'banz')
    filename = '../results/spreadShapAndBanz'+str(numNodes)+str(cutoffFraction)
    f = open(filename,'a')
    nodesRankedByShap = rankNodes(shapMCGame5)       
    nodesRankedByBanz = rankNodes(banzMCGame5)
    print >> f, "nodes ranked by Shap" 
    print >> f, nodesRankedByShap
    print "nodes Ranked By Shap", nodesRankedByShap
    print >> f, "nodes ranked by Banz"
    print >> f, nodesRankedByBanz
    print "nodes ranked by Banz", nodesRankedByBanz
    for i in range(selectedNodes):     
        print >> f, computeSpread(numNodes, adj, nodesRankedByShap[0:i],  Wcutoff, numSteps)

        print >> f, computeSpread(numNodes, adj, nodesRankedByBanz[0:i], Wcutoff, numSteps)
    print time.time() - start_time
    f.close()

def execute(handler, Graph):
    if Graph in ['twitter.txt', 'DirectedEdgeWeighted4000node']:
        numNodes, numEdges,  adj =  readDirectedGraph(Graph, True)
    else:
        numNodes, numEdges,  adj =  readUndirectedGraph(Graph, True)
    if handler is 1:
        maxIterations = 5000 
        numRuns = 20
        avgMarginalError = np.zeros(maxIterations/5)
        maxMarginalError = np.zeros(maxIterations/5)
        minMarginalError = np.zeros(maxIterations/5) 
        minMarginalError.fill(100)
        for i in range(numRuns): 
            marginalError, timeForSteps = compareBanzGame1(numNodes, adj, maxIterations)
            avgMarginalError += np.array(marginalError)
            maxMarginalError = np.maximum(maxMarginalError, marginalError)
            minMarginalError = np.minimum(minMarginalError, marginalError)  
        print avgMarginalError/numRuns , maxMarginalError, minMarginalError, timeForSteps
        filename = '../results/compareBanzGame1_marginalErr'+str(numNodes)    
        f = open(filename, 'a')
        print >> f, " ".join([str(x) for x in avgMarginalError/numRuns] )
        print >> f, " ".join([str(x) for x in maxMarginalError] )
        print >> f, " ".join([str(x) for x in minMarginalError] )
        print >> f, " ".join([str(x) for x in timeForSteps] )
        f.close
    
    if handler is 2:    
        intercepts = [2, 0, 0]
        weights = [0, 0.5, 0.75]
        for type in range(len(intercepts)):
          maxIterations = 5000 
          numRuns =20
          avgMarginalError = np.zeros(maxIterations/5)
          maxMarginalError = np.zeros(maxIterations/5)
          minMarginalError = np.zeros(maxIterations/5) 
          minMarginalError.fill(100)
          for i in range(numRuns): 
              print i
              marginalError, timeForSteps = compareBanzGame2(numNodes, adj, getThreshold(numNodes, adj, intercepts[type], weights[type]), maxIterations)
              avgMarginalError += np.array(marginalError)
              maxMarginalError = np.maximum(maxMarginalError, marginalError)
              minMarginalError = np.minimum(minMarginalError, marginalError)  
          print avgMarginalError/numRuns , maxMarginalError, minMarginalError, timeForSteps
          filename = '../results/compareBanzGame2_marginalErr'+ str(numNodes)+str(type)    
          f = open(filename, 'a')
          print >> f, " ".join([str(x) for x in avgMarginalError/numRuns] )
          print >> f, " ".join([str(x) for x in maxMarginalError] )
          print >> f, " ".join([str(x) for x in minMarginalError] )
          print >> f, " ".join([str(x) for x in timeForSteps] )
          f.close


    if handler is 3:
        avgdist = 6.91
        weights = avgdist*np.array([0.125, 0.25])
        for type in range(len(weights)):
          maxIterations = 5000 
          numRuns = 3
          avgMarginalError = np.zeros(maxIterations/5)
          maxMarginalError = np.zeros(maxIterations/5)
          minMarginalError = np.zeros(maxIterations/5) 
          minMarginalError.fill(100)
          distances,D = getDcutoff(numNodes, adj, weights[type])
          for i in range(numRuns): 
              marginalError, timeForSteps = compareBanzGame3(numNodes, adj, D, maxIterations)
              avgMarginalError += np.array(marginalError)
              maxMarginalError = np.maximum(maxMarginalError, marginalError)
              minMarginalError = np.minimum(minMarginalError, marginalError)  
          print avgMarginalError/numRuns , maxMarginalError, minMarginalError, timeForSteps
          filename = '../results/compareBanzGame3_marginalErr'+ str(numNodes)+str(type)    
          f = open(filename, 'a')
          print >> f, " ".join([str(x) for x in avgMarginalError/numRuns] )
          print >> f, " ".join([str(x) for x in maxMarginalError] )
          print >> f, " ".join([str(x) for x in minMarginalError] )
          print >> f, " ".join([str(x) for x in timeForSteps] )
          f.close

    if handler is 4:
        start_time = time.time()
        distances, D = getDcutoff(numNodes, adj, 'inf')
        print "time taken for distances, D is",time.time() - start_time  
        start_time = time.time()
        maxIterations = 30000         
        for f in [f1, f2]:
            compareBanzGame4(numNodes, adj,  f, distances, D, maxIterations)
            print "function is completed in time ", time.time() - start_time

    if handler is 5:
        maxIterations = 500000
        Wcutoff =[0.1]
        for w in Wcutoff:
            numTimeSteps = 100
            selectedNodes = int (0.05*numNodes)
            compareShapAndBanz(numNodes, adj, maxIterations, getWcutoff(numNodes, adj, 0, w).tolist(), numTimeSteps, selectedNodes,w)
    if handler is 7:
        start_time = time.time()
        dist,D = dijkstraDistances(numNodes, adj, 100000)
        for i in range(numNodes):
            for j in range(numNodes):
                print dist[i][j],
            print    
        for i in range(numNodes):
            for j in D[i]:
                print D[i][j],
            print          
        print time.time() - start.time()
    if handler is 9:
            maxIterations = 1200 
            numRuns = 2
            avgMarginalError = np.zeros(maxIterations/5)
            maxMarginalError = np.zeros(maxIterations/5)
            minMarginalError = np.zeros(maxIterations/5) 
            minMarginalError.fill(100)
            for i in range(numRuns): 
                marginalError, timeForSteps = compareShapGame1(numNodes, adj, maxIterations)
                avgMarginalError += np.array(marginalError)
                maxMarginalError = np.maximum(maxMarginalError, marginalError)
                minMarginalError = np.minimum(minMarginalError, marginalError)  
            print avgMarginalError/numRuns , maxMarginalError, minMarginalError, timeForSteps   
        
            filename = '../results/compareShapGame1_marginalErr'    
            f = open(filename, 'a')
            print >> f, " ".join([str(x) for x in avgMarginalError/numRuns] )
            print >> f, " ".join([str(x) for x in maxMarginalError] )
            print >> f, " ".join([str(x) for x in minMarginalError] )
            print >> f, " ".join([str(x) for x in timeForSteps] )
            f.close
    if handler is 10:
        d_cutoff = 0.4
        distances, D = getDcutoff(numNodes, adj, d_cutoff)
        maxIterations = 2 
        banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, f1, distances, D, 'banz' )
        
        filename = '../results/compareBanzGame4_MC_marginalErr'+ str(numNodes)+str('f1')+str(maxIterations)+str(d_cutoff)    
        f = open(filename, 'a')
        print >> f, banzMCGame4
        print >> f, banzMCForIter
        print >> f, timeMCForIter
        f.close
    if handler is 11:
        d_cutoff = 0.8
        distances, D = getDcutoff(numNodes, adj, d_cutoff)
        maxIterations = 2 
        banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, f1, distances, D, 'banz' )
        
        filename = '../results/compareBanzGame4_MC_marginalErr'+ str(numNodes)+str('f1')+str(maxIterations)+str(d_cutoff)    
        f = open(filename, 'a')
        print >> f, banzMCGame4
        print >> f, banzMCForIter
        print >> f, timeMCForIter
        f.close
    if handler is 12:
        files =['spreadShap49410.1','spreadBanz49410.1']
        for input_file in files:
            with open(input_file, 'r') as file:
                seeds = file.read()
                print seeds
                seeds = [int (i) for i in seeds.split(',')]    
                Wcutoff = 0.1
                numSteps = 100
                filename = "computeSpread"+str(input_file)+str(numNodes)+str(Wcutoff)
                f = open(filename, 'a')
                for i in range(len(seeds)):
                    print >> f, computeSpread(numNodes, adj, seeds[0:i], getWcutoff(numNodes, adj, 0, Wcutoff), numSteps)
                f.close
    if handler is 13:
        filename = '../results/DirectedEdgeWeightedGraph'
        f = open(filename,'a')
        for v in adj:
            for u in adj[v]:
                print v, u, 1.0/float(len(adj[v]))
        f.close        

if __name__ =="__main__":
    np.set_printoptions(threshold=np.nan) 
    G1 = '../data/power.txt'
    G2 = '../data/astro-ph.txt'
    G3 = '../data/4node.txt'
    G4 = '../data/8node.txt'
    G5 = '../data/twitter.txt' 
    G6 = '../data/hep_data8k'
    G7 = '../data/DirectedEdgeWeighted4000node'
    G8 = '../data/hep.txt1'
    G9 = '../data/phy.txt'
    execute(1, G4)
