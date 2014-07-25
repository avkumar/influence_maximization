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
from dijkstra_improvable import dijkstraDistances
from exactMethods import computeBanzExactGame3
from mcMethods import computeMCGame3
from exactMethods import computeBanzExactGame4
from mcMethods import computeMCGame4
#from exactMethods import computeBanzExactGame2


def readUndirectedGraph(filename, isWeighted) :
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

#    print adj        
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

#    print adj        
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
  #  return 1.0

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
    filename = 'compareBanzGame4_Err'+str(numNodes)+str(maxIterations)+str(payoffFunction)    
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
    filename = 'compareBanzGame4_marginalErr'+ str(numNodes)+str(payoffFunction)
    f = open(filename, 'a')
    print >> f, marginalError
    print >> f, timeMCForIter 
    f.close


def compareShapAndBanz(numNodes, adj, maxIterations, Wcutoff, numSteps, selectedNodes, cutoffFraction):
    start_time = time.time()
    shapMCGame5, shapMCForIter, timeMCForIter = computeMCGame5(numNodes, adj, maxIterations, Wcutoff, 'shap')
    filename = 'compareShapAndBanz_Err'+str(numNodes)+str(cutoffFraction)    
    f = open(filename, 'a')
    marginalError = []
    for iter in range(len(shapMCForIter)):
        marginalError.append(max(np.absolute(np.array(shapMCGame5) -np.array(shapMCForIter[iter]))/ (np.array(shapMCGame5) + np.array([0.00001]))))
    print marginalError     
    print >> f, marginalError
    print >> f, timeMCForIter
    print shapMCGame5
    print shapMCGame5
    print "time for ShapMCgame5 is", time.time() - start_time
    start_time = time.time()
    banzMCGame5, banzMCForIter, timeMCForIter =  computeMCGame5(numNodes, adj, maxIterations,  Wcutoff , 'banz')
    print banzMCGame5
    print "time for BanzMCgame5 is", time.time() - start_time
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame5) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame5) + np.array([0.00001]))))
    print marginalError     
     
    print >> f, marginalError
    print >> f, timeMCForIter 
    f.close
    filename = 'spreadShapAndBanz'+str(numNodes)+str(cutoffFraction)
    f = open(filename,'a')
    print >> f, shapMCGame5 
    print >> f, banzMCGame5

    nodesRankedByShap = rankNodes(shapMCGame5)       
    nodesRankedByBanz = rankNodes(banzMCGame5)
    print >> f, "nodes ranked by Shap" 
    print >> f, nodesRankedByShap
    print >> f, "nodes ranked by Banz"
    print >> f, nodesRankedByBanz
    for i in range(selectedNodes):     
        
        print >> f, computeSpread(numNodes, adj, nodesRankedByShap[0:i],  Wcutoff, numSteps)

        print >> f, computeSpread(numNodes, adj, nodesRankedByBanz[0:i], Wcutoff, numSteps)
    print time.time() - start_time
    f.close()

def execute(handler, Graph):
    if Graph in ['AdobeCaptivate_1kFollowersGraph.list', 'DirectedEdgeWeighted4000node']:
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
        filename = 'compareBanzGame1_marginalErr'+str(numNodes)    
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
          filename = 'compareBanzGame2_marginalErr'+ str(numNodes)+str(type)    
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
          filename = 'compareBanzGame3_marginalErr'+ str(numNodes)+str(type)    
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
      #      filename = 'compareBanzGame4_marginalErr'+ str(numNodes)+str(type)    
       #     f = open(filename, 'a')
        #    banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, f1, distances, D, 'banz' )
#            print  banzMCGame4, banzMCForIter, timeMCForIter

    if handler is 5:
        maxIterations = 100000
        Wcutoff =[0.1,  0.25, 0.5]
        for w in Wcutoff:
            numTimeSteps = 100
            selectedNodes = int (0.05*numNodes)
            compareShapAndBanz(numNodes, adj, maxIterations, getWcutoff(numNodes, adj, 0, w).tolist(), numTimeSteps, selectedNodes,w)
            
    if handler is 6:
#        print "executing 5"
         print getDcutoff(numNodes, adj, 2)
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
        
            filename = 'compareShapGame1_marginalErr'    
            f = open(filename, 'a')
            print >> f, " ".join([str(x) for x in avgMarginalError/numRuns] )
            print >> f, " ".join([str(x) for x in maxMarginalError] )
            print >> f, " ".join([str(x) for x in minMarginalError] )
            print >> f, " ".join([str(x) for x in timeForSteps] )
            f.close
    if handler is 10:
        d_cutoff = 0.4
        distances, D = getDcutoff(numNodes, adj, d_cutoff)
#        print f1(2)
        maxIterations = 2 
        banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, f1, distances, D, 'banz' )
        
        filename = 'compareBanzGame4_MC_marginalErr'+ str(numNodes)+str('f1')+str(maxIterations)+str(d_cutoff)    
        f = open(filename, 'a')
        print >> f, banzMCGame4
        print >> f, banzMCForIter
        print >> f, timeMCForIter
        f.close
    if handler is 11:
        d_cutoff = 0.8
        distances, D = getDcutoff(numNodes, adj, d_cutoff)
#        print f1(2)
        maxIterations = 2 
        banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, maxIterations, f1, distances, D, 'banz' )
        
        filename = 'compareBanzGame4_MC_marginalErr'+ str(numNodes)+str('f1')+str(maxIterations)+str(d_cutoff)    
        f = open(filename, 'a')
        print >> f, banzMCGame4
        print >> f, banzMCForIter
        print >> f, timeMCForIter
        f.close
    if handler is 12:
#        input_file = 'power_res_seeds'
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
        filename = 'DirectedEdgeWeightedGraph'
        f = open(filename,'a')
        for v in adj:
            for u in adj[v]:
                print v, u, 1.0/float(len(adj[v]))
        f.close        

if __name__ =="__main__":
    np.set_printoptions(threshold=np.nan) 
    G1 = 'power.txt'
    G2 = 'astro-ph.txt'
    G3 = '4node.txt'
    G4 = '8node.txt'
    G5 = 'AdobeCaptivate_1kFollowersGraph.list' 
    G6 = 'hep_data8k'
    G7 = 'DirectedEdgeWeighted4000node'
    G8 = 'hep.txt1'
    execute(5, G8)
    
#    execute(9, G4)

##   print numNodes, numEdges, adj

#    print "shapExactGame1 is ", shapExactGame1 
#    start_banzExactGame1 = time.time() 
#    banzExactGame1 = computeBanzExactGame1(numNodes, adj)
#    print time.time() - start_banzExactGame1
#    print "banzExactGame1 is", banzExactGame1
##    nodesRankedByShap = rankNodes(shapExactGame1)       
##    nodesRankedByBanz = rankNodes(banzExactGame1)
##    shapMCGame1 = computeMCGame1(numNodes, adj, maxIteration, 'shap')
##    print "shapMCGame1 is", shapMCGame1
#    start_banzMCGame1 = time.time()
#    banzMCGame1, banzMCForIter, timeMCForIter =  computeMCGame1(numNodes, adj, maxIteration, 'banz')
#    print time.time() - start_banzMCGame1
#    print "banMCGame1 is", banzMCGame1
##    print "banzMCForIter is", banzMCForIter
#    marginalError = []
#    for iter in range(len(banzMCForIter)):
#        marginalError.appen(max(np.absolute(np.array(banzMCGame1) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame1) + np.array([0.00001]))))
#    print marginalError     
                      
 
 #    shapExactGame2 = computeShapExactGame2(numNodes, adj, numThreshold) 
#    print "shapExactGame2 is ", shapExactGame2 
  #  numThreshold = [1] * numNodes
   # banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshold)
    #print "banzExactGame2 is", banzExactGame2
#    shapMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshold, 'shap')
#    print "shapMCGame2 is", shapMCGame2
#    banzMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'banz')
#    print "banMCGame2 is", banzMCGame2
#    shapNodes = [2564,335,1027,1232,231,256,913,18,3272,118,2393,37,40,347,1078,2868,1489,1354,2879,5508,1165,469,1836,984,2339,236,5503,80,346,560,409,456,1,2567,2394,4674,4741,4287,1797,2909,5027,333,871,493,1190,4073,516,3143,97,61,1484,2294,211,1558,401,1619,5701,6199,149,832,752,1710,3692,792,327,94,600,1909,306,931,481,556,6846,1476,3507,436,1975,62,4970,1369,357,987,1203,2165,3865,4719,511,711,232,3260,1202,1976,235,265,847,43,1605,748,1137,1533,3180,952,3028,6202]
 #   banzNodes =  [2564,231,37,3272,913,516,1558,1078,1232,335,2294,256,2879,1484,2868,1836,871,401,556,2394,18,4741,1165,1203,1489,984,5995,2393,1027,1354,409,236,1044,118,3507,3671,3860,347,1,5159,3143,7918,493,1190,2492,1996,560,5508,346,5943,4670,792,456,436,5227,8033,97,2228,61,602,4165,1369,3965,715,5027,152,2669,1202,4970,518,952,4073,3273,2360,2007,987,344,1024,40,1619,3046,2600,1124,497,2165,1932,511,8000,1975,3765,247,5079,861,5878,3457,4674,5558,8657,80,3308,1179,1255,5007,2567]
#    shapMCGame5 = computeMCGame5(numNodes, adj, 12, Wcutoff, 'shap')
#    print "shapMCGame5 is", shapMCGame5
#    banzMCGame5 = computeMCGame5(numNodes, adj, 3, Wcutoff, 'banz')
#    print "banzMCGame5 is", banzMCGame5
#    numThreshhold = [2]*numNodes
#
#    shapExactGame2 = computeShapExactGame2(numNodes, adj, numThreshhold) 
#    print "shapExactGame2 is ", shapExactGame2 
#    banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshhold)
#    print "banzExactGame2 is", banzExactGame2
#    shapMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'shap')
#    print "shapMCGame2 is", shapMCGame2
#    banzMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'banz')
#
#    shapExactBruteGame5 =  computeShapBruteForceGame5(numNodes, adj, Wcutoff)
#    print shapExactBruteGame5
#    banzExactBrute =  computeBanzBruteForceGame5(numNodes, adj, Wcutoff)
#    banzExactBruteGame5 = computeBanzBruteForceGame5(numNodes, adj, Wcutoff)


