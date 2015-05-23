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
#from exactMethods import computeBanzExactGame2


def readGraph(filename, isWeighted) :
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
            if not adj.__contains__(node2):
                adj[node2] = {}
            adj[node2][node1] = float(weight)

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
    return distances, D


def f1(x):
  #  return 1.0/(1 + x)
    return 1.0

def f2(x):
    return 1.0/pow(2, x)
        


def compareBanzGame1(numNodes, adj):
    start_time = time.time()
    banzExactGame1 = computeBanzExactGame1(numNodes, adj) 
    print "banzExactGame1 is ", banzExactGame1 
    start_banzMCGame1 = time.time() 
    print "computation time for banzExactGame1 is ", start_banzMCGame1 - start_time 
    banzMCGame1, banzMCForIter, timeMCForIter =  computeMCGame1(numNodes, adj, 12000, 'banz')
    print "banzMCGame1 is ", banzMCGame1
    print "computation time for banzMCGame1 is", time.time() - start_banzMCGame1
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame1) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame1) + np.array([0.00001]))))
    filename = 'compareBanzGame1_marginalErr'    
    f = open(filename, 'a')
    print >> f, marginalError
    f.close
    print marginalError     
    


def compareBanzGame2(numNodes,adj,numThreshold):
    start_time = time.time()
    banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshold) 
    print "banzExactGame2 is ", banzExactGame2 
    start_banzMCGame2 = time.time()
    print "computation time for banzExactGame2 is ", start_banzMCGame2 - start_time 
    banzMCGame2, banzMCForIter, timeMCForIter =  computeMCGame2(numNodes, adj, 120000, numThreshold, 'banz')
    print "banzMCGame2 is", banzMCGame2
    print "computation time for banzMCGame2 is", time.time() - start_banzMCGame2
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame2) - np.array(banzMCForIter[iter]))/ (np.array(banzMCGame2) + np.array([0.00001]))))
    print marginalError     
    filename = 'compareBanzGame2_marginalErr'
    f = open(filename, 'a')
    print >> f, marginalError
    f.close



def compareBanzGame3(numNodes, adj, D):
    start_time = time.time()
    banzExactGame3 = computeBanzExactGame3(numNodes, adj, D) 
    print "banzExactGame3 is ", banzExactGame3 
    start_banzMCGame3 = time.time()
    print "computation time for banzExactGame3 is ", start_banzMCGame3 - start_time 
    banzMCGame3, banzMCForIter, timeMCForIter =  computeMCGame3(numNodes, adj, 120000, D, 'banz')
    print "banzMCGame3 is", banzMCGame3
    print "computation time for banzMCGame3 is", time.time() - start_banzMCGame3
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame3) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame3) + np.array([0.00001]))))
    print marginalError     
    filename = 'compareBanzGame3_marginalErr'
    f = open(filename, 'a')
    print >> f, marginalError
    f.close


def compareBanzGame4(numNodes, adj, f, distances, D):
    start_time = time.time()
    print "D" ,D
    print "distances", distances
    banzExactGame4 = computeBanzExactGame4(numNodes, adj,f ,distances, D) 
    print "banzExactGame4 is ", banzExactGame4 
    start_banzMCGame4 = time.time()
    print "computation time for banzExactGame4 is ", start_banzMCGame4 - start_time 
    banzMCGame4, banzMCForIter, timeMCForIter =  computeMCGame4(numNodes, adj, 12, D, 'banz')
    print "banzMCGame4 is", banzMCGame4
    print "computation time for banzMCGame4 is", time.time() - start_banzMCGame4
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame4) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame4) + np.array([0.00001]))))
    print marginalError     
    filename = 'compareBanzGame4_marginalErr'
    f = open(filename, 'a')
    print >> f, marginalError
    f.close


def compareShapAndBanz(numNodes, adj, maxIterations, Wcutoff, numSteps, selectedNodes):
    start_time = time.time()
    shapMCGame5, shapMCForIter, timeMCForIter = computeMCGame5(numNodes, adj, maxIterations, Wcutoff, 'shap')
    filename = 'compareShapAndBanz_mErr'    
    f = open(filename, 'a')
    marginalError = []
    for iter in range(len(shapMCForIter)):
        marginalError.append(max(np.absolute(np.array(shapMCGame5) -np.array(shapMCForIter[iter]))/ (np.array(shapMCGame5) + np.array([0.00001]))))
    print marginalError     
    print >> f, marginalError, timeMCForIter
    print shapMCGame5
    print time.time() - start_time
    banzMCGame5, banzMCForIter, timeMCForIter =  computeMCGame5(numNodes, adj, maxIterations,  Wcutoff , 'banz')
    print banzMCGame5
    marginalError = []
    for iter in range(len(banzMCForIter)):
        marginalError.append(max(np.absolute(np.array(banzMCGame5) -np.array(banzMCForIter[iter]))/ (np.array(banzMCGame5) + np.array([0.00001]))))
    print marginalError     
     
    print >> f, marginalError
    f.close
    nodesRankedByShap = rankNodes(shapMCGame5)       
    nodesRankedByBanz = rankNodes(banzMCGame5)
    for i in range(selectedNodes):     
        
        print computeSpread(numNodes, adj, nodesRankedByShap[0:5*i],  Wcutoff, numSteps)

        print computeSpread(numNodes, adj, nodesRankedByBanz[0:5*i], Wcutoff, numSteps)
    print time.time() - start_time


def execute(handler, Graph):
    numNodes, numEdges,  adj =  readGraph(Graph, True)
    if handler is 1:
        compareBanzGame1(numNodes, adj)
    if handler is 2:    
        compareBanzGame2(numNodes, adj, getThreshold(numNodes, adj, 2, 0))
        compareBanzGame2(numNodes, adj, getThreshold(numNodes, adj, 0, 0.5))
        compareBanzGame2(numNodes, adj, getThreshold(numNodes, adj, 0, 0.75))
    if handler is 3:
        distances, D = getDcutoff(numNodes, adj, 'inf')
        avgDis = computeAvgDist(numNodes, distances)
        print avgDis

        distances,D = getDcutoff(numNodes, adj, 3)
        print distances, D
        compareBanzGame3(numNodes, adj, D)
#        distances,D = getDcutoff(numNodes, adj, avgDis*0.25)
#        compareBanzGame3(numNodes, adj, D)
#        distances,D = getDcutoff(numNodes, adj, avgDis/8)
#        compareBanzGame3(numNodes, adj, D)
    if handler is 4:
        distances, D = getDcutoff(numNodes, adj, 1)
                
        compareBanzGame4(numNodes, adj, f1, distances, D)

    if handler is 5:
        compareShapAndBanz(numNodes, adj, 300, getWcutoff(numNodes, adj, 0, 0).tolist(), 1000, 500)
    if handler is 6:
#        print "executing 5"
         print getDcutoff(numNodes, adj, 2)
    if handler is 7:
        start_time = time.time()
        dijkstraDistances(numNodes, adj, 100000)
#        for i in range(numNodes):
#            for j in range(numNodes):
#                print distances[i][j],
#            print    
#        for i in range(numNodes):
#            for j in D[i]:
#                print D[i][j],
#            print          
#        print time.time() - start_time()
if __name__ =="__main__":
    G1 = 'power.txt'
    G2 = 'astro-ph.txt'
    G3 = '4node.txt'
    G4 = '8node.txt'
    
    execute(7, G2)
#    execute(3, G4)

##   print numNodes, numEdges, adj
#    shapExactGame1 = computeShapExactGame1(numNodes, adj) 
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


