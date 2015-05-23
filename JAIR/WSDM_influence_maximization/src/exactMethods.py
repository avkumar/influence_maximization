'''
This file contains all the methods that computes exact Shapley values and exact Banzhaf  values for the games discussed in 
"Banzhaf index based approach for influence maximization"
'''



import math
import itertools 
from itertools import permutations
from itertools import chain, combinations
import numpy as np
from operator import mul    
from fractions import Fraction

def nCk(n,k): 
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )



def computeBanzExactGame1(numNodes, adj):
    banzExact = [0]*numNodes
    for i in range(numNodes):
        banzExact[i] += (1.0 / pow(2.0,  len(adj[i]) ) )
        for u in adj[i]:
            banzExact[i] += ( 1.0 / pow( 2.0, len( adj[u] ))  )
    return banzExact


def computeShapExactGame1(numNodes, adj):
    shapExact = [0]*numNodes
    for i in range(0, numNodes):
        if adj.__contains__(i):
            shapExact[i] += (1.0 / ( 1 +  len(adj[i]))  )
            for u in adj[i]:
                if adj[i].__contains__(u):
                    shapExact[i] += ( 1.0 / (1 +  len( adj[u] ))  )
    return shapExact


def computeShapExactGame2(numNodes, adj, numThreshold):
    shapExact = [0]*numNodes
    for i in range(numNodes):
        shapExact[i] = min(1.0,  (( float) (numThreshold[i])) /(1.0+ (float) (len(adj[i]))) )
        for u in adj[i]:
            degree = float (len(adj[u]))
            shapExact[i] = max(0.0, ( degree - (float) (numThreshold[u]) + 1.0)/	( degree*(1.0+degree))  )
    return shapExact
    					 

def computeBanzExactGame2(numNodes, adj, numThreshold):
    banzExact = [0]*numNodes
    for vert in range(numNodes):
        if adj.__contains__(vert):
            for r in range(int(numThreshold[vert])):
                banzExact[vert] +=     nCk ( len(adj[vert]), r) * (1.0 / pow(2.0,  len(adj[vert]) )  )
            for u in adj[vert]:
                banzExact[vert] +=     nCk( len(adj[u]) - 1, int (numThreshold[u] ) - 1 ) * (1.0 / pow(2.0,  len(adj[u]) )  )
    return banzExact
 

def computeShapExactGame3(numNodes, adj, D):
    shapExact = [0]*numNodes
    for vert in range(numNodes):
        for u in D[i]:
            extDegree = len(D[u]) - 1
            shapExact[vert] +=  (1.0 / (1.0 + extDegree)) 
    return shapExact

def computeBanzExactGame3(numNodes, adj, D):
    banzExact = [0]*numNodes
    for vert in range(numNodes):
        for u in D[vert]:
            extDegree = len(D[u]) - 1 
            if extDegree < 200: 
                banzExact[vert] +=  (1.0 / pow(2.0, extDegree)) 
    return banzExact
            
            
def computeBanzExactGame4(numNodes, adj, f, distances, D):
   banzExact = [0]*numNodes
   for i in range(numNodes):
        sum = 0.0 
        currSV = 0.0 
        index = numNodes - 1
        prevDistance = -1.0
        prevSV = -1.0
        wart = 0.0
        while index > 0:
            print D[i][index]
            print  distances[i][D[i][index]]
            if index < 100:                   # To prevent from out of precision error
                wart = f(distances[i][D[i][index]])/pow(2.0, index)
            else:
                wart = 0
            if distances[i][D[i][index]] == prevDistance:
                currSV = prevSV
            else:
                currSv = wart - sum
            banzExact[D[i][index]] += currSV
            if index < 100:
                sum += (wart/pow(2.0, index - 1))
            prevDistance = distances[i][D[i][index]]
            prevSV = currSV
            index = index - 1
        banzExact[i] += (f(0) - sum)    
    
   return banzExact


#12 node graph 


def computeFactorial(n):
    factorial = 1
    for i in range(1, n+1):
        factorial = factorial * i
    return factorial   


factorials_12 = [0,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600]

def computeShapBruteForceGame5(numNodes, adj, Wcutoff):
    factorial =  computeFactorial(numNodes)
    nodes = [i for i in range(numNodes)]
    shapExact = [0]*numNodes

    for permutation in itertools.permutations(nodes):
        Counted = [False] * numNodes
        Weights = [0] * numNodes
        for i in range(numNodes):
            vert = permutation[i]
            for u in adj[vert]:
                weight = adj[vert][u]
                if not Counted[u]:
                    Weights[u] += weight
                    if Weights >= Wcutoff[u] :
                        Counted[u]= True
                        shapExact[vert] += 1

            if not Counted[vert]:
                shapExact[vert] += 1
                Counted[vert] = True

    shapExact = [(shapExact[i] / float(factorial )) for i in range(numNodes)]
    
    return shapExact


def computeBanzBruteForceGame5(numNodes, adj, Wcutoff):
    Nodes = [i for i in range(numNodes)]
    numCoalitions = 0
    BanzExact = [0]*numNodes
    for coalition in chain.from_iterable(combinations(Nodes, r) for r in range(len(Nodes)+1)): 
        Counted = [False]*numNodes
        Weights = [0]*numNodes
        for vert in coalition:
            Counted[vert] = True
            for u in adj[vert]:
                Weights[u] += adj[u][vert]
                if Weights[u] >= Wcutoff[u]:
                    Counted[u] = True
                    

        for vert in range(numNodes):
            if vert not in coalition:
                if not Counted[vert]:
                    BanzExact[vert] += 1
                if adj.__contains__( vert ): 
                    for u in adj[vert]:
                        if not Counted[u] and Weights[u] + adj[vert][u] >= Wcutoff[u]:
                            BanzExact[vert] += 1
        numCoalitions += 1
                       
    for i in range(numNodes):
        BanzExact[i] = BanzExact[i] * 2/ float(numCoalitions)                     
    return  BanzExact
                     









