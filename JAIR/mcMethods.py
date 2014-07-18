import random
import numpy as np


def game1MCShapstep(numNodes, adj, shuffling, shapMC):
    Counted =  [False]*numNodes
    for i in range(numNodes):
        vert = shuffling[i]
        for u in adj[vert]:
            if not Counted[u]:
                shapMC[vert] += 1
                Counted[u] = True
        if not Counted[vert]:
            shapMC[vert] += 1
            Counted[vert] = True
#    print shapMC        
    return shapMC



def game1MCBanzstep(numNodes, adj, coalition,  banzMC):
        
    Nodes = [i for i in range(numNodes)]
    Counted=[False]*numNodes
    for vert in coalition:
        Counted[vert] = True
        if adj.__contains__(vert):
            for u in adj[vert]:
                Counted[u] = True

    nonCoalition = set(Nodes) - set(coalition)
    for vert in nonCoalition:
          if not Counted[vert]:
              banzMC[vert] += 1
          if adj.__contains__(vert):
              for u in adj[vert]:
                  if not Counted[u]:
                      banzMC[vert] += 1
                                
                  
    return banzMC


        


def game1MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game1MCShapstep(numNodes, adj, copyOfNodes, indexMC)
    return shapMC


def game1MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game1MCBanzstep(numNodes, adj, copyOfNodes, indexMC)
    return banzMC


def game2MCShapstep (numNodes, adj, shuffling, numThresholdEdges, shapMC):
    Counted = [0]*False
    Edges = [0]*False

    for i in range(numNodes):
        vert = shuffling[i]
        for u in adj[i]:
            if not Counted[u]:
                Edges[u] += 1
                if Edges[u] >= numThresholdEdges[u]:
                    shapMC[vert] += 1
                    Counted[u] = True
        if not Counted[vert]:
             shapMC[vert] += 1
             Counted[vert] = True
        shapMC[vert] += 1
    return shapMC         


def game2MCBanzstep (numNodes,  adj, coalition , numThresholdEdges,  banzMC):
    Nodes = [i for i in range(numNodes)]
    Edges = [0]*numNodes
    Counted=[False]*numNodes
    for vert in coalition:
        Counted[vert] = True
        if adj.__contains__(vert):
            for u in adj[vert]:
                Edges[u] += 1
                if Edges[u] >= numThresholdEdges[u]:    
                    Counted[u] = True


    nonCoalition = set(Nodes) - set(coalition)
    for vert in nonCoalition:
          if not Counted[vert]:
              banzMC[vert] += 1
          if adj.__contains__(vert):
              for u in adj[vert]:
                  if not Counted[u]:
                     if  Edges[u] + 1 == numThresholdEdges[u] :
                         banzMC[vert] += 1
                                
    return banzMC
            

def game2MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game2MCShapstep(numNodes, adj, copyOfNodes, numThresholdEdges, indexMC)
    return shapMC


def game2MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game2MCBanzstep(numNodes, adj, copyOfNodes, numThresholdEdges, indexMC)
    return banzMC




def game3MCShapstep(numNodes, adj, shuffling, influence_dist):
    Counted = [0] * False
    Edges = [0] * False
    for i in range(numNodes):
        vert = shuffling[i]
        for u in D[vert]:
            if not Counted[u]:
                shapMC[vert] += 1
                Counted[u] = True
                
        if not Counted[vert]:
            shapMC[vert] += 1
            Counted[vert] = True
            

def game3MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game3MCShapstep(numNodes, adj, copyOfNodes, influence_dist, indexMC)
    return shapMC


def game3MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game3MCBanzstep(numNodes, adj, copyOfNodes, influence_dist, indexMC)
    return banzMC



def game4MCShapstep(numNodes, adj, shuffling, payOff_function, distances):
    Dist =  [1000]*numNodes
    for i in range(numNodes):
        vert = shuffling[i]
        for j in range(numNodes):
            if distance[vert][j] < Dist[j]:
                shapMC[vert] +=( payOff_function(distances[vert][j] ) - payOff_function(Dist[j]) )
                Dist[j] = distances[vert][j]

def game4MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game4MCShapstep(numNodes, adj, copyOfNodes, payOff_function, distances, indexMC)
    return shapMC


def game4MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game4MCBanzstep(numNodes, adj, copyOfNodes, payOff_function, distances, indexMC)
    return banzMC




def game5MCBanzstep (numNodes,  adj, coalition , Wcutoff,  banzMC):
    Nodes = [i for i in range(numNodes)] 
    Counted = [False]*numNodes
    Weights = [0]*numNodes
    for vert in coalition:
        Counted[vert] = True
        if adj.__contains__(vert):
            for u in adj[vert]:
                Weights[u] += adj[u][vert]
                if Weights[u] >= Wcutoff[u]:
                    Counted[u] = True
                    
    nonCoalition =  set(Nodes) - set(coalition) 
    for vert in nonCoalition :
            if not Counted[vert]:
                banzMC[vert] += 1
            if adj.__contains__( vert ): 
                for u in adj[vert]:
                    if not Counted[u] and Weights[u] + adj[vert][u] >= Wcutoff[u]:
                       banzMC[vert] += 1

    return banzMC



def game5MCShapstep (numNodes,  adj, shuffling , Wcutoff, shapMC):
   Weights = [0] * numNodes
   Counted = [0] * numNodes
   for i in range(numNodes):
        vert = shuffling[i]
        newV = 0
        if adj.__contains__(vert):
            for u in adj[vert]:
                weight = adj[vert][u]
                if not Counted[u]:
                    Weights[u]+= weight
                    if(Weights[u] >= Wcutoff[u]):
                        newV = newV +1
                        Counted[u] = True
        if not Counted[vert]:
            newV = newV+ 1
            Counted[vert] = True
        shapMC[vert] += newV
   return shapMC





def game5MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game5MCShapstep(numNodes, adj, copyOfNodes, Wcutoff, indexMC)
    return shapMC


def game5MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game5MCBanzstep(numNodes, adj, copyOfNodes, Wcutoff, indexMC)
    return banzMC


def computeShapMC(numNodes, adj, maxIteration, stepfunc, numThresholdEdges, Wcutoff, influence_dist,  payOff_function, distances):


    shapMC = [0]* numNodes    
    myIter = 1
    shapMCForIter = np.zeros( shape=(maxIteration/5, numNodes))
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        random.shuffle(copyOfNodes)

        shapMC =  stepfunc(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, shapMC)
    	myIter = myIter + 1

        if myIter % 5 is 0:
                shapMCForIter[myIter/5 - 1] = np.array(shapMC) / float(myIter)

              
    for i in range(numNodes):    
        shapMC[i] = shapMC[i] / float(maxIteration)
    return shapMC, shapMCForIter 



def computeBanzMC(numNodes, adj, maxIteration, stepfunc,  numThresholdEdges, Wcutoff, influence_dist,  payOff_function, distances):

    banzMC = [0]* numNodes    
    myIter = 1
    banzMCForIter = np.zeros( shape=(maxIteration/5, numNodes))
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        coalition = []
        for j in copyOfNodes:
            if (random.randint(0,1)):
                 coalition.append(j)
        banzMC =  stepfunc(numNodes, adj, coalition, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, banzMC)
        myIter = myIter + 1
        if myIter % 5 is 0:
                banzMCForIter[myIter/5 - 1] = np.array(banzMC) * 2 / float(myIter)

              
    for i in range(numNodes):    
        banzMC[i] = banzMC[i] * 2 / float(maxIteration)
    return banzMC, banzMCForIter 




'''
/**
 *
 * n     - the size of the graph
 * adj   - the adjencent matrix of the graph
 * maxIteretion - the maximum number of MC iteration
 * banzExact - the exact SV to compute error
 * k - parameter for Game 2
 * D - the sequence of vertices limited by d_cutoff parameter from Game 3
 * f - function from Game 4
 * distances - distances between each pair of vertices computed with Dijkstra algorithm
 *
 */
'''

def computeMCGame1(numNodes, adj, maxIteration, powerIndex ):
    if (powerIndex == 'shap'):
        return(computeShapMC(numNodes, adj, maxIteration, game1MCShapstepAdapter,  0, 0, 0, 0, 0))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, game1MCBanzstepAdapter, 0, 0, 0 , 0, 0))

        



def computeMCGame2(numNodes, adj, maxIteration, numThresholdEdges,  powerIndex ):
    if (powerIndex == 'shap'):
        return(computeShapMC(numNodes, adj, maxIteration, game2MCShapstepAdapter, numThresholdEdges, 0,  0, 0, 0))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, game2MCBanzstepAdapter, numThresholdEdges, 0,  0 , 0, 0))




def computeMCGame3(numNodes, adj, maxIteration, Dcutoff, powerIndex ):
    if (powerIndex == 'shap'):
        return(computeShapMC(numNodes, adj, maxIteration, game3MCShapstepAdapter,  0, 0, Dcutoff, 0, 0))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, game3MCBanzstepAdapter, 0, 0, Dcutoff, 0 ,  0))





def computeMCGame4(numNodes, adj, maxIteration, f, powerIndex ):
    if (powerIndex == 'shap'):
        return(computeShapMC(numNodes, adj, maxIteration, game4MCShapstepAdapter,  0,  0, 0, f, distances))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, game4MCBanzstepAdapter, 0,  0 , 0, f, distances))




def computeMCGame5(numNodes, adj, maxIteration, Wcutoff, powerIndex ):
    if (powerIndex == 'shap'):
        return(computeShapMC(numNodes, adj, maxIteration, game5MCShapstepAdapter,  0,  Wcutoff, 0, 0, 0))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, game5MCBanzstepAdapter, 0,  Wcutoff , 0, 0, 0 ))


















