import random



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
    numCounted = [0]*numNodes
    for vert in coalition:
        numCounted[vert] += 1
        for u in adj[vert]:
            numCounted[u] += 1 
            
    for vert in range(numNodes):
        if vert in coalition:
            if numCounted[vert] <= 1:
                banzMC[vert] += 1
            for u in adj[vert]:
                if numCounted[u] <= 1:
                    banzMC[vert] += 1
        else:
            if numCounted[vert] == 0:
                 banzMC[vert] += 1
            for u in adj[vert]:
                 if numCounted[u] is 0:
                     banzMC[vert] += 1

    return banzMC


        


def game1MCShapstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game1MCShapstep(numNodes, adj, copyOfNodes, indexMC)
    return shapMC


def game1MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    banzMC = game1MCBanzstep(numNodes, adj, copyOfNodes, indexMC)
    return banzMC

'''
void game2MCstep(int n, vector< pair<short,double> >* adj, short* shuffling , int *k){
	int newV = 0;
	short u = 0;

	for (short i=0; i<n; i++){
		Counted[i] = false;
		Edges[i] = 0;
	}

	 for (short i=0; i<n; i++) {
		short vert = shuffling[i];
		newV = 0;
		for (unsigned short j=0; j<adj[vert].size(); j++) {
			u = adj[vert][j].first;
			if ( !Counted[u] ) {
				Edges[u]++;
				if(Edges[u] >= k[u]) {
					newV++;
					Counted[u] = true;
				}
			}
		}
		if (!Counted[vert]){
			newV++;
			Counted[vert] = true;
		}
		shapMC[vert] += newV;
	}
}

void game2MCstepAdapter(int n, vector< pair<short,double> > *adj, short * shuffling,
					    int *k, double *Wcutoff, vector<short> *D, double (*f)(double), vector<double>* waga){
	game2MCstep(n, adj, shuffling, k);
}
'''

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

    def findBanzMC(vert, coalition, Edges):

        banzMC_vert_val = 0
        if adj.__contains__(vert):
            for u in adj[vert]:
                if coalition:
                     if u not in coalition:
                         Edges[u] += 1
                         if Edges[u] == numThresholdEdges[u]:
                             banzMC_vert_val += 1
            if Edges[vert] < numThresholdEdges[vert]:
                banzMC_vert_val += 1

        return banzMC_vert_val

    def computeEdges(coalition):  
        Edges = [0] * numNodes  
        for vert in coalition:
            for u in adj[vert]:
                Edges[u] += 1
        return Edges

    for vert in range(numNodes):
        if vert in coalition:
            coalition_cpy = list(coalition)
            coalition_cpy.remove(vert)
            Edges = computeEdges(coalition_cpy)
            banzMC[vert] = findBanzMC(vert, coalition_cpy, Edges)
        else:
            Edges = computeEdges(coalition)
            banzMC[vert] = findBanzMC(vert, coalition, Edges)
     

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
        for u in adj[vert]:
            Weights[u] += adj[u][vert]
            if Weights[u] >= Wcutoff[u]:
                Counted[u] = True
                

    for vert in range(numNodes):
        if vert not in coalition:
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
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        random.shuffle(copyOfNodes)

        shapMC =  stepfunc(numNodes, adj, copyOfNodes, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, shapMC)
    	myIter = myIter + 1

    for i in range(numNodes):    
        shapMC[i] = shapMC[i] / float(maxIteration)
    return shapMC 



def computeBanzMC(numNodes, adj, maxIteration, stepfunc,  numThresholdEdges, Wcutoff, influence_dist,  payOff_function, distances):

    banzMC = [0]* numNodes    
    myIter = 1
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        coalition = []
        for j in copyOfNodes:
            if (random.randint(0,1)):
                 coalition.append(j)
#        print coalition         
        banzMC =  stepfunc(numNodes, adj, coalition, numThresholdEdges, Wcutoff, influence_dist, payOff_function, distances, banzMC)
        myIter = myIter + 1

    for i in range(numNodes):    
        banzMC[i] = banzMC[i] * 2 / float(maxIteration)
    return banzMC 




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

















'''


void computeMCGame1(int n, vector< pair<short,double> > *adj, int maxIteration, double *shapExact){
	computeShapMC(n, adj, maxIteration, shapExact, game1MCstepAdapter,0,0,0,0,0);
}
void computeMCGame2(int n, vector< pair<short,double> >* adj, int maxIteration,
					int *k, double* shapExact){
	computeShapMC(n, adj, maxIteration, shapExact, game2MCstepAdapter,k,0,0,0,0);
}
void computeMCGame3(int n, vector< pair<short,double> >* adj, int maxIteration,
					vector<short> *D, double* shapExact){
	computeShapMC(n, adj, maxIteration, shapExact, game3MCstepAdapter,0,0,D,0,0);
}
void computeMCGame4(int n, vector< pair<short,double> >* adj, int maxIteration,
				 	double (*f)(double), vector<double>* distances, double* shapExact){
	computeShapMC(n, adj, maxIteration, shapExact, game4MCstepAdapter,0,0,0,(*f),distances);
}
void computeMCGame5(int n, vector< pair<short,double> >* adj, int maxIteration,
					double *Wcutoff, double* shapExact){
	computeShapMC(n, adj, maxIteration, shapExact, game5MCstepAdapter,0,Wcutoff,0,0,0);
}


'''


