import random
'''

short shuffling[MAX_NODES];
/* in MC methods (games 1-3) we count which vertices have already contributed
 * to coalition */
bool Counted[MAX_NODES];
/* in game 2 we count how many edges for each vertex is connected with coalition */
int Edges[MAX_NODES];
/* distance from coalition to each vertex for game 4 */
double Dist[MAX_NODES];

/* sum of weights from coalition to each vertex for game 5 */
double Weights[MAX_NODES];

double shapMC[MAX_NODES];

/* Time stamps */
//LARGE_INTEGER start, stop, proc_freq;

/* step in game 1 */

void game1MCstep(int n, vector< pair<short,double> > *adj, short* shuffling){
	int newV = 0;
	short u = 0;

	for (short i=0; i<n; i++)
		Counted[i] = false;

	 for (short i=0; i<n; i++) {
		short vert = shuffling[i];
		newV = 0;
		for (unsigned short j=0; j<adj[vert].size(); j++) {
			u = adj[vert][j].first;
			if ( !Counted[u] ) {
				newV++;
				Counted[u] = true;
			}
		}
		if (!Counted[vert]){
			newV++;
			Counted[vert] = true;
		}
		shapMC[vert] += newV;
	}
}



void game1MCstepAdapter(int n, vector< pair<short,double> > *adj, short * shuffling,
					    int *k, double *Wcutoff, vector<short> *D,double (*f)(double), vector<double>* waga){
	game1MCstep(n, adj, shuffling);
}
'''

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
    Counted = [False]*numNodes
    for vert in range(numNodes):
        for otherVert in coalition:
            if otherVert != vert:
                for u in adj[otherVert]:
                    Counted[u] = True
                Counted[otherVert] = True
        for u in adj[vert]:
            if not Counted[u]:
                banzMC[vert] += 1
                Counted[u] = True
        if not Counted[vert]:
            banzMC[vert] += 1
            Counted[vert] = True
#    print shapMC        
    return banzMC



def game1MCShapstepAdapter(numNodes, adj, copyOfNodes, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
    shapMC = game1MCShapstep(numNodes, adj, copyOfNodes, indexMC)
    return shapMC



def game1MCBanzstepAdapter(numNodes, adj, copyOfNodes, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances, indexMC):
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


void game3MCstep(int n, vector< pair<short,double> >* adj, short* shuffling , vector<short> *D){
	int newV = 0;
	short u = 0;

	for (short i=0; i<n; i++){
		Counted[i] = false;
		Edges[i] = 0;
	}

	 for (short i=0; i<n; i++) {
		short vert = shuffling[i];
		newV = 0;
		for (unsigned short j=0; j<D[vert].size(); j++) {
			u = D[vert][j];
			if ( !Counted[u] ) {
				newV++;
				Counted[u] = true;

			}
		}
		if (!Counted[vert]){
			newV++;
			Counted[vert] = true;
		}
		shapMC[vert] += newV;
	}
}

void game3MCstepAdapter(int n, vector< pair<short,double> > *adj, short * shuffling,
					    int *k, double *Wcutoff, vector<short> *D, double (*f)(double), vector<double>* waga){
	game3MCstep(n, adj, shuffling, D);
}


void game4MCstep(int n, vector< pair<short,double> >* adj, short * shuffling,
				 double (*f)(double), vector<double>* distances){

	for (short i=0; i<n; i++)
		Dist[i] = INFI;

	 for (short i=0; i<n; i++) {
		 short vert = shuffling[i];

		 for (unsigned short j=0; j<n; j++){
			if ( distances[vert][j] < Dist[j] ) {
				shapMC[vert] += ((*f)(distances[vert][j]) - (*f)(Dist[j]));
					Dist[j] = distances[vert][j];
				}
		 }
	}
}

void game4MCstepAdapter(int n, vector< pair<short,double> > *adj, short * shuffling,
					    int *k, double *Wcutoff, vector<short> *D, double (*f)(double), vector<double>* waga){
	game4MCstep(n, adj, shuffling, (*f), waga);
}

void game5MCstep(int n, vector< pair<short,double> >* adj, short* shuffling , double *Wcutoff){
	int newV = 0;
	short u = 0;
	double weight = 0;

	for (short i=0; i<n; i++){
		Counted[i] = false;
		Weights[i] = 0;
	}

	 for (short i=0; i<n; i++) {
		short vert = shuffling[i];
		newV = 0;
		for (unsigned short j=0; j<adj[vert].size(); j++) {
			u = adj[vert][j].first;
			weight = adj[vert][j].second;
			if ( !Counted[u] ) {
				Weights[u]+= weight;
				if(Weights[u] >= Wcutoff[u]) {
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

void game5MCstepAdapter(int n, vector< pair<short,double> > *adj, short * shuffling,
					    int *k, double *Wcutoff, vector<short> *D, double (*f)(double), vector<double>* waga){
	game5MCstep(n, adj, shuffling, Wcutoff);
}

/*	This function randomly shuffle table. Each permutation is equally probable  */
void shuffle(short* tab, int n) {
	for (short i=0; i<n; i++) {
		short r = i + (rand() % (n-i)); // Random remaining position.
		short temp = tab[i]; tab[i] = tab[r]; tab[r] = temp;
	}
}


void coalition(short* tab, int n) {
    for (short i = 0; i < n; i ++) {
        if(rand() % 2) {
             tab[i] = -1;
        }      
        else 
             tab[i] = i;
    }         
}



void computeShapMC(int n, vector< pair<short,double> > *adj, int maxIteration, double *shapExact,
			void (*step)(int , vector< pair<short,double> >* ,  short*, int*, double*,  vector<short>*,
					double (*f)(double), vector<double>* ),int *k, double *Wcutoff, vector<short> *D,
					double (*f)(double),vector<double>* distances) {
	 for (short i=0; i<n; i++) {
		 shuffling[i] = i;
		 shapMC[i] = 0;
	 }


	int myIter = 1;
	while (myIter <= maxIteration){
		shuffle(shuffling,n);
	    (*step)(n,adj,shuffling,k,Wcutoff,D,(*f),distances);

		if (myIter % 5 == 0) {
			double precision=1.0;
			for (short i=0; i<n; i++) {
				double shap = (shapMC[i])/((double)(myIter));
				precision = min(precision,min(shapExact[i],shap)/max(shapExact[i],shap)) ;
			}
		}
	myIter++;
	}

	 for (short i=0; i<n; i++)
		 shapMC[i] /= maxIteration;

}
'''
def computeShapMC(numNodes, adj, maxIteration,  numThreshholdEdges, Wcutoff, influence_dist,  payOff_function, distances):


#        game1MCstepAdapter(numNodes, adj, copyOfNodes, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances):
    shapMC = [0]* numNodes    
    myIter = 1
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        random.shuffle(copyOfNodes)

        shapMC =  game1MCShapstepAdapter(numNodes, adj, copyOfNodes, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances, shapMC)
    	myIter = myIter + 1

    print shapMC
    for i in range(numNodes):    
        shapMC[i] = shapMC[i] / float(maxIteration)
    print shapMC     
    return shapMC 



def computeBanzMC(numNodes, adj, maxIteration,  numThreshholdEdges, Wcutoff, influence_dist,  payOff_function, distances):


#        game1MCstepAdapter(numNodes, adj, copyOfNodes, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances):
    banzMC = [0]* numNodes    
    myIter = 1
    while myIter <= maxIteration:
        copyOfNodes = range(numNodes)
        coalition = []
        for j in copyOfNodes:
            if (random.randint(0,1)):
                 coalition.append(j)
        banzMC =  game1MCBanzstepAdapter(numNodes, adj, coalition, numThreshholdEdges, Wcutoff, influence_dist, payOff_function, distances, banzMC)
        myIter = myIter + 1

    print banzMC
    for i in range(numNodes):    
        banzMC[i] = banzMC[i] / float(maxIteration)
    print banzMC     
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
        return(computeShapMC(numNodes, adj, maxIteration,  0, 0, 0, 0, 0))
    if (powerIndex == 'banz'):
        return (computeBanzMC(numNodes, adj, maxIteration, 0, 0, 0 , 0, 0))

        

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


