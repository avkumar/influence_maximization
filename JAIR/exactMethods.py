import math
import itertools 
from itertools import permutations
from itertools import chain, combinations

from operator import mul    # or mul=lambda x,y:x*y
from fractions import Fraction

def nCk(n,k): 
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )



print "in exactMethods.py"
def computeBanzExactGame1(numNodes, adj):
    banzExact = []
    print "in computeExactgame1"
    print numNodes, adj[0], adj[1]
    for i in range(0, numNodes):
        banzExact.append(1.0 / pow(2.0,  len(adj[i]) ) )
        if i == 1:
            print banzExact[1]
        for u in adj[i]:
            banzExact[i] += ( 1.0 / pow( 2.0, len( adj[u] ))  )
            if i == 1:
                print (1.0 / pow(2, len(adj[u]))) 
    print len(banzExact)    
    return banzExact

print "in exactMethods.py"
def computeShapExactGame1(numNodes, adj):
    shapExact = [0]*numNodes
    print "in computeExactgame1"
    print numNodes, adj[0], adj[1]
    for i in range(0, numNodes):
        if adj.__contains__(i):
            shapExact[i] += (1.0 / ( 1 +  len(adj[i]))  )
            for u in adj[i]:
                if adj[i].__contains__(u):
                    shapExact[i] += ( 1.0 / (1 +  len( adj[u] ))  )
    print len(shapExact)    
    return shapExact




def computeShapExactGame2(numNodes, adj, numThreshold):
    shapExact = [0]*numNodes
    for i in range(numNodes):
        shapExact[i] = min(1.0,  (( float) (numThreshold[i])) /(1.0+ (float) (len(adj[i]))) )
        for u in adj[i]:
            degree = float (len(adj[u]))
            shapExact[i] = max(0.0, ( degree - (float) (numThreshold[u]) + 1.0)/	( degree*(1.0+degree))  )
    return shapExact
    					 

def computeShapExactGame2(numNodes, adj, numThreshold):
    banzExact = [0]*numNodes
    for vert in range(numNodes):
        if adj.__contains__(vert):
            for r in range(numThreshold[vert]):
                banzExact[vert] +=     nCk ( len(adj[vert]), r) * (1.0 / pow(2.0,  len(adj[vert]) )  )
        for u in adj[vert]:
            banzExact[vert] +=     nCk( len(adj[u]) - 1, numThreshold[u] - 1) * (1.0 / pow(2.0,  len(adj[u]) )  )
    return banzExact
 

'''
void computeExactGame3(int n, vector< pair<short,double> >* adj, vector<short> *D){

	for (unsigned short i=0; i<n; i++) {
		 for (unsigned short j=0; j<D[i].size(); j++) {
			 double extDegree = D[D[i][j]].size() - 1;
			 shapExact[i] += (1.0 / (1.0 + extDegree));
		 }
	}
}

'''
def computeShapExactGame3(numNodes, adj, D):
    for i in range(numNodes):
        for u in D[i]:
            extDegree = len(D[u]) - 1
            shapExact +=  (1.0 / (1.0 + extDegree)) 


'''


void computeExactGame4(int n, vector< pair<short,double> >* adj, double (*f)(double), vector<double> *distances, vector<short> *D) {

	for (short i=1; i<=n; i++)
		shapExact[i] = 0;

	for (short i=1; i<=n; i++) {

		long double  sum=0.0;

		double currSV = 0.0;
		short index = n-1;
		
		double prevDistance = -1.0;
		double prevSV = -1.0;

		double wart= 0.0;
		while (index > 0) {
			wart = ( (*f)(distances[i][D[i][index]])/(1.0+ ((double)index) ) );
				if (distances[i][D[i][index]] == prevDistance) {
					currSV = prevSV;
				} else {			
					currSV = ( wart - sum);
				}

	 			shapExact[D[i][index]] += currSV;
	 			sum += ( wart/((double)index));
				
				prevDistance = distances[i][D[i][index]];
				prevSV = currSV;

	 			index--;
	 		}
		shapExact[i] +=  ((*f)(0.0) - sum);
	}

}


double alfa[MAX_NODES];
double beta[MAX_NODES];

void computeApproximateGame5(int n, vector< pair<short,double> >* adj, double *Wcutoff) {

	double weight;
	short u  = 0;

	for (unsigned short i=0; i<n; i++) {
		 alfa[i] = 0;
		 beta[i] = 0;
		 for (unsigned short j=0; j<adj[i].size(); j++) {
			 weight = adj[i][j].second;
			 alfa[i] += weight;
			 beta[i] += (weight*weight);
		 }
	}


	for (unsigned short i=0; i<n; i++) {
		shapApprox[i] = 0;
		int deg = adj[i].size();
		for (unsigned short m=0; m<=deg; m++){
			double my = (m*alfa[i])/(deg);
			double sigma = ( m*(deg-m)*(beta[i]-( (alfa[i]*alfa[i])/deg )) ) / (deg*(deg - 1));
			double p = (1+ erf( (Wcutoff[i]-my)/(sqrt(2.0*sigma))))/2.0;
			shapApprox[i] += p/(1+deg);
		}

		for (unsigned short j=0; j<adj[i].size(); j++) {
			double p = 0;
			u = adj[i][j].first;
			weight = adj[i][j].second;
			deg = adj[u].size();
			int deg_bis = deg - 1;
			double alfa_bis = alfa[u] - weight;
			double beta_bis = beta[u] - (weight*weight);

			for (unsigned short m=0; m<=deg_bis; m++){
				double my = (m*alfa_bis)/(deg_bis);
				double sigma = ( m*(deg_bis-m)*(beta_bis-( (alfa_bis*alfa_bis)/deg_bis )) ) / (deg_bis*(deg_bis - 1));
				double z = (erf( (Wcutoff[u]-my)/sqrt(2.0*sigma) ) - erf( (Wcutoff[u]-weight-my)/sqrt(2.0*sigma) ) )/2.0;
				p += ((deg-m)*z)/(deg*(deg+1));
			}
			shapApprox[i] += p;
		}

	}
}

short Nodes[MAX_NODES];
/* sum of weights from coalition to each vertex for game 5 */
double Weights2[MAX_NODES];
bool Counted2[MAX_NODES];

long factorials[] = {0,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600};

void computeBruteForceGame5(int n, vector< pair<short,double> >* adj, double *Wcutoff) {
/*
	double weight;
	short u  = 0;

	for (short i=0; i<n; i++){
		Nodes[i]=i;
	}
	do {
		for (short i=0; i<n; i++) {
			Counted2[i] = false;
			Weights2[i] = 0;
		}

		 for (short i=0; i<n; i++) {
			short vert = Nodes[i];
			double newV = 0;
			for (unsigned short j=0; j<adj[vert].size(); j++) {
				u = adj[vert][j].first;
				weight = adj[vert][j].second;
				if ( !Counted2[u] ) {
					Weights2[u]+= weight;
					if(Weights2[u] >= Wcutoff[u]) {
						newV++;
						Counted2[u] = true;
					}
				}
			}
			if (!Counted2[vert]){
				newV++;
				Counted2[vert] = true;
			}
			shapExact[vert] += newV;
		}

	} while ( next_permutation (Nodes,Nodes+n) );

	for (short i=0; i<n; i++){
		shapExact[i] /= factorials[n];
	}
*/
}

'''
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
                     









