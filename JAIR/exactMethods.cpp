/**
 *
 * The implementation of five games from paper "Efficient Computation of the Shapley Value for
 * Game-Theoretic Network Centrality" and also brute force method for computing Shapley value
 *
 * THIS SOURCE CODE IS SUPPLIED "AS IS" WITHOUT WAR-
 * RANTY OF ANY KIND, AND ITS AUTHOR AND THE JOURNAL OF
 * ARTIFICIAL INTELLIGENCE RESEARCH (JAIR) AND JAIR'S PUB-
 * LISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES
 * INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PUR-
 * POSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER
 * ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS
 * SOURCE CODE, AND NEITHER THE AUTHOR NOR JAIR, NOR JAIR'S
 * PUBLISHER AND DISTRIBUTORS, WILL BE LIABLE FOR DAM-
 * AGES OF ANY KIND RESULTING FROM ITS USE. Without limiting
 * the generality of the foregoing, neither the author, nor JAIR, nor JAIR's
 * publisher and distributors, warrant that the Source Code will be error-free,
 * will operate without interruption, or will meet the needs of the user.
 *
 *
 */

#include <vector>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "global.h"

using namespace std;


double shapExact[MAX_NODES];
double shapApprox[MAX_NODES];

void computeExactGame1(int n, vector< pair<short,double> >* adj) {
	short u = 0;
	for (unsigned short i=1; i<=n; i++) {
		shapExact[i] = (1.0 /  (1 +  ( adj[i].size())) );

		 for (unsigned short j=0; j<adj[i].size(); j++) {
			 u = adj[i][j].first;
			 shapExact[i] += (1.0 / (1.0 + adj[u].size()));
		 }
	}
}

void computeExactGame2(int n, vector< pair<short,double> >* adj, int *k){
	short u = 0;

	for (unsigned short i=1; i<=n; i++) {
		shapExact[i] = min(1.0,  ((double)k[i])/(1.0+ (double)adj[i].size()) );

		 for (unsigned short j=0; j<adj[i].size(); j++) {
			 u = adj[i][j].first;
			 double degree = (double)adj[u].size();

			 shapExact[i] += max(0.0, ( degree - (double)k[u] + 1.0)/
					 	( degree*(1.0+degree))  );
		 }
	 }
}

void computeExactGame3(int n, vector< pair<short,double> >* adj, vector<short> *D){

	for (unsigned short i=0; i<n; i++) {
		 for (unsigned short j=0; j<D[i].size(); j++) {
			 double extDegree = D[D[i][j]].size() - 1;
			 shapExact[i] += (1.0 / (1.0 + extDegree));
		 }
	}
}

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



