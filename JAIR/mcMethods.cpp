/*
 *
 * In this file we implement five Monte Carlo methods for computing Shapley value in games
 * defined in paper "Efficient Computation of the Shapley Value for Game-Theoretic Network Centrality"
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
//#include <windows.h>

#include "global.h"

using namespace std;

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


/**
 *
 * n     - the size of the graph
 * adj   - the adjencent matrix of the graph
 * maxIteretion - the maximum number of MC iteration
 * shapExact - the exact SV to compute error
 * k - parameter for Game 2
 * D - the sequence of vertices limited by d_cutoff parameter from Game 3
 * f - function from Game 4
 * distances - distances between each pair of vertices computed with Dijkstra algorithm
 *
 */

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





