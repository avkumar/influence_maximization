/*
 *
 * The header of the file where we implement five Monte Carlo methods for computing Shapley value in games
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

#ifndef MCMETHODS_H_
#define MCMETHODS_H_

#include <vector>
#include <cstdlib>
using namespace std;

void computeMCGame1(int n, vector< pair<short,double> >* adj, int maxIteration,
					double* shapExact);

void computeMCGame2(int n, vector< pair<short,double> >* adj, int maxIteration,
					int *k, double* shapExact);

void computeMCGame3(int n, vector< pair<short,double> >* adj, int maxIteration,
					vector<short> *D, double* shapExact);

void computeMCGame4(int n, vector< pair<short,double> >* adj, int maxIteration,
				 	double (*f)(double), vector<double>* distances, double* shapExact);

void computeMCGame5(int n, vector< pair<short,double> >* adj, int maxIteration,
					double *Wcutoff, double* shapExact);

/* The approximation of Shapley value */
extern double shapMC[];



#endif /* MCMETHODS_H_ */
