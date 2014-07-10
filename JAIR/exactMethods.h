/**
 *
 * The header file of implementation of five games from paper "Efficient Computation of the Shapley Value for
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

#ifndef EXACTMETHODS_H_
#define EXACTMETHODS_H_

#include <vector>
#include <cstdlib>
using namespace std;

void computeExactGame1(int n, vector< pair<short,double> >* adj);
void computeExactGame2(int n, vector< pair<short,double> >* adj, int *k);
void computeExactGame3(int n, vector< pair<short,double> >* adj, vector<short> *D);
void computeExactGame4(int n, vector< pair<short,double> >* adj, double (*f)(double), vector<double> *distances, vector<short> *D);



void computeApproximateGame5(int n, vector< pair<short,double> >* adj, double *Wcutoff);
void computeBruteForceGame5(int n, vector< pair<short,double> >* adj, double *Wcutoff);

/* the exact SV  */
extern double shapExact[];

/* the approx SV  */
extern double shapApprox[];


#endif /* EXACTMETHODS_H_ */
