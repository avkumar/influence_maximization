/*
 *
 * Dijkstra algorithm
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

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include <vector>
#include <cstdlib>
using namespace std;


/* Distances counted by Dijkstra algorithm */
extern vector<double> distances[];

/* Collection of nodes ordered by the distance to source */
extern vector<short> D[];

/* Compute all distances in graph */
void dijkstraDistances(int n, vector< pair<short,double> >* adj, double maxDist = -1);

/* Compute all distances in graph from given source*/
void dijkstra(int n, vector< pair<short,double> >* adj, short source, double maxDist = -1);

/* Allocate space required for storing all distances O(n^2) */
void allocateDijkstraSpace(int n);

#endif /* DIJKSTRA_H_ */
