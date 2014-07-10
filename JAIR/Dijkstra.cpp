/*
 *
 * In this file we implement Dijkstra algorithm
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
#include <set>

#include "global.h"
#include "Dijkstra.h"

using namespace std;

 vector<double> distances[MAX_NODES];

int SOURCE=0;

struct cmp
{
    // answer if distance to "a" from SOURCE is lower than distance to "b"
    bool operator() (const short &a, const short &b)
    {
        if (distances[SOURCE][a] < distances[SOURCE][b]) return true;
        if (distances[SOURCE][a] > distances[SOURCE][b]) return false;
        return a<b;
    }
};

set<short,cmp> heap;
vector<short> D[MAX_NODES];



void dijkstra(int n, vector< pair<short,double> >* adj, short source, double maxDist)
{
	short v, u;
	double c;
	SOURCE = source;

	distances[SOURCE][SOURCE] = 0.0;
	heap.clear();

    for (short i=0; i<n; i++) heap.insert(i);


    while( !heap.empty() )
    {
        u = *(heap.begin());
        heap.erase(heap.begin());

        /* if we reach max distance we skip further exploration */
        if (maxDist != -1 && distances[SOURCE][u] > maxDist)
        	break;

        if (distances[SOURCE][u] != INFI) {
			for (unsigned short i=0; i<adj[u].size(); i++)
			{
				v = adj[u][i].first;
				c = adj[u][i].second;
				if ( distances[SOURCE][u] + c < distances[SOURCE][v]){
					heap.erase(v);
					distances[SOURCE][v] = distances[SOURCE][u] + c ;
					heap.insert(v);
				}
			}
        }

		D[SOURCE].push_back(u);
    }
}
/**
 * This function compute distance between all pairs of nodes
 */
void dijkstraDistances(int n, vector< pair<short,double> >* adj, double maxDist){
    for (int i=0; i<n; i++) {
    	dijkstra(n, adj, i, maxDist);
    }
}


void allocateDijkstraSpace(int n){
    for (short i=0; i<n; i++) {
    	distances[i].resize(n, INFI);
    	D[i].reserve(n);

    }
}
