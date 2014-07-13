
/**
 *
 * The example of usage our methods
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

#include <iostream>
#include <fstream>
//#include <windows.h>

#include "global.h"
#include "Dijkstra.h"
#include "mcMethods.h"
#include "exactMethods.h"

#include <cmath>

using namespace std;


/* number of nodes */
short n;
/* number of edges */
int m;
/* Graph */
vector< pair<short,double> > adj[MAX_NODES];




void readGraph(char *name, bool weighted){

	ifstream infile;
	infile.open (name);
    short node1,node2;
    double weight ;
    infile >> n >> m;
    cout<<n<<" "<<m<<endl;
    for (int i=0; i<m; i++)
    {
    	weight = 1.0;
    	if (weighted) {
    		infile >> node1 >> node2 >> weight;
    	} else {
    		infile >> node1 >> node2;
    	}
        adj[node1].push_back( pair<short,double>(node2,weight) );
        adj[node2].push_back( pair<short,double>(node1,weight) );
    }
}


void readRandomGraph(int size){
    double weight ;
    n = size;
    for (int i=0; i<n; i++)
    	adj[i].clear();

    for (int i=0; i<n; i++)
    {
    	for (int j=0; j<n; j++) {
    		weight = (double)rand()/(double)RAND_MAX;
            adj[i].push_back( pair<short,double>(j,weight) );
            adj[j].push_back( pair<short,double>(i,weight) );
    	}
    }
}


double countAvgDist(){
    double avgDis=0;
    int count = 0;
    for (int i=0; i<n; i++) {
    	 for (int j=0; j<n; j++) {
    		 if (distances[i][j] != INFI) {
    			 avgDis +=distances[i][j];
    			 count++;
    		 }
    	 }
    }
    avgDis /= (double)count;
    return avgDis;
}



int main(void)
{

    srand(time(0));
    int maxIteration = 0;
    cout.precision(15);

    ofstream myfile ("example.txt");
	maxIteration = 3200;
	readGraph(G3,true);
    computeExactGame1(n, adj);
//    computeMCGame1(n, adj, maxIteration, shapExact);
    double Wcutoff[n+1];
    for (short i = 1; i < n + 1; i ++){
        Wcutoff[i] = 1;
    }    
    computeMCGame5(n, adj, maxIteration, Wcutoff, shapMC);
    if (myfile.is_open())
    {
      myfile << "This is a line.\n";
   	for (short i=1; i<n+1; i++) {
		myfile<<i<<" "<<shapExact[i]<<" "<<shapMC[i]<<endl;
	}

     myfile << "This is another line.\n";
      myfile.close();
    }


    return 0;
}

