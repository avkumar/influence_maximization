#! /usr/bin/env python

import exactMethods
import mcMethods
from exactMethods import computeExactGame1
from mcMethods import computeMCGame1
'''

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
'''


def readGraph(filename, isWeighted) :
    weight = 1.0
    adj= {}
    with open(filename) as f:
        numNodes, numEdges = f.readline().split(',')
        numNodes , numEdges = int (numNodes), int (numEdges)
        for i in range(0,int (numEdges)):
            if (isWeighted):
                node1, node2, weight = f.readline().split(',')
            else :
                node1, node2 = f.readline().split(',')
            node1, node2, weight = int (node1), int (node2), float (weight)    

            if not adj.__contains__(node1):
                adj[node1] = {}
            adj[node1][node2] = float(weight)
            if not adj.__contains__(node2):
                adj[node2] = {}
            adj[node2][node1] = float(weight)

#    print adj        
    return (numNodes, numEdges, adj)

                         
'''

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
'''
'''
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
'''

'''
int main(void)
{

    srand(time(0));
    int maxIteration = 0;
    cout.precision(15);

    ofstream myfile ("example.txt");
	maxIteration = 3200;
	readGraph(G1,true);
    computeExactGame1(n, adj);
    computeMCGame1(n, adj, maxIteration, shapExact);
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

'''
if __name__ =="__main__":
    maxIteration = 0
    maxIteration = 12
    G1 = 'power.txt'
    G2 = 'astro-ph.txt'
    numNodes, numEdges,  adj =  readGraph(G1, True)
#    print numNodes, numEdges, adj
    shapExact1 = computeExactGame1(numNodes, adj)
    shapComputeMCGame1 = computeMCGame1(numNodes, adj, maxIteration)
    print shapExact1
    print shapComputeMCGame1
