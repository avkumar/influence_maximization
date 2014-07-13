#! /usr/bin/env python
import time
import exactMethods
import mcMethods
from exactMethods import computeShapExactGame1
from mcMethods import computeMCGame1
from exactMethods import computeBanzExactGame1
from mcMethods import computeMCGame5
from exactMethods import computeBanzBruteForceGame5
from exactMethods import computeShapBruteForceGame5
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


def readRandomGraph(size):
    adj = {}
    for i in range(size):
        adj[i] = {}
        for j in range(i + 1, size):
            if not adj.__contains__(j):
                adj[j] = {}
            weight = random.randint(0,5)
            adj[i][j] = weight 
            adj[j][i] = weight
            
    return adj       
            
                        

def countAvgDist(dist ):
    return 1 


def getWcutoff (numNodes):
    Wcutoff = [1]*numNodes
    return Wcutoff


if __name__ =="__main__":
    start_time = time.time()
    maxIteration = 0
    maxIteration = 1200
    G1 = 'power.txt'
    G2 = 'astro-ph.txt'
    G3 = '4node.txt'
    numNodes, numEdges,  adj =  readGraph(G1, True)
#    print numNodes, numEdges, adj
#    shapExactGame1 = computeShapExactGame1(numNodes, adj) 
#    print "shapExactGame1 is ", shapExactGame1 
#    banzExactGame1 = computeBanzExactGame1(numNodes, adj)
#    print "banzExactGame1 is", banzExactGame1
#    shapMCGame1 = computeMCGame1(numNodes, adj, maxIteration, 'shap')
#    print "shapMCGame1 is", shapMCGame1
#    banzMCGame1 = computeMCGame1(numNodes, adj, maxIteration, 'banz')
#    print "banMCGame1 is", banzMCGame1
    Wcutoff = getWcutoff(numNodes)
#    shapMCGame5 = computeMCGame5(numNodes, adj, 12, Wcutoff, 'shap')
#    print "shapMCGame5 is", shapMCGame5
#    banzMCGame5 = computeMCGame5(numNodes, adj, 3, Wcutoff, 'banz')
#    print "banzMCGame5 is", banzMCGame5
#    numThreshhold = [2]*numNodes
#
#    shapExactGame2 = computeShapExactGame2(numNodes, adj, numThreshhold) 
#    print "shapExactGame2 is ", shapExactGame2 
#    banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshhold)
#    print "banzExactGame2 is", banzExactGame2
#    shapMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'shap')
#    print "shapMCGame2 is", shapMCGame2
#    banzMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'banz')
#
#    shapExactBruteGame5 =  computeShapBruteForceGame5(numNodes, adj, Wcutoff)
#    print shapExactBruteGame5
#    banzExactBrute =  computeBanzBruteForceGame5(numNodes, adj, Wcutoff)
#    banzExactBruteGame5 = computeBanzBruteForceGame5(numNodes, adj, Wcutoff)
#    print banzExactBruteGame5
    shapMCGame5 = computeMCGame5(numNodes, adj, maxIteration, Wcutoff, 'shap')
    print shapMCGame5
    banzMCGame5 = computeMCGame5(numNodes, adj, maxIteration, Wcutoff , 'banz')
    print banzMCGame5

    print("--- %s seconds ---" % time.time() - start_time)
#
