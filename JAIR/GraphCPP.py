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
from computeInfluence import computeSpread
#from exactMethods import computeBanzExactGame2
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
    Wcutoff = [0]*numNodes
    return Wcutoff

def rankNodes(powerIndices):
    idx = sorted(xrange(len(powerIndices)), key=powerIndices.__getitem__)
    idx.reverse()
    return idx

    

if __name__ =="__main__":
    start_time = time.time()
    maxIteration = 0
    maxIteration = 120
    G1 = 'power.txt'
    G2 = 'astro-ph.txt'
    G3 = '4node.txt'
    numNodes, numEdges,  adj =  readGraph(G1, True)
#    print numNodes, numEdges, adj
    shapExactGame1 = computeShapExactGame1(numNodes, adj) 
    print "shapExactGame1 is ", shapExactGame1 
    banzExactGame1 = computeBanzExactGame1(numNodes, adj)
    print "banzExactGame1 is", banzExactGame1
    nodesRankedByShap = rankNodes(shapExactGame1)       
    nodesRankedByBanz = rankNodes(banzExactGame1)
#    shapMCGame1 = computeMCGame1(numNodes, adj, maxIteration, 'shap')
#    print "shapMCGame1 is", shapMCGame1
#    banzMCGame1 = computeMCGame1(numNodes, adj, maxIteration, 'banz')
#    print "banMCGame1 is", banzMCGame1
 #    shapExactGame2 = computeShapExactGame2(numNodes, adj, numThreshold) 
#    print "shapExactGame2 is ", shapExactGame2 
  #  numThreshold = [1] * numNodes
   # banzExactGame2 = computeBanzExactGame2(numNodes, adj, numThreshold)
    #print "banzExactGame2 is", banzExactGame2
#    shapMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshold, 'shap')
#    print "shapMCGame2 is", shapMCGame2
#    banzMCGame2 = computeMCGame2(numNodes, adj, maxIteration, numThreshhold, 'banz')
#    print "banMCGame2 is", banzMCGame2
    Wcutoff = getWcutoff(numNodes)
    shapNodes = [2564,335,1027,1232,231,256,913,18,3272,118,2393,37,40,347,1078,2868,1489,1354,2879,5508,1165,469,1836,984,2339,236,5503,80,346,560,409,456,1,2567,2394,4674,4741,4287,1797,2909,5027,333,871,493,1190,4073,516,3143,97,61,1484,2294,211,1558,401,1619,5701,6199,149,832,752,1710,3692,792,327,94,600,1909,306,931,481,556,6846,1476,3507,436,1975,62,4970,1369,357,987,1203,2165,3865,4719,511,711,232,3260,1202,1976,235,265,847,43,1605,748,1137,1533,3180,952,3028,6202]
    banzNodes =  [2564,231,37,3272,913,516,1558,1078,1232,335,2294,256,2879,1484,2868,1836,871,401,556,2394,18,4741,1165,1203,1489,984,5995,2393,1027,1354,409,236,1044,118,3507,3671,3860,347,1,5159,3143,7918,493,1190,2492,1996,560,5508,346,5943,4670,792,456,436,5227,8033,97,2228,61,602,4165,1369,3965,715,5027,152,2669,1202,4970,518,952,4073,3273,2360,2007,987,344,1024,40,1619,3046,2600,1124,497,2165,1932,511,8000,1975,3765,247,5079,861,5878,3457,4674,5558,8657,80,3308,1179,1255,5007,2567]
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

#    shapMCGame5 = computeMCGame5(numNodes, adj, 12000, Wcutoff, 'shap')
#    print shapMCGame5
#    print time.time() - start_time
#    banzMCGame5 = computeMCGame5(numNodes, adj, 12000, Wcutoff , 'banz')
#    print banzMCGame5
#    print len(shapNodes),  len(banzNodes)
    for i in range(20):
        
        print computeSpread(numNodes, adj, nodesRankedByShap[0:5*i], Wcutoff, 30)

        print computeSpread(numNodes, adj, nodesRankedByBanz[0:5*i], Wcutoff, 30)
    print time.time() - start_time
#
