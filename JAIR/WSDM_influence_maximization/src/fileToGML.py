with open ('../data/4node.txt') as f:
    out = open('4node.gml', 'a')
    content = f.read().splitlines()
    numNodes = int (content[0].split(',')[0])
    numEdges = int (content[0].split(',')[1])	
    print >> out,  "graph"
    print >> out,  "[" 
    for i in range(0, numNodes):
	print >> out,  "node"
	print >> out,  "["
	print >> out,  "id", i
	print >> out,  "]"
	

    for i in range(1, len(content)):
	print >> out,   "edge"
	print >> out,  "["
	print >> out,   "source",content[i].split(',')[0]
	print >> out,  "target",content[i].split(',')[1]
	print >> out,   "value",content[i].split(',')[2]
	print >> out,  "]"

    print >> out,  "]" 
