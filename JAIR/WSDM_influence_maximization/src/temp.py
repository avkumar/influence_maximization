import random
import time
start_time = time.time()   
copyOfNodes = [i for i in xrange(16000)]
for iter in range(10000):
    coalition = []
    for j in  copyOfNodes:
        if (random.randint(0,1)):
             coalition.append(j)

print time.time() - start_time
