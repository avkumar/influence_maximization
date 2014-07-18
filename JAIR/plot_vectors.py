import matplotlib.pyplot as plt
import numpy as np


ExactrunTime = 0.0119438171387
MCrunTime=115 
# runtime is 12 milli seconds or 11900 microseconds divided among length of y i.e. 2400. so unit is 50 microseconds:: i.e. 
with open('compareBanzGame1_marginalErr') as f:
    read_data = f.read()
    split = int (len(read_data)*0.1)
    read_data = read_data.split(',')[ :split]
    read_data = [read_data[i].strip(']') for i in range(len(read_data))] 

    read_data = [read_data[i].strip('[') for i in range(len(read_data) -2)] 
    x = np.array([float (read_data[i]) for i in range(len(read_data) - 2)])
    y = np.array([0.05 * i for i in range(len(x))])
    y1 = [ExactrunTime for i in range(len(x))] 

    plt.plot(x, y, linestyle=':', label = 'Monte Carlo Solution')
    plt.plot(x, y1, linestyle='--', color='r', label = 'Our Own Solution')
#    plt.plot(x,2*y, marker = 'o', linestyle = '--', color='r', label='banzha_ Nodes')
#    plt.title('Marginal Weight_Cutoff')
    plt.legend(loc='lower right')
    plt.savefig('shapleyAndBanzhaf', rmat = 'eps', dpi=1000)
