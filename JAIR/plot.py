with open ('out_w1.0', 'r') as f:
    read_data = f.read()
    read_data1 = read_data.split('\n')[9:]
    print len(read_data1)
    print read_data1[0], read_data1[1]

    shapley_Nodes = [int (read_data1[2*i].split(',')[299].strip(']')) for i in range(99)] 
    banzhaf_Nodes = [int (read_data1[2*i+1].split(',')[299].strip(']')) for i in range(99)] 
    nodes = [5*x for x in range(99)]
    print shapley_Nodes

