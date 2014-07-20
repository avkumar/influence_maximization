with open ('out_banz_shap_0weight', 'r') as f:
    read_data = f.read()
    read_data1 = read_data.split('\n')[8:]
    print len(read_data1)
    print read_data1[0], read_data1[1]

    shapley_Nodes = [int (read_data1[2*i].split(',')[499].strip(']')) for i in range(499)] 
    banzhaf_Nodes = [int (read_data1[2*i+1].split(',')[499].strip(']')) for i in range(499)] 
    nodes = [5*x for x in range(499)]
    print shapley_Nodes, banzhaf_Nodes, nodes

