import pickle
if __name__ == "__main__": 
 
  measures=[                                                                                      
  'eigen',                                                                                        
  'degree', 
  'indegree', 
  'betweenness',                                                                                  
  'closeness', 
  'outdegree'] 
  files = ['power.pkl', 'hep_data.pkl', 'astro-ph.pkl'] 
  for file in files: 
    eigen = pickle.load( open( file, "rb" ) )                                                     
#    for measure in measureValues: 
#        print measure                                                                            
                

