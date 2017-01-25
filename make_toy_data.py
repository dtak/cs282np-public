# Creates some toy data: either random A or bars using a fixed number
# of features and Gaussian noise.  

# Imports 
import numpy as np
import numpy.random as npr 
import matplotlib.pyplot as plt 
import pdb 

# The little blocks from G&G, too easy for most applications but
# useful for debugging algorithms
def generate_gg_blocks(): 
    A = np.zeros( [ 4 , 36 ] ) 
    A[0, [ 1 , 6 , 7 , 8 , 13 ] ] = 1
    A[1, [ 3 , 4 , 5 , 9 , 11 , 15 , 16 , 17  ] ] = 1 
    A[2, [ 18 , 24 , 25 , 30 , 31 , 32 ] ] = 1
    A[3, [ 21 , 22 , 23 , 28 , 34 ] ] = 1 
    return A 
    
# Alternative: Random A.  Note that one should probably do something
# in between that's more structured but less toy for exploring!
def generate_random_A( feature_count , dim_count , sigma_a=10 ):
    A = sigma_a * npr.normal( 0 , 1 , [ feature_count , dim_count ] ) 
    return A 
    
# Generate the data
def generate_data( data_count , data_type , sigma=0.1 , sigma_a=5 ):
    # data_type can be 'gg' for the basic blocks, 'finite-random' for
    # a random A but a finite Z (same size as gg), and
    # 'infinite-random' for a simulated Z from the IBP and a random A
    dim_count = 36
    if data_type == 'gg': 
        A = generate_gg_blocks()
        feature_count = A.shape[0]
        dim_count = A.shape[1]
        Z = npr.binomial( 1 , .25 , [ data_count , feature_count ] )
    if data_type == 'finite-random':
        feature_count = 4
        A = generate_random_A( feature_count , dim_count )
        Z = npr.binomial( 1 , .25 , [ data_count , feature_count ] )
    if data_type == 'infinite-random':
        from simulate import sample_Z_restaurant
        Z , feature_count_by_data = sample_Z_restaurant( data_count , 2 )
        feature_count = Z.shape[1]
        A = generate_random_A( feature_count , dim_count ) 
        
    # add noise to the data and return 
    data_set = np.dot( Z , A ) + sigma * npr.normal( 0 , 1 , [ data_count , dim_count ] ) 
    return data_set , Z , A 
    
# Just to see what the features look like 
if __name__ == "__main__":
    data_count = 50
    data_set , Z , A = generate_data( data_count , 'infinite-random' ) 

    plt.figure(1)
    plt.subplot( 1 , 3 , 1 )
    plt.imshow( data_set , cmap = 'gray' , interpolation='none' )
    plt.subplot( 1 , 3 , 2 ) 
    plt.imshow( Z , cmap = 'gray' , interpolation='none' )
    plt.subplot( 1 , 3 , 3 ) 
    plt.imshow( A , cmap = 'gray' , interpolation='none' ) 
    plt.show( block=False ) 

    # Explore anything else that you want to 
    pdb.set_trace() 
