# A skeleton file if you find it helpful 

# Imports 
import numpy as np
import numpy.random as npr 
import scipy.special as sps 
from scipy.misc import logsumexp 
from copy import copy
import pdb

# --- Helper Functions --- # 
# Collapsed Likelihood 
def cllikelihood( data_set , Z , sigma_a=5 , sigma_n=.1 ):

    return ll 

# Uncollapsed Likelihood 
def ullikelihood( data_set , Z , A , sigma_n=.1 ):

    return ll 

# Compute how many times each history occurs [Kh]
def history_count_set( Z ):

    return count_set 

# Prior on Z
def lpriorZ( Z , alpha ):

    return lp 
    
# Prior on A
def lpriorA( A , sigma_a=5 ):

    return lp

# Mean function for A 
def mean_A( data_set , Z , sigma_a=5 , sigma_n=.1 ):

    return A

# --- Resample Functions --- # 
# Resample A 
def resample_A( data_set , Z , sigma_a=5 , sigma_n=.1 ): 

    return A

# Sample next stick length in the stick-breaking representation
def sample_pi_next( pi_min , alpha , data_count ):

    return pi_next 

# --- Samplers --- # 
# Slice sampler from Teh, Gorur, and Ghahramani, fully uncollapsed 
def slice_sampler( data_set , alpha , sigma_a=5 , sigma_n=.1 , iter_count=35 , init_Z=None ):
    data_count = data_set.shape[0]
    dim_count = data_set.shape[1] 
    ll_set = np.zeros( [ iter_count ] )
    lp_set = np.zeros( [ iter_count ] ) 
    Z_set = list()
    A_set = list() 
    
    # Initialize the variables 
    
    # MCMC loop 
    for mcmc_iter in range( iter_count ):
        
        # Sampling existing pi
        
        # Sampling slice_var
            
        # Extending the matrix
        
        # Sampling existing Z

        # Sampling existing A 

        # Compute likelihoods and store 
        ll_set[ mcmc_iter ] = ullikelihood( data_set , Z , A , sigma_n ) 
        lp_set[ mcmc_iter ] = lpriorA( A , sigma_a ) + lpriorZ( Z , alpha ) 
        A_set.append( A ); Z_set.append( Z ) 

        # print
        print mcmc_iter , Z.shape[1] , ll_set[ mcmc_iter ] , lp_set[ mcmc_iter ] 

    # return 
    return Z_set , A_set , ll_set , lp_set 
    
# The uncollapsed LG model. In a more real setting, one would want to
# additionally sample/optimize the hyper-parameters!  
def ugibbs_sampler( data_set , alpha , sigma_a=5 , sigma_n=.1 , iter_count=35 , init_Z=None ):
    data_count = data_set.shape[0]
    dim_count = data_set.shape[1] 
    ll_set = np.zeros( [ iter_count ] )
    lp_set = np.zeros( [ iter_count ] ) 
    Z_set = list()
    A_set = list() 
    
    # Initialize Z randomly (explore how different initializations matter)
        
    # MCMC loop 
    for mcmc_iter in range( iter_count ):
        
        # Sampling existing A 
        
        # Sampling existing Z 

        # Loop through instantiated Z 
                            
        # Remove any unused
       
        # Consider adding new features - decide whether it should be
        # collapsed or uncollapsed.
        
        # Compute likelihood and prior 
        ll_set[ mcmc_iter ] = ullikelihood( data_set , Z , A , sigma_n ) 
        lp_set[ mcmc_iter ] = lpriorA( A , sigma_a ) + lpriorZ( Z , alpha ) 
        A_set.append( A ); Z_set.append( Z ) 

        # print
        print mcmc_iter , Z.shape[1] , ll_set[ mcmc_iter ] , lp_set[ mcmc_iter ] 

    # return 
    return Z_set , A_set , ll_set , lp_set 

            
# The collapsed LG model from G&G.  In a more real setting, one would
# want to additionally sample/optimize the hyper-parameters!
def cgibbs_sampler( data_set , alpha , sigma_a=5 , sigma_n=.1 , iter_count=35 , init_Z=None ):  
    data_count = data_set.shape[0] 
    ll_set = np.zeros( [ iter_count ] )
    lp_set = np.zeros( [ iter_count ] ) 
    Z_set = list()
    A_set = list() 
    
    # Initialize Z randomly (explore how different initializations matter)
        
    # MCMC loop 
    for mcmc_iter in range( iter_count ):
        for data_index in range( data_count ):

            # Consider adding new features
            
            # Loop through instantiated Z 
                           
        # Remove any unused
        
        # Compute likelihood and also the mean value of A, just so we
        # can visualize it later
        ll_set[ mcmc_iter ] = cllikelihood( data_set , Z , sigma_a , sigma_n )
        lp_set[ mcmc_iter ] = lpriorZ( Z , alpha ) 
        A = mean_A( data_set , Z , sigma_a , sigma_n )
        A_set.append( A ); Z_set.append( Z )

        # print
        print mcmc_iter , Z.shape[1] , ll_set[ mcmc_iter ] , lp_set[ mcmc_iter ] 
        
    # return 
    return Z_set , A_set , ll_set , lp_set 
        
    
