# Imports 
from make_toy_data import generate_data 
from mcmc import cgibbs_sampler , ugibbs_sampler , slice_sampler 
from variational import run_vi 
import matplotlib.pyplot as plt 
import numpy as np 
import pdb 

# ---- Set up ---- #
# Generate data from the toy G&G blocks (easy to visualize) 
data_count = 50
data_type = 'gg' 
data_set , Z , A = generate_data( data_count , data_type )

# What kind of inference?
inference_algorithm = 'cgibbs'
# inference_algorithm = 'ugibbs' 
# inference_algorithm = 'slice'
# inference_algorithm = 'vi' 

# Inference Parameters (note, these do not necessarily match the true parameters!) 
alpha = 2.0  

# ---- Inference ---- #
if inference_algorithm is 'cgibbs': 
    Z_set , A_set , ll_set , lp_set = cgibbs_sampler( data_set , alpha )
    inference_type = 'mcmc' 
if inference_algorithm is 'ugibbs': 
    Z_set , A_set , ll_set , lp_set = ugibbs_sampler( data_set , alpha )
    inference_type = 'mcmc' 
if inference_algorithm is 'slice': 
    Z_set , A_set , ll_set , lp_set = slice_sampler( data_set , alpha )
    inference_type = 'mcmc'
if inference_algorithm is 'vi': 
    nu_set , phi_set , Phi_set , tau_set , elbo_set = run_vi( data_set , alpha )
    inference_type = 'vi'

# ---- A few plots ---- #
# Plot the true and the last sample of features 
plt.figure(1)
plt.subplot( 1 , 3 , 1 )
plt.imshow( A , cmap = 'gray' , interpolation='none' )
plt.xlabel( 'True' )
plt.subplot( 1 , 3 , 2 )
if inference_type is 'mcmc': 
    plt.imshow( A_set[-1] , cmap = 'gray' , interpolation='none' )
    plt.xlabel( 'Final Sample' ) 
else:
    plt.imshow( phi_set[-1] , cmap = 'gray' , interpolation='none' )
    plt.xlabel( 'Final Posterior Mean' ) 

# Plot either the log joint or the elbo 
plt.subplot( 1 , 3 , 3 )
if inference_type is 'mcmc': 
    plt.plot( ll_set )
    plt.xlabel( 'Log-likelihood' ) 
else: 
    plt.plot( elbo_set )
    plt.xlabel( 'Evidence Lower Bound' ) 

# If we are using the 'gg' model, we can also look at the blocks more
# like blocks, to get a better visual 
if data_type is 'gg':
    plt.figure(2)
    if inference_type is 'mcmc': 
        my_A = A_set[-1] 
    else:
        my_A = phi_set[-1] 
    feature_max = np.max( ( 4 , my_A.shape[0] ) )
    for feature_index in range( feature_max ):
        if feature_index < 4: 
            plt.subplot( 2 , feature_max , 1 + feature_index )
            plt.imshow( np.reshape( A[feature_index,:] , (6,6) ) ,
                        cmap = 'gray' , interpolation='none' )
        if feature_index < my_A.shape[0]: 
            plt.subplot( 2 , feature_max , 1 + feature_max + feature_index )
            plt.imshow( np.reshape( my_A[feature_index,:] , (6,6) ) ,
                        cmap = 'gray' , interpolation='none' )
        
# Show plots 
plt.show( block=False )

# Anything else to explore?
pdb.set_trace() 
