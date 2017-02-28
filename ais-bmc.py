# For exploring the bidirectional MC from Grosse, Ancha, and Roy... we
# will stick to the finite case for simplicity for now... 
import numpy as np
import numpy.random as npr
from scipy.misc import logsumexp
import matplotlib.pyplot as plt 
import pdb 

from make_toy_data import generate_data, generate_gg_blocks
from mcmc import cllikelihood, lpriorZ 
from simulate import sample_Z_weaklimit

# ---- Set up ---- #
# Generate data from the toy G&G blocks (easy to visualize).  Since
# this is simulated data, we can easily create more!
data_count = 5
data_type = 'gg'
data_set , Z , A = generate_data( data_count , data_type , sigma_n = 0.1 )

# Choose an alpha and sigma 
sigma_a = 5
sigma_n = 0.1
alpha = 2.0 

# Choose number of beta-steps
beta_count = 100 
beta_step = 1.0 / beta_count 

# ---- Helper Functions ---- #
def resample_Z_weighted( Z , data_set , alpha , sigma_n , beta ):
    data_count , feature_count = Z.shape 
    for feature_index in range( feature_count ):
        for data_index in range( data_count ): 
            Z[ data_index , feature_index ] = 1.0                   
            cll1 = cllikelihood( data_set , Z , sigma_a , sigma_n )
            Z[ data_index , feature_index ] = 0.0                  
            cll0 = cllikelihood( data_set , Z , sigma_a , sigma_n )
            mk = np.sum( Z[ : , feature_index ] ) - Z[ data_index , feature_index ]
            p_z1 = ( mk + alpha / feature_count ) / ( data_count + alpha / feature_count )
            lp_z1 = np.log( p_z1 ) + beta * cll1 
            lp_z0 = np.log( 1 - p_z1 ) + beta * cll0
            if np.log( npr.random() ) < ( lp_z1 - logsumexp( np.array( [ lp_z0 , lp_z1 ] ) ) ):
                Z[ data_index , feature_index ] = 1.0
            else: 
                Z[ data_index , feature_index ] = 0.0
    return Z 

# ---- AIS Forward ---- #
k_forward = 5
Z_set = list() 
lower_bound_set = np.zeros( ( k_forward ) )
for k_index in range( k_forward ):
    print "On forward index " , k_index 
    Z = sample_Z_weaklimit( data_count , truncation=4 , alpha=alpha )
    lweight = 0; beta = 0.0 
    ll_set = np.zeros( ( beta_count ) ) # for debugging 
    for beta_index in range( beta_count ):
        beta = beta + beta_step
        cll = cllikelihood( data_set , Z , sigma_a , sigma_n )
        lweight = lweight + beta_step * cll 
        Z = resample_Z_weighted( Z , data_set , alpha , sigma_n , beta )
        ll_set[ beta_index ] = cll

    # save things we may wish to inspect later...     
    lower_bound_set[ k_index ] = lweight 
    Z_set.append( Z ) 
        
# ---- AIS Backward ---- # 
# Generate a Z and a data set 
Z = sample_Z_weaklimit( data_count , truncation=4 , alpha=alpha )
A = generate_gg_blocks()
dim_count = A.shape[1] 
new_data_set = np.dot( Z , A ) + sigma_n * npr.normal( 0 , 1 , [ data_count , dim_count ] )

# Run it backwards
print "On backward" 
lweight = 0; beta = 1.0 
ll_set = np.zeros( ( beta_count ) ) # for debugging 
for beta_index in range( beta_count ):
    beta = beta - beta_step
    cll = cllikelihood( new_data_set , Z , sigma_a , sigma_n )
    lweight = lweight - beta_step * cll 
    Z = resample_Z_weighted( Z , new_data_set , alpha , sigma_n , beta )
upper_bound = lweight 
        
# ---- Explore the Results ---- #
print "Summary" 
print upper_bound
print lower_bound_set
print upper_bound - lower_bound_set 
pdb.set_trace() 

 








