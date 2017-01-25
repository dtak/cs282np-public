# This is a very basic file for variational inference, taken from the
# tech report.  
import numpy as np
import numpy.random as npr 
import scipy.special as sps
import pdb 

# Compute the elbo 
def compute_elbo( data_set , alpha , sigma_a , sigma_n , phi , Phi , nu , tau ):

    return elbo     

# Compute q and Elogstick variables for a given feature index
def compute_q_Elogstick( tau , feature_index ):

    # return
    return qk , Elogstick


# Run the VI  
def run_vi( data_set , alpha , sigma_a=5 , sigma_n=.1 , iter_count=35 , feature_count=10 ):
    data_count = data_set.shape[0]
    dim_count = data_set.shape[1] 
    elbo_set = np.zeros( [ iter_count ] )
    nu_set = list()   # nu are the varitional parameters on Z 
    phi_set = list()  # phi mean param of A 
    Phi_set = list()  # Phi cov param of A, per feat -> same for all dims 
    tau_set = list()  # tau are the variational parameters on the stick betas
    
    # Initialize objects 

    # Optimization loop 
    for vi_iter in range( iter_count ):

        # Update Phi and phi
            
        # Get the intermediate variables
            
        # Update tau, nu
        for feature_index in range( int( feature_count ) ):
            
            # update nu_k

            # update tau

        # Compute the ELBO
        elbo = compute_elbo( data_set , alpha , sigma_a , sigma_n , \
                             phi , Phi , nu , tau )

        # Store things and report  
        elbo_set[ vi_iter ] = elbo 
        nu_set.append( nu )
        phi_set.append( phi )
        Phi_set.append( Phi ) 
        tau_set.append( tau )
        print vi_iter , elbo 

    # return
    return nu_set , phi_set , Phi_set , tau_set , elbo_set


