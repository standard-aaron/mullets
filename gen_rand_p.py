from __future__ import division 

import numpy as np 
import itertools

def generate_phi(prior):
	"""
	Generates pairwise phi function from a Dirichlet distribution

	input:
		- prior - list w/ entries (A_00, A_01, A_10, A_00)
	
	returns:
		- phi - length 4 list that sums to one 
				(prob 00, prob 01, prob 10, prob 00) 

	Note: Higher P_ij means that the 
	"""
	return np.random.dirichlet(prior)
 
def calc_joint_prob(phi_fns, x, Z):
	"""
	Add documentation 
	"""
	L = len(phi_fns) + 1
	assert len(x) == L
	
	if Z is False:  
		prob = 1
	else: 
	    prob = 1 / Z
	for i in range(L - 1):
		tup = (x[i], x[i+1])
		prob *= phi_fns[i][tup]
	return prob 

def get_normalizer(phi_fns):
	L = len(phi_fns) + 1
	all_tups = list(itertools.product([0, 1], repeat=L)) # All binary permutations 
	total_prob = 0 
	for tup in all_tups:
		total_prob += calc_joint_prob(phi_fns, tup, False) # Get unnormalized prob 
	return total_prob

def make_dist(phi_fns, normalize):
	"""
	Add documentation 
	"""
	if normalize:
		normalizer = get_normalizer(phi_fns) # Might need to delete if L too large... 
	else:
		normalizer = False 

	return lambda x: calc_joint_prob(phi_fns, x, normalizer) # e.g. x = (0, 1, 0, .., 1) length l 

# Generate P from a prior distribution
def generate_P(L, prior_mat, normalize=False):
	"""
	inputs:
		- L:  number of Loci 
		- prior_mat: L-1 x 4 numpy array 
					 prior_mat[i, :] = prior for pairwise phi function F_(i, i+1)

	returns:
		- P: realized prob distribution over {0, 1}^L 
	"""
	phi_fns = []
	for i in range(L - 1):
		prior_i = prior_mat[i, :]
		dist = generate_phi(prior_i)
		dist.shape = (2, 2)
		phi_fns.append(dist)

	return make_dist(phi_fns, normalize)

# Random Example
if __name__ == '__main__': 
	
	L = 3
	prior = np.array([1, 2, 3, 4])
	prior_mat = np.zeros((3, 4))
	for i in range(L):
		prior_mat[i, :] = prior

	P_unorm = generate_P(L, prior_mat, False) # Unormalized prob

	# Whats the probability of seeing (0, 1, 1) under P?
	prob_unorm = P_unorm((0, 1, 1)) 

	P_norm = generate_P(L, prior_mat, True) # True prob distribution 
	prob_norm = P_norm((0, 1, 1))
