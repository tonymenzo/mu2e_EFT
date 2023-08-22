import sys
import mu2e_matching
import numpy as np

# Set the EFT scale
scale = 100 # GeV

for _ in range(100):
	# Give initial conditions for Wilson coefficients as a python dictionary:
	dict1 = {'C61u':1./scale**2}
	
	# The allowed keys depend on the DM type (Dirac, Majorana, ... ) 
	# and the number of active quark flavors. They can be printed via
	
	# 3-flavor, Dirac:
	print('Allowed keys in the three-flavor theory:\n')
	print(mu2e_matching.WC_3f({}).wc_name_list)
	
	# Initialize an instance of the 3-flavor Wilson coefficient class. 
	wc3f = mu2e_matching.WC_3f(dict1)
	print('Parton-level Wilson coefficients: ', wc3f.coeff_dict)
	
	# The main method is to output the NR coefficients. 
	p1 = np.array([wc3f.input_dict['mmu'], 1., 1., 1.])
	q = np.array([0., 1.3, 1.2, 1.2])
	
	wc3f_nucleon = wc3f.match_nucleon_basis(p1,q)
	print('Nucleon-level Wilson coefficients: ', wc3f_nucleon)
	
	# We can create the relevant yaml file and run the decay rate computation
	element     = 'Al'
	isotope     = 0
	interaction = 'bw'
	oscb        = 0
	isochar     = 'proton'
	plots       = 'none'
	mL          = 0.1
	
	wc3f.write_yaml(p1, q, element, isotope, interaction, oscb, isochar, plots, mL, compute_rate = True)