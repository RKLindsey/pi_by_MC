# Pi by Monte Carlo
# Example code by R. K. Lindsey (2024)

# Run this script with: time python3 pi_by_MC.py <nsamples>

import sys
from numpy.random import uniform as unirand

nsamples = int(sys.argv[1]) # Number of samples to use for calculation
ncircle  = 0                # Number of samples that fall within the circle

for i in range(nsamples):
    
    # Generate coordinates for the random point in a square
    # Assume the square is bounded in x = [0,1], y = [0,1]
    
    x = unirand()
    y = unirand()
    
    # Check whether point falls within circle. If so, increment ncircle
    
    if pow(x,2.0)+pow(y,2.0) <= 1:
        ncircle += 1
        
# Estimate value of pi based on results
print( 4.0*float(ncircle)/float(nsamples))

 


