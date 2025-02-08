# This is a program that generates a random boolean function and checks if it is degenerate.
import sys
import os

# Add the parent directory of the examples folder to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import numpy as np
from boolnetanalyzer import bn_tools as bnt
from boolnetanalyzer import randomize

# Generate a random boolean function and count how many are degenerate
# As n increases, the probability of generating a degenerate function decreases to 0
# This is because the number of functions is 2^2^n
# The number of degenerate functions is at most n.2^2^(n-1)
n = 4
count = 0
for i in range(10000):
    
    function = randomize.booleanfunction(num_nodes=n,seed = i)
    #print(f'Function: {function}')
    #print(f'Degenerated: {bnt.is_degenerated(function)}')
    if bnt.is_degenerated(function):
        count += 1

print(f"Probability of generating a degenerate function: {count/10000}")

