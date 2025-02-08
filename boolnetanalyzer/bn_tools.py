# This file contains tools to analyze boolean networks
import numpy as np
import networkx as nx

import numpy as np
import itertools
import os
import sys
import matplotlib.pyplot as plt
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, script_path)

from boolnetanalyzer.load import BooleanNet
import boolnetanalyzer.randomize as rn
import boolnetanalyzer.utils as ut
import time
import json 
# load the environmental variable from the .env file
from dotenv import load_dotenv
load_dotenv()

# Define the path to your Julia file
# Check if the environment variable for julia is set to True
def check_julia_installed():
    try:
        import juliacall
        return True
    except ImportError:
        return False
if check_julia_installed():    
    if os.getenv("JULIA") == "True":
        from juliacall import Main as jl
        julia_package_path = os.path.join(script_path,"boolnetanalyzer")
        julia_package_path = os.path.join(julia_package_path,"cubicalhomology")
        jl.eval('using Pkg')
        jl.eval('Pkg.activate(julia_package_path)')
        jl.eval('Pkg.update(julia_package_path)')
        jl.eval("Pkg.build(julia_package_path)")
        jl.eval('Pkg.status()')
        julia_file_path = os.path.join(julia_package_path,"src")
        julia_file_path = os.path.join(julia_file_path,"cubicalhomology.jl")
        # Check that the environment has been activated
        jl.include(julia_file_path)
        # Use the package within the Julia environment
        jl.eval("using cubicalhomology")
        jl.eval("using SparseArrays")
        jl.eval("using LinearAlgebra")

class Statespace():
    def __init__(self, attractor_states, basin_sizes,num_nodes):
        self.attractors = attractor_states
        self.basin_sizes = basin_sizes
        self.num_attractors = len(attractor_states)
        self.update_dict = {}
        self.attr_dict = {}
        self.num_nodes = num_nodes


    def get_robustness_exact(self):
        '''computes the proportion of neighbors in the Boolean hypercube
        who, following synchronous update, transition to the same attractor,
        
        Note: attr_dict must contain a fully sampled state space'''
        
        attractors = self.attractors
        basin_sizes = self.basin_sizes
        N = self.num_nodes 
        n_attractors = len(attractors)
        
        if n_attractors==1:
            return 1
        b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
        count_of_neighbors_who_transition_to_same_attractor = 0
        for xbin,x in enumerate(self.attr_dict):
            for i in range(N):
                if x[i]==0:
                    ybin = xbin + b_for_bin2dec[i]
                else:
                    continue
                if self.attr_dict[xbin] == self.attr_dict[ybin]:
                    count_of_neighbors_who_transition_to_same_attractor += 1
        return count_of_neighbors_who_transition_to_same_attractor/2**(N-1)/N

def essential_variables(function,regulators):
    """
    Identify the essential variables in a Boolean function.
    The essential variables are the variables that are necessary to determine the output of the function.
    We identify the essential variables by checking if the function output changes when the variable is flipped.

    Parameters:
    ----------
        function (list or np.array): 
            Boolean function represented as a list or numpy array.
            The function is represented as a vector of 2^n elements, where n is the number of variables.
            The function at index i is the output of the function for the binary representation of i.
        regulators (list):
            List of regulators in the network.

    Returns:
    -------
        essential_vars (list): 
            A list of essential variables in the Boolean function.

    Examples:
    --------
        >>> essential_variables(np.array([0, 1, 1, 0]),['A','B'])
        ['A', 'B']
        >>> essential_variables(np.array([0, 1, 1, 1]),['A','B'])
        ['B']
        >>> essential_variables(np.array([1, 1, 1, 1]),['A','B'])
        []
    """
    
    if len(function)==0:
        return regulators
    num_nodes = len(regulators)

    # Initialize a list to store the essential variables
    essential_vars = []

    # Loop over each variable index to see if the function output changes when the variable is flipped
    for var_index in range(num_nodes):
        # Compute the step size for the current variable index
        step_size = 2**(num_nodes-1-var_index)
        
        # Create a mask to identify which entries depend on the current variable
        # If n=3, the mask for the first variable is [0,0,0,0,1,1,1,1]. 
        # The mask for the second variable is [0,0,1,1,0,0,1,1]. 
        # The mask for the third variable is [0,1,0,1,0,1,0,1] and so on.
        variable_mask = (np.arange(2**num_nodes) % (2**(num_nodes-var_index))) // step_size
        
        # Check if the Boolean function output changes when the variable is flipped
        depends_on_variable = np.any(function[variable_mask == 0] != function[variable_mask == 1])
        if depends_on_variable:
            essential_vars.append(var_index)
            
    return [regulators[i] for i in essential_vars]

def essential_truth_tables(truth_tables,regulators):
    """
    Identify the essential variables in a Boolean function.
    The essential variables are the variables that are necessary to determine the output of the function.
    We identify the essential variables by checking if the function output changes when the variable is flipped.
    
    
    
    Parameters:
    ----------
        truth_tables (list or np.array): 
            Boolean function represented as a list or numpy array.
            The function is represented as a vector of 2^n elements, where n is the number of variables.
            The function at index i is the output of the function for the binary representation of i.
        regulators (list):
            List of regulators in the network.
            
    Returns:
    -------
        essential_vars (list): 
            A list of essential variables in the Boolean function.
            
    Examples:
    --------
        >>> essential_truth_tables([np.array([0, 1, 1, 0]),np.array([0, 1, 1, 1])],['A','B'])
        [['A', 'B'], ['B']]
        >>> essential_truth_tables([np.array([0, 1, 1, 1]),np.array([1, 1, 1, 1])],['A','B'])
        [['B'], []]
        >>> essential_truth_tables([np.array([1, 1, 1, 1]),np.array([1, 1, 1, 1])],['A','B'])
        [[], []]
    """
    ess_vars = [essential_variables(np.array(truth_table), regulators) for truth_table, regulators in zip(truth_tables, regulators)]
    ess_truth_tables = []
    for i in range(len(ess_vars)):
        if len(ess_vars[i])==0:
            ess_truth_tables.append([truth_tables[i][0]])
        else:
            degree = len(ess_vars[i])
            X = list(itertools.product([0, 1], repeat = degree))
            f = np.array([],dtype=int)
            for j in range(2**degree):
                        x = X[j]
                        # first set the essential variables to their values
                        variable_values = {var: val for var, val in zip(ess_vars[i], x)}
                        # set the non essential variables to 0
                        variable_values.update({var: 0 for var in regulators[i] if var not in ess_vars[i]})
                        state = [variable_values[var] for var in regulators[i]]
                        f = np.append(f,truth_tables[i][int(''.join(map(str, state)), 2)])
            # convert f to list and append to ess_truth_tables
            ess_truth_tables.append(f.tolist())
    return ess_truth_tables, ess_vars
def is_degenerated(function):
    """
    Check if a Boolean function is degenerated (does not depend on at least one variable).
    We check this by verifying if the function depends on each variable.
    In particular, we check if the function changes its output when the variable is flipped.

    Parameters:
    ----------
        function (list or np.array): 
            Boolean function represented as a list or numpy array.
            The function is represented as a vector of 2^n elements, where n is the number of variables.
            The function at index i is the output of the function for the binary representation of i.


    Returns:
    -------
        bool: 
            True if the function is degenerated, False otherwise.

    Examples:
    --------
        >>> is_degenerated(np.array([0, 1, 1, 0]))
        False
        >>> is_degenerated(np.array([0, 1, 1, 1]))
        False
        >>> is_degenerated(np.array([1, 1, 1, 1]))
        True
    """
    num_nodes = int(np.log2(len(function)))

    # Loop over each variable index to see if the function depends on it

    for var_index in range(num_nodes):
        # Compute the step size for the current variable index
        step_size = 2**(num_nodes-1-var_index)
        
        # Create a mask to identify which entries depend on the current variable
        # If n=3, the mask for the first variable is [0,0,0,0,1,1,1,1]. 
        # The mask for the second variable is [0,0,1,1,0,0,1,1]. 
        # The mask for the third variable is [0,1,0,1,0,1,0,1] and so on.
        # If the Boolean function does not depend on the current variable, return True
        variable_mask = (np.arange(2**num_nodes) % (2**(num_nodes-var_index))) // step_size
        
        # Check if the Boolean function depends on the current variable
        depends_on_variable = np.any(function[variable_mask == 0] != function[variable_mask == 1])
        
        
        if not depends_on_variable:
            return True

            
    return False

def num_attractors(boolean_network, max_steps=1000, num_simulations=100, initial_sample_points=[], state_transition_graph=False, exact=True,output_file="output.txt"):
    """
    Calculate the number of attractors in a Boolean network.
    The algorithm uses a random sampling approach to identify attractors in the network.
    In particular, we start with a random initial state and simulate the network until we reach either a known attractor basin or a new attractor state.
    We repeat this process for a number of simulations to identify all attractors in the network.

    Parameters:
    ----------
        boolean_network (class boolean network): 
            Boolean network represented as a list or numpy array.
            The network is represented as a vector of 2^n elements, where n is the number of variables.
            The network at index i is the output of the network for the binary representation of i.
        max_steps (int, optional): 
            The maximum number of steps to simulate the network. Default is 1000.
        num_simulations (int, optional): 
            The number of simulations to run. Default is 100.
        initial_sample_points (list, optional): 
            A list of initial sample points to start the simulations from. If not provided, points will be chosen randomly. Default is an empty list.
        exact (bool, optional):
            If True, the function will consider all possible initial states. Default is True.
        state_transition_graph (bool, optional):
            If True, the function will output the state transition graph. Default is False.
        output_file (str, optional):
            The name of the output file to save the state transition graph. Default is "output.txt".
        
    
    Returns:
    -------
        attractors (list):
            A list of attractors in the network. Each attractor is a list of states.
        basin_sizes (list):
            A list of basin sizes for each attractor.
        attr_dict (dict):   
            A dictionary mapping each state to its corresponding attractor index.
        update_dict (dict):
            A dictionary mapping each state to its next state.
    
    Examples:
    --------
        >>> attractors,  basin_sizes,  attr_dict, update_dict = num_attractors(
        rn.boolnet(5, seed=9), initial_sample_points=[i for i in range(25)])
        >>> # Print the results
        >>> print(f"Attractors: {attractors}")
        >>> print(f"Number of attractors: {len(attractors)}")
        >>> print(f"Basin sizes: {basin_sizes}")
        >>> print(f"Number of indices considered: {len(attr_dict)}")
        >>> print(f"Attractor index dict: {attr_dict}")
        >>> print(f"Update dict: {update_dict}")
        Attractors: [[11, 22], [27, 30, 15, 21]]
        Number of attractors: 2
        Basin sizes: [7, 20]
        Number of attractor index dict entries: 27
        Attractor index dict: {0: 0, 11: 0, 22: 0, 1: 1, 19: 1, 27: 1, 30: 1, 15: 1, 21: 1, 2: 1, 3: 1, 4: 1, 5: 1, 17: 1, 6: 1, 7: 1, 8: 1, 14: 1, 9: 0, 18: 0, 10: 1, 12: 1, 13: 1, 16: 0, 20: 0, 23: 1, 24: 1}
        Update dict: {0: 11, 11: 22, 22: 11, 1: 19, 19: 27, 27: 30, 30: 15, 15: 21, 21: 27, 2: 3, 3: 19, 4: 1, 5: 17, 17: 19, 6: 1, 7: 17, 8: 14, 14: 5, 9: 18, 18: 11, 10: 6, 12: 5, 13: 17, 16: 11, 20: 11, 23: 19, 24: 14}
        
    """
    
    
    if isinstance(boolean_network, BooleanNet):
        # Identify the format of the input network
        if boolean_network.input_type == "rtt-format":
                        
            num_nodes = boolean_network.num_nodes
            
            # Initialize dictionaries and lists for tracking attractors and basin sizes
            
            attractors = []
            basin_sizes = []
            attractor_index_dict = {}
            update_dict = {}
            
            queue = []
            if initial_sample_points:
                num_simulations = len(initial_sample_points)

            if exact:
                initial_sample_points = [i for i in range(2**num_nodes)]
                num_simulations = len(initial_sample_points)
            # Sample initial points and update the network state
            # We use two variables to track the state vector
            # state_vector: The current state of the network as a list of binary values
            # state_dec: The current state of the network as a decimal value. For example, [1, 0, 1] would be 5.
            # The reason for this redundancy is two fold:
            # The binary format is needed to update the state vector
            # The decimal format allows us to efficiently store the state transition graph and the attractor index dictionary
            
            for step in range(num_simulations):
                np.random.seed(step)
    
                if not initial_sample_points:
                    # Randomly initialize the state vector if no initial points are provided
                    state_vector = np.random.randint(2, size=num_nodes)
                    state_dec = ut.bin2dec(state_vector)
                    
                else:
                    # Use provided initial sample points
                    
                    state_dec = initial_sample_points[step]
                    state_vector = ut.dec2bin(state_dec, num_nodes)
                
                # Initialize the queue with the current state
                # The queue is used to track the state vector till we reach an attractor basin or a new attractor
                # We use decimal values to save space
                queue = [state_dec]
                
                # Check if the current state is part of a known attractor basin
                # If not, continue the search till we reach an attractor basin or a new attractor
                # Else we ignore that state and move to the next state
                
                
                if state_dec not in attractor_index_dict:
                    # Update the state vector till we reach an attractor basin or a new attractor
                    while len(queue) < max_steps:
                        # Update the state vector using the boolean network
                        next_state_vector = boolean_network.update(state_vector)
                        next_state_dec = ut.bin2dec(next_state_vector)
                        
                        # Optionally update the state transition graph
                        if state_dec not in update_dict:
                            update_dict[state_dec] = next_state_dec
                        
                        state_vector = next_state_vector
                        
                        try:
                            # Check if the current state is part of a known attractor basin
                            
                            attractor_index = attractor_index_dict[next_state_dec]
                            
                            
                            # All elements in the queue are part of the same attractor
                            # Update attractor index dictionary for all states in the queue
                            basin_sizes[attractor_index] += len(queue)
                            attractor_index_dict.update(dict(zip(queue, [attractor_index] * len(queue))))
                            break
                        except KeyError:
                            try:
                                # Check if the current state is part of a new attractor
                                attractor_start_index = queue.index(next_state_dec)
                                
                                # New attractor identified
                                # Update attractor index dictionary for all states in the queue
                                attractor_index_dict.update(dict(zip(queue, [len(attractors)] * len(queue))))
                                attractors.append(queue[attractor_start_index:])
                                basin_sizes.append(len(queue))
                                break
                            except ValueError:
                                # Continue the search if no attractor is found
                                pass
                        
                        queue.append(next_state_dec)
                        state_dec = next_state_dec

            # Save the state transition graph if required
            if state_transition_graph:
                with open(output_file, "w") as f:
                    json.dump(update_dict, f)

            # Return the results
            return attractors, basin_sizes, attractor_index_dict, update_dict



        

    else:
        # If the boolean network is not of type list or numpy array, raise an error
        raise ValueError("Invalid input type for boolean_network. Expected list or numpy array.")
    return 1


def canalizing_variables(boolean_func, variables=[]):
    '''
    Determines if a Boolean function is canalizing.
    
    A Boolean function f(x_1,...,x_n) is canalizing if it is canalizing in at least one variable.
    A Boolean function f(x_1,...,x_n) is canalizing in x_i if fixing x_i to a specific value (0 or 1) 
    results in a constant function value for all other variables.

    Input:
        boolean_func: A Boolean function represented as a vector. This is the output of the truth table 
                      (a list or array of length 2^n, where n is the number of input variables).
        num_vars (optional): The number of variables in the Boolean function. If not provided, it is 
                             inferred from the length of `boolean_func`.

    Output:
        True if the Boolean function is canalizing, False otherwise.

    ..note:: The function has serious issues with the current implementation. It is not working as expected.

    Examples:

        >>>is_canalizing([1, 0, 0, 1])
        False
        >>>is_canalizing([1, 1, 1, 1])
        True
        >>>is_canalizing([1, 0, 1, 1, 1, 0, 1, 1])
        True
    '''
    
    # Convert boolean_func to a numpy array if it is a list
    if isinstance(boolean_func, list):
        boolean_func = np.array(boolean_func)
        
    # Calculate the number of variables if not provided
    if len(variables) == 0:
        num_vars = int(np.log2(len(boolean_func)))
        variables = [i for i in range(num_vars)]
    # Compute half of the length of the boolean function vector (2^(num_vars-1))
    half_length = 2**(num_vars - 1)
    
    # Generate all possible input combinations for `num_vars` variables
    input_combinations = np.array(list(itertools.product([0, 1], repeat=num_vars))).T
    
    # Compute the complement of each input combination (i.e., 1 - x_i)
    complement_combinations = 1 - input_combinations
    
    # Stack the input combinations and their complements vertically
    # For n=2, the stacked_combinations matrix will look like:
    # [[0 0 1 1]
    #  [0 1 0 1]]  # Input combinations

    # [[1 1 0 0]
    #  [1 0 1 0]]   # Complement combinations
    # The stacked_combinations matrix has 2*n rows and 2^n columns
    # [[0 0 1 1]
    #  [0 1 0 1]
    #  [1 1 0 0]
    #  [1 0 1 0]]
    stacked_combinations = np.vstack((input_combinations, complement_combinations))
    
    # Calculate the dot product of the stacked combinations matrix with the boolean function
    # This gives us the values of the function when each variable is fixed to 0 or 1
    result_matrix = np.dot(stacked_combinations, boolean_func)
    
    
    # Check if the function is canalizing by looking for any result that equals half_length or 0
    canalizing_vars = []
    for i in range(num_vars):
        if np.any(result_matrix[i] == half_length) or np.any(result_matrix[i] == 0):
            #append the canalizing variable to the list with the canalizing input
            canalizing_vars.append((num_vars-i-1,1))
    for i in range(num_vars, 2*num_vars):
        if np.any(result_matrix[i] == half_length) or np.any(result_matrix[i] == 0):
            #append the canalizing variable to the list with the canalizing input
            canalizing_vars.append((2*num_vars-i-1,0))
            
    # if np.any(result_matrix == half_length) or np.any(result_matrix == 0):
    #     return True
    return [(variables[i],j) for i,j in canalizing_vars]

def get_robustness_exact(net,left_side_of_truth_table=[]):
    '''computes the proportion of neighbors in the Boolean hypercube
    who, following synchronous update, transition to the same attractor,
    
    Note: attr_dict must contain a fully sampled state space'''
    
    attractors, basin_sizes, attractor_dict, update_dict = num_attractors(net)
    N = net.num_nodes 
    n_attractors = len(attractors)
    
    if left_side_of_truth_table==[]: #to reduce run time, this should be calculated once and then passed as argument
        left_side_of_truth_table = list(itertools.product([0, 1], repeat = N))    
    if n_attractors==1:
        return 1
    b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    count_of_neighbors_who_transition_to_same_attractor = 0
    for xbin,x in enumerate(left_side_of_truth_table):
        for i in range(N):
            if x[i]==0:
                ybin = xbin + b_for_bin2dec[i]
            else:
                continue
            if attractor_dict[xbin] == attractor_dict[ybin]:
                count_of_neighbors_who_transition_to_same_attractor += 1
    return count_of_neighbors_who_transition_to_same_attractor/2**(N-1)/N


def get_robustness_and_attractors_simulation(net,number_different_IC = 500):
    ''' 
    Computes the number of attractors as well as their robustness. It does this
    sampling the attractor landscape as well as the proportion of neighbors in 
    the Boolean hypercube who, following synchronous update, transition to the
    same attractor.
    
    Inputs: 
        boolean_network (class boolean network): 
            Boolean network represented as a list or numpy array.
            The network is represented as a vector of 2^n elements, where n is the number of variables.
            The network at index i is the output of the network for the binary representation of i.
        number_different_IC:
            An integer which specifies the number of different initial conditions used
            for probing the state space for attractors.
    Outputs:
        attractors: a list where each entry is represents an attractor 
        len(attractors): the number of attractors
        
            
    '''

    dictF = dict()
    attractors = []
    basin_sizes = []
    attractor_dict = dict()
    N = net.num_nodes
    
    if net.input_type == 'rtt-format':
        F = net.truth_tables
        I = net.regulators    
        degrees = list(map(len,I))
    
        b_for_bin2decs = [np.array([2**i for i in range(NN)])[::-1] for NN in range(max(degrees)+1)]
        b_for_bin2dec = np.array([2**i for i in range(N)])[::-1]
    
        robustness_approximation = 0
    
        for i in range(number_different_IC):
            index_attractors = []
            for j in range(2):
                if j==0:
                    x = np.random.randint(2, size = N)
                    xbin = np.dot(x,b_for_bin2dec)
                    x_old = x.copy()
                else:
                    x = x_old
                    random_flipped_bit = np.random.choice(N)
                    x[random_flipped_bit] = 1 - x[random_flipped_bit]
                    xbin = np.dot(x,b_for_bin2dec)                
                queue = [xbin]
                while True: #as long as we haven't reached an attractor state, keep updating
                    try:
                        fxbin = dictF[xbin]
                    except KeyError:
                        fx = []
                        for jj in range(N):
                            #fx.append(F[i][sum([x[I[i]][degrees[i]-j-1]*b_for_bin2dec[j] for j in range(degrees[i])])])
                            #fx.append(F[i][sum([x[I[i]][j]*b_for_bin2dec[j] for j in range(degrees[i])])])
                            #fx.append(F[i][np.dot(x[I[i]] , b_for_bin2dec[N-degrees[i]:])])
                            fx.append(F[jj][np.dot(x[I[jj]] , b_for_bin2decs[degrees[jj]])])
                            
                        #fx = update(F,I,N,x)
                        fxbin = np.dot(fx,b_for_bin2dec)#sum([fx[N-i-1]*b_for_bin2dec[i] for i in range(N)])
                        dictF.update({xbin:fxbin})
                        x = np.array(fx) #works only because if we don't know fx now we also won't know F[fx] 
                    try: # check if this state has a known attractor
                        index_attr = attractor_dict[fxbin] #returns a KeyError if we haven't identified fxbin as an attractor before
                        basin_sizes[index_attr] += 1
                        attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                        break
                    except KeyError:
                        try: #check if fxbin is part of a new attractor
                            index = queue.index(fxbin) #returns a ValueError if fxbin is not part of queue
                            #new attractor
                            index_attr = len(attractors)
                            attractor_dict.update( list(zip( queue , [index_attr]*len(queue) )) )
                            attractors.append(queue[index:])
                            basin_sizes.append(1)
                            break
                        except ValueError:
                            pass
                    queue.append(fxbin)
                    xbin = fxbin
                index_attractors.append(index_attr)
            if index_attractors[0] == index_attractors[1]:
                robustness_approximation += 1
                            
        return (attractors, len(attractors), basin_sizes, robustness_approximation/number_different_IC)
    
def cubical_homology_poset(vertices,edges,simplex_dim,codimension,cubical_homology_dim):
    """
    This function calculates the cubical homology of a graph given the vertices and edges of the graph.
    This requires Julia to be installed and the JULIA environment variable to be set to True.
    Parameters:
    ----------
        vertices (list): 
            List of vertices in the graph.
        edges (list): 
            List of edges in the graph.
        simplex_dim (int): 
            The dimension of the simplices we are considering
        codimension (int): 
            The codimension of the simplices that are shared by the simplices.

    Returns:
    -------
        int: 
            The betti number of the graph.

    Examples:   
    --------
        >>> vertices = [0,1,2,3,4,5,6,7]
        >>> edges = [(0,1),(1,2),(2,3),(3,4),(4,0),(5,0),(5,1),(6,1),(6,2),(7,2),(7,3),(7,4),(7,5),(7,6)]
        >>> cubical_homology_poset(vertices,edges,2,0)
        1
    ..note:: 
        If the graph has cycles then the function will throw StackOverflowError.
    
    
    """
    
    vertices_jl = jl.Array(vertices)
    edges_jl = jl.Array(edges)
    edges_jl=jl.cubicalhomology.make_transitive(edges_jl)
    poset = jl.cubicalhomology.create_poset(vertices_jl,edges_jl)
    

    graph = jl.cubicalhomology.poset_to_graph(poset,simplex_dim,codimension) 
    #log.write(f"Preprocessing the graph  with {len(graph.vertices)} vertices and {len(graph.edges)} edges")

    graph = jl.cubicalhomology.relabel_vertices(graph)
    post_graph = jl.cubicalhomology.preprocess_graph(graph)
    # save the graph to a file
    G = nx.Graph()
    G.add_nodes_from(post_graph.vertices)
    G.add_edges_from(post_graph.edges)
    plt.figure()
    plt.title("Graph")
    nx.draw(G)

    plt.savefig("graph.png")
    log.write("Image of the graph stored at graph.png")
    #input("Press Enter to continue...")
    #log.write(f"Calculating cubical homology in dimension {cubical_homology_dim} for graph with {len(graph.vertices)} vertices and {len(graph.edges)} edges")
    return jl.cubicalhomology.cubicalHomology(post_graph,cubical_homology_dim)

def order_complex(vertices,edges):
    """
    This function orders the simplicial complex given the vertices and edges of the graph.
    Parameters:
    ----------
        vertices (list): 
            List of vertices in the graph.
        edges (list): 
            List of edges in the graph.
    
    Returns:
    -------
        Julia list: 
            The simplicial complex ordered by dimension.
    
    Examples:
    --------
        >>> vertices = [0,1,2,3,4,5,6,7]
        >>> edges = [(0,1),(1,2),(2,3),(3,4),(4,0),(5,0),(5,1),(6,1),(6,2),(7,2),(7,3),(7,4),(7,5),(7,6)]
        >>> order_complex(vertices,edges)
    
    ..note::
        The object will throw a StackOverflowError if the graph has cycles. Furthermore, the function will return a Julia object.
    """
    vertices_jl = jl.Array(vertices)
    edges_jl = jl.Array(edges)
    edges_jl=jl.cubicalhomology.make_transitive(edges_jl)
    poset = jl.cubicalhomology.create_poset(vertices_jl,edges_jl)
    return jl.cubicalhomology.poset_to_simplicial(poset)

def a_theory_simplicial(simplicial_complex,cubical_homology_dim):
    """
    This function calculates the cubical homology of a simplicial complex.
    Parameters:
    ----------
        simplicial_complex (Julia list): 
            The simplicial complex.
        cubical_homology_dim (int): 
            The dimension of the cubical homology to be computed.
    
    Returns:
    -------
        int: 
            The betti number of the simplicial complex.
    
    Examples:
    --------
        >>> simplicial_complex = [[0,1,2],[1,2,3],[2,3,4],[3,4,0],[4,0,1],[5,0,1],[5,1,2],[6,1,2],[6,2,3],[7,2,3],[7,3,4],[7,4,0],[7,5,0],[7,6,1]]
        >>> a_theory_simplicial(simplicial_complex,2)
        1
    """
    graph = jl.cubicalhomology.simplicial_to_graph(simplicial_complex)
    graph = jl.cubicalhomology.relabel_vertices(graph)
    post_graph = jl.cubicalhomology.preprocess_graph(graph)
    return graph,post_graph,jl.cubicalhomology.cubicalHomology(post_graph,cubical_homology_dim)
    
def simplicial_to_graph(simplicial_complex,simplex_dim,codimension):
    """
    This function converts a simplicial complex to a graph.
    Parameters:
    ----------
        simplicial_complex (Julia list): 
            The simplicial complex.
        simplex_dim (int): 
            The dimension of the simplices we are considering
        codimension (int): 
            The codimension of the simplices that are shared by the simplices.
    
    Returns:
    -------
        Julia list: 
            The graph.

    Examples:
    --------
        >>> simplicial_complex = [[0,1,2],[1,2,3],[2,3,4],[3,4,0],[4,0,1],[5,0,1],[5,1,2],[6,1,2],[6,2,3],[7,2,3],[7,3,4],[7,4,0],[7,5,0],[7,6,1]]
        >>> simplicial_to_graph(simplicial_complex,2,1)
    
    
    """
    
    graph = jl.cubicalhomology.simplicial_to_graph(simplicial_complex,simplex_dim,codimension)
    graph = jl.cubicalhomology.relabel_vertices(graph)
    log.write(f"Preprocessing the graph  with {len(graph.vertices)} vertices and {len(graph.edges)} edges")
    graph = jl.cubicalhomology.preprocess_graph(graph)
    return graph

def wiring_diagram(network):
    """
    Parameters:
    ----------
        network (BooleanNet): 
            A Boolean network object.
        

    Returns:
    -------
        list: 
            The wiring diagram

    Examples:   
    --------
        >>> num_nodes = 5
        >>> wiring_diagram_list= [(i,i+1) for i in range(num_nodes-1)]+ [(num_nodes//2,0)] + [(num_nodes-1,num_nodes//2)]
        >>> print(wiring_diagram_list)
        [(0, 1), (1, 2), (2, 3), (3, 4), (2, 0), (4, 2)]
        
        >>> boolean_network = rn.boolnet(num_nodes=num_nodes,seed = 0, edges_wiring_diagram = wiring_diagram_list)   
        >>> print(wiring_diagram(boolean_network))
        [(2, 0), (0, 1), (1, 2), (4, 2), (2, 3), (3, 4)]
    
    """
    if network.input_type == "text-format":
        wiring_diagram = [(network.variables.index(regulator),i)  for i, regulators in enumerate(network.regulators) for regulator in regulators]


    if network.input_type == "rtt-format":
        wiring_diagram = [(regulator, i) for i, regulators in enumerate(network.regulators) for regulator in regulators]

    return wiring_diagram
def strongly_connected_components(wiring_diagram):
    """
    Computes the DAG of strongly connected components of a wiring diagram.
    Parameters:
    ----------
        wiring_diagram (list): 
            The wiring diagram of the network.
            
    Returns:
    -------
        compoents (list): 
            The DAG of strongly connected components.
        DAG (list): 
            The DAG of strongly connected components.
    Examples:
    --------
        >>> wiring_diagram = [(1, 0), (0, 1), (0, 2)]
        >>> nodes,edges = strongly_connected_components(wiring_diagram)
        >>> print(nodes)
        [(0, {'members': {2}}), (1, {'members': {0, 1}})]
        >>> print(edges)
        [(1, 0)]
    
    """
    

    # Create a directed graph
    G = nx.DiGraph()

    # Add edges to the graph (example graph)
    G.add_edges_from(wiring_diagram)
    
    directed_acyclic_graph = nx.condensation(G)
    
    return directed_acyclic_graph.nodes.data(),directed_acyclic_graph.edges
def network_scc(network):
    """
    Computes the DAG of strongly connected components of a network
    Parameters:
    ----------
        network (list): 
            A Boolean Network.
            
    Returns:
    -------
        components (list): 
            The DAG of strongly connected components.
        DAG (list): 
            The DAG of strongly connected components.
    Examples:
    --------
        >>> wiring_diagram = [(1, 0), (0, 1), (0, 2)]
        >>> boolean_network = rn.boolnet(num_nodes=3,seed = 0, edges_wiring_diagram = wiring_diagram_list)
        >>> nodes,edges = network_scc(boolean_network)
        >>> print(nodes)
        [(0, {'members': {2}}), (1, {'members': {0, 1}})]
        >>> print(edges)
        [(1, 0)]
    
    """
    return strongly_connected_components(wiring_diagram(network))

def q_values(adj_matrix):
    if adj_matrix.size==0:
        return [0]
    adj_matrix_transpose = adj_matrix.T
    shared_face_matrix = np.dot(adj_matrix, adj_matrix_transpose)-1
    if shared_face_matrix.size==0:
        return [0]    
    q_values = []
    dim_simpl_complx = np.max(np.sum(adj_matrix,axis=1))
    for q in range(0,int(dim_simpl_complx)):
        q_adj_matrix = [
            (i, j)
            for i in range(shared_face_matrix.shape[0])
            for j in range(shared_face_matrix.shape[1])
            if shared_face_matrix[i,j] >= q
        ]
        if q_adj_matrix==[]:
            q_values.append(0)
            continue
        #print(adj_matrix)
        G = nx.Graph()
        G.add_edges_from(q_adj_matrix)
        # count the number of connected components
        connected_components = nx.number_connected_components(G)
        q_values.append(connected_components)
        #print(connected_components)
    return q_values


def eccentricity(adj_matrix):
    if adj_matrix.size==0:
        return [0]
    adj_matrix_transpose = adj_matrix.T
    shared_face_matrix = np.dot(adj_matrix, adj_matrix_transpose)
    
    # Create a mask that zeros out the diagonal elements
    mask = np.ones(shared_face_matrix.shape, dtype=bool)
    np.fill_diagonal(mask, 0)

    # Apply the mask to zero out diagonal elements
    masked_matrix = np.where(mask, shared_face_matrix, 0)
    #print(masked_matrix)
    # Find the maximum value in the masked matrix
    max_shared_faces = np.max(masked_matrix,axis=0)
    


    
    eccentricity = [np.inf if max_shared_face == 0 else (shared_face_matrix[i, i] - max_shared_face) / max_shared_face for i, max_shared_face in enumerate(max_shared_faces)]

    return eccentricity
    

if __name__ == "__main__":

    adj_matrix = np.array([[0, 1, 1, 0, 0, 0],[1, 0, 1, 1, 0, 0],[1, 1, 0, 0, 1, 0],[0, 1, 0, 0, 1, 1],[0, 0, 1, 1, 0, 1],[0, 0, 0, 1, 1, 0]])
    print(eccentricity(adj_matrix))
    adj_matrix = np.array([[1,0,1,1],[0,1,1,1],[0,0,1,0],[0,0,0,1]])
    print(eccentricity(adj_matrix))
    input("Press Enter to continue...")
    # Usage examples
    # Essential variables
    
    print(essential_truth_tables([np.array([1, 1, 1, 0,1,1,1,0]),np.array([1, 0, 1, 0])],[['A','B','c'],['A','B']]))
    # >>> ([[1, 1, 1, 0], [1, 0]], [['B', 'c'], ['B']])
    print(essential_variables(np.array([1, 0, 1, 0]),['A','B']))
    # >>> ['B']
    # Degenerated function
    print(is_degenerated(np.array([1, 1, 1, 1])))
    # >>> True
    attractors,  basin_sizes,  attr_dict, update_dict = num_attractors(
        rn.boolnet(5, seed=9), initial_sample_points=[i for i in range(25)])
    # Print the results
    print(f"Attractors: {attractors}")
    print(f"Number of attractors: {len(attractors)}")
    print(f"Basin sizes: {basin_sizes}")
    print(f"Number of indices considered: {len(attr_dict)}")
    print(f"Attractor index dict: {attr_dict}")
    print(f"Update dict: {update_dict}")
    # >>> Attractors: [[11, 22], [27, 30, 15, 21]]
    # >>> Number of attractors: 2
    # >>> Basin sizes: [7, 20]
    # >>> Number of attractor index dict entries: 27
    # >>> Attractor index dict: {0: 0, 11: 0, 22: 0, 1: 1, 19: 1, 27: 1, 30: 1, 15: 1, 21: 1, 2: 1, 3: 1, 4: 1, 5: 1, 17: 1, 6: 1, 7: 1, 8: 1, 14: 1, 9: 0, 18: 0, 10: 1, 12: 1, 13: 1, 16: 0, 20: 0, 23: 1, 24: 1}
    # >>> Update dict: {0: 11, 11: 22, 22: 11, 1: 19, 19: 27, 27: 30, 30: 15, 15: 21, 21: 27, 2: 3, 3: 19, 4: 1, 5: 17, 17: 19, 6: 1, 7: 17, 8: 14, 14: 5, 9: 18, 18: 11, 10: 6, 12: 5, 13: 17, 16: 11, 20: 11, 23: 19, 24: 14}
    print(canalizing_variables([1, 0, 1, 1, 1, 0, 1, 1]))
    # >>> 
    input("Press Enter to continue...")
    # Generate attractors and related data
    start_time = time.time()
    # attractors,  basin_sizes,  attr_dict, update_dict = num_attractors(
    #     rn.boolnet(20, seed=9), initial_sample_points=[i for i in range(2**20)], state_transition_graph=True, output_file="output.txt"
    # )
    # stop_time = time.time()
    # print(f"Time taken: {stop_time-start_time}")
    # # Print the results
    # print(f"Attractors: {attractors}")
    # print(f"Number of attractors: {len(attractors)}")
    # print(f"Basin sizes: {basin_sizes}")
    # print(f"Number of indices considered: {len(attr_dict)}")
    # Code to test GPU implementation
    #attractors,  basin_sizes,  attr_dict, update_dict = num_attractors_gpu(rn.boolnet(10, seed=9), initial_sample_points=[i for i in range(1526)])
    #print(is_canalizing([1, 0, 1, 1, 1, 0, 1, 1]))
    
    print(essential_variables(np.array([0, 1, 1, 0]),[0,2]))
    print(essential_variables(np.array([0, 1, 1, 1]),[0,2]))
    print(essential_variables(np.array([1, 1, 1, 1]),[0,2]))
    
    # add_numbers = Main.eval('add_numbers')
    # # Call specific Julia functions
    # result_add = Main.add_numbers(3,5)
    print("Testing A-theory functions")
    # s = jl.eval("SparseVector(10, [1, 3, 7], [1.0, -2.0, 3.0])")

    # # Retrieve the sparse vector 's' from Julia
    
    # print(s)
    # result = jl.cubicalhomology.nCube(3)
    # print(result)
    # print("############################################")
    # result = jl.cubicalhomology.hyperOctahedrial(3)
    # print(result)
    # print("############################################")
    # print(jl.cubicalhomology.generate_reversals(3))
    # print("############################################")
    # print(jl.cubicalhomology.add_numbers(3,5))
    # print("############################################")
    # print(jl.cubicalhomology.permutations([i for i in range(3)]))
    # print("############################################")
    # print(jl.cubicalhomology.sign_of_permutation([1,2,0]))
    # print("############################################")
    # for g in jl.cubicalhomology.hyperOctahedrial(3):
    #     print(jl.cubicalhomology.calculateImageCube(jl.cubicalhomology.nCube(3),g))
    # #print(jl.cubicalhomology.permuteCoords(vertex, gpElet))
    # v=[0,1,2,3,4, # regular C5 vertices
    # 5,6,7,8,9] # extra vertices
    # e = [(0,1),(1,2),(2,3),(3,4),(4,0), # regular edges
    #     (5,0),(5,1),(6,1),(6,2),(7,2),(7,3),(8,3),(8,4),(9,4),(9,0)] # extra edges 

    # e = jl.cubicalhomology.makeSymmetricReflexive(e)
    # C5_star = jl.cubicalhomology.graph(v,e)
    # C5_star = jl.cubicalhomology.relabel_vertices(C5_star)
    # print("############################################")
    # print(C5_star)
    # print("############################################")
    # print(jl.cubicalhomology.cubicalHomology(jl.cubicalhomology.preprocess_graph(C5_star),1))
    # print("############################################")
    # before_import = set(dir(jl))
    

    # # Define the Poset struct and create an instance within a Julia module
    
    # after_import = set(dir(jl))
    # imported_items = after_import - before_import
    # print(imported_items)
    # print(jl.cubicalhomology.add_numbers(3,5))
    vertices = [1, 2, 3, 4,5,6,7]
    relations = [(1, 2), (1,3),(1,4),(2, 5),(2,6),(3,7),(3,6),(4,5),(4,7)]
    print(cubical_homology_poset(vertices,relations,2,1,1))
    #relations=jl.cubicalhomology.make_transitive(relations)
    #print(jl.cubicalhomology.create_poset(vertices,relations))

    # # Access the Poset instance from the Julia module
    # # Use module name to access variables and functions defined in it
    # poset_obj = jl.eval("MyModule.P")

    # # Access fields of the Poset instance
    # vertices = jl.eval("MyModule.P.vertices")
    # relations = jl.eval("MyModule.P.relations")

    # print("Vertices:", vertices)
    # print("Relations:", relations)


    # # Access fields of the Poset instance
    # vertices = jl.eval("MyModule.P.vertices")
    # relations = jl.eval("MyModule.P.relations")

    # print("Vertices:", vertices)
    # print("Relations:", relations)
    # jl.eval("""
    # mutable struct Poset
    #     vertices::Vector{Int}
    #     relations::Vector{Tuple{Int, Int}}
    # end
            
    # P = Poset([1, 2, 3, 4], [(1, 2), (2, 3)])
    # """)
    # after_import = set(dir(jl))
    # imported_items = after_import - before_import
    # print(imported_items)

    # # Define your data in Python
    # vertices = [1, 2, 3, 4]
    # relations = [(1, 2), (2, 3)]

    # # Convert Python lists to Julia arrays
    # vertices_jl = jl.eval(f"Vector({vertices})")
    # relations_jl = jl.eval(f"Vector({relations})")
    # print(vertices_jl)
    # print(relations_jl)
    # # Create the Poset object in Julia
    # jl.eval(f"poset_obj = Poset({vertices_jl}, {relations_jl})")

    # poset_obj = jl.eval("poset_obj")

    # # Print vertices and relations
    # print("Vertices:", jl.eval("poset_obj.vertices"))
    # print("Relations:", jl.eval("poset_obj.relations"))


    # Now call the Julia function
    # print("############################################")
    # poset = jl.cubicalhomology.create_poset(vertices,relations)
    # print(poset)
    # #print(jl.cubicalhomology.poset_to_simplicial(jl.cubicalhomology.create_poset(vertices,relations)))
    # simplicial = jl.cubicalhomology.poset_to_graph(poset,2,0)
    # simplicial = jl.cubicalhomology.relabel_vertices(simplicial)
    # simplicial = jl.cubicalhomology.preprocess_graph(simplicial)
    # # graph = jl.cubicalhomology.simplicial_to_graph(simplicial,2,1)
    # # graph = jl.cubicalhomology.relabel_vertices(graph)
    # # print(jl.cubicalhomology.cubicalHomology(jl.cubicalhomology.preprocess_graph(graph),1))
    # print(jl.cubicalhomology.cubicalHomology(simplicial,1))
    num_nodes = 3
    wiring_diagram_list= [(i,i+1) for i in range(num_nodes-1)]+ [(num_nodes//2,0)] + [(num_nodes-1,num_nodes//2)]
    print(wiring_diagram_list)
    wiring_diagram_list = [(0,1),(1,0),(0,2)]
    boolean_network = rn.boolnet(num_nodes=num_nodes,seed = 0, edges_wiring_diagram = wiring_diagram_list)
    input_data = [
            ("A", "A and B"),
            ("B", "not A or C")
        ]

    net1= BooleanNet(input_type="text-format")
    print(net1.input_type)
    net1.load(input_data)
    print(wiring_diagram(boolean_network))
    nodes,edges = strongly_connected_components(wiring_diagram_list)
    print(nodes)
    print(edges)
    print(wiring_diagram(net1))
    print(network_scc(boolean_network))
    print(network_scc(net1))
