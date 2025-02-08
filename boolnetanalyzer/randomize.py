
# This file contains the different kinds of network generator
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import random
import sympy
import pandas as pd
from collections import defaultdict
from matplotlib import colors as mcolors

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, script_path)
#import BooleanNet Class from load.py
from boolnetanalyzer.load import BooleanNet
import boolnetanalyzer.bn_tools as bnt
import boolnetanalyzer.utils as utils
import boolnetanalyzer.load as load
def boolnet(num_nodes, indegree_distribution='uniform', dist_parameter=None, indegree_vector=None, 
            STRONGLY_CONNECTED=True, NO_SELF_REGULATION=True, edges_wiring_diagram=None, bias=0.5, seed=None):
    """
    Generates a random boolean network. The default is a strongly connected network with uniform in-degree distribution and no self loops.
    

    Parameters:
    ----------
        num_nodes (int): 
            Number of nodes in the network.

        indegree_distribution (str): 
            Distribution type for in-degrees (`constant`, `dirac`, `delta`, `uniform`, `poisson`).
            The default is `constant`. Note that `dirac`, `delta` and `constant` are equivalent and we allow this for convenience.
            Default is `uniform`.


        dist_parameter (int/float): 
            Parameter for the in-degree distribution.
            For `constant`/`dirac`/`delta`, this is the in-degree of every node.
            For `poisson` distribution, it is the parameter `lambda`.
            For `uniform` distribution, it is the maximum in-degree.
            Default is set to num_nodes/num_nodes - 1 depending on the `NO_SELF_REGULATION` flag.


        indegree_vector(list/np.array, optional): 
            User-defined in-degree vector.
            This is the vector of in-degrees for each node in the network. Default is None.

        STRONGLY_CONNECTED (bool): 
            Ensure the network is strongly connected.
            Default is True.

        NO_SELF_REGULATION (bool): 
            If True, no node will have a self-regulatory edge. Default is True.

        edges_wiring_diagram (list, optional): 
            Predefined edge list. Default is None.
            The edge list is a list of tuples (source, target).

        bias (float):
             Bias for the Boolean functions i.e number of ones vs zeroes. This is the probability of having a 1 in the Boolean function. Default is 0.5.

    Returns:
    --------
        BooleanNet: A custom class representing the Boolean network. The network will be in reduced truth tables format.

    Examples:
    --------
        >>> from booleantools import randomize as rn
        >>> net = rn.boolnet(10)
        >>> print(net.truth_tables,net.regulators)
    """
    if seed is not None:
        np.random.seed(seed)
    if dist_parameter is None:
        dist_parameter = num_nodes-int(NO_SELF_REGULATION)
    if edges_wiring_diagram is None:
        if utils.generate_indegrees(num_nodes, indegree_distribution, dist_parameter, indegree_vector, NO_SELF_REGULATION) is None:
            return
        else:
            while True:
                indegree_vector_processed = utils.generate_indegrees(num_nodes, indegree_distribution, dist_parameter, indegree_vector, NO_SELF_REGULATION)
                
                edges_wiring_diagram = edge_list(num_nodes, indegree_vector_processed, NO_SELF_REGULATION)
                if STRONGLY_CONNECTED:
                    G = nx.from_edgelist(edges_wiring_diagram, create_using=nx.MultiDiGraph())
                    strongly_connected = nx.is_strongly_connected(G)
                    if not strongly_connected:
                        continue
                break
    else:
        try:
            assert len(set(np.array(edges_wiring_diagram).flatten())) == num_nodes
        except AssertionError:
            print("number of nodes provided in edges_wiring_diagram != N")
            return
        indegree_vector_processed = np.zeros(num_nodes, dtype=int)
        for target in np.array(edges_wiring_diagram)[:, 1]:
            indegree_vector_processed[target] += 1

    truth_tables = [nondegeneratefunction(indegree_vector_processed[i], bias) for i in range(num_nodes)]
    regulators = [[] for _ in range(num_nodes)]
    for edge in edges_wiring_diagram:
        regulators[edge[1]].append(edge[0])
    for i in range(num_nodes):
        regulators[i] = np.sort(regulators[i])
    net=BooleanNet(input_type="rtt-format",num_nodes = num_nodes)
    net.load((truth_tables, regulators))  
    return net


def edge_list(num_nodes, indegree_vector, NO_SELF_REGULATION=True, NO_SINK=False, seed=None):
    """
    Generates a random edge list for a network.

    Parameters:
    ----------
        num_nodes (int): 
            Number of nodes in the network.

        indegree_vector (np.array):    
            Vector of in-degrees.

        NO_SELF_REGULATION (bool): 
            Prevent self-regulation.

        NO_SINK (bool): 
            Ensure there are no sinks in the network. In other words, every node has an outdegree > 0.

    Returns:
    --------
        list: 
            Edge list representing the network.
            The edge list is a list of tuples (source, target).

    Examples:
    --------
        >>> from booleantools import randomize as rn
        >>> edges = rn.edge_list(5, np.array([2, 2, 2, 2, 2]), NO_SELF_REGULATION=True)
        >>> print(edges)
        [(1, 0), (3, 0), (2, 1), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (3, 4), (1, 4)]

        >>> edges = rn.edge_list(5, np.array([2, 2, 2, 2, 2]), NO_SELF_REGULATION=False, NO_SINK=True)
        >>> print(edges)
        [(1, 0), (3, 0), (2, 1), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (3, 4), (1, 4)]
    """
    if seed is not None:
        np.random.seed(seed)
        
    if not NO_SINK:
        edge_list = []
        for i in range(num_nodes):
            if NO_SELF_REGULATION:
                # exclude self-regulation i.e i -> i
                indices = np.random.choice(np.append(np.arange(i), np.arange(i+1, num_nodes)), indegree_vector[i], replace=False)
            else:
                indices = np.random.choice(np.arange(num_nodes), indegree_vector[i], replace=False)
            edge_list.extend(list(zip(indices, i * np.ones(indegree_vector[i], dtype=int))))
    else:
        
        edge_list = []
        outdegree = np.zeros(num_nodes, dtype=int)
        sum_indegree = sum(indegree_vector)  # total number of regulations

        for i in range(num_nodes):
            if NO_SELF_REGULATION:
                # Choose indices excluding self (node i)
                indices = np.random.choice(np.append(np.arange(i), np.arange(i+1, num_nodes)), indegree_vector[i], replace=False)
            else:
                # Choose indices potentially including self (node i)
                indices = np.random.choice(np.arange(num_nodes), indegree_vector[i], replace=False)
            
            # Update outdegree for chosen indices
            outdegree[indices] += 1
            # Add edges to edge list
            edge_list.extend(list(zip(indices, i * np.ones(indegree_vector[i], dtype=int))))

        # Ensure no node has outdegree of 0
        while min(outdegree) == 0:
            index_sink = np.where(outdegree == 0)[0][0]  # Node with outdegree 0
            index_edge = np.random.randint(sum_indegree)  # Random edge index

            if NO_SELF_REGULATION:
                # Ensure the chosen edge doesn't point to index_sink
                while edge_list[index_edge][1] == index_sink:
                    index_edge = np.random.randint(sum_indegree)
            
            # Update outdegrees
            outdegree[index_sink] += 1
            outdegree[edge_list[index_edge][0]] -= 1
            
            # Modify edge to have index_sink as new source node
            edge_list[index_edge] = (index_sink, edge_list[index_edge][1])
        
    return edge_list



def nondegeneratefunction(num_nodes, probability_one=0.5, seed=None):
    """
    Generates a random non-degenerated Boolean function.

    Parameters:
    ----------
        num_nodes (int): 
            Number of variables.
        probability_one (float): 
            Probability of having a 1 in the Boolean function.

    Returns:
    -------
        np.array: 
            Boolean function as a vector.

    Examples:
    --------
        >>> from booleantools import randomize as rn
        >>> f = rn.nondegeneratefunction(3)
    """
    if seed is not None:
        np.random.seed(seed)
    while True:
        f = np.array(np.random.random(2**num_nodes) < probability_one, dtype=int)
        
        if not bnt.is_degenerated(f):
            return f
        # Change the seed to get a different function
        if seed is not None:
            seed += 1
            np.random.seed(seed)
   
def linearboolnet(num_nodes,wiring_diagram):
    """
    Generates a linear Boolean network from a given wiring diagram.

    Parameters:
    ----------
        num_nodes (int): 
            Number of nodes in the network.

        wiring_diagram (list): 
            Wiring diagram of the network.
            The wiring diagram is a list of tuples (source, target).

    Returns:
    -------
        BooleanNet: 
            A custom class representing the Boolean network. The network will be in reduced truth tables format.

    Examples:
    --------
        >>> from booleantools import randomize as rn
        >>> edges = rn.edge_list(5, np.array([2, 2, 2, 2, 2]), NO_SELF_REGULATION=True)
        >>> net = rn.linearboolnet(5, edges)
        >>> print(net.truth_tables,net.regulators)
    """
    regulators = [[] for _ in range(num_nodes)]
    for edge in wiring_diagram:
        regulators[edge[1]].append(edge[0])
    truth_tables= []
    for i in range(num_nodes):
        #truth_table = np.zeros(2**len(regulators[i]), dtype=int)
        # Generate the left hand side of the truth table
        possible_values = list(itertools.product([0, 1], repeat=len(regulators[i])))
        truth_table = [np.sum(values)%2 for values in possible_values]
        truth_tables.append(truth_table)
    #print(truth_tables,regulators)
    net = BooleanNet(input_type="rtt-format",num_nodes=num_nodes)
    net.load((truth_tables, regulators))  
    return net
def booleanfunction(num_nodes, probability_one=0.5, seed=None):
    """
    Generates a random Boolean function.

    Parameters:
    ----------
        num_nodes (int): 
            Number of variables.
        probability_one (float): 
            Probability of having a 1 in the Boolean function.
        seed (int):
            Seed for the random number

    Returns:
    -------
        np.array: 
            Boolean function as a vector.

    Examples:
    --------
        >>> from booleantools import randomize as rn
        >>> f = rn.booleanfunction(3)
    """
    if seed is not None:
        np.random.seed(seed)
    return np.array(np.random.random(2**num_nodes) < probability_one, dtype=int)
if __name__ == "__main__":

    
    f = booleanfunction(3, seed=2)
    print(f)
    # >>> [1 1 0 1 1 1 1 0]
    f = nondegeneratefunction(3,seed=2)
    print(f)
    # >>> [1 1 0 1 1 1 1 0]
    
    net = boolnet(4,seed=1) # The seed forces the same network to be generated # Use this only when you need to debug
    print(net.truth_tables,net.regulators)
    
    # >>> [array([1, 0, 0, 0]), array([1, 0]), array([0, 1]), array([0, 0, 1, 0])] [array([2, 3]), array([3]), array([0]), array([0, 1])]
    
    N = 5
    indegree_vector = np.array([2, 2, 2, 2, 2])
    print(edge_list(N, indegree_vector, NO_SELF_REGULATION = False, NO_SINK=True, seed=1))
    # >>> [(2, 0), (1, 0), (0, 1), (2, 1), (4, 2), (3, 2), (3, 3), (0, 3), (3, 4), (0, 4)]

    num_nodes = 5
    edges = [(1, 0), (3, 0), (2, 1), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (3, 4), (1, 4),(4,0)]
    net = linearboolnet(num_nodes, edges)
    print(net.truth_tables,net.regulators)

    
    