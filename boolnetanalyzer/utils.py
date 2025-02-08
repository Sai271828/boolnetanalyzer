# This is a utils file that will conatins helper functions
import numpy as np
def generate_indegrees(num_nodes, distribution="constant", dist_parameter=2, indegree_vector=None, NO_SELF_REGULATION=False, seed=None):
    """
    Generate a degree distribution for a network of a given size. 
    The function generates a vector of in-degrees for each node in the network.
    The function can generate a constant, uniform, Poisson, Dirac or user-defined degree.
    The function can also generate a degree distribution with or without self-regulation.
    Self-regulation is defined as a node having an edge to itself.
    The seed parameter is used to ensure reproducibility of the generated degree distribution.

    Parameters
    ----------
        num_nodes : int
            The number of nodes in the network.
        distribution : str
            The type of degree distribution to generate. Options are `constant`, `uniform`, `poisson`, `dirac`, and `delta`. Note that `dirac`, `delta` and `constant` are equivalent.
            
        dist_parameter : int or float
            The parameter of the degree distribution. For `constant`/`dirac`/`delta` this is the in-degree of every node. 
            For `poisson` distribution it is the parameter `lambda`.
            For `uniform` distribution it is the maximum indegree.

        indegree_vector : list or numpy array. 
            This is the vector of in-degrees for each node in the network. Default is None.

        NO_SELF_REGULATION : bool. 
            If True, no node will have a self-regulatory edge. Default is False.
        seed : int
            Seed for reproducibility. Default is None.

    Returns
    -------
        np.array
            An array of integers representing the degree of each node in the network.

    Examples
    --------    
    >>> generate_indegrees(10, 'poisson', 2)
    [1 3 2 1 4 1 1 2 1 1]
    >>> generate_indegrees(10, 'uniform', 3)
    [1 2 1 2 1 1 3 3 3 1]
    >>> generate_indegrees(10, 'constant', 2)
    [2 2 2 2 2 2 2 2 2 2]
    >>> generate_indegrees(10,indegree_vector=[1, 2, 1, 2, 1, 1, 3, 3, 3, 1])
    [1 2 1 2 1 1 3 3 3 1]
    """
    if indegree_vector is not None:
        if type(indegree_vector) in [list, np.array]:
            try:
                assert type(indegree_vector) in [list, np.array]
                assert np.all([type(el) in [int, np.int64] for el in indegree_vector])
                assert len(indegree_vector) == num_nodes
                assert min(indegree_vector) >= 1
                assert max(indegree_vector) <= num_nodes
                if type(indegree_vector) == list:    
                    indegree_vector = np.array(indegree_vector)     
                return indegree_vector
            except AssertionError:
                print('Error: A vector indegree_vector was submitted.\nEnsure that indegree_vector is an N-dimensional vector where each element of indegree_vector is an integer between 1 and N.')
                return
        else:
            print("Here")
            print('Error: indegree_vector must be a list or numpy array.')
            return
    else:
        if seed is not None:
            np.random.seed(seed)

        if distribution in ['constant', 'dirac', 'delta']:
            try:
                assert type(dist_parameter) in [int, np.int64]
                assert dist_parameter >= 1
                assert dist_parameter <= num_nodes
                return  np.ones(num_nodes, dtype=int) * dist_parameter
            except AssertionError:
                print('Error: n must be a single integer between 1 and N when using a constant degree distribution.')
                return
        elif distribution == 'uniform':
            try:
                assert type(dist_parameter) in [int, np.int64]
                assert dist_parameter >= 1
                assert dist_parameter <= num_nodes
                return 1 + np.random.randint(dist_parameter, size=num_nodes)
            except AssertionError:
                print('Error: n must be a single integer between 1 and N representing the upper bound of a uniform degree distribution (lower bound==1).')
                return
        elif distribution == 'poisson':
            try:
                assert type(dist_parameter) in [int, np.int64, float, np.float64]
                assert dist_parameter >= 1
                #assert dist_parameter <= num_nodes
                indegree_vector = np.random.poisson(lam = dist_parameter , size = num_nodes)
                # Ensure that no node has an in-degree of 0 or more than N-1
                indegree_vector[indegree_vector==0] = 1
                indegree_vector[indegree_vector>num_nodes-int(NO_SELF_REGULATION)] = num_nodes-int(NO_SELF_REGULATION)
                return indegree_vector
            except AssertionError:
                print('Error: n must be a single number > 0 representing the Poisson parameter.')
                return
        else:
            print('None of the predefined indegree distributions were chosen.\nTo use a user-defined in-degree vector, use the input indegree_vector to submit an N-dimensional vector where each element of n must be between 1 and N')
            return


def bin2dec(binary_vector):
        """
        Convert a binary vector to an integer.

        Parameters:
        binary_vector (list): List containing binary digits (0 or 1).

        Returns:
        int: Integer value converted from the binary vector.
        """
        binary_string = ''.join(str(bit) for bit in binary_vector)
        return int(binary_string, 2)


def dec2bin(integer_value, num_bits):
    """
    Convert an integer to a binary vector.

    Parameters:
    integer_value (int): Integer value to be converted.
    num_bits (int): Number of bits in the binary representation.

    Returns:
    list: List containing binary digits (0 or 1).
    """
    binary_string = bin(integer_value)[2:].zfill(num_bits)
    return [int(bit) for bit in binary_string]
# main function to test the indegree_distribution function
if __name__ == "__main__":
    num_nodes = 10
    distribution = 'constant'
    dist_parameter = 2
    indegree_vector = generate_indegrees(10, 'uniform', 3, seed = 0)
    print(indegree_vector)
    # >>> [1 2 1 2 2 3 1 3 1 1]
 