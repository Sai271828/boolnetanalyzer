# This is used for loading boolean networks
import numpy as np
import re
import pandas as pd
import warnings
# load the environmental variable from the .env file
from dotenv import load_dotenv
import os
import sys
load_dotenv()
script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, script_path)

# Function to parse text file
def parse_text_file(file_path, delimiter="=", char_not="NOT", char_and="AND", char_or="OR"):
    """
    Parses a text file containing node-formula pairs, processes the formulas by replacing specified logical operators 
    (e.g., "NOT", "AND", "OR") with their Python equivalents, and returns a list of tuples with the node and updated formula.

    Parameters:
    ----------`
        file_path (str): The path to the text file that needs to be parsed.
        delimiter (str, optional): The character used to separate the node and the formula in each line of the file.
                                   Defaults to "=".
        char_not (str, optional): The string representing logical negation in the file. It will be replaced with "not" in Python.
                                  Defaults to "NOT".
        char_and (str, optional): The string representing logical conjunction in the file. It will be replaced with "and" in Python.
                                  Defaults to "AND".
        char_or (str, optional): The string representing logical disjunction in the file. It will be replaced with "or" in Python.
                                 Defaults to "OR".

    Returns:
    -------
        list: A list of tuples where each tuple contains a node (str) and its corresponding updated formula (str).
    
    Example:
    --------
        Input file:
        ```
        node1 = NOT A AND B
        node2 = A OR B
        ```

        Output:
        [('node1', 'not A and B'), ('node2', 'A or B')]
    """
    
    update_formulae = []

    with open(file_path, 'r') as file:
        data = file.readlines()
        for line in data:
            line = line.strip()
            if line:
                node, formula = line.split(delimiter)
                node = node.strip()
                formula = formula.strip()
                formula = formula.replace(char_not, "not")
                formula = formula.replace(char_and, "and")
                formula = formula.replace(char_or, "or")
                update_formulae.append((node, formula))
    
    return update_formulae


# Define the class BooleanNet
class BooleanNet():
    def __init__(self, input_type,num_nodes=1):
        """
        Initializes the class instance with a specified number of nodes and an input format type.
        
        Parameters
        ----------
        num_nodes : int
            The number of nodes (variables) in the network.
        
        input_type : str
            The format of the input data. It must be one of the following:
            - "NCF-format": Format specific to NCF.
            - "rtt-format": Format that uses truth tables.
            - "text-format": Format that uses text-based formulas.
        
        Raises
        ------
        ValueError
            If the provided `input_type` is not one of the accepted formats, a ValueError is raised.
        """
        
        # Store the number of nodes in the network.
        self.num_nodes = num_nodes
        # Store the format type for input data (e.g., NCF, rtt, or text format).
        self.input_type = input_type

        # Validate that the input type is one of the supported formats.
        if self.input_type != "NCF-format" and self.input_type != "rtt-format" and self.input_type != "text-format":
            raise ValueError("Invalid input type")
        
        # If the input type is "NCF-format", no specific initialization is needed yet.
        if self.input_type == "NCF-format":
            pass  # Future setup for NCF format can go here.
        
        # If the input type is "rtt-format", initialize data structures for truth tables and regulators.
        elif self.input_type == "rtt-format":
            self.truth_tables = []  # List to hold the reduced truth tables for each node.
            self.regulators = []  # List to hold the regulators for each node.
        
        # If the input type is "text-format", initialize data structures for formulas, variables, and regulators.
        elif self.input_type == "text-format":
            self.update_formulae = []  # List to store the update formulas for each node.
            self.variables = []  # List to store the variable names.
            self.regulators = []  # List to store the regulators for each node.
            self.original_variables = {}  # Dictionary to map original variable names.
        
        # If somehow an invalid input type is provided, raise a ValueError (extra check for safety).
        else:
            raise ValueError("Invalid input type")

    def replace_in_variables(self, symbol, replacement):
        """
        Many biological networks have variable names that contain special characters like semi-colons, forward slashes, etc.
        This causes many issues with Python code. This function replaces the special characters in the variable names with a placeholder text.

        Parameters:
        ----------
        symbol (str): The special character to be replaced in the variable names.
        replacement (str): The placeholder text to replace the special character with.

        Returns:
        --------
        None

        Modifies internal attributes:
        - `original_variables`: Dictionary mapping original variable names to cleaned names.
        - `variables`: List of variable names with the special character replaced by the placeholder text.
        - `regulators`: List of lists of regulator names with the special character replaced by the placeholder text.
        - `update_formulae`: List of tuples where each tuple consists of a variable name and its corresponding formula, with the special character replaced by the placeholder text.

        """
        # Replace in original_variables
        for variable in self.original_variables:
            self.original_variables[variable] = self.original_variables[variable].replace(symbol,replacement)

        # Replace in variables
        self.variables = [variable.replace(symbol, replacement) for variable in self.variables]

        # Replace in regulators
        self.regulators = [[variable.replace(symbol, replacement) for variable in regulators] for regulators in self.regulators]

        # Replace in update_formulae
        self.update_formulae = [(variable.replace(symbol, replacement), formula.replace(symbol, replacement)) for variable, formula in self.update_formulae]
        
         

        
    ## Define the clean variables function
    ## In biological networks, the variables are usually represented with very wierd names such APC/C, APC/C:cdc20, etc. 
    ## This can cause a lot of complications in the code. So, we need to clean the variables and store a dictionary that remebers the original variable names
    def clean_variables(self):
        
        """
        Clean and standardize variable names and update formulae in the Boolean network.
        The method is specifically designed to handle variables and formulae in the "text-format".

        This method performs the following operations:
            - Replaces semi-colons (;) in variable names and update formulae with the placeholder text "_sc_".
            - Replaces forward slashes (/) in variable names and update formulae with the placeholder text "_fs_".
            - Replaces dashes (-) in variable names and update formulae with the placeholder text "_dash_".
            - Replaces plus signs (+) in variable names and update formulae with the placeholder text "_plus_".
            - Replaces periods (.) in variable names and update formulae with the placeholder text "_dot_".
            - Updates the `original_variables` dictionary to reflect these changes.

        

        Parameters:
        -----------
            None

        Returns:
        --------
            None

        Modifies internal attributes:
            - `original_variables`: Dictionary mapping original variable names to cleaned names.
            - `variables`: List of variable names with semi-colons and forward slashes replaced by placeholders.
            - `regulators`: List of lists of regulator names with semi-colons and forward slashes replaced by placeholders.
            - `update_formulae`: List of tuples where each tuple consists of a variable name and its corresponding formula, with semi-colons and forward slashes replaced by placeholders.

        Examples:
        ---------
        >>> network = BooleanNet("text-format",3)
        >>> network.variables = ['var1;', 'var2/']
        >>> network.regulators = [['var1;'], ['var1;', 'var2/']]
        >>> network.update_formulae = [('var1;', 'var1; AND var2/')]
        >>> network.clean_variables()
        >>> print(network.variables)
        ['var1_sc_', 'var2_fs_']
        >>> print(network.regulators)
        [['var1_sc_'], ['var1_sc_', 'var2_fs_']]
        >>> print(network.update_formulae)
        [('var1_sc_', 'var1_sc_ AND var2_fs_')]

        """
        if self.input_type == "text-format":
            # replace the semi-colon in the variables and regulators with the text _sc_
            self.replace_in_variables(";", "_sc_")
            # replace the forward slash in the variables and regulators with the text _fs_
            self.replace_in_variables("/", "_fs_")
            # replace the - in the variables and regulators with the text _dash_
            self.replace_in_variables("-", "_dash_")
            # replace the + in the variables and regulators with the text _plus_
            self.replace_in_variables("+", "_plus_")
            # replace the . in the variables and regulators with the text _dot_
            self.replace_in_variables(".", "_dot_")
            
        
        
    ## Define the load function
    def load(self, input_data):
        """
        Loads the input data into the class instance based on the specified input format type.

        This method determines the appropriate *internal* loading function (`load_ncf`, `load_rtt`, or `load_text`) 
        based on the format type set during initialization. The format must match the input type provided 
        when the class was initialized, otherwise a `ValueError` will be raised.

        Parameters
        ----------
        input_data : varies
            The input data to load into the class. The structure of the input depends on the format type:

            - **"NCF-format"**: Data structured according to the NCF format. Currently, this format is **not implemented** in the class.

            - **"rtt-format"**: Data structured using **reduced truth tables**.
                - This is a tuple consisting of two lists:
                    1. The **first list** contains the **reduced truth tables** for each node. Each entry in this list represents the truth table for the corresponding node.
                    2. The **second list** contains the **regulators** for each node. Each entry represents the nodes that regulate the output of the corresponding node.
                - Example:
                    
                    ```python
                    ([[1, 1, 0, 0], [0, 1], [0, 0, 0, 0]], [[0, 2], [0], [1, 2]])
                    ```
                    
                - In this case, the **truth tables** are stored in the variable `self.truth_tables` and the **regulators** in the variable `self.regulators`.

            - **"text-format"**: Data structured as **text-based formulas**.
                - This is a list of tuples where each tuple consists of:
                    1. The **first string** representing the node (variable) name.
                    2. The **second string** representing the **update formula** associated with that node.
                - Example:
                    
                    ```python
                    [("x1", "x2 and not x1"), ("x2", "x1 or x3")]
                    ```
                
                - In this case, the **variables** are stored in `self.variables`, the **regulators** are stored in `self.regulators` and the **update formulas** in `self.update_formulae`.
                .. note::

                    The constants found in the update formulas but not initially listed as nodes are added to both `self.variables` and `self.update_formulae`.

        Raises
        ------
        ValueError
            If the input format type is not one of the expected types ("NCF-format", "rtt-format", or "text-format"),
            a `ValueError` is raised.

        Examples:
        ---------
        ```python
        >>> network = BooleanNet("rtt-format", 3)
        >>> network.load(([[1, 1, 0, 0], [0, 1], [0, 0, 0, 0]], [[0, 2], [0], [1, 2]]))
        >>> print(network.truth_tables), print(network.regulators)
        [[1, 1, 0, 0], [0, 1], [0, 0, 0, 0]]
        [[0, 2], [0], [1, 2]]
        >>> network = BooleanNet("text-format", 2)
        >>> network.load([("A", "A and B"), ("B", "not A")])
        >>> print(network.update_formulae)
        [('A', 'A and B'), ('B', 'not A')]
        """
                
        if self.input_type == "NCF-format":
            self.load_ncf(input_data)
        elif self.input_type == "rtt-format":
            self.load_rtt(input_data)
        elif self.input_type == "text-format":
            self.load_text(input_data)
        else:
            raise ValueError("Invalid input type")

        
    ## Define the update function
    def update(self,state):
        """
        Update the state of the network.
        
        Parameters:
        ----------
        state (list): Current state of the network.
        
        Returns:
        ------- 
        list: Updated state of the network.
        
        Examples:
        --------
        
        >>> network = BooleanNet("rtt-format",3)
        >>> network.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> state = [1, 0, 1]
        >>> network.update(state)
        [0, 0, 0]
        
        """
        # update the state of the network
        if self.input_type == "NCF-format":
            return self.update_ncf(state)
        elif self.input_type == "rtt-format":
            return self.update_rtt(state)
        elif self.input_type == "text-format":
            return self.update_text(state)
        else:
            raise ValueError("Invalid input type")
        pass
    
        ## Define the add function
    def add(self, network2):
        """
        Add two boolean networks together.

        Parameters:
        ----------
        network2 (BooleanNet): The boolean network to be added to the current network.

        Output:
        -------
        The current network is updated with the contents of the second network merged into it.

        Examples:
        ---------
        >>> network1 = BooleanNet("rtt-format",3)
        >>> network1.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> network2 = BooleanNet("rtt-format",2)
        >>> network2.load(([[1,0],[1,1]],[[0],[0]]))
        >>> network1.add(network2)
        >>> print(network1.truth_tables), print(network1.regulators)
        [[1, 1, 0, 0], [0, 1], [0, 0, 0, 0], [1, 0], [1, 1]]
        [[0, 2], [0], [1, 2], [3], [3]]
        """
        if self.input_type == "NCF-format":
            self.add_ncf(network2)
        elif self.input_type == "rtt-format":
            self.add_rtt(network2)
        elif self.input_type == "text-format":
            self.add_text(network2)
        else:
            raise ValueError("Invalid input type")


    ## Define the add_node function
    def add_node(self, node):
        """
        Add a new node to the boolean network.

        Parameters:
        ----------
        node (int): The node to be added to the network.

        Output:
        -------
        Updates the network by adding the specified node.

        Examples:
        ---------
        >>> network = BooleanNet("rtt-format",3)
        >>> network.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> network.add_node(4)
        >>> print(network.variables)
        [0, 1, 2, 3, 4]
        """
        if self.input_type == "NCF-format":
            self.add_node_ncf(node)
        elif self.input_type == "rtt-format":
            self.add_node_rtt(node)
        elif self.input_type == "text-format":
            self.add_node_text(node)
        else:
            raise ValueError("Invalid input type")


    ## Define the add_edge function
    def add_edge(self, edge):
        """
        Add a new edge to the boolean network.

        Parameters:
        ----------
        edge (tuple): A tuple where the first element is the node and the second element is a list of input nodes.

        Output:
        -------
        Updates the network by adding the specified edge.

        Examples:
        ---------
        >>> network = BooleanNet("rtt-format",3)
        >>> network.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> network.add_edge((3, [1, 2]))
        >>> print(network.regulators)
        [[0, 2], [0], [1, 2], [1, 2]]
        """
        if self.input_type == "NCF-format":
            self.add_edge_ncf(edge)
        elif self.input_type == "rtt-format":
            self.add_edge_rtt(edge)
        elif self.input_type == "text-format":
            self.add_edge_text(edge)
        else:
            raise ValueError("Invalid input type")


    ## Define the change_format function
    def change_format(self, new_format):
        """
        Change the format of the boolean network.

        Parameters:
        ----------
        new_format (str): The new format to convert the network to. Options are "NCF-format", "rtt-format", and "text-format".

        Output:
        -------
        Updates the network to the new format.

        .. note::
            Currently only changing from "text-format" to "rtt-format" is implemented. 
            In this case, the method works well as long as the degree of the nodes is less than or equal to 15. 
            If the degree is greater than 15, a warning is issued. To modify this, change the max_degree parameter in the `change_format_text` method.

        Examples:
        ---------
        >>> network = BooleanNet("rtt-format",3)
        >>> network.change_format("text-format")
        >>> print(network.input_type)
        text-format
        """
        if self.input_type == "NCF-format":
            self.change_format_ncf(new_format)
        elif self.input_type == "rtt-format":
            self.change_format_rtt(new_format)
        elif self.input_type == "text-format":
            self.change_format_text(new_format)
        else:
            raise ValueError("Invalid input type")

    ## Define the essential_regulators function
    def essential_regulators(self):
        """
        Identify the essential regulators in the network.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Modifies internal attribute:
            - `regulators`: A list of essential regulators for each node in the network.

        .. note::
            Code to correspondingly modify the truth tables needs to be added.
            

        Examples:
        ---------
        >>> network = BooleanNet("rtt-format",3)
        >>> network.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> print(network.essential_regulators())
        [[0, 2], [0], [1, 2]]
        """
        import booleantools.bn_tools as bnt
        self.regulators = [bnt.essential_variables(np.array(truth_table), regulators) for truth_table, regulators in zip(self.truth_tables, self.regulators)]
        ## Need to edit the truth tables as well!
    

        


##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
    ## Functions for the NCF format
    ## These have not been implemented yet
    def load_ncf(self, input_data):
        pass

    def update_ncf(self,state):
        pass

    def add_ncf(self, network2):
        pass

    def add_node_ncf(self, node):
        pass
    
    def add_edge_ncf(self, edge):
        pass
    
    def change_format_ncf(self, new_format):
        pass

##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
    ## Functions for the text format
    def load_text(self, input_data):
        """
        Loads a network in text format, where the network is represented as a list of tuples. Each tuple consists of a node 
        and its corresponding update formula. The method also calculates the network variables and their regulators.

        Parameters
        ----------
        input_data : list of tuple of str
            A list where each tuple contains two strings:
            - The first string represents the node (variable) name.
            - The second string represents the update formula associated with that node.

        Raises
        ------
        ValueError
            If `input_data` is not a list of tuples, a ValueError is raised.

        Notes
        -----
            - After loading the data, the method computes two main properties:
                - `self.variables`: A list of all node names.
                - `self.regulators`: A list of the regulators for each variable, derived from the update formulas by ignoring logical 
                operators ('and', 'or', 'not').
            - The method also handles any constants found in the update formulas but not initially listed as nodes, 
            adding them to both `self.variables` and `self.update_formulae`.
            - It replaces semicolons and forwardslashes in variable names with '_sc_' via `self.clean_variables()`.
        
        Example
        -------
        input_data = [
            ("node1", "A and B"),
            ("node2", "not A or C")
        ]
        This method would extract the variables and regulators for these nodes.
        """
        
        if isinstance(input_data, list):
            if all(isinstance(item, tuple) for item in input_data):
                self.update_formulae = input_data
                # Calculate variables (node names)
                self.variables = [item[0] for item in input_data]
                # Calculate regulators from update formulas
                self.regulators = []
                logical_operators = {'and', 'or', 'not'}
                for item in input_data:
                    # Extract variables from the update formula, excluding logical operators
                    right_variables = re.findall(r'[^()\s]+', item[1].strip())
                    self.regulators.append(list(set([element for element in right_variables if element not in logical_operators])))
                
                for regulators in self.regulators:
                    for variable in regulators:
                        if variable not in self.variables:
                            self.variables.append(variable)
                            self.regulators.append([variable])
                            self.update_formulae.append((variable, variable))
                for variable in self.variables:
                    self.original_variables[variable] = variable
                
                # Replace special characters in variable names and regulators
                self.clean_variables()
                self.num_nodes = len(self.variables)
            else:
                raise ValueError("Invalid input format. Input data should be a list of tuples")
        pass


    def update_text(self,state):
        pass

    def add_text(self, network2):
        pass

    def add_node_text(self, node):
        pass

    def add_edge_text(self, edge):
        pass

    def change_format_text(self, new_format,max_degree=15):
        if new_format == "NCF-format":
            pass  
        elif new_format == "text-format":
            print("The network is already in text format")
        elif new_format == "rtt-format":
            import itertools
            ## Convert the text format to rtt format
            truth_tables = []
            regulators = []
            for i in range(len(self.update_formulae)):
                degree = len(self.regulators[i])
                regulator_indices = [self.variables.index(regulator) for regulator in self.regulators[i]]
                regulators.append(regulator_indices)
                f = np.array([],dtype=int)
                if degree<=15:
                    X = list(itertools.product([0, 1], repeat = degree))
                    for j in range(2**degree):
                        x = X[j]
                        
                        variable_values = {var: val for var, val in zip(self.regulators[i], x)}
                        
                        f = np.append(f,eval(self.update_formulae[i][1], {}, variable_values) % 2)
                else:
                    warnings.warn(f"The degree is {degree} which exceeds the maximum degree {max_degree}. Continuing but the truth table may not be accurate.", UserWarning)
                        
                truth_tables.append(f)
            self.truth_tables = truth_tables
            self.input_type = "rtt-format"
            self.regulators = regulators
            self.num_nodes = len(self.variables)
            
                

        else:
            raise ValueError("Invalid input type")
##############################################################################################################################################################
##############################################################################################################################################################
##############################################################################################################################################################
    ## Functions for the rtt format
    def load_rtt(self, input_data):
        if isinstance(input_data, tuple):
            if all(isinstance(item, list) for item in input_data):
                self.truth_tables = input_data[0]
                self.regulators = input_data[1]
            else:
                raise ValueError("Invalid input format. Input data should be a tuple of two lists")
            if len(self.truth_tables) != self.num_nodes:
                raise ValueError("The number of rows in truth_tables should be equal to the number of nodes ")
            if len(self.regulators) != self.num_nodes:
                raise ValueError("The number of rows in regulators should be equal to the number of nodes ")
            
            # Perform all checks in a single loop
            for f_row, i_row in zip(self.truth_tables, self.regulators):
                # Check that the length of each row in F matches 2**n where n is the length of the corresponding row in I
                if len(f_row) != 2**len(i_row):
                    raise ValueError("The length of each row in truth_tables should be 2**n where n is the length of the corresponding row in regulators")
                
                # Check that every row of F contains only 0 or 1
                if not all(element in [0, 1] for element in f_row):
                    raise ValueError("Each row in truth_tables should contain only 0 or 1 elements")
                
                # Check that every row of I contains elements from 0 to n-1 where n is the number of nodes
                if not all(element in range(self.num_nodes) for element in i_row):
                    raise ValueError("Each row in regulators should contain elements ranging from 0 to n-1 where n is the number of nodes")
                
                # Check that every row of I has unique elements
                if len(i_row) != len(set(i_row)):
                    raise ValueError("Each row in regulators should contain unique elements")

        else:
            raise ValueError("Invalid input format. Input data should be a tuple of two lists")
        # check that the input data is of the correct shape
    
    
    
    
    def update_rtt(self, state):
        """
        Update the state of the network in rtt format.

        Parameters:
        state (list): Current state of the network.

        Returns:
        list: Updated state of the network.
        """
        updated_state = [
            self.truth_tables[node][int(''.join(str(state[index]) for index in input_nodes), 2)]
            for node, input_nodes in enumerate(self.regulators)
        ]
        
        return updated_state
    
    
    
    
    def add_rtt(self, network2):
        """
        Add two boolean networks in rtt format.

        Parameters:
        network2 (BooleanNet): The network to be added to the current network.

        Output:
        The current network with the second network added to it inplace.
        """
        # Shift the indices in regulator networks2 by the number of nodes in the rttrst network
        I2 = [[index + self.num_nodes for index in input_nodes] for input_nodes in network2.regulators]
        self.truth_tables.extend(network2.truth_tables)
        self.regulators.extend(I2)
        self.num_nodes += network2.num_nodes

    

    
    def add_node_rtt(self, node):
        
        self.num_nodes += 1
        self.truth_tables.append([])
        self.regulators.append([])
        
    

    
    def add_edge_rtt(self, edge):
        """
        This function does not work as expected. It needs to be fixed.
        """
        
        node, input_node = edge
        self.truth_tables[node].append(0)
        self.regulators[node].append(input_node)
        """
        fix this code"""
        pass

    def add_canalizing_edge_rtt(self, edge, canalizing_value):
        """
        Add a canalizing edge to the boolean network in rtt format.
        Canalizing value is a tuple (a,b).
        Suppose that the old function is f(x_1, x_2, x_3), then the new function will be f(x_1, x_2, x_3) = b if x_1 = a and f(x_1, x_2, x_3)+b otherwise.

        Parameters:
        edge (tuple): A tuple where the first element is the node and the second element is the input node.
        canalizing_value (tuple): A tuple (a,b) where a is the canalizing value and b is the value of the function when the input node is equal to a.

        Output:
        Updates the network by adding the specified canalizing edge.

        Examples:
        ---------
        >>> network = BooleanNet("rtt-format",3)
        >>> network.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
        >>> network.add_canalizing_edge_rtt((3, 1), (0, 1))
        >>> print(network.truth_tables)
        [[1, 1, 0, 0], [0, 1], [0, 0, 0, 0], [1, 0, 1, 1]]"""

        node, regulator = edge
        a, b = canalizing_value

        # Update the regulators
        # Append the input_node to the front of the regulators list
        self.regulators[node].insert(0,regulator)
        self.truth_tables[node] = [(x + b) % 2 for x in self.truth_tables[node]]


        if a == 0:
            
            new_truth_table = [b] *  2**(len(self.regulators[node])-1)
            new_truth_table.extend(self.truth_tables[node])
            self.truth_tables[node] = new_truth_table

        else:
            new_truth_table = [b] *  2**(len(self.regulators[node])-1)
            self.truth_tables[node].extend(new_truth_table)
        # Calculate the new truth table


    def add_linear_edge_rtt(self, edge):

        node, regulator = edge
        # Update the regulators
        # Append the input_node to the front of the regulators list
        self.regulators[node].insert(0,regulator)
        # Update the truth table
        # Update the truth table by flipping each value (mod 2 operation) 
        new_truth_table = [(x + 1) % 2 for x in self.truth_tables[node]] 
        # Extend the existing truth table with the new truth table self
        self.truth_tables[node].extend(new_truth_table)

    

        
    def change_format_rtt(self, new_format):
        pass

    def compose_rtt(self,network2):
        """
        Compose two boolean networks in rtt format.

        Parameters:
        network2 (BooleanNet): The network to be composed with the current network.

        Output:
        The composed network. The network1 formulae will be substituted into the network2 formulae.
        """
        # Calculate the regulators of the composition.
        # If 0,1 regulate node 2 in network1 and 0,2 regulate node 3 in network 1. 
        # Further suppose that 2 and 3 in network2 regulate 0 in network 2.
        # Then 0,1,2 regulate 0 in the composed network.
        import itertools
        total_regulators = []
        composite_truth_tables = []
        for regulators in network2.regulators:
            all_regulators = []
            for regulator in regulators:
                all_regulators.extend(self.regulators[regulator])
            all_regulators = list(set(all_regulators))
            total_regulators.append(all_regulators)
        
        # Calculate the truth tables of the composition
        for node in range(self.num_nodes):
            degree = len(total_regulators[node])
            X = list(itertools.product([0, 1], repeat = degree))
            f = np.array([],dtype=int)
            if degree == 0:
                f = np.append(f, network2.truth_tables[node][0])
            else:
                for j in range(2**degree):
                    x = X[j]
                    variable_values = {var: val for var, val in zip(total_regulators[node], x)}
                    # Calculate the values of the regulators of network 2
                    # Using list comprehensions for concise code
                    
                    networks2_state = [
                        self.truth_tables[regulator][
                            int(''.join(map(str, [variable_values[var] for var in self.regulators[regulator]])), 2) 
                            if ''.join(map(str, [variable_values[var] for var in self.regulators[regulator]])) != '' 
                            else 0
                        ]
                        for regulator in network2.regulators[node]
                    ]




                    # Appending to f using np.append
                    f = np.append(f, network2.truth_tables[node][int(''.join(map(str, networks2_state)), 2)])
            composite_truth_tables.append(f)

        return composite_truth_tables, total_regulators




        
if __name__=="__main__":

    # net = BooleanNet(input_type="rtt-format",num_nodes=3)
    # print("Boolean network loaded successfully")
    # # >>> Boolean network loaded successfully
    # print("Number of nodes: ", net.num_nodes)
    # # >>> Number of nodes: 3
    # print("Input type: ", net.input_type)
    # # >>> Input type: rtt-format

    # net.load(([[0,0,1,1],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
    # net1 =  net
    # truth_tables, regulators = net.compose_rtt(net1)
    # print(truth_tables)
    # print(regulators)
    # import booleantools.bn_tools as bnt
    # print(bnt.essential_truth_tables(truth_tables,regulators))
    # net2 = BooleanNet(input_type="rtt-format",num_nodes=3)
    # net2.load(bnt.essential_truth_tables(truth_tables,regulators))
    # print(net2.truth_tables)
    # print(net2.regulators)
    # net2.compose_rtt(net)
    # print(net2.truth_tables)
    # print(net2.regulators)
    # net2.compose_rtt(net)
    # print(net2.truth_tables)
    # print(net2.regulators)
    # input("Press Enter to continue...")
    # # Load the boolean network
    # # Loading a boolean network in text format
    # net = BooleanNet(input_type="text-format")
    # print("Boolean network initialized successfully")
    # # >>> Boolean network initialized successfully
    # print("Number of nodes: ", net.num_nodes) # Should print the default number of nodes. Will be updated to the correct number after loading the network.
    # # >>> Number of nodes: 1
    # print("Input type: ", net.input_type)
    # # >>> Input type: text-format
    # input_data = [
    #         ("node1", "A and B"),
    #         ("node2", "not A or C")
    #     ]
    
    # net.load(input_data)
    # # Print the network
    # print(net.update_formulae)
    # # >>> [('node1', 'A and B'), ('node2', 'not A or C')]
    # print(net.variables)
    # # >>> ['node1', 'node2', 'A', 'B', 'C']
    # print(net.regulators)
    # # >>> [['A', 'B'], ['A', 'C'], ['A'], ['B'], ['C']]
    # print(net.num_nodes)
    # # >>> 5
    # print(net.original_variables)
    # # >>> {'node1': 'node1', 'node2': 'node2', 'A': 'A', 'B': 'B', 'C': 'C'}
    # # Change the format of the network to rtt format
    # net.change_format("rtt-format")
    # print(net.truth_tables)
    # # >>> [array([0, 0, 0, 1]), array([1, 1, 0, 1]), array([0, 1]), array([0, 1]), array([0, 1])]
    # print(net.regulators)
    # # >>> [[2, 3], [2, 4], [2], [3], [4]]
    # print(net.num_nodes)
    # # >>> 5
    # print(net.input_type)
    # # >>> rtt-format
    # # print the reduced truth tables
    # print(net.truth_tables)
    # # print the regulators for each node
    # print(net.regulators)

    # # Loading a boolean network in rtt format
    # net = BooleanNet(input_type="rtt-format",num_nodes=3)
    # print("Boolean network loaded successfully")
    # # >>> Boolean network loaded successfully
    # print("Number of nodes: ", net.num_nodes)
    # # >>> Number of nodes: 3
    # print("Input type: ", net.input_type)
    # # >>> Input type: rtt-format

    # net.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
    # print(net.truth_tables)
    # # >>> [[1, 1, 0, 0], [0, 1], [0, 0, 0, 0]]
    # print(net.regulators)
    # # >>> [[0,2],[0],[1,2]]
    
    # net1 = BooleanNet(input_type="rtt-format",num_nodes=3)
    # net1.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
    # print("Boolean network loaded successfully")
    # # >>> Boolean network loaded successfully
    # net.add(net1)
    # print("Boolean network added successfully")
    # # >>> Boolean network added successfully
    # print(net.truth_tables), print(net.regulators)
    # # >>> [[1, 1, 0, 0], [0, 1], [0, 0, 0, 0], [1, 1, 0, 0], [0, 1], [0, 0, 0, 0]]
    # # >>> [[0, 2], [0], [1, 2], [3, 5], [3], [4, 5]]

    # # Loading a textfile from the cell collective database
    # # We will load one of the large networks from the cell collective database
    # # The indegree of several nodes is greater than 15. So, we will get a warning.
    # # When we print the truth tables, we will see that the truth tables are empty for nodes with indegree greater than 15.
    # net = BooleanNet( input_type="text-format")
    # file_path = os.path.join(os.path.dirname(__file__),"..", "update_rules_cell_collective", "ErbB (1-4) Receptor Signaling_23637902.txt")
    # net.load(parse_text_file(file_path=file_path))
    # net.change_format("rtt-format")
    # for formula,table in zip(net.update_formulae,net.truth_tables):
    #         print(formula)
    #         print(table)
    #         # >>> ('ErbB3_Y1180', '( EGFR_ErbB3 )  or ( ErbB2_ErbB3 )  or ( ErbB3_ErbB4 )')
    #         # >>> [0 1 1 1 1 1 1 1]
    #         # >>> ('EGFR_EGFR_EGF_PM', '( ( EGFR_EGFR_EGF_PM  ) and not ( EGFR_EGFR_EGF_CCP  ) )  or ( ( ( EGFR_Free and ( ( ( EGF ) )  and ( ( not EGFR_T654 ) ) )     ) and not ( EGFR_EGFR_EGF_CCP  )  ) and not ( ErbB2_Free  ) )')
    #         # >>> [0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 0]
    #         if len(table)==0:
    #             # >>> ('PI3K', '( ( Gbg_i  ) and not ( PI3K_I  ) )  or ( ( Fak  ) and not ( PI3K_I  ) )  or ( ( E_cadherin  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1257  ) and not ( PI3K_I  ) )  or ( ( Ras  ) and not ( PI3K_I  ) )  or ( ( Src and ( ( ( Cbl_RTK ) ) )     ) and not ( PI3K_I  ) )  or ( ( EGFR_Y920  ) and not ( PI3K_I  ) )  or ( ( ErbB4_Y1056  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1203_05  ) and not ( PI3K_I  ) )  or ( ( Gab1  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1270  ) and not ( PI3K_I  ) )  or ( ( EGFR_Y845  ) and not ( PI3K_I  ) )  or ( ( Crk  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1178  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1035  ) and not ( PI3K_I  ) )  or ( ( ErbB3_Y1241  ) and not ( PI3K_I  ) )')
    #             # >>> []
    #             # This is because the indegree of the node is greater than 15
    #             input("Press Enter to continue...")
    # print(net.regulators)
    # >>> [[216, 100, 42, 143, 111, 160], [61], [61, 149, 224, 88, 174, 2], [84, 210, 3, 162, 172, 176, 169], [78, 216, 104, 128, 180, 95, 93, 4, 81, 67, 50], [27, 82, 23], [6, 61, 2, 199], [7, 216, 2], [4, 93], [155, 82, 31, 18], [3], [136, 11, 57, 212], [179], [216, 191, 0, 131, 156, 35], [61, 179], [221, 62, 144, 23, 46, 31], [216, 3, 96, 109, 61, 56, 30, 2], [82, 17, 23, 27, 76], [201, 173, 52, 225], [41], [58], [113, 74, 142, 80, 21, 211, 203, 73], [182, 92, 37, 68, 34, 161], [226, 227, 166, 191, 228, 229, 23, 52, 173, 225], [221, 62, 144, 23, 46, 31], [220], [26, 63, 172, 107, 85], [226, 173, 177, 52, 27, 225], [108, 216, 113, 106, 171, 93, 50, 29], [131, 80], [224, 149, 88, 174, 30], [201, 191, 227, 166, 173, 52, 225], [216, 42, 61, 218, 30, 2], [33, 171, 87], [34, 22, 222], [221, 204, 31, 13, 74, 144, 62, 142, 23, 34, 46, 203, 73, 222], [94], [161], [182, 36, 161], [4, 93], [141, 40, 207], [106, 3, 61, 58, 93, 4, 6, 2], [216, 42, 59, 149, 7], [221, 204, 62, 74, 144, 127, 142, 23, 34, 222, 46, 203, 73, 31], [105, 17, 76], [216], [191, 227, 166, 173, 52, 225], [59, 91, 230, 3], [231, 179, 188], [167, 49, 232, 187, 95, 14, 91], [216, 114, 0, 138, 93, 4, 43, 91, 50], [171], [18, 52, 46, 27, 233], [43, 72], [171, 91, 50, 179], [221, 204, 62, 74, 144, 162, 142, 23, 34, 24, 222, 46, 218, 160, 203, 73, 31], [152, 215, 110], [11, 61], [78, 128, 58, 125, 81, 57, 91], [216, 1, 59, 172, 7, 230], [216, 171, 150, 158, 95, 99, 10, 2, 212], [224, 61, 149, 174, 30, 2], [62, 191, 234, 166, 173, 203], [8, 208, 209, 39, 195, 90, 192, 194], [27, 82, 23], [216, 47, 0, 42, 83, 4, 43, 2], [3], [169, 179], [221, 119, 62], [27, 82, 23], [70, 112, 91, 95, 172, 179, 169], [140], [185], [204, 113, 106, 80, 21, 211, 203, 2, 73], [113, 74, 106, 142, 222], [221, 62, 122, 75, 92, 32, 196], [17, 235], [127, 83, 61], [162, 2, 3], [79, 172, 118, 98], [216, 148, 94, 21, 211, 2, 213, 29], [216, 42], [201, 173, 52, 225], [216, 42, 127, 47], [236, 140], [8, 208, 209, 39, 195, 90, 115, 192, 72, 194], [136, 60], [33, 216], [49, 3, 58, 219, 6], [216], [135, 4, 93], [87, 91, 40, 172], [75], [77, 104, 59, 110, 93, 180, 95, 4, 230, 50, 16], [61, 222], [167, 49, 95], [12, 179], [155, 82, 31, 18], [132, 79, 56, 98], [219, 99, 124], [62, 74, 15, 97, 34, 203, 222, 18, 142, 9, 5, 24, 221, 204, 82, 23, 149, 69, 64, 216, 144, 42, 218, 178, 46, 73, 31], [216, 204, 210, 74, 144, 221, 3, 62, 142, 23, 34, 19, 222, 46, 169, 73, 31, 203], [127], [75], [48], [23, 27, 228, 82], [25, 200, 126, 67], [171, 117, 40, 26], [116, 125], [216, 49, 96, 109, 91, 40], [210, 4, 179], [27, 82, 23], [181, 2, 55], [221, 62, 0, 206, 2], [216, 12, 1, 114, 61, 30, 172, 87], [140], [237, 214, 93, 125, 4, 67, 220], [137, 172, 121, 107, 85], [132, 79, 118], [221, 62, 13, 142, 203], [216, 152], [8, 208, 123, 90, 115, 72, 50], [103], [236], [238, 150, 124, 134], [216, 108, 204, 74, 144, 221, 62, 142, 23, 34, 58, 222, 46, 203, 73, 31], [65], [216, 42, 127, 47], [216, 32, 42, 128], [216, 42, 95, 152, 120, 99, 91, 212], [36], [13], [216, 114, 132, 171, 56, 98, 10, 91, 50], [27, 82, 23], [167, 60], [239], [86, 136, 150, 240], [8, 208, 209, 123, 39, 90, 115, 192, 72, 194, 50], [221, 62, 74, 144, 236, 142, 23, 34, 222, 46, 203, 73, 31], [182, 36, 94], [135], [59, 230, 183, 95], [221, 13, 74, 80, 142, 75, 211, 68, 2, 213], [221, 62, 144, 23, 46, 31], [167, 191, 124, 166, 136, 169], [221, 204, 62, 74, 144, 145, 0, 142, 61, 23, 34, 152, 222, 46, 203, 73, 31], [60, 187], [61, 75], [2], [236, 224, 58, 180, 4, 241], [154, 134, 146, 86], [155, 82, 31, 18], [145], [46, 27, 18], [60, 124], [201, 173, 225], [221, 62, 144, 23, 46, 31], [216, 221, 166, 127, 59], [158, 187, 183], [46, 27, 18], [221, 62, 144, 23, 46, 31], [36, 130], [163], [96, 171, 61, 169, 2], [34, 139], [216, 13, 4, 180], [221, 204, 62, 144, 166, 23, 242, 46], [167, 243, 150, 134], [216, 106, 7, 2, 168], [70, 244], [46, 27, 18], [107, 171, 186, 172], [221, 204, 31, 74, 62, 144, 142, 23, 34, 172, 46, 203, 73, 222], [82, 173, 228, 23, 245, 27], [216, 42, 6, 19, 20, 174, 7], [224, 149, 88, 174, 30, 175], [3, 87, 176], [27, 229, 44], [46, 27, 18], [169], [61, 205, 30], [181, 91, 219, 95], [36, 161], [183, 158, 187, 14, 91], [21, 113], [117, 140], [209, 132, 186, 115, 172, 192, 195, 50], [232, 146, 150, 187], [167, 3, 124, 136, 61, 56, 187, 87, 150], [1, 3, 132, 20, 91, 50], [46, 27, 18], [3], [93, 135, 145, 140], [61, 75], [4, 26, 93], [100, 0, 71, 56, 58, 140, 93, 4, 53], [56], [216, 221], [46, 27, 18], [2], [200, 1, 20, 54, 179, 189, 50], [201, 155, 82, 18, 246, 31], [46, 27, 18], [62, 80, 75, 211, 68, 203, 2, 73, 213], [147, 204, 184, 166, 61, 103, 75, 21, 73], [1], [206, 131, 68, 2, 213], [171, 91, 150], [84], [221, 204, 62, 74, 209, 216, 0, 142, 224, 23, 34, 222, 46, 203, 73, 31, 55], [210, 101, 181, 2, 55], [221, 62, 61, 131, 68, 2, 213, 29], [136, 11, 212], [221, 119, 62, 2], [216, 116, 169], [13, 0, 171, 83, 61, 43, 2], [216, 49, 183, 42, 138, 127, 143, 187, 15, 218, 129, 150, 66], [46, 27, 18], [221, 62, 144, 23, 152, 218, 46, 31], [124, 219, 99, 14, 181], [108, 204, 62, 74, 221, 216, 144, 138, 142, 23, 34, 222, 46, 203, 73, 31], [221, 227, 191, 166, 142, 173], [147, 74, 184, 61, 103, 75, 21, 34, 222], [46, 27, 18], [216, 217, 13, 145, 190, 42, 83, 56, 202, 95, 116, 247, 157, 198, 153, 197, 170, 151], [225], [226], [227], [228], [229], [230], [231], [232], [233], [234], [235], [236], [237], [238], [239], [240], [241], [242], [243], [244], [245], [246], [247]]
    print("====================================================================================================")   
    print("====================================================================================================")
    # testing add linear edfe
    net = BooleanNet(input_type="rtt-format",num_nodes=3)
    net.load(([[1,1,0,0],[0,1],[0,0,0,0]],[[0,2],[0],[1,2]]))
    net.add_linear_edge_rtt((0, 1))
    print(net.truth_tables)
    print(net.regulators)
    
    