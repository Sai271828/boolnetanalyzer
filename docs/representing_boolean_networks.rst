Boolean Network Representations
===============================

This package provides several subclasses of Boolean Networks, each offering more efficient representations than the general case. Below, we outline the supported representations:

Random Boolean Networks
-----------------------
For random Boolean networks, we offer two key representations:

1. **Formula-Based Representation (FBR):**
    - This representation uses a logical formula to define the network.
    - The formula consists of a conjunction of disjunctions of literals.
    - The formula is stored as a string.
    
    **Example:**
    ::
    
        f(x1, x2) = ((x1 or x2) and (not x1 or x2), (x1 or not x2))
    For the exact format, see the 

2. **Reduced Truth Table (RTT):**
    - This representation uses the truth table of the network.
    - The truth table is stored as a list of 0's and 1's.
    - Each integer represents the output of the network for a particular input.
    - Each truth table also comes with a list of integers representing the regulators of each node.
    
    **Example:**
    ::
    
        f(x1, x2) = (x1 and x2, x1)

    is represented as:

    ::

        [[0, 0, 0, 1], [0, 1]]
        [[0,1],[0]]

    In this case, the first list contains two lists each representing the output of the first node when the inputs are 00, 01, 10, 11 respectively, and the outputs of the second node when the inputs are 0 and 1 respectively.
    The second list contains the regualtors of each node.
    The second list also contains the order of regulators for each node.
    For instance [0,1] is different from [1,0] as the order of regulators is different.
    In the first case, the first regulator is x1 and the second regulator is x2.
    Therefore the truth table values are x1=0, x2=0, f(x1,x2)=0, x1=0, x2=1, f(x1,x2)=0, x1=1, x2=0, f(x1,x2)=0, x1=1, x2=1, f(x1,x2)=1.
    In particular, x1 is the most significant bit and x2 is the least significant bit.
    In the second case, the first regulator is x2 and the second regulator is x1.
    Therefore the truth table values are x2=0, x1=0, f(x1,x2)=0, x2=0, x1=1, f(x1,x2)=0, x2=1, x1=0, f(x1,x2)=0, x2=1, x1=1, f(x1,x2)=1.
    So, in this case, x2 is the most significant bit and x1 is the least significant bit.
.. note::

   Currently, support for Nested Canalizing Functions (NCF) is unavailable.



Nested Canalizing Functions (NCF)
---------------------------------
This representation is based on the concept of canalizing functions:

- A function is *canalizing* if it has a specific input that can determine the output regardless of other inputs.
- A function is *nested canalizing* if it has canalizing inputs that are themselves canalizing.

The NCF representation is stored as a list of integers, where each integer represents the canalizing input of the corresponding node.

**Example:**
::
    
    f(x1, x2) = (x1 and x2, x1 or x2)

is represented as:

::

    [[0, 0], [0, 0]] for x1 and x2
    [[1, 0], [1, 0]] for x1 or x2
