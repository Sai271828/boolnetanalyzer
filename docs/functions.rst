Package Overview
================

This package is designed to analyze boolean networks. It provides tools to load boolean networks from different sources, generate random boolean networks, and analyze the networks.
The package is organized into the following modules:

* :mod:`boolnetanalyzer.utils` - Utility functions

* :mod:`boolnetanalyzer.randomize` - Functions to generate random boolean networks

* :mod:`boolnetanalyzer.bn_tools` - Functions to analyze boolean networks

* :mod:`boolnetanalyzer.load` - Functions to load boolean networks from files


Loading Boolean Networks
================================
We implement a Boolean Network class to represent boolean networks. 
The class includes methods to load boolean networks from different sources, such as NCF files, RTT files, and text files. 

Boolean Network Class
=====================


The :class:`BooleanNet` class represents a boolean network with the following attributes:

- **num_nodes** (int):
  The number of nodes in the network.

- **input_type** (str):
  The type of input, such as "rtt", "text", "ncf", etc.

The main function of the class is to load Boolean networks from different sources.
.. automodule:: boolnetanalyzer.load
    :members:
.. autofunction:: boolnetanalyzer.load.BooleanNet.load

Below we describe the class and other methods in detail.

.. autoclass:: boolnetanalyzer.load.BooleanNet
   :members: update, add, add_node, add_edge, change_format
   :exclude-members: __init__, add_ncf, add_rtt, add_text, add_edge_ncf, add_edge_rtt, add_edge_text, add_node_ncf, add_node_rtt, add_node_text, load_ncf, load_rtt, load_text, update_ncf, update_rtt, update_text,change_format_ncf, change_format_rtt, change_format_text  

Generating Random Boolean Networks
============================

.. automodule:: boolnetanalyzer.randomize
    :members:
    :undoc-members:
    :private-members:
    :special-members:
    :inherited-members:
    :show-inheritance:

Tools for Analyzing Boolean Networks
=====================

.. automodule:: boolnetanalyzer.bn_tools
    :members:
    :undoc-members:
    :private-members:
    :special-members:
    :inherited-members:
    :show-inheritance:

