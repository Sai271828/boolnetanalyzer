module cubicalhomology
using SparseArrays; # for sparse matrix and vecotr handling 
using LinearAlgebra; # used for rank computations 
using Base.Threads;

# Import the original .jl file
include("a_theory.jl")

# Export the main functions so they can be used externally
export nCube,add_numbers,hyperOctahedrial,create_poset

end

