# structures for graphs and posets
using Base.Threads
# Define Graph Structure
mutable struct graph
    # araray of vertex labels
    vertices::Array 

    # array of edges in the graph. Needs to be simple graph so it should be symmetric and reflexive 
    edges::Array 
end   

# Define poset Structure
mutable struct poset
    # the underlying set of the poset
    vertices::Array

    # the relations, (a,b) means a < b should be transitive, the function make_transitive can do this if just have the generators
    relations::Array 
end   

# Define simplicial copmplex Structure
mutable struct simplicialComplex
    # underlying vertex set
    vertices::Array
    # the set of simplices
    simplices::Array
end   

#=======================================================================================
Functions for generating the equivalence classes according to the hyperOctahedrial group
========================================================================================#

#= function to generate the n-cube
    Inputs: 
        1) a natural number n

    Outputs: 
        1) n-cube as an array of two arrays, first being the vertices, second being the faces
    
    Notes: 
        - kind of messy but works 
=#
function nCube(n::Int)

    # handling for zero cube
    if n==0
        return [[(0)],[]]
    end

    map=[[[0],[1]],[[[[0]],[[1]]]]]

    if n==1
        return map
    end

    
    # iteratively define higher dimensional cubes 

    for i=2:n
        new_map = [[],[]]
        p2 = deepcopy(map)

        map4 = deepcopy(map)

        for x = 1:length(map[1])
            push!(map4[1][x],0) # origonal becomes 0 face, new becomes 1 face

            push!(p2[1][x],1)
        end

        map = deepcopy(map4)
        
        map3=deepcopy(map)

        for x = 1:length(map[2]) # loop thru pairs of faces
            for y = 1:length(map[2][x][1]) # loop thru vertices in faces
                push!(map3[2][x][1][y],0)
                push!(map3[2][x][2][y],0)

                push!(p2[2][x][1][y],1)
                push!(p2[2][x][2][y],1)
                
            end
        end
        
        map=deepcopy(map3)
        

        new_map[1]=vcat(map[1],p2[1])
        for x = 1:length(map[2])
            push!(new_map[2],[vcat(map[2][x][1],p2[2][x][1]),vcat(map[2][x][2],p2[2][x][2])])
        end

        push!(new_map[2],[deepcopy(map[1]),deepcopy(p2[1])])

        map = deepcopy(new_map)

    

    end

    return map
end
 
#= function to generate the reversals as lists of -1 and 1
    Inputs: 
        1) a natural number n, representing the dimension you want to generate (length of list)

    Outputs: 
        1) An array of all reversals for that dimension
    
    Notes: 
        - 1 does nothing to a coordinate, -1 changes form 1 to 0 and 0 to 1 
=#
 function generate_reversals(n::Int)
    if n == 0
        return [[]] # empty array for 0th dimension
    else
        shorter_lists = generate_reversals(n - 1)
        return vcat([[1; l] for l in shorter_lists], [[-1; l] for l in shorter_lists])
    end
end

#= function to generate all permutations of a list (used for symmetries)
    Inputs: 
        1) a list 

    Outputs: 
        1) a list of all permutations 
    
    Notes: 
        
=#
 function permutations(lst)
    if length(lst) == 0
        return [[]]
    else
        result = []
        for i in 1:length(lst)
            # Take the current element
            current = lst[i]
            # Remaining list without the current element
            remaining = [lst[j] for j in 1:length(lst) if j != i]
            # Recursively find permutations of the remaining list
            for perm in permutations(remaining)
                push!(result, [current; perm])
            end
        end
        return result
    end
end

#= function to determine the sign of a permutation
    Inputs: 
        1) a permutation

    Outputs: 
        1) its sign 
    
    Notes: 
        
=#
 function sign_of_permutation(perm)
    n = length(perm)
    inversions = 0
    for i in 1:(n-1)
        for j in (i+1):n
            if perm[i] > perm[j]
                inversions += 1
            end
        end
    end
    return (-1)^inversions
end

#= function to generate all permutations of a list (used for symmetries)
    Inputs: 
        1) an integer n

    Outputs: 
        1) a list of all elements of the n-th hyperoctahedrial group, written as [symmetry, reversal, sign] 
    
    Notes: 
        
=#
 function hyperOctahedrial(n::Int)
    grp = []

    list = 1:n 
    # Generate all permutations
    symmetric_Gp = permutations(list) # lists represent image of ordered list
    reversals = generate_reversals(n)
    for g in symmetric_Gp
        for r in reversals

            sgn = prod(r)*sign_of_permutation(g)

            push!(grp,[g,r,sgn])
        end
    end

    return grp
end

#= function to permute the coordinates of a vertex in the n-cube according to a group element in the hyperOctahedrial group
    Inputs: 
        1) vertex in coordinate representation
        2) the group element

    Outputs: 
        1) the permuted vertex
    
    Notes: 
        
=#
 function permuteCoords(vertex, gpElet)
    
    pm1=[]
    for n = 1:length(vertex)
        push!(pm1,vertex[gpElet[1][n]])
    end

    for n = 1:length(vertex)
        pm1[n] = Int(mod(pm1[n] + (gpElet[2][n]+3)/2 , 2))
    end
    return pm1
end

#= Function that calculates the image of an n-cube under a group element
    Inputs: 
        1) then n-cube
        2) the group element

    Outputs: 
        1) the image
    
    Notes: 
        
=#
 function calculateImageCube(cube,gpElet)
    cube2=deepcopy(cube)

    if !(length(cube2[1][1]) == length(gpElet[1]))
        println("incompatible sizes")
        return 
    end

    for n=1:length(cube2[1])
        cube2[1][n] = permuteCoords(cube2[1][n],gpElet)
    end

    for n=1:length(cube2[2])
        for m = 1:length(cube2[2][n][1])
            cube2[2][n][1][m] = permuteCoords(cube2[2][n][1][m], gpElet)
            cube2[2][n][2][m] = permuteCoords(cube2[2][n][2][m], gpElet)
        end
    end

    return cube2

end

#= Function to take a cube and its image under a transformation, and give an array representing the map from the cube to its image
    Inputs: 
        1) the n-cube
        2) the image

    Outputs: 
        1) an array representing the map
    
    Notes: transformation[i] = vertex in image that cube.vertices[i] gets sent to
        
=#
 function transformationCoords(cube,image)
    transformation = []
    for n = 1:length(image[1])
        m = findfirst(==(image[1][n]), cube[1]) # coord of image in origonal cube
        push!(transformation,m)
    end
    return transformation
end

#= Function to generate the list of the image maps, how the faces get swapped, and the sign of the permutation for the entire n-th hyperoctahedrial group
    Inputs: 
        1) a natural number n

    Outputs: 
        1) list of maps as an array of arrays
    
    Notes: 
        - each inner array is an array of three arrays, first is the map, second is how the faces are mapped, and third is the sign of the permutation
        
=#
function generateEQClist(n)
    EQClist = []

    cube = nCube(n)
    grp = hyperOctahedrial(n)

    for g in grp # if parralelized, need first thing in list to be identity
        transformation = transformationCoords(cube,calculateImageCube(cube, g))
        push!(EQClist,[transformation,g[1],g[3]])
    end
    
    return EQClist
end

#= Function that generates the equivalence class of a map
    Inputs: 
        1) the map as an array
        2) a list of degenerate coordinates
        3) the list of maps from the generateEQClsit function

    Outputs: 
        1) the equivalence class of the map
    
    Notes: 
        
=#
function generateEQC(map,degens,EQClist)
    EQC = [] # output array of maps
    for el in EQClist
        newDegen = [el[2][Int(d)] for d in degens] # compute the new degenerate coordinates
        cmap = [map[im] for im in el[1]] # new map under transformation
        push!(EQC,[cmap,newDegen,el[3]])
    end        
    return EQC
end

#= Function that creates a list of the faces of a cube in terms of position in the origonal map
    Inputs: 
        1) a natural number n

    Outputs: 
        1) a list of the faces 
    
    Notes: 
        - faces are paied [A,B], where A is negative face and B is positive
        
=#
function faceList(n)
    cube = nCube(n) # initialize n cube
    facesList=[] # output array of faces
    for f in cube[2]
        cface1=[] # negative face
        cface2=[] # positive face
        for i in f[1]
            m = findfirst(==(i), cube[1])
            push!(cface1,m)
        end

        for i in f[2]
            m = findfirst(==(i), cube[1])
            push!(cface2,m)
        end
        push!(facesList,[cface1,cface2]) # add the pair to the list
    end
    return facesList
end

#= Function to compute the faces of a map
    Inputs: 
        1) the map as an array
        2) the facelist given from the faceList function
    Outputs: 
        1) the faces of the map, paired as [negativeFace, positiveFace]
    
    Notes: 
        
=#
function faces(map,facesList)
    face=[]
    for f in facesList
        m1=[map[i] for i in f[1]] # negative face
        m2=[map[j] for j in f[2]] # positive face
        push!(face,[m1,m2])
    end
    return face
end

#===================================================================#
# Functions to compute the homology of a Graph 
#===================================================================#


#= Function to generate a dictionary of neighborhoods of a graph
    Inputs: 
        1) the graph G

    Outputs: 
        1) the neighborhood dictionary of the graph
    
    Notes: 
        - dictionary of the form a => A, where A is the set of vertices a is connected to
=#
function get_nhood_dict(G)
    # initialize nhood dict for checking maps
    nhoodDict = Dict{Any,Array{Any}}()
    for v in G.vertices
        cnhood = []
        for e in G.edges
            if e[1]==v
                push!(cnhood,e[2]) # add neighbors
            end
        end
        nhoodDict[v]=cnhood
    end
    return nhoodDict
end

#= Function to determine if a pair of n-1 cubes forms an n cube by making f the n-th negative face and g the n-th positive face
    Inputs: 
        1) the first n-1 cube
        2) the second n-1 cube
        3) the neighborhood dictionary of the graph

    Outputs: 
        1) true/false value of if they form an n cube
    
    Notes: 
        
=#
 function isPairNcube(f,g,nhoodDict)
    map=true
    
    for i in 1:length(f[1]) # since maps are pairs [A,B,C], where A is the map, B is the degenerate coordinates, and C is the sign relative to the EQC rep
        if !(g[1][i] in nhoodDict[f[1][i]]) #check that pairs are connected
            map=false
            break
        end
    end
    return map
end

#= Function that generates n-cubes of a graph
    Inputs: 
        1) an array of the n-1 cubes
        2) the graph G
        3) the equivalence class list for the n-th dimension
        4) the neighborhood dictionary of the graph

    Outputs: 
        1) the n cubes
    
    Notes: 
        
=#
function graphMaps(A,G,EQClist,nhoodDict) # A = array of the n-1 cubes, G = graph
    
    # check if generating zero cubes
    if length(A)==0
        C=[]
        for i in G.vertices
            push!(C,[[[i],[],1]]) # last coord is automorphism sign wrt first cube in class
        end
        return C
    end

    B = Channel{Vector{Any}}(length(A)*length(A)*length(A[1])) # output channel of maps
    
    @threads for i = 1:length(A) 
        for j = i:length(A)
            f=A[i]
            g=A[j]
            
            
            for t=1:length(f)
                if isPairNcube(f[1],g[t],nhoodDict) 
                    h = vcat(f[1][1],g[t][1])
                    degens = []
                    for d in f[1][2]
                        if d in g[t][2]
                            push!(degens,d)
                        end
                    end

                    if f[1]==g[t]
                        push!(degens,log2(length(h)))
                    end

                    eq = generateEQC(h,degens,EQClist)

                    
                    put!(B,eq)
                    
                end
            end
        end 
    end
    close(B) # close the channel
    cubes = collect(B) # collect the cubes
    return(cubes)
end

#= Function that removes degenerate maps from a list
    Inputs: 
        1) the maps to check

    Outputs: 
        1) the set of non degenerate maps

    Notes: 
        
=#
function nonDegeneracies(maps)
    B = Channel{Vector{Any}}(length(maps))
    @threads for i in maps
        f=i[1] # check representitive
        if length(f[2])==0 #check if degenerate in any direction
            put!(B,i)
        end
    end
    close(B)
    nonDegen=collect(B)
    return nonDegen
end

#= Function that checks if a map is semi degenerate, i.e. if x~-x
    Inputs: 
        1) the map
        2) the n-th equivalence class list

    Outputs: 
        1) a true/false value of if its semi degenerate
    
    Notes: 
        
=#
function checkSemiDegen(map,EQClist)

    for g in EQClist
        if g[3]==-1 # check all maps with sign -1 in EQC
            if map == [map[im] for im in g[1]] # check that they're equal
                return true
            end
        end
    end
    return false
    
end

#= Function to remove semi degenerate maps
    Inputs: 
        1) an array of maps (here they are equivalence classes)

    Outputs: 
        1) the non semi degenerate maps
    
    Notes: 
        
=#
function semiNonDegen(maps)
    
    B=Channel{Vector{Any}}(length(maps))
    @threads for i in maps
        semiDegen=false
        for j in i
            if i[1][1]==j[1] && prod(i[1][3])*prod(j[3]) == -1 # check if semi degenerate
                semiDegen=true
                break
            end
        end
        if !(semiDegen)
            put!(B,i) # if not semi degenerate, add to list
        end
    end
    close(B)
    semiNonDegen = collect(B)
    return semiNonDegen
end

#= Function to check if a map f appears in the equivalence class g
    Inputs: 
        1) the map f
        2) the equivalence class g

    Outputs: 
        1) a true/false value
    
    Notes: 
        
=#
function is_related(f,g) # f = map, g = equivalence class

    if !(sort(f[1])==sort(g[1][1])) # check if they have same vertices
        return false
    end
    
    for h in g
        if h[1]==f[1] # check if equal 
            return true
        end
    end

    return false
end

#= Function to remove duplicate equivalence classes
    Inputs: 
        1) the array of equivalence calsses of maps

    Outputs: 
        1) the array of unique equivalence classes
    
    Notes: 
        
=#
function remove_duplicates(maps) 
    B=Channel{Vector{Any}}(length(maps))
    @threads for i = 1:length(maps)
        k=maps[i]
        f = k[1]
        not_dupe = true
        for j = i+1:length(maps)
            g = maps[j]
            if is_related(f,g)
                not_dupe = false
                break
            end
        end

        if not_dupe
            put!(B,k)
        end
    end

    close(B)
    C = collect(B)
    return C
end

#= Function to generate a dictionary of coordinates of maps
    Inputs: 
        1) array of equivalence classes

    Outputs: 
        1) the coordinate dictionary
    
    Notes: 
        
=#
function coordDict(lowerCubes)
    cdict = Dict{}()
    for i = 1:length(lowerCubes)
        for j = 1:length(lowerCubes[i])
            cdict[lowerCubes[i][j][1]] = [i,lowerCubes[i][j][3]]
        end
    end

    return cdict 
end

#= Function converts an array of sparse vectors to a matrix with the vectors as column vectors
    Inputs: 
        1) the array of vectors

    Outputs: 
        1) the matrix
    
    Notes: 
        
=#
function sparse_col_concat(vectors)
    # Determine the dimensions of the resulting matrix
    n_rows = length(vectors[1])
    n_cols = length(vectors)
    
    # Initialize arrays to store the row indices, column indices, and values of the non-zero elements
    row_indices = Int[]
    col_indices = Int[]
    values = eltype(vectors[1])[]
    
    # Iterate through each vector and extract its non-zero elements
    for (j, vec) in enumerate(vectors)
        for i in 1:length(vec)
            if vec[i] != 0
                push!(row_indices, i)
                push!(col_indices, j)
                push!(values, vec[i])
            end
        end
    end
    
    # Construct the sparse matrix using the collected non-zero elements
    return sparse(row_indices, col_indices, values, n_rows, n_cols)
end

#= Function computes the boundary of a map and stores it as a sparse vector
    Inputs: 
        1) the map as an array
        2) the lower non degenerate cubes
        3) the n-dimensional face list

    Outputs: 
        1) the boundary as a coordinate vector, stored sparsly
    
    Notes: 
        
=#
function boundarySum(map,lowerNonDegen,faceList)
    
    image = spzeros(length(lowerNonDegen))

    cdict = coordDict(lowerNonDegen)
    
    fc = faces(map[1],faceList) # generate faces

    for i = 1:length(fc) # loop through faces
        # negative face
        key = fc[i][1]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            image[ind] += (-1)^(i+1)*sgn
        end
        
    
        # positive face
        key = fc[i][2]
        if key in keys(cdict)
            indSgn=cdict[key]
            ind = indSgn[1]
            sgn = indSgn[2]

            image[ind] += (-1)^(i)*sgn
        end
    end

    return image
end

#= Function to calculate image of boundary map del_n: L_n -> L_n-1
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        
=#
function calculateImage(nonDegens,lowerNonDegen,faceList)
    B = Channel{SparseVector{Float64, Int64}}(length(nonDegens))
    @threads for map in nonDegens
        im = boundarySum(map[1],lowerNonDegen,faceList) # take first rep of classes
        put!(B,im)        
    end 
    close(B)
    image = collect(B)
    return image
end

#= function to generate the matrix representing the n-th boundary map, stored sparsly
    Inputs: 
        1) the array of n dimensional maps (EQCs of maps)
        2) the array of n-1 dimensional maps (EQCs of maps)
        3) the n-dimensional face list

    Outputs: 
        1) the matrix of the n-th boundary map stored as a sparse matrix
    
    Notes: 
        
=#
function boundaryMapMatrixSparse(nonDegens,lowerNonDegens, faceList)  # for now cant do sparse, fix
    im = calculateImage(nonDegens, lowerNonDegens,faceList)

    if length(im)==0
        return []
    end

    m=sparse_col_concat(im)
    return m
end

#= Function to generate the matrix of the n+1 boundary map directly without storing any intermediate variables
    Inputs: 
        1) an array of n dimensional EQCs
        2) an array of n dimensional non degenerate EQCs
        3) n+1 face list
        4) the nhood dictionary of the graph
        5) the n+1 dimensional EQClist

    Outputs: 
        1) an array of sparse coordinate vectors 
    
    Notes: 
        - right now using the conjecture, if false remove the break line 
        
=#
function graphMapsMatrix(A,lowerNonDegen,facesList,nhoodDict,EQClist) # A = array of the n-1 cubes, G = graph
    
    B = Channel{SparseVector{Float64, Int64}}(length(A)*length(A)) # channel of coordinate vectors

    cdict = coordDict(lowerNonDegen)

    lenVec = length(lowerNonDegen)

    @threads for x = 1:length(A)-1  
 
        for y = x+1:length(A)
            
            h=A[x]
            k=A[y]
            
            f=h[1]
            
                for g in k
                    if isPairNcube(f,g,nhoodDict) # check if forms an n+1 cube
                        
                        # check first to see if degenerate
                        degen=checkSemiDegen(vcat(f[1],g[1]),EQClist)
                        
                        if !(degen) 

                            image = spzeros(Int,lenVec) # generate coordinate vector
                            
                            fc = faces(vcat(f[1],g[1]),facesList) # generate the faces
                            
                            for i = 1:length(fc) # loop through the faces
                                
                                # negative face
                                key = fc[i][1]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    image[ind] += (-1)^(i+1)*sgn
                                end
                                
                            
                                # positive face
                                key = fc[i][2]
                                if key in keys(cdict)
                                    indSgn=cdict[key]
                                    ind = indSgn[1]
                                    sgn = indSgn[2]

                                    image[ind] += (-1)^(i)*sgn
                                end
                                
                            end
                             
                            put!(B,image)
                             
                        end
                        
                        break # only combine once per equivalence class, remove if conj false
                    end
                end
            
        end 
        
    end
    
    close(B)
    img = collect(B)
  
    if length(img)==0
        return []
    end
    
    m=sparse_col_concat(img)

    row_indices, col_indices, values = findnz(m)
    m=sparse(row_indices, col_indices, Float64.(values), size(m)...) # convert to float so rank function can handle

    return m
end

#= Function to generate the n-th homology of a graph
    Inputs: 
        1) the graph G
        2) the dimension to compute
    Outputs: 
        1) the n-th homology group
    
    Notes: 
        
=#
function cubicalHomology(G::graph,n::Int)
    
    #Error handling
    #verify graph has number entries
    if !all(x -> x isa Number, G.vertices)
        error("Error: all vertices of the graph must be numbers. Try running the relabel_vertices function on the graph before computing homology") 
    end

    #verify positive dimension
    if n < 0
        error("Error: dimension must be at least 0") 
    end


    #computing the homology
    #Initialize dicts of EQC and faces
    EQCdict = Dict{Int, Array{Any}}()

    faceDict = Dict{Int, Array{Any}}()

    #fill the dictionaries
    for i = 1:n+1
        EQCdict[i] = generateEQClist(i)
    end

    for i = 0:n+1
        faceDict[i] = faceList(i)
    end
    
    # initialize nhood dict for checking maps
    nhoodDict = get_nhood_dict(G)

    # Initialize map the dictionary
    maps = Dict{Int, Array{Any}}()  

    # Populate the dictionary
    maps[-1]= [[[nothing,[],1]]] # for 0-th homology groups
    maps[0] = graphMaps([],G,[],[])
    for i=1:n

        mi=graphMaps(maps[i-1],G,EQCdict[i],nhoodDict)
        
        maps[i]=mi
    
        maps[i]=remove_duplicates(maps[i])
        
    end

    # remove degeneracies and semi degeneracies
    nonDegen = nonDegeneracies(maps[n])

    lowerNonDegen = nonDegeneracies(maps[n-1])
    
    nonDegen = semiNonDegen(nonDegen)

    lowerNonDegen = semiNonDegen(lowerNonDegen)

    # n-th boundary map 
    delN = boundaryMapMatrixSparse(nonDegen, lowerNonDegen, faceDict[n])

    # determine the matrix for the n+1-th boundary map
    delN1 = graphMapsMatrix(maps[n],nonDegen,faceDict[n+1],nhoodDict,EQCdict[n+1])
    
    dimIM=rank(delN1)

    # handling for different cases
    if delN==[]
        dimKer = 0
    else
        dimKer=size(delN,2)-rank(delN)
    end
    
    # homology dimension
    return dimKer-dimIM
    
end

#===================================================================#
# Functions for preprocessing graphs
#===================================================================#

#= Function that removes a vertex from a graph and nhood dict
    Inputs: 
        1) the vertex to remove
        2) the graph G 
        3) the nhood dict of the graph
    Outputs: 
        1) the new graph
        2) the new nhood dict
    
    Notes: 
        
=#
function remove_vertex(v,G,nhoodDict)

    # remove from graph
    filter!(x -> x != v, G.vertices)
    filter!(t -> v ∉ t, G.edges)    

    pop!(nhoodDict,v)

    for key in keys(nhoodDict)
        filter!(x -> x != v, nhoodDict[key])
    end

    return G, nhoodDict
end

#= Function to detect and removes a vertex of degree n,  if none found returns the graph along with a false flag
    Inputs: 
        1) the graph G
        2) the nhood dict of the graph
    Outputs: 
        1) the graph with one vertex removed (if found)
        2) the nhood dict with one vertex removed (if found)
        3) true/false value if vertex was removed
    
    Notes: 
        
=#
function remove_vert_deg_n(G,nhoodDict)

    for key in keys(nhoodDict) 
        for v in filter(x -> x != key, nhoodDict[key]) # dont consider loops

            # check if v is connected to everything the key is
            conn = true
            for w in filter(x -> x != key, nhoodDict[key]) # still dont consider loops bc already know key ~ v
                if !(w in nhoodDict[v])
                    conn = false
                    break
                end
            end

            if conn == true
                G, nhooodDict = remove_vertex(key,G,nhoodDict)
                return G, nhoodDict, true
            end
        end
    end

    return G, nhoodDict, false # no case found
end

#= Function to preproccess a graph, making homology computations quicker
    Inputs: 
        1) the graph G
    Outputs: 
        1) the optimal graph to compute on given results about preprocessing
    
    Notes: 
        
=#
function preprocess_graph(G::graph)
    nhoodDict = get_nhood_dict(G)
    flg = true

    while flg
        G, nhoodDict, flg = remove_vert_deg_n(G,nhoodDict)
    end

    return G
end


#===================================================================#
# Functions relevant to making inputting graphs easier
#===================================================================#

#= Function to make a set symmetric and reflexive, useful for edge sets
    Inputs: 
        1) the set A
    Outputs: 
        1) the symmetric and reflexiove version of A
    
    Notes: 
        
=#
function makeSymmetricReflexive(A)
    B=copy(A)
    for a in A
        if !( (a[2],a[1]) in B)
            push!(B,(a[2],a[1]) )
        end
        if !( (a[1],a[1]) in B)
            push!(B,(a[1],a[1]) )
        end
        if !( (a[2],a[2]) in B)
            push!(B,(a[2],a[2]) )
        end
    end
    return B
end

#= Function to quotient a graph by relating two vertices
    Inputs: 
        1) the graph G
        2) the first vertex
        3) the second vertex
    Outputs: 
        1) the graph with v1 = v2
    
    Notes: 
        
=#
 function quotient_vertices(G::graph,v1::Any,v2::Any)
    new_vert=[]
    new_edges=[]

    for v in G.vertices
        if !(v==v2)
            push!(new_vert,v)
        end
    end

    for i in G.edges
        if !(i[1]==v2 || i[2]==v2)
            push!(new_edges,i)
        end
        if i[1]==v2 && !(i[2]==v2) && !(i[2]==v1)
            push!(new_edges,(v1,i[2]))
        end
        if !(i[1]==v2) && i[2]==v2 && !(i[1]==v1)
            push!(new_edges,(i[1],v1))
        end
    end
    return graph(new_vert,new_edges)
end

#= Function to generate the suspension of a graph
    Inputs: 
        1) the graph G
        2) the length of suspension
    Outputs: 
        1) the suspension of the graph
    
    Notes: 
        
=#
 function suspension(G::graph,n::Int)
    n=n-1
    vertices=[]
    edges=[]
    for i = 1:n
        for j in G.vertices
            push!(vertices,(j,i))
        end
    end
    for i in vertices
        for j in vertices
            if (i[1]==j[1]) & (abs(i[2]-j[2])==1)
                push!(edges,(i,j))
                push!(edges,(j,i))
            end
            if (i[2]==j[2]) & ((i[1],j[1]) in G.edges)
                push!(edges,(i,j))
                push!(edges,(j,i))
            end
            if i[2]==n
                push!(edges,(i,'s'))
                push!(edges,('s',i))
            end
            if i[2]==1
                push!(edges,(i,'n'))
                push!(edges,('n',i))
            end
        end 
    end
    push!(vertices,'n')
    push!(vertices,'s') 

    push!(edges,('s','s')) 
    push!(edges,('n','n')) 

    susp=graph(vertices,unique!(edges))
    return susp
end

#= Function to generate the cartesian product of an array of arrays
    Inputs: 
        1) the array of arrays
    Outputs: 
        1) the cartesian product
    
    Notes: 
        
=#
 function cartesian_product(arrays) 
    # Base case: if arrays is empty, return an array with an empty tuple
    if isempty(arrays)
        return [()]
    end
    
    # Get the first array and the rest of the arrays
    first_array, rest_arrays = arrays[1], arrays[2:end]
    
    # Recursive call to get the Cartesian product of the rest of the arrays
    rest_product = cartesian_product(rest_arrays)
    
    # Combine each element of the first array with each tuple from the rest product
    result = [(x, rest...) for x in first_array for rest in rest_product]
    
    return result
end

#= Function to compute the box product of a list of graphs
    Inputs: 
        1) the list of graphs
    Outputs: 
        1) the graph that is a box product of all graphs in the list
    
    Notes: 
        
=#
function multBoxProd(graphs) 

    vert=cartesian_product([G.vertices for G in graphs])
    edge=[]
    
    for v in vert
        
        for w in vert
            tot = 0 # check if its an edge in the product
            
            for i = 1:length(v)
                
                if !(v[i] == w[i]) && ( (v[i],w[i]) in graphs[i].edges) #check if theye connected in g_i
                    tot += 1
                elseif !((v[i],w[i]) in graphs[i].edges)
                    tot += 2
                end
                
            end
            if tot <2
                push!(edge,(v,w))
            end
        end
    end

    return(graph(vert,edge))

end

#= Function to generate the complete graph on n vertices
    Inputs: 
        1) a number n
    Outputs: 
        1) the complete graph on n vertices
    
    Notes: 
        
=#
function completeGraph(n)
    v=[]
    e=[]
    for i = 1:n
        push!(v,i)
    end

    for w1 in v
        for w2 in v
            push!(e,(w1,w2))
        end
    end

    return graph(v,e)
end

#= Function to relabel the vertices of a graph as integers
    Inputs: 
        1) the graph G
    Outputs: 
        1) the relabeled graph
    
    Notes: 
        
=#
function relabel_vertices(G::graph)
    new_vert=[]
    
    transformDict = Dict{Any,Any}()
    for i in 1:length(G.vertices)
        push!(new_vert,i)
        transformDict[G.vertices[i]] = i
    end
    B = Channel{Any}(length(G.edges))
    @threads for e in G.edges
        put!(B,(transformDict[e[1]],transformDict[e[2]]))
    end
    close(B)
    new_edge=collect(B)

    return graph(new_vert,new_edge)
end

#============================================================================#
# FUNCTIONS FOR POSETS
#============================================================================#


#= Function to take a poset and returns a dictionary of X: vertices X is < than
    Inputs: 
        1) the poset P
    Outputs: 
        1) the dictionary
    
    Notes: 
        
=#
function rel_dict(P) 
    rel=Dict{}()

    for v in P.vertices
        rel[v] = []
    end

    for x in P.relations
        push!(rel[x[1]],x[2])
    end
    return rel
end

#= Function to generate the set of chains of varying lengths of a poset
    Inputs: 
        1) the poset P
    Outputs: 
        1) the set of chains
    
    Notes: 
        
=#
function posetChains(P)
    chains = []
    rel=rel_dict(P)

    #All chains starting at each vertex
    # Function to perform DFS and find all chains
    function dfs(v, path, chains)
        push!(path, v)               # Add the current vertex to the current path
        push!(chains, copy(path))    # Add a copy of the current path to chains
        for neighbor in rel[v]  # For each neighbor of the current vertex
            dfs(neighbor, path, chains) # Recursively call dfs on the neighbor
        end
        pop!(path)                   # Backtrack by removing the current vertex from the path
    end

    # Find all chains starting from each vertex
    for v in P.vertices
        dfs(v, [], chains)
    end

    return chains
end

function posetChains(P, k=3)
    chains = []
    rel = rel_dict(P)

    # Function to perform DFS and find all chains up to length k
    function dfs(v, path, chains, depth)
        push!(path, v) 
        if depth ≤ k 
            push!(chains, copy(path))
        end
        if depth < k 
            for neighbor in rel[v]  
                dfs(neighbor, path, chains, depth + 1)
            end
        end
        pop!(path)                   
    end

    # Find all chains starting from each vertex
    for v in P.vertices
        dfs(v, [], chains, 1)
    end

    return chains
end

#= Function to convert a poset to a simplicial complex
    Inputs: 
        1) the poset P
    Outputs: 
        1) simplicial complex of the poset
    
    Notes: 
        
=#
function poset_to_simplicial(P)
    return simplicialComplex(P.vertices,posetChains(P))
end

#= Function to determine the dimension of intersection of two sets
    Inputs: 
        1) the first list
        2) the second list
    Outputs: 
        1) the dimension of intersection
    
    Notes: 
        
=#
function intersectionDim(list1, list2)
    set1 = Set(list1)
    set2 = Set(list2)
    intersection = intersect(set1, set2)
    return length(intersection)
end

#= Function to convert a simplicial complex to a graph
    Inputs: 
        1) the simplicial complex K
        2) the dimension of simplices to be vertices
        3) the codimension of a shared face to be an edge
    Outputs: 
        1) the resulting graph
    
    Notes: 
        
=#
function simplicial_to_graph(K,n,k) # m = dim you want vertices to be, n = dim of shared face for edge
    vert=[]
    edge=[]

    # add vertices
    for x in K.simplices 
        if length(x)==n+1
            push!(vert,x)
        end
    end

    for x in vert
        for y in vert
            if intersectionDim(x,y) >= n-k+1
                push!(edge,(x,y))
            end
        end
    end

    return graph(vert,edge)
end

#= Function to convert a poset to a graph
    Inputs: 
        1) the poset
        2) the dimension of simplices to be vertices
        3) the codimension of a shared face to be an edge
    Outputs: 
        1) the resulting graph
    
    Notes: 
        - k=0 means shared n-1 dimension face becomes edge
        - k=1 means n-2 ...
=#
function poset_to_graph(P,n,k)
    K=poset_to_simplicial(P)
    return simplicial_to_graph(K,n,k)
end

#= Function to make a relation transitive
    Inputs: 
        1) the relation
    Outputs: 
        1) the transitive relation
    
    Notes: 
        
=#
function make_transitive(relation)
    # Extract all unique elements from the pairs to form the set
    elements = unique(vcat([a for (a, b) in relation], [b for (a, b) in relation]))

    # Create a dictionary to map elements to indices
    index_map = Dict(element => i for (i, element) in enumerate(elements))

    # Initialize the adjacency matrix
    n = length(elements)
    adj_matrix = falses(n, n)

    # Populate the adjacency matrix based on the relation
    for (a, b) in relation
        adj_matrix[index_map[a], index_map[b]] = true
    end

    # Apply the Floyd-Warshall algorithm to compute the transitive closure
    for k in 1:n
        for i in 1:n
            for j in 1:n
                adj_matrix[i, j] = adj_matrix[i, j] || (adj_matrix[i, k] && adj_matrix[k, j])
            end
        end
    end

    # Convert the adjacency matrix back to a set of pairs
    transitive_relation = []
    for i in 1:n
        for j in 1:n
            if adj_matrix[i, j]
                push!(transitive_relation, (elements[i], elements[j]))
            end
        end
    end

    return transitive_relation
end

function add_numbers(a::Int, b::Int)
        return a + b
    end

function create_poset(vertices::Array, relations::Array)
        return poset(vertices, relations)
    end
