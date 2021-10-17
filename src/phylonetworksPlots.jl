# auxiliary functions:
# calculate plotting coordinates of nodes & edges for PhyloNetworks objects

"""
    getEdgeNodeCoordinates(net, useEdgeLength::Bool, useSimpleHybridLines::Bool)

Calculate coordinates for plotting later with Gadfly or RCall.

Actually modifies some (minor) attributes of the network,
as it calls `directEdges!`, `preorder!` and `cladewiseorder!`.

output: tuple with the following elements, in which the order of
nodes corresponds to the order in `net.node`, and the order of
edges corresponds to that in `net.edge` (filtered to minor edges as needed).

1. `edge_xB`: x coordinate for the Beginning and ...
2. `edge_xE`: ...  End of each edge, in the same order as in `net.edge`
3. `edge_yB`: y coordinate for edges, Begin ...
4. `edge_yE`: ... and End.
   Each major edge is drawn as a single horizontal line. Minor hybrid edges are
   drawn as: a single diagonal segment if `useSimpleHybridLines` is true,
   or as 2 connected segments otherwise: one horizontal (whose length on the
   x axis can be used to represent the edge length), and the other diagonal to
   connect the horizontal segment to the child node.
   edge_* contains the coordinates for the horizontal segment only, which is
   reduced to a single point (Begin = End) when using "SimpleHybridLines".
   minoredge_* (see below) contains information for the diagonal segment.
   Agreed, edge_yB = edge_yE always (relic from before v0.3:
   no minoredge output back then, and simple diagonal lines only)
5. `node_x`: x and ...
6. `node_y`: ... y coordinate at the middle of the vertical bar that represents a node.
   The (or each) parent edge of the node connects to this middle point,
   but the node itself is drawn as a vertical bar connected to all it children edges.
   order: same as in `net.node`
7. `node_yB`: y coordinates of the Beginning and the ...
8. `node_yE`: ... End of the vertical bar representing each node.
   The x coordinate (Begin & End) of the end points of the vertical bar is
   the same as that of the mid-point, given by `node_x`.
9. `minoredge_xB`: x coordinate for the Beginning and ...
10. `minoredge_xE`: ... End of the diagonal segment of each minor hybrid edge,
   in the same order as in `filter(e -> !e.isMajor, net.edge)`.
11. `minoredge_yB`: y coordinate for the beginning and ...
12. `minoredge_yE`: ... end of the diagonal segment of each minor hybrid edge.
13-16. `xmin`, `xmax`, `ymin`, `ymax`: ranges for the x and y axes.
"""
function getEdgeNodeCoordinates(net::HybridNetwork, useEdgeLength::Bool, useSimpleHybridLines::Bool)
    try
        directEdges!(net)   # to update isChild1
    catch e
        if isa(e, PhyloNetworks.RootMismatch)
            e = PhyloNetworks.RootMismatch( e.msg * "\nPlease change the root, perhaps using rootatnode! or rootatedge!")
        end
        rethrow(e)
    end
    preorder!(net)       # to update net.nodes_changed: true pre-ordering

    # determine y for each node = y of its parent edge: post-order traversal
    # also [yB,yE] for each internal node: range of y's of all children nodes
    # y max is the numTaxa + number of minor edges
    ymin = 1.0;
    ymax = net.numTaxa
    if !useSimpleHybridLines
        ymax += sum(!e.isMajor for e in net.edge)
    end

    node_y  = zeros(Float64, net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    node_yB = zeros(Float64,net.numNodes) # min (B=begin) and max (E=end)
    node_yE = zeros(Float64,net.numNodes) #   of at children's nodes
    edge_yB = zeros(Float64,net.numEdges) # yE of edge = y of child node
    # set node_y of leaves: follow cladewise order
    # also sets edge_yB of minor hybrid edges
    nexty = ymax # first tips at the top, last at bottom
    cladewise_queue = copy(net.node[net.root].edge) # the child edges of root
    # print("queued the root's children's indices: "); @show queue
    while !isempty(cladewise_queue)
        cur_edge = pop!(cladewise_queue); # deliberate choice over shift! for cladewise order
        # increment spacing and add to node_y if leaf
        if getChild(cur_edge).leaf
            node_y[findfirst(x->x===getChild(cur_edge), net.node)] = nexty
            nexty -= 1
        end

        # only for new hybrid lines:
        # increment spacing and add to edge_yB if parent edge is minor
        if !cur_edge.isMajor && !useSimpleHybridLines
            edge_yB[findfirst(x->x===cur_edge, net.edge)] = nexty
            nexty -= 1
        end

        # push children edges if this is a major edge:
        if cur_edge.isMajor
            for e in getChild(cur_edge).edge
                if getParent(e) === getChild(cur_edge)# don't go backwards
                    push!(cladewise_queue, e)
                end
            end
        end
    end

    # set node_y of internal nodes: follow post-order
    for i=length(net.node):-1:1
        nn = net.nodes_changed[i]
        !nn.leaf || continue # previous loop took care of leaves
        ni = findfirst(x -> x===nn, net.node)
        node_yB[ni]=ymax; node_yE[ni]=ymin;
        minor_yB  = ymax; minor_yE  = ymin;
        nomajorchild=useSimpleHybridLines # only use this var if using simple hybrid lines
        for e in nn.edge
            if nn == PhyloNetworks.getParent(e) # if e = child of node
                if useSimpleHybridLines
                    # old simple hybrid lines
                    if e.isMajor || nomajorchild
                        cc = PhyloNetworks.getChild(e)
                        yy = node_y[findfirst(x -> x===cc, net.node)]
                        yy!==nothing || error("oops, child $(cc.number) has not been visited before node $(nn.number).")
                    end
                    if e.isMajor
                        nomajorchild = false # we found a child edge that is a major edge
                        node_yB[ni] = min(node_yB[ni], yy)
                        node_yE[ni] = max(node_yE[ni], yy)
                    elseif nomajorchild # e is minor edge, and no major found so far
                        minor_yB = min(minor_yB, yy)
                        minor_yE = max(minor_yE, yy)
                    end
                else
                    # new pretty hybrid lines
                    if e.isMajor
                        cc = PhyloNetworks.getChild(e)
                        child_y = node_y[findfirst(x -> x===cc, net.node)]
                        child_y!==nothing || error("oops, child $(cc.number) has not been visited before node $(nn.number).")
                    else
                        child_y = edge_yB[findfirst(x->x===e, net.edge)]
                    end
                    node_yB[ni] = min(node_yB[ni], child_y)
                    node_yE[ni] = max(node_yE[ni], child_y)
                end
            end
        end
        if nomajorchild # if all children edges are minor hybrid edges. If so: not level 1, not tree-child
            if minor_yB == minor_yE # one single child. jitter by 0.1 to make the plot readable
                minor_yB += (minor_yB < (ymax+ymin)/2 ? 0.1 : -0.1)
                minor_yE = minor_yB
            end
            node_yB[ni] = minor_yB
            node_yE[ni] = minor_yE
        end
        node_y[ni] = (node_yB[ni]+node_yE[ni])/2
        if nomajorchild #since the minor edges are leaving from the center of the node's y pos.
            node_yB[ni] = node_y[ni]
            node_yE[ni] = node_y[ni]
        end
    end

    # setting branch lengths for plotting
    elenCalculate = !useEdgeLength
    if (useEdgeLength)
        allBLmissing = true; nonBLmissing = true;
        for e in net.edge
            if (nonBLmissing && e.length==-1.0) nonBLmissing=false; end
            if (allBLmissing && e.length!=-1.0) allBLmissing=false; end
        end
        if (allBLmissing)
            println("All edge lengths are missing, won't be used for plotting.")
            elenCalculate = true
        end
        if (!nonBLmissing && !allBLmissing) # not all, but some are missing
            @warn "At least one non-missing edge length: plotting any missing length as 1.0"
        end
    end
    elen = Float64[] # edge lengths to be used for plotting. same order as net.edge.
    if (elenCalculate)
        # setting elen such that the age of each node = 1 + age of oldest child
        # (including minor hybrid edges): need true post-ordering.
        # calculating node ages first, elen will be calculated later.
        elen     = zeros(Float64,net.numEdges)
        node_age = zeros(Float64,net.numNodes)
        for i=length(net.node):-1:1 # post-order traversal
            if (net.nodes_changed[i].leaf) continue; end
            ni = findfirst(x -> x===net.nodes_changed[i], net.node)
            for e in net.nodes_changed[i].edge # loop over children only
                if net.nodes_changed[i] == (e.isChild1 ? e.node[2] : e.node[1])
                    node_age[ni] = max(node_age[ni], 1 +
                     node_age[findfirst(x -> x=== PhyloNetworks.getChild(e), net.node)])
                end
            end
        end
    else
        for e in net.edge
            push!(elen, (e.length==-1.0 ? 1.0 : e.length))
        end
    end

    # determine xB,xE for each edge: pre-order traversal, uses branch lengths
    # then x and yB,yE for each node: x=xE of parent edge
    xmin = 1.0; xmax=xmin
    node_x  = zeros(Float64,net.numNodes) # order: in net.nodes, *!not in nodes_changed!*
    edge_xB = zeros(Float64,net.numEdges) # min (B=begin) and max (E=end)
    edge_xE = zeros(Float64,net.numEdges) # xE-xB = edge length
    node_x[net.root] = xmin # root node: x=xmin=0
    for i=2:length(net.node)              # true pre-order, skipping the root (i=1)
        ni = findfirst(x -> x===net.nodes_changed[i], net.node)
        ei = nothing # index of major parent edge of current node
        for e in net.nodes_changed[i].edge
            if (e.isMajor && net.nodes_changed[i] == e.node[e.isChild1 ? 1 : 2]) # major parent edge
                ei = findfirst(x -> x===e, net.edge)
                break
            end
        end
        ei !== nothing || error("oops, could not find major parent edge of node number $ni.")
        edge_yB[ei] = node_y[ni]
        pni = findfirst(x -> x===PhyloNetworks.getParent(net.edge[ei]), net.node) # parent node index
        edge_xB[ei] = node_x[pni]
        if elenCalculate
            elen[ei] = node_age[pni] - node_age[ni]
        end
        edge_xE[ei] = edge_xB[ei] + elen[ei]
        node_x[ni] = edge_xE[ei]
    end
    edge_yE = copy(edge_yB) # true for tree and major edges

    # coordinates of the diagonal lines that connect hybrid edges with their targets
    minoredge_xB = Float64[]
    minoredge_xE = Float64[]
    minoredge_yB = Float64[]
    minoredge_yE = Float64[]

    for i=1:net.numEdges
        if (!net.edge[i].isMajor) # minor hybrid edges
            # indices of child and parent nodes
            cni = findfirst(x -> x===PhyloNetworks.getChild( net.edge[i]), net.node)
            pni = findfirst(x -> x===PhyloNetworks.getParent(net.edge[i]), net.node)

            edge_xB[i] = node_x[pni]
            edge_xE[i] = useSimpleHybridLines ? edge_xB[i] : (useEdgeLength ? edge_xB[i] + elen[i] : node_x[cni])

            if useSimpleHybridLines
                edge_yB[i] = node_y[pni]
            end
            edge_yE[i] = edge_yB[i]

            push!(minoredge_xB, edge_xE[i])
            push!(minoredge_yB, edge_yE[i])
            push!(minoredge_xE, node_x[cni])
            push!(minoredge_yE, node_y[cni])
            #@show i; @show net.edge[i]; @show pni; @show net.node[pni]; @show cni; @show net.node[cni]
        end
    end

    xmax = max(xmax, edge_xE...)

    #@show node_x;  @show node_yB; @show node_y;  @show node_yE
    #@show edge_xB; @show edge_xE; @show edge_yB; @show edge_yE
    return edge_xB, edge_xE, edge_yB, edge_yE,
           node_x, node_y, node_yB, node_yE,
           minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE,
           xmin, xmax, ymin, ymax
end


"""
    checkNodeDataFrame(net, nodeLabel)

Check data frame for node annotations:
- check that the data has at least 2 columns (if it has any)
- check that the first column has integers (to serve as node numbers)
- remove rows with no node numbers
- warning if some node numbers in the data are not in the network.
"""
function checkNodeDataFrame(net::HybridNetwork, nodeLabel::DataFrame)
    labelnodes = size(nodeLabel,1)>0
    if (labelnodes && (size(nodeLabel,2)<2 ||
            !(nonmissingtype(eltype(nodeLabel[!,1])) <: Integer)))
        @warn "nodeLabel should have 2+ columns, the first one giving the node numbers (Integer)"
        labelnodes = false
    end
    if labelnodes # remove rows with no node number, check if at least one row remains
        filter!(row->!ismissing(row[1]), nodeLabel)
        labelnodes = size(nodeLabel,1)>0
    end
    if labelnodes
      tmp = setdiff(nodeLabel[!,1], [n.number for n in net.node])
      if length(tmp)>0
        msg = "Some node numbers in the nodeLabel data frame are not found in the network:\n"
        for a in tmp msg *= string(" ",a); end
        @warn msg
      end
    end
    return(labelnodes, nodeLabel)
end

"""
    prepareNodeDataFrame(net, nodeLabel::DataFrame,
        showNodeNumber::Bool, showIntNodeLabel::Bool, labelnodes::Bool,
        node_x, node_y)

Make data frame for node annotation. `node_*` should be Float64 vectors.
`nodeLabel` should have columns as required by [`checkNodeDataFrame`](@ref)

Columns of output data frame:
- x, y: coordinates on the plots (from `node_*`)
- name: node name
- num: node number
- lab: node label
- lea: is leaf?
"""
function prepareNodeDataFrame(net::HybridNetwork, nodeLabel::DataFrame,
        showNodeNumber::Bool, showIntNodeLabel::Bool, labelnodes::Bool,
        node_x::Array{Float64,1}, node_y::Array{Float64,1})
    nrows = (showNodeNumber || showIntNodeLabel || labelnodes ? net.numNodes : net.numTaxa)
    ndf = DataFrame(:name => Vector{String}(undef,nrows),
        :num => Vector{String}(undef,nrows), :lab => Vector{String}(undef,nrows),
        :lea => Vector{Bool}(  undef,nrows), :x => Vector{Float64}( undef,nrows),
        :y => Vector{Float64}( undef,nrows), copycols=false)
    j=1
    for i=1:net.numNodes
    if net.node[i].leaf  || showNodeNumber || showIntNodeLabel || labelnodes
        ndf[j,:name] = net.node[i].name
        ndf[j,:num] = string(net.node[i].number)
        if labelnodes
          jn = findfirst(isequal(net.node[i].number), nodeLabel[!,1])
          ndf[j,:lab] = (jn===nothing || ismissing(nodeLabel[jn,2]) ? "" :  # node label not in table or missing
            (nonmissingtype(eltype(nodeLabel[!,2])) <: AbstractFloat ?
              @sprintf("%0.3g",nodeLabel[jn,2]) : string(nodeLabel[jn,2])))
        end
        ndf[j,:lea] = net.node[i].leaf # use this later to remove H? labels
        ndf[j,:y] = node_y[i]
        ndf[j,:x] = node_x[i]
        j += 1
    end
    end
    # @show ndf
    return(ndf)
end

"""
    prepareEdgeDataFrame(net, edgeLabel::DataFrame, mainTree::Bool,
        edge_xB, edge_xE, edge_yB, edge_yE,
        minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE)

Check data frame for edge annotation.
`edge_*`: Float64 vectors giving the coordinates for the beginning and end of edges.
Return data frame with columns
- x, y: coordinates on the plots
- len: node name
- gam: gamma (inheritance value)
- num: node number
- lab: node label
- hyb: is hybrid?
- min: is minor?
"""
function prepareEdgeDataFrame(net::HybridNetwork, edgeLabel::DataFrame, mainTree::Bool,
        edge_xB::Array{Float64,1}, edge_xE::Array{Float64,1},
        edge_yB::Array{Float64,1}, edge_yE::Array{Float64,1},
        minoredge_xB::Array{Float64,1}, minoredge_xE::Array{Float64,1},
        minoredge_yB::Array{Float64,1}, minoredge_yE::Array{Float64,1})
    nrows = net.numEdges - (mainTree ? net.numHybrids : 0)
    edf = DataFrame(:len => Vector{String}(undef,nrows),
        :gam => Vector{String}(undef,nrows), :num => Vector{String}(undef,nrows),
        :lab => Vector{String}(undef,nrows), :hyb => Vector{Bool}(undef,nrows),
        :min => Vector{Bool}(  undef,nrows), :x => Vector{Float64}(undef,nrows),
        :y  => Vector{Float64}(undef,nrows), copycols=false)
    labeledges = size(edgeLabel,1)>0
    if (labeledges && (size(edgeLabel,2)<2 ||
            !(nonmissingtype(eltype(edgeLabel[!,1])) <: Integer)))
        @warn "edgeLabel should have 2+ columns, the first one giving the edge numbers (Integer)"
        labeledges = false
    end
    if labeledges # remove rows with no edge number and check if at least one remains
        filter!(row->!ismissing(row[1]), edgeLabel)
        labeledges = size(edgeLabel,1)>0
    end
    if labeledges
      tmp = setdiff(edgeLabel[!,1], [e.number for e in net.edge])
      if length(tmp)>0
        msg = "Some edge numbers in the edgeLabel data frame are not found in the network:\n"
        for a in tmp msg *= string(" ",a); end
        @warn msg
      end
    end
    j=1
    for i = 1:length(net.edge)
        if (!mainTree || !net.edge[i].hybrid || net.edge[i].isMajor)
            edf[j,:len] = (net.edge[i].length==-1.0 ? "" : @sprintf("%0.3g",net.edge[i].length))
            # @sprintf("%c=%0.3g",'γ',net.edge[i].length)
            edf[j,:gam] = (net.edge[i].gamma==-1.0  ? "" : @sprintf("%0.3g",net.edge[i].gamma))
            edf[j,:num] = string(net.edge[i].number)
            if labeledges
              je = findfirst(isequal(net.edge[i].number), edgeLabel[!,1])
              edf[j,:lab] = (je===nothing || ismissing(edgeLabel[je,2]) ? "" :  # edge label not found in table
                (nonmissingtype(eltype(edgeLabel[!,2])) <: AbstractFloat ?
                  @sprintf("%0.3g",edgeLabel[je,2]) : string(edgeLabel[je,2])))
            end
            edf[j,:hyb] = net.edge[i].hybrid
            edf[j,:min] = !net.edge[i].isMajor
            minorIndex = 1;
            if net.edge[i].isMajor
                edf[j,:y] = (edge_yB[i] + edge_yE[i])/2
                edf[j,:x] = (edge_xB[i] + edge_xE[i])/2
            else
                edf[j,:y] = (minoredge_yB[minorIndex] + minoredge_yE[minorIndex])/2
                edf[j,:x] = (minoredge_xB[minorIndex] + minoredge_xE[minorIndex])/2
                minorIndex += 1
            end
            j += 1
        end
    end
    # @show edf
    return labeledges, edf
end

