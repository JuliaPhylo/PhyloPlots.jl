indexin_net(n::PhyloNetworks.Node, net::HybridNetwork) = findfirst(x -> x === n, net.node)
indexin_net(e::PhyloNetworks.Edge, net::HybridNetwork) = findfirst(x -> x === e, net.edge)

"""
    edgenode_coordinates(
        net::HybridNetwork,
        useedgelength::Bool,
        style::Symbol,
        majorcurved::Bool,
        preorder::Bool=true,
    )

Calculate coordinates of edges segments and node midpoints & segments,
that can be used later for plotting.

Actually modifies some (minor) attributes of the network, as it calls
`directedges!` and `preorder!`, unless with argument `preorder=false`.

output: tuple with the following elements, in which the order of
nodes corresponds to the order in `net.node`, and the order of
edges corresponds to that in `net.edge` (filtered to minor edges as needed).

1. `edge_xB`: x coordinate for the Beginning and ...
2. `edge_xE`: ...  End of each edge, in the same order as in `net.edge`
3. `edge_yB`: y coordinate for edges, Begin ...
4. `edge_yE`: ... and End.
   * Each major edge is drawn as a horizontal line if `majorcurved` is false,
     and curved otherwise. This choice affects the beginning `yB` of hybrid edges;
     and the `node_yB/E` interval of node segments, which covers the y positions
     of a node's children drawn as straight line only.
   * Minor hybrid edges are drawn as:
     + a single diagonal segment (straight or curved) in the `:fulltree` style,
     + or as 2 connected segments otherwise: one horizontal whose length on
       the x axis can be used to represent the edge length, and the other
       diagonal (straight or curved) to connect the horizontal segment
       to the child node.

   `edge_*` contains the coordinates for the horizontal segment only, which is
   reduced to a single point (Begin = End) with the `:fulltree` style.
   `minoredge_*` (see below) contains information for the diagonal segment.
   Agreed, `edge_yB` = `edge_yE` always (relic from before v0.3:
   no minoredge output back then, and simple diagonal lines only)
5. `node_x`: x and ...
6. `node_y`: ... y coordinate along the vertical bar that represents a node.
   The (or each) parent edge of the node connects to this central point,
   but the node itself is drawn as a vertical bar connected to all it children edges.
   order: same as in `net.node`
7. `node_yB`: y coordinates of the Beginning and the ...
8. `node_yE`: ... End of the vertical bar representing each node.
   The x coordinate (Begin & End) of the end points of the vertical bar is
   the same as that of the mid-point, given by `node_x`.
9. `minoredge_xB`: x coordinate for the Beginning and ...
10. `minoredge_xE`: ... End of the diagonal segment of each minor hybrid edge,
   in the same order as in `filter(e -> !e.ismajor, net.edge)`.
11. `minoredge_yB`: y coordinate for the beginning and ...
12. `minoredge_yE`: ... end of the diagonal segment of each minor hybrid edge.
13-16. `xmin`, `xmax`, `ymin`, `ymax`: ranges for the x and y axes.
"""
function edgenode_coordinates(
    net::HybridNetwork,
    useedgelength::Bool,
    style::Symbol,
    majorcurved::Bool,
    preorder::Bool=true,
)
    usedirecthybridline = style != :fulltree
    if preorder
      try
        directedges!(net)   # to update ischild1
      catch e
        if isa(e, PhyloNetworks.RootMismatch)
            e = PhyloNetworks.RootMismatch( e.msg * "\nPlease change the root, perhaps using rootatnode! or rootatedge!")
        end
        rethrow(e)
      end
      preorder!(net)       # to update net.vec_node: true pre-ordering
    end

    # determine y for each node = y of its parent edge: post-order traversal
    # also [yB,yE] for each internal node: range of y's of all children nodes
    # y max is the numtaxa + number of minor edges
    ymin = 1.0;
    ymax = net.numtaxa
    if style == :fulltree # reserve space for "corner" nodes: from minor edges
        ymax += sum(!e.ismajor for e in net.edge)
    elseif style == :lsatree
        # add number of internal nodes that are leaves in the LSA tree
        for n in net.node
            n.leaf && continue
            isleaf = true
            for e in n.edge
                getchild(e) === n && continue
                if !e.hybrid
                    isleaf = false
                    break
                end
            end
            if isleaf
                ymax += 1
            end
        end
    end

    node_y  = zeros(Float64, net.numnodes) # order: in net.nodes, *!not in vec_node!*
    node_yB = zeros(Float64,net.numnodes) # min (B=begin) and max (E=end)
    node_yE = zeros(Float64,net.numnodes) #   of at children's nodes
    edge_yB = zeros(Float64,net.numedges) # yE of edge = y of child node
    edge_yE = Vector{Float64}(undef, net.numedges)
    node_w = [(nn.leaf ? 1 : 0) for nn in net.node] # weight: number of descendant tips

    # set node_y of leaves: follow cladewise order along a spanning tree
    # also sets edge_yB of minor hybrid edges
    cladewise_node2children =  prepare_cladewiseorder(net, style)
    leafcoordinates_cladewiseorder!(
        (node_y, node_yB, node_yE, edge_yB, edge_yE), Ref(ymax), # modified
        cladewise_node2children, net, style==:fulltree)
    # fixit: check that it work below. then delete the commented-out block below
    #=
    nexty = ymax # first tips at the top, last at bottom
    cladewise_queue = copy(getroot(net).edge) # the child edges of root
    # print("queued the root's children's indices: "); @show queue
    while !isempty(cladewise_queue)
        cur_edge = pop!(cladewise_queue); # deliberate choice over shift! for cladewise order
        # increment spacing and add to node_y if leaf
        cur_child = getchild(cur_edge)
        ni = findfirst(x->x===cur_child, net.node)
        if cur_child.leaf # later: instead check if the node/index is not a key in the dictionary
            ni = indexin_net(cur_child, net)
            node_y[ni]  = nexty
            node_yB[ni] = nexty
            node_yE[ni] = nexty
            nexty -= 1
        end

        # only for new hybrid lines:
        if cur_edge.hybrid
            if cur_edge.ismajor && style == :lsatree
                # fixit: do *not* do this, unless cur_child is a leaf in the LSA tree
                node_y[ni]  = nexty
                nexty -= 1
            elseif !cur_edge.ismajor && !usedirecthybridline
            # increment spacing and add to edge_yB if parent edge is minor
                ei = indexin_net(cur_edge, net)
                edge_yB[ei] = nexty
                edge_yE[ei] = nexty
                nexty -= 1
            end
        end

        # push the appropriate children edges to the "queue"
        if cur_edge.ismajor
            for e in cur_child.edge
                if getparent(e) === cur_child # don't go backwards
                    if style != :lsatree || !e.hybrid
                        push!(cladewise_queue, e)
                    end
                end
            end
        end
       if style == :lsatree
        # to follow the LSA tree: push the major parent edge of h when we visit lsa(h).
        # 1. loop over hybrid2lsa: h_ni = hybrid node index, lsa_ni = its LSA node index
        if haskey(lsa2hybrid, ni)
            for h_ni in lsa2hybrid[ni]
                push!(cladewise_queue, getparentedge(net.node[h_ni]))
            end
        end
       end
    end
    =#
    @show node_y # only leaves should be non-zero

    # set node_y of internal nodes: follow post-order
    for i in length(net.node):-1:1
        nn = net.vec_node[i]
        !nn.leaf || continue # previous loop took care of leaves
        ni = indexin_net(nn, net)
        node_yB[ni]=ymax; node_yE[ni]=ymin
        minor_yB  = ymax; minor_yE  = ymin;
        nomajorchild = usedirecthybridline # only use this var if using simple hybrid lines
        for e in nn.edge
            if nn == getparent(e) # if e = child of node
                if usedirecthybridline # unbroken one-segment hybrid lines
                    if e.ismajor || nomajorchild
                        cc = getchild(e)
                        ci = indexin_net(cc, net)
                        cy = node_y[ci]
                        cy != 0 || error("child $(cc.number) has not been visited before node $(nn.number).")
                    end
                    if e.ismajor
                        nomajorchild = false # we found a child edge that is a major edge
                        if !majorcurved || !e.hybrid
                            node_yB[ni] = min(node_yB[ni], cy)
                            node_yE[ni] = max(node_yE[ni], cy)
                        end
                        node_y[ni] += node_w[ci] * cy # running average; was initialized at 0
                        node_w[ni] += node_w[ci]
                    elseif nomajorchild # e is minor edge, and no major found so far
                        minor_yB = min(minor_yB, cy)
                        minor_yE = max(minor_yE, cy)
                    end
                else
                    if e.ismajor
                        cc = getchild(e)
                        ci = indexin_net(cc, net)
                        cy = node_y[ci]
                        cy != 0 || error("child $(cc.number) has not been visited before node $(nn.number).")
                        node_y[ni] += node_w[ci] * cy
                        node_w[ni] += node_w[ci]
                    else
                        cy = edge_yB[indexin_net(e, net)]
                        node_y[ni] += cy
                        node_w[ni] += 1
                    end
                    if !majorcurved || !e.ismajor || !e.hybrid
                        node_yB[ni] = min(node_yB[ni], cy)
                        node_yE[ni] = max(node_yE[ni], cy)
                    end
                end
            end
        end
        if nomajorchild # children edges are all minor hybrids
            if minor_yB == minor_yE # one single child. jitter by 0.1 to make the plot readable
                minor_yB += (minor_yB < (ymax+ymin)/2 ? 0.1 : -0.1)
                minor_yE = minor_yB
            end
            node_yB[ni] = minor_yB
            node_yE[ni] = minor_yE
            node_y[ni]  = (minor_yB + minor_yE)/2
        else
            # node_y[ni] = (node_yB[ni]+node_yE[ni])/2 ## v2.1.0 and earlier
            node_y[ni] /= node_w[ni] # weight > 0 necessarily if !nomajorchild
        end
        if majorcurved
            node_yB[ni] = min(node_yB[ni], node_y[ni])
            node_yE[ni] = max(node_yE[ni], node_y[ni])
        end
        if nomajorchild #since the minor edges are leaving from the center of the node's y pos.
            node_yB[ni] = node_y[ni]
            node_yE[ni] = node_y[ni]
        end
    end

    # setting branch lengths for plotting
    elenCalculate = !useedgelength
    if useedgelength
        allBLmissing = true; nonBLmissing = true;
        for e in net.edge
            if (nonBLmissing && e.length==-1.0) nonBLmissing=false; end
            if (allBLmissing && e.length!=-1.0) allBLmissing=false; end
        end
        if allBLmissing
            println("All edge lengths are missing, won't be used for plotting.")
            elenCalculate = true
        end
        if !nonBLmissing && !allBLmissing # not all, but some are missing
            @warn "At least one non-missing edge length: plotting any missing length as 1.0"
        end
    end
    elen = Float64[] # edge lengths to be used for plotting. same order as net.edge.
    if elenCalculate
        # setting elen such that the age of each node = 1 + age of oldest child
        # (including minor hybrid edges): need true post-ordering.
        # calculating node ages first, elen will be calculated later.
        elen     = zeros(Float64,net.numedges)
        node_age = zeros(Float64,net.numnodes)
        for nn in Iterators.reverse(net.vec_node) # post-order traversal
            nn.leaf && continue
            ni = indexin_net(nn, net)
            for e in nn.edge # loop over children only
                nn == getparent(e) || continue
                node_age[ni] = max(node_age[ni],
                    1 + node_age[indexin_net(getchild(e), net)])
            end
        end
    else
        for e in net.edge
            push!(elen, (e.length==-1.0 ? 1.0 : e.length))
        end
    end

    # determine xB,xE for each edge: pre-order traversal, uses branch lengths
    # then x and yB,yE for each node: x=xE of parent edge
    xmin = 0.0; ## was 1.0 in v2.1.0 and earlier
    xmax = xmin
    node_x  = zeros(Float64,net.numnodes) # order: in net.nodes, *!not in vec_node!*
    edge_xB = zeros(Float64,net.numedges) # min (B=begin) and max (E=end)
    edge_xE = zeros(Float64,net.numedges) # xE-xB = edge length
    node_x[net.rooti] = xmin # root node: x=xmin=0
    for i in 2:length(net.node)           # true pre-order, skipping the root (i=1)
        ni = indexin_net(net.vec_node[i], net)
        ee = getparentedge(net.vec_node[i])
        ei = indexin_net(ee, net) # index of major parent edge of current node
        pni = indexin_net(getparent(ee), net) # parent node index
        edge_yE[ei] = node_y[ni]
        edge_yB[ei] = (majorcurved && ee.hybrid ? node_y[pni] : edge_yE[ei] )
        edge_xB[ei] = node_x[pni]
        if elenCalculate
            elen[ei] = node_age[pni] - node_age[ni]
        end
        edge_xE[ei] = edge_xB[ei] + elen[ei]
        node_x[ni] = edge_xE[ei]
    end

    # coordinates of the diagonal lines that connect hybrid edges with their targets
    minoredge_xB = Float64[]
    minoredge_xE = Float64[]
    minoredge_yB = Float64[]
    minoredge_yE = Float64[]

    for (i, e) in enumerate(net.edge) # minor hybrid edges: arrow (& trivial segment)
        e.ismajor && continue # skip major edges
        cni = indexin_net(getchild(e), net) # indices of child and parent nodes
        pni = indexin_net(getparent(e), net)
        edge_xB[i] = node_x[pni]
        edge_xE[i] = (usedirecthybridline ? edge_xB[i] :
                        (useedgelength ? edge_xB[i] + elen[i] : node_x[cni]))
        if usedirecthybridline # trivial horizontal segment: reduced to a point
            edge_yB[i] = node_y[pni]
            edge_yE[i] = edge_yB[i]
        end
        x0 = edge_xE[i];  y0 = edge_yE[i]
        x2 = node_x[cni]; y2 = node_y[cni]
        push!(minoredge_xB, x0); push!(minoredge_yB, y0)
        push!(minoredge_xE, x2); push!(minoredge_yE, y2)
        #@show i; @show net.edge[i]; @show pni; @show net.node[pni]; @show cni; @show net.node[cni]
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
    quadraticbezier_control(x0, y0, x2, y2, bend, rtol=1e-10)

Coordinates of point P1 to serve as middle control point for a quadratic
Bézier curve between the anchor points P0 = (x0,y0) and P2 = (x2,y2).
- P1 = (x0,y2) if x0 ≤ x2, which is most frequent
  (i.e. when not using branch lengths, or if the network is time-consistent).
- P1 = (x2,y0) if x0 > x2, which looks like the edge goes back in time.

If P1≈P0 (P0 → P2 is horizontal) or if P1≈P2 (P0 → P2 is vertical),
then the Bézier curve is almost straight. In this case:
- if P1≈P2 (vertical) or if `bend` is not 0, `(missing, missing)` is returned
- if `bend==0` and if P1≈P2 (the default Bézier curve would be horizontal),
  then `((x0+x2)/2, y2 + bend)` is returned to force a bend.
"""
function quadraticbezier_control(
    x0::Real,
    y0::Real,
    x2::Real,
    y2::Real,
    bend::Real,
    rtol=1e-10,
)
    dy = abs(y2 - y0) # d(P1,P0)
    dx = abs(x2 - x0) # d(P1,P2)
    if dy < rtol * dx
        if bend==0
            return (missing,missing)
        else
            return ((x0+x2)/2, y2 - bend)
        end
    end
    if dx < rtol * dy
        return (missing,missing)
    end
    return (x0>x2 ? (x2,y0) : (x0,y2))
end

"""
    quadraticbezier_midpoint(x0, y0, x2, y2, bend)

Midpoint around the Bézier curve from [`quadraticbezier_control`](@ref),
with `x = (x0+x2)/2`. In most cases, the Bézier control point is (x0,y2)
and the midpoint is at time `t=1/sqrt(2)`.
It is at `t=1/2` if a bent is forced to avoid a straight horizontal curve.
"""
function quadraticbezier_midpoint(x0, y0, x2, y2, bend)
    xmid = (x0 + x2)/2
    x1,y1 = quadraticbezier_control(x0,y0, x2,y2, bend)
    ymid = (ismissing(x1) ? (y0 + y2)/2 :
        (x1==x0 ? # then take Bézier at t =   1/sqrt(2)
            0.085786437626905*y0 + 0.4142135623730951*y1 + y2/2 :
        (y1==y0 ? # then take Bézier at t = 1-1/sqrt(2)
            0.085786437626905*y2 + 0.4142135623730951*y1 + y0/2 :
            y0/4 + y1/2 + y2/4))) # take Bézier at t=1/2
    return (xmid, ymid)
end

"""
    check_nodedataframe(net, nodelabel)

Check data frame for node annotations:
- check that the data has at least 2 columns (if it has any)
- check that the first column has integers (to serve as node numbers)
- remove rows with no node numbers
- warning if some node numbers in the data are not in the network.
"""
function check_nodedataframe(
    net::HybridNetwork,
    nodelabel::DataFrame
)
    labelnodes = size(nodelabel,1)>0
    if (labelnodes && (size(nodelabel,2)<2 ||
            !(nonmissingtype(eltype(nodelabel[!,1])) <: Integer)))
        @warn "nodelabel should have 2+ columns, the first one giving the node numbers (Integer)"
        labelnodes = false
    end
    if labelnodes # remove rows with no node number, check if at least one row remains
        nodelabel = filter(row->!ismissing(row[1]), nodelabel) # creates a copy: don't modify the user's df
        labelnodes = size(nodelabel,1)>0
    end
    if labelnodes
      tmp = setdiff(nodelabel[!,1], [n.number for n in net.node])
      if length(tmp)>0
        msg = "Some node numbers in the nodelabel data frame are not found in the network:\n"
        for a in tmp msg *= string(" ",a); end
        @warn msg
      end
    end
    return(labelnodes, nodelabel)
end

"""
    prepare_nodedataframe(net, nodelabel::DataFrame,
        shownodenumber::Bool, shownodename::Bool, labelnodes::Bool,
        node_x, node_y)

Make data frame for node annotation. `node_*` should be Float64 vectors.
`nodelabel` should have columns as required by [`check_nodedataframe`](@ref).
`shownodename` is to show the name of internal nodes.
Leaf names are always included.

Columns of output data frame:
- x, y: coordinates on the plots (from `node_*`)
- name: node name
- num: node number
- lab: node label
- lea: is leaf?
"""
function prepare_nodedataframe(
    net::HybridNetwork,
    nodelabel::AbstractDataFrame,
    shownodenumber::Bool,
    shownodename::Bool,
    labelnodes::Bool,
    node_x::Array{Float64,1},
    node_y::Array{Float64,1}
)
    nrows = (shownodenumber || shownodename || labelnodes ? net.numnodes : net.numtaxa)
    ndf = DataFrame(:name => Vector{String}(undef,nrows),
        :num => Vector{String}(undef,nrows), :lab => fill(""::String,nrows),
        :lea => Vector{Bool}(  undef,nrows), :x => Vector{Float64}( undef,nrows),
        :y => Vector{Float64}( undef,nrows), copycols=false)
    j=1
    for i=1:net.numnodes
    if net.node[i].leaf  || shownodenumber || shownodename || labelnodes
        ndf[j,:name] = net.node[i].name
        ndf[j,:num] = string(net.node[i].number)
        if labelnodes
          jn = findfirst(isequal(net.node[i].number), nodelabel[!,1])
          ndf[j,:lab] = (jn===nothing || ismissing(nodelabel[jn,2]) ? "" :  # node label not in table or missing
            (nonmissingtype(eltype(nodelabel[!,2])) <: AbstractFloat ?
              @sprintf("%0.3g",nodelabel[jn,2]) : string(nodelabel[jn,2])))
        end
        ndf[j,:lea] = net.node[i].leaf
        ndf[j,:y] = node_y[i]
        ndf[j,:x] = node_x[i]
        j += 1
    end
    end
    # @show ndf
    return(ndf)
end

"""
    prepare_edgedataframe(net, edgelabel::DataFrame, style::Symbol,
        edge_xB, edge_xE, edge_yB, edge_yE,
        minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE,
        curved, bend)

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
function prepare_edgedataframe(
    net::HybridNetwork,
    edgelabel::AbstractDataFrame,
    style::Symbol,
    edge_xB::Array{Float64,1},
    edge_xE::Array{Float64,1},
    edge_yB::Array{Float64,1},
    edge_yE::Array{Float64,1},
    minoredge_xB::Array{Float64,1},
    minoredge_xE::Array{Float64,1},
    minoredge_yB::Array{Float64,1},
    minoredge_yE::Array{Float64,1},
    curved::Symbol,
    bend::Real,
)
    nrows = net.numedges
    edf = DataFrame(:len => Vector{String}(undef,nrows),
        :gam => Vector{String}(undef,nrows), :num => Vector{String}(undef,nrows),
        :lab => fill(""::String,nrows),      :hyb => Vector{Bool}(undef,nrows),
        :min => Vector{Bool}(  undef,nrows), :x => Vector{Float64}(undef,nrows),
        :y  => Vector{Float64}(undef,nrows), copycols=false)
    labeledges = size(edgelabel,1)>0
    if (labeledges && (size(edgelabel,2)<2 ||
            !(nonmissingtype(eltype(edgelabel[!,1])) <: Integer)))
        @warn "edgelabel should have 2+ columns, the first one giving the edge numbers (Integer)"
        labeledges = false
    end
    if labeledges # remove rows with no edge number and check if at least one remains
        edgelabel = filter(row->!ismissing(row[1]), edgelabel) # creates a copy: don't modify the user's df
        labeledges = size(edgelabel,1)>0
    end
    if labeledges
      tmp = setdiff(edgelabel[!,1], [e.number for e in net.edge])
      if length(tmp)>0
        msg = "Some edge numbers in the edgelabel data frame are not found in the network:\n"
        for a in tmp msg *= string(" ",a); end
        @warn msg
      end
    end
    j=1   # index of row in edf
    imh=1 # index of minor hybrid edge in filter(ee -> !ee.ismajor, net.edge)
    for i = 1:length(net.edge)
        ee = net.edge[i]
        edf[j,:len] = (ee.length==-1.0 ? "" : @sprintf("%0.3g",ee.length))
        # @sprintf("%c=%0.3g",'γ',ee.length)
        edf[j,:gam] = (ee.gamma==-1.0  ? "" : @sprintf("%0.3g",ee.gamma))
        edf[j,:num] = string(ee.number)
        if labeledges
            je = findfirst(isequal(ee.number), edgelabel[!,1])
            edf[j,:lab] = (je===nothing || ismissing(edgelabel[je,2]) ? "" :  # edge label not found in table
            (nonmissingtype(eltype(edgelabel[!,2])) <: AbstractFloat ?
                @sprintf("%0.3g",edgelabel[je,2]) : string(edgelabel[je,2])))
        end
        edf[j,:hyb] = ee.hybrid
        edf[j,:min] = !ee.ismajor
        if ee.ismajor || style != :majortree # use first segment
            x0,y0, x2,y2 = (edge_xB[i], edge_yB[i], edge_xE[i], edge_yE[i])
        else # minor edge, use second segment (arrow)
            x0,y0, x2,y2 = (minoredge_xB[imh], minoredge_yB[imh],
                            minoredge_xE[imh], minoredge_yE[imh])
            imh += 1
        end
        if curved==:none || !ee.hybrid || (ee.ismajor && curved != :both) ||
                (!ee.ismajor && style != :majortree)
            edf[j,:x] = (x0 + x2)/2
            edf[j,:y] = (y0 + y2)/2
        else # mid-point depends on the Bézier control point
            edf[j,:x], edf[j,:y] = quadraticbezier_midpoint(x0,y0, x2,y2, (ee.ismajor ? 0 : bend))
        end
        j += 1
    end
    # @show edf
    return labeledges, edf
end

"""
    prepare_LSAtreetraversal(net::HybridNetwork)

Dictionaries `(hybrid2lsa, lsa2hybrid)`, containing indices of
hybrid nodes and, for each hybrid, the index of the LSA of its parents.
Indices are in `net.node`.
- `hybrid2lsa` maps each hybrid to its parents' LSA.
- `lsa2hybrid` maps each node that is the LSA of some hybrid to a
  vector listing these hybrids.

Also, the network is modified as follows:
- for each hybrid node `n` in `net`, `n.prev` stores its LSA node.
- `net.vec_int1` ... fixit

**Warning**: assume that `net` is already preordered, that is,
with its nodes listed in a preorder in `net.vec_node`
"""
function prepare_LSAtreetraversal(net::HybridNetwork)
    lsaM = leaststableancestor_matrix(net, false)[:all] # do not preorder again
    # get LSA(parents(n)) for each hybrid node n
    hybrid2lsa = Dict{Int,Int}()
    for (ni,nn) in enumerate(net.node)
        nn.hybrid || continue
        # switch to indices in net.vec_node
        lsa_i = findfirst(x->x===nn, net.vec_node)
        for e in nn.edge # loop over parents of nn only
            getchild(e) === nn || continue
            pi_inlsaM = findfirst(x->x===getparent(e), net.vec_node)
            newlsa = lsaM[lsa_i, pi_inlsaM]
            lsa_i = findfirst(x->x===newlsa, net.vec_node)
        end
        # back to indices in net.node
        lsa_node = net.vec_node[lsa_i]
        nn.prev = lsa_node
        push!(hybrid2lsa, ni => indexin_net(lsa_node, net))
    end
    lsa2hybrid = Dict{Int, Vector{Int}}()
    for (h_ni, lsa_ni) in hybrid2lsa
        if haskey(lsa2hybrid, lsa_ni)
            push!(lsa2hybrid[lsa_ni], h_ni)
        else
            lsa2hybrid[lsa_ni] = [h_ni]
        end
    end
    return (hybrid2lsa, lsa2hybrid)
end

"""
    prepare_cladewiseorder(net::HybridNetwork, style::Symbol)

Dictionary `d` mapping each internal node index to the vector of children indices,
where 'internal' and 'children' correspond to the tree defined by `net`
and the `style`; and where indices are in `net.node`.
- `:majortree` style: tree obtained by removing all minor hybrid edges in `net`
- `:fulltree` style: tree obtained by creating a new leaf for each minor hybrid edge,
  whose name is that of the hybrid node. A hybrid with `p` parents is then represented
  as an internal node (the hybrid node itself) along the major hybrid edge, and by
  `p-1` fake leaves (or "corners" in the plot) sharing the same name.
- `:lsatree`: tree obtained by removing all hybrid edges, then connecting each
  original hybrid node to the LSA of its parents in `net`.

For each node with index `ni`, `d[ni]` is the following tuple:
- (vector of children indices, `true`) for the major tree & LSA tree styles
- (vector of children indices, vector of children type) for the full tree style
  where the children type is a boolean: true if child via a major edge,
  false if child via a minor edge, that is, a corner / fake leaf.
  For such children, the index is that of the minor hybrid child *edge*,
  and this index is in `net.edge`.

**Warning**: assume that `net` is already preordered, that is,
with its nodes listed in a preorder in `net.vec_node`
"""
function prepare_cladewiseorder(net::HybridNetwork, style::Symbol)
    lsa2hybrid = (style == :lsatree ? prepare_LSAtreetraversal(net)[2] : nothing)
    fulltree = (style == :fulltree)
    childType = Tuple{Int,Bool}
    node2childvec = Dict{Int,Vector{childType}}()
    cladewise_stack = copy(getroot(net).edge) # the child edges of root
    if style == :lsatree && haskey(lsa2hybrid, net.rooti)
        for h_ni in lsa2hybrid[net.rooti]
            push!(cladewise_stack, getparentedge(net.node[h_ni]))
        end
    end
    while !isempty(cladewise_stack)
        cur_edge = pop!(cladewise_stack); # deliberate choice over shift! for cladewise order
        nn = getchild(cur_edge)
        ni = findfirst(x->x===nn, net.node)
        pn = (style == :lsatree && nn.hybrid ? nn.prev : getparent(cur_edge)) # parent in *tree*
        pni = indexin_net(pn, net)  # pni = parent node index
        # add ni to the list of pni's children
        if cur_edge.ismajor
            push!(get!(node2childvec, pni, childType[]), (ni, true))
        elseif fulltree
            ei = indexin_net(cur_edge, net)
            push!(get!(node2childvec, pni, childType[]), (ei, false))
        end
        # push the appropriate children edges to the stack
        if cur_edge.ismajor || fulltree
            for e in nn.edge
                if getparent(e) === nn # don't go backwards
                    if style != :lsatree || !e.hybrid
                        push!(cladewise_stack, e)
                    end
                end
            end
        end
        if style == :lsatree  && haskey(lsa2hybrid, ni)
            # nn = lsa(h) for some hybrid h: push major parent edge of h
            for h_ni in lsa2hybrid[ni]
                push!(cladewise_stack, getparentedge(net.node[h_ni]))
            end
       end
    end
    return node2childvec
end

"""
    leafcoordinates_cladewiseorder!(node_edge_y, fixit..., style)

Set the node y-axis coordinates of leaves, and the edge y-axis coordinates of
minor diagonal segments (placement of "corners", or fake leaves, for
minor hybrid edges) under the full tree style.

The first 2 arguments are modified.

fixit
"""
function leafcoordinates_cladewiseorder!(
    node_edge_y,
    nexty::Base.RefValue,
    cladewisedict::Dict,
    net::HybridNetwork,
    fulltree::Bool
)
    leafcoordinates_cladewiseorder!(node_edge_y, nexty, net.rooti, cladewisedict, net, fulltree)
end
function leafcoordinates_cladewiseorder!(
    node_edge_y,
    nexty::Base.RefValue,
    pni::Int, # node index in net.node, and in y vectors
    cladewisedict::Dict,
    net::HybridNetwork,
    fulltree::Bool,
)
    node_y, node_yB, node_yE, edge_yB, edge_yE = node_edge_y
    if haskey(cladewisedict, pni) # parent node: not a leaf in traversal tree
        for (ni,ismajor) in cladewisedict[pni]
            if ismajor
                leafcoordinates_cladewiseorder!(node_edge_y, nexty, ni, cladewisedict, net, fulltree)
            else # fake leaf, corner edge: stop recursion, assign next y
                @info "parent index: $pni, minor child edge index: $ni"
                edge_yB[ni] = nexty[]
                edge_yE[ni] = nexty[]
                nexty[] -= 1
            end
        end
    else # leaf in the traversal tree: stop recursion; assign the next y
        @info "node index: $pni, nexty = $(nexty[])"
        node_y[ pni] = nexty[]
        node_yB[pni] = nexty[]
        node_yE[pni] = nexty[]
        nexty[] -= 1
    end
end
