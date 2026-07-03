indexin_net(n::PhyloNetworks.Node, net::HybridNetwork) = findfirst(x -> x === n, net.node)
indexin_net(e::PhyloNetworks.Edge, net::HybridNetwork) = findfirst(x -> x === e, net.edge)

"""
    edgenode_coordinates(
        net::HybridNetwork,
        useedgelength::Bool,
        usedirecthybridline::Bool,
        curved::Symbol=:none,
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
   * Each major edge is drawn as a horizontal line by default. But when
     `majorcurved` is true, major hybrid edges are drawn curved instead.
   * Minor hybrid edges are drawn as:
     + a single diagonal segment (or curve) if `usedirecthybridline` is true,
     + or as 2 connected segments otherwise: one horizontal (whose length on
       the x axis can be used to represent the edge length), and the other
       diagonal (or curved) to connect the horizontal segment to the child node.

   `edge_*` contains the coordinates for the horizontal segment only, which is
   reduced to a single point (Begin = End) when `usedirecthybridline` is true.
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
    usedirecthybridline::Bool,
    curved::Symbol=:none,
    bend::Real=0.3,
    preorder::Bool=true,
)
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

    majorcurved = curved==:both

    # determine y for each node = y of its parent edge: post-order traversal
    # also [yB,yE] for each internal node: range of y's of all children nodes
    # y max is the numtaxa + number of minor edges
    ymin = 1.0;
    ymax = net.numtaxa
    if !usedirecthybridline
        ymax += sum(!e.ismajor for e in net.edge)
    end

    node_y  = zeros(Float64, net.numnodes) # order: in net.nodes, *!not in vec_node!*
    node_yB = zeros(Float64,net.numnodes) # min (B=begin) and max (E=end)
    node_yE = zeros(Float64,net.numnodes) #   of at children's nodes
    edge_yB = zeros(Float64,net.numedges) # yE of edge = y of child node
    edge_yE = Vector{Float64}(undef, net.numedges)
    # set node_y of leaves: follow cladewise order
    # also sets edge_yB of minor hybrid edges
    node_w  = zeros(Int, net.numnodes) # weight: number of descendant tips
    nexty = ymax # first tips at the top, last at bottom
    cladewise_queue = copy(getroot(net).edge) # the child edges of root
    # print("queued the root's children's indices: "); @show queue
    while !isempty(cladewise_queue)
        cur_edge = pop!(cladewise_queue); # deliberate choice over shift! for cladewise order
        # increment spacing and add to node_y if leaf
        cur_child = getchild(cur_edge)
        if cur_child.leaf
            ni = indexin_net(cur_child, net)
            node_y[ni]  = nexty
            node_yB[ni] = nexty
            node_yE[ni] = nexty
            nexty -= 1
            node_w[ni] = 1
        end

        # only for new hybrid lines:
        # increment spacing and add to edge_yB if parent edge is minor
        if !cur_edge.ismajor && !usedirecthybridline
            ei = indexin_net(cur_edge, net)
            edge_yB[ei] = nexty
            edge_yE[ei] = nexty
            nexty -= 1
        end

        # push children edges if this is a major edge:
        if cur_edge.ismajor
            for e in cur_child.edge
                if getparent(e) === cur_child # don't go backwards
                    push!(cladewise_queue, e)
                end
            end
        end
    end

    # set node_y of internal nodes: follow post-order
    for i=length(net.node):-1:1
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
        for i=length(net.node):-1:1 # post-order traversal
            if net.vec_node[i].leaf continue; end
            ni = indexin_net(net.vec_node[i], net)
            for e in net.vec_node[i].edge # loop over children only
                if net.vec_node[i] == (e.ischild1 ? e.node[2] : e.node[1])
                    node_age[ni] = max(node_age[ni], 1 +
                     node_age[indexin_net(getchild(e), net)])
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
    minoredge_xC = Union{Missing,Float64}[] # control point for Bézier curve
    minoredge_yB = Float64[]
    minoredge_yE = Float64[]
    minoredge_yC = Union{Missing,Float64}[] # missing for straight line

    verticalhybridlines = !usedirecthybridline && !useedgelength
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
        if verticalhybridlines
            push!(minoredge_xC, missing)
            push!(minoredge_yC, missing)
        else
            cx, cy = quadraticbezier_control(x0, y0, x2, y2, bend)
            push!(minoredge_xC, cx)
            push!(minoredge_yC, cy)
        end
        #@show i; @show net.edge[i]; @show pni; @show net.node[pni]; @show cni; @show net.node[cni]
    end

    xmax = max(xmax, edge_xE...)

    #@show node_x;  @show node_yB; @show node_y;  @show node_yE
    #@show edge_xB; @show edge_xE; @show edge_yB; @show edge_yE
    @show minoredge_xC; @show minoredge_yC
    return edge_xB, edge_xE, edge_yB, edge_yE,
           node_x, node_y, node_yB, node_yE,
           minoredge_xB, minoredge_xE, minoredge_xC,
           minoredge_yB, minoredge_yE, minoredge_yC,
           xmin, xmax, ymin, ymax
end


"""
    quadraticbezier_control(x0, y0, x2, y2, bend, rtol=1e-10)

Coordinates P1 = (x0,y2), to serve as middle control point for a quadratic
Bézier curve between the anchor points P0 = (x0,y0) and P2 = (x2,y2).

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
            return ((x0+x2)/2, y2 + bend)
        end
    end
    if dx < rtol * dy
        return (missing,missing)
    end
    return (x0,y2)
end

const _HybridSegment = NTuple{4, Float64}

"""`true` when two chords share the same start and end."""
@inline function _same_chord(seg_a::_HybridSegment, seg_b::_HybridSegment)
    x0, y0, x2, y2 = seg_a
    px0, py0, px2, py2 = seg_b
    return abs(px0 - x0) < 1e-10 && abs(py0 - y0) < 1e-10 &&
           abs(px2 - x2) < 1e-10 && abs(py2 - y2) < 1e-10
end

"""
    _segs_overlap(seg_a, seg_b; tol=1e-8)

Return `true` when `seg_a` and `seg_b` are collinear (both endpoints of `seg_b`
lie within perpendicular distance `tol` of the line through `seg_a`) and their
projections onto the dominant axis overlap. Detects both full and partial overlap
between major and minor hybrid edge segments.
"""
function _segs_overlap(seg_a::_HybridSegment, seg_b::_HybridSegment; tol::Float64=1e-8)
    x0a, y0a, x2a, y2a = seg_a
    x0b, y0b, x2b, y2b = seg_b
    dxa = x2a - x0a;  dya = y2a - y0a
    len_a = sqrt(dxa^2 + dya^2)
    len_b = sqrt((x2b - x0b)^2 + (y2b - y0b)^2)
    (len_a < 1e-6 || len_b < 1e-6) && return false
    # Both endpoints of seg_b must lie on the line through seg_a
    abs(dxa * (y0b - y0a) - dya * (x0b - x0a)) > tol * len_a && return false
    abs(dxa * (y2b - y0a) - dya * (x2b - x0a)) > tol * len_a && return false
    # Collinear: check range overlap on the dominant axis
    if abs(dxa) >= abs(dya)
        lo_a, hi_a = minmax(x0a, x2a)
        lo_b, hi_b = minmax(x0b, x2b)
        return lo_a <= hi_b + tol && lo_b <= hi_a + tol
    else
        lo_a, hi_a = minmax(y0a, y2a)
        lo_b, hi_b = minmax(y0b, y2b)
        return lo_a <= hi_b + tol && lo_b <= hi_a + tol
    end
end

"""Signed side of point `(qx, qy)` relative to the chord direction."""
@inline function _chord_side(dx::Float64, dy::Float64,
                             mx::Float64, my::Float64, qx::Float64, qy::Float64)
    return dx * (qy - my) - dy * (qx - mx)
end

"""
    _hybrid_bow_sign(seg; offset, partner=nothing)

Choose `force_bow_sign` for a hybrid Bézier chord `seg = (x0, y0, x2, y2)`.
When the chord is the same as the partner chord, fans them in opposite directions.
Otherwise bows left (away from right-side taxa).
"""
function _hybrid_bow_sign(seg::_HybridSegment; offset::Float64,
                          partner::Union{Nothing, _HybridSegment}=nothing)
    x0, y0, x2, y2 = seg
    if partner !== nothing && _same_chord(seg, partner)
        px0, py0, px2, py2 = partner
        pcx, pcy = quadraticbezier_control(px0, py0, px2, py2, 0.3)
        dx, dy = x2 - x0, y2 - y0
        mx, my = (x0 + x2) / 2, (y0 + y2) / 2
        side = _chord_side(dx, dy, mx, my, pcx, pcy)
        return side >= 0 ? -1.0 : 1.0
    end
    return _bow_left_sign(x0, y0, x2, y2, offset)
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
    minoredge_yE::Array{Float64,1}
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
        edf[j,:y] = (edge_yB[i] + edge_yE[i])/2
        edf[j,:x] = (edge_xB[i] + edge_xE[i])/2
        if style == :majortree && !ee.ismajor
            edf[j,:y] = (minoredge_yB[imh] + minoredge_yE[imh])/2
            edf[j,:x] = (minoredge_xB[imh] + minoredge_xE[imh])/2
            imh += 1
        end
        j += 1
    end
    # @show edf
    return labeledges, edf
end
