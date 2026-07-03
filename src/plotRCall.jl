"""
    plot(net::HybridNetwork; ...)

Plot a network with edges going from left to right and taxa (leaves) placed on
the right, using R graphics. Optional arguments are listed below.

## lines forming the network:

- `useedgelength = false`: if true, the tree edges and major hybrid edges are
  drawn proportionally to their length. Minor hybrid edges are not, however.
  Note that edge lengths in coalescent units may scale very poorly with time.
- `style = :fulltree`: symbol indicating the style of the diagram
  * `:majortree` will simply draw minor edges onto the major tree.
  * `:fulltree` will draw minor edges as their own branches in the tree,
    in the same style used by [icytree](https://icytree.org). This is
    useful for overlapping or confusing networks.
- `arrowlen`: the length of the arrow tips in the full tree style.
  The default is 0.1 if `style = :fulltree`,
  and 0 if `style = :majortree` (making the arrows appear as segments).
- `minorlinetype`: type of lines used for minor edges, represented by arrows.
  Default is "solid" under the major-tree style, and "longdash" under the
  full tree style.
- `edgewidth=1`: width of horizontal (not diagonal) edges. To vary them,
  use a dictionary to map the number of each edge to its desired width.
- `curved = :none`: curvature for hybrid edges (`:none`, `:minor`, `:both`).
- `bend = 0.3`: y-bend for a downward curvature of minor hybrid edges that
  would otherwise be horizontal, to avoid overlap with other edges.
  Only used when `curved` is requested.
- `xlim`, `ylim`: array of 2 values, to determine the axes limits.

## tip annotations:

- `showtiplabel = true`: if true, taxon labels (names) are shown.
- `tipoffset = 0`: to offset tip labels.
- `tipcex = 1`: character expansion for tip and internal node names.

## nodes & edges annotations:

- `shownodelabel = false`: if true, internal nodes are labelled with their names.
  Useful for hybrid nodes, which have labels such as "H1".
- `shownodenumber = false`: if true, nodes are labelled with the number used internally.
- `showedgenumber = false`: if true, edges are labelled with the number used internally.
- `showedgelength = false`: if true, edges are labelled with their length (above).
- `showgamma = false`: if true, hybrid edges are labelled with their heritability (below).
- `edgelabel = DataFrame()`: dataframe with two columns: the first with edge numbers,
  the second with labels (like bootstrap values) to annotate edges. empty by default.
- `nodelabel = DataFrame()`: dataframe with two columns: the first with node numbers,
  the second with labels (like bootstrap values for hybrid relationships)
  to annotate nodes. empty by default.
- `nodecex = 1` and `edgecex = 1`: character expansion for labels in the
  `nodelabel` and `edgelabel` data frames.
- `nodelabeladj = 1` and `edgelabeladj = [.5,0]`: position adjustment to place
   the labels from `edgelabel` and `nodelabel` data frames.

## colors:

- `edgecolor = "black"`: color for tree edges if a single color is given. To
   vary edge colors, use a dictionary to map the number of each edge to its
   desired color. In that case, unmapped edges are given the `defaultedgecolor`.
- `majorhybridedgecolor = "deepskyblue4"`: color for major hybrid edges.  
  `minorhybridedgecolor = "deepskyblue"`: color for minor hybrid edges.  
   Both are ignored if `edgecolor` is a dictionary.
- `defaultedgecolor`: default is "black" if `edgecolor` is a dictionary, or
  `edgecolor` otherwise. Used to draw the segments representing each node.
- `edgenumbercolor = "grey"`: color for edge numbers.
- `edgelabelcolor = "black"`: color for labels in the `edgelabel` data frame.
- `nodelabelcolor = "black"`: color for labels in the `nodelabel` data frame.

Output the following named tuple, that can be used for downstream annotations
with RCall:

```
(xmin, xmax, ymin, ymax,
 node_x,     node_y,     node_y_lo, node_y_hi,
 edge_x_lo,  edge_x_hi,  edge_y_lo,  edge_y_hi,
 arrow_x_lo, arrow_x_hi, arrow_y_lo, arrow_y_hi,
 node_data, edge_data)
```

1. `:xmin`: minimum x value of the plot
2. `:xmax`: maximum x value of the plot
3. `:ymin`: minimum y value of the plot
4. `:ymax`: maximum y value of the plot
5. `:node_x`: x values of the nodes in net.node in their respective order
6. `:node_y`: y values of the nodes
7. `:node_y_lo`: y value of the beginning of the vertical bar representing the
    clade at each node
8. `:node_y_hi`: y value of the end of the vertical bar
9. `:edge_x_lo`: x value of the beginning of the edges in `net.edge` in their
    respective order. An edge (or branch) in the network is represented by 2
    segments if it is a minor edge under the `:fulltree` style.
    `:edge_*` give the coordinates of the first (horizontal) segment for these
    minor edges, and the coordinates of the unique (horizontal) segment for
    major edges (tree of major hybrid edges).
10. `:edge_x_hi`: x value of the end of the edges
11. `:edge_y_lo`: y value of the beginning of the edges
12. `:edge_y_hi`: y value of the end of the edges
13. `:arrow_x_lo`: x value for the beginning of the arrows, for minor hybrid
     edges, listed in the same order as in `filter(e -> !e.ismajor, net.edge)`.
14. `:arrow_x_hi`: same, but for the end of arrows
15. `:arrow_y_lo`: y values for the beginning of arrows
16. `:arrow_y_hi`: same, but for the end of arrows
17. `:node_data`: node data frame: see section [Adding labels](@ref) for more
18. `:edge_data`: edge data frame

Note that `plot` actually modifies some (minor) attributes of the network,
as it calls `PhyloNetworks.directedges!` and `PhyloNetworks.preorder!`
(unless with option `preorder = false`, which is not recommended).

If hybrid edges cross tree and major edges, you may choose to rotate some tree
edges to eliminate crossing edges, using
[`PhyloNetworks.rotate!`](https://juliaphylo.github.io/PhyloNetworks.jl/dev/lib/public/#PhyloNetworks.rotate!-Tuple%7BHybridNetwork,%20Integer%7D).

**Alternative**: a tree or network can be exported with [`sexp`](@ref)
and then displayed with R's "plot" and all its options.
"""
function plot(
    net::HybridNetwork;
    useedgelength::Bool=false,
    showtiplabel::Bool=true,
    shownodenumber::Bool=false,
    showedgelength::Bool=false,
    showgamma::Bool=false,
    edgecolor = "black",
    majorhybridedgecolor::AbstractString="deepskyblue4",
    minorhybridedgecolor::AbstractString="deepskyblue",
    defaultedgecolor = nothing,
    showedgenumber::Bool=false,
    shownodelabel::Bool=false,
    edgelabel::AbstractDataFrame=DataFrame(),
    nodelabel::AbstractDataFrame=DataFrame(),
    xlim = nothing,
    ylim = nothing,
    tipoffset = 0,
    tipcex = 1,
    nodecex = 1,
    edgecex = 1,
    style::Symbol=:fulltree,
    arrowlen::Real=(style==:majortree ? 0 : 0.1),
    minorlinetype = nothing,
    edgewidth = 1,
    edgenumbercolor = "grey", # don't limit the type because R accepts many types
    edgelabelcolor = "black", # and these colors are used as is
    nodelabelcolor = "black",
    edgelabeladj = [.5,0],
    nodelabeladj = 1,
    curved::Symbol = :none,
    bend::Real = 0.3,
    preorder::Bool=true,
)
    if getroot(net).leaf
        @warn """The network is rooted at a leaf: the plot won't look good.
            Try rooting the network on the edge adjacent to that leaf, with
            rootonedge!(network_name, $(getroot(net).edge[1].number))"""
    end
    (edge_xB, edge_xE, edge_yB, edge_yE, node_x, node_y, node_yB, node_yE,
     hybridedge_xB, hybridedge_xE, hybridedge_xC,
     hybridedge_yB, hybridedge_yE, hybridedge_yC,
     xmin, xmax, ymin, ymax) = edgenode_coordinates(
        net, useedgelength, style==:majortree, curved, bend, preorder)
    nedges = length(net.edge)
    nminor = length(hybridedge_xB)
    labelnodes, nodelabel = check_nodedataframe(net, nodelabel)
    ndf = prepare_nodedataframe(net, nodelabel, shownodenumber,
            shownodelabel, labelnodes, node_x, node_y)
    if showtiplabel || shownodenumber || shownodelabel || labelnodes
        expfac = 0.1  # force 10% more space to show tip/node/root name
        expfacy = 0.5 # additive expansion for y axis
        xmin -= (xmax-xmin)*expfac
        xmax += (xmax-xmin)*expfac
        ymin -= expfacy
        ymax += expfacy
    end
    xmax += tipoffset
    if !isnothing(xlim)
        length(xlim) == 2 ||
          error("xlim needs to contain 2 values: lower and upper limits. defaults: [$xmin,$xmax]")
        xmin=xlim[1]; xmax=xlim[2]
    end
    if !isnothing(ylim)
        length(ylim) == 2 ||
          error("ylim needs to contain 2 values: lower and upper limits. defaults: [$ymin,$ymax]")
        ymin=ylim[1]; ymax=ylim[2]
    end
    leaves = [n.leaf for n in net.node]
    if isa(edgecolor, AbstractDict) # then ignore {maj|min}orhybridedgecolor
      defaultedgecolor = (isnothing(defaultedgecolor) ? "black" : string(defaultedgecolor) )
      eCol = Vector{String}(undef,nedges)
      hybmincol_vec = Vector{String}()
      for (ie,ee) in enumerate(net.edge)
        ec = string(get(edgecolor, ee.number, defaultedgecolor))
        eCol[ie] = ec
        if !ee.ismajor
          push!(hybmincol_vec, ec)
        end
      end
    else
      if isnothing(defaultedgecolor)
        defaultedgecolor = edgecolor
      end
      eCol = fill(edgecolor, nedges)
      eCol[ [ e.hybrid  for e in net.edge] ] .= majorhybridedgecolor
      eCol[ [!e.ismajor for e in net.edge] ] .= minorhybridedgecolor
      hybmincol_vec = Iterators.repeated(minorhybridedgecolor, nminor)
    end

    if isa(edgewidth, Number)
      edgewidth_vec = Iterators.repeated(edgewidth, nedges)
      hybridedgewidth_vec = Iterators.repeated(edgewidth, nminor)
    elseif isa(edgewidth, AbstractDict)
      ewtype = valtype(edgewidth)
      ewtype <: Number || error("edgewidth should be numerical")
      edgewidth_vec = Vector{ewtype}(undef,nedges)
      hybridedgewidth_vec = Vector{ewtype}()
      for (ie,ee) in enumerate(net.edge)
        # fill in edgewidth vector, with default 1 for non-listed edges
        ew = (haskey(edgewidth, ee.number) ? edgewidth[ee.number] : one(ewtype))
        edgewidth_vec[ie] = ew
        if !ee.ismajor
          push!(hybridedgewidth_vec, ew)
        end
      end
    end
    # this makes the arrows dashed if :fulltree is used
    if isnothing(minorlinetype)
        minorlinetype = (style==:majortree ? "solid" : "longdash")
    end

    if !(style in [:fulltree, :majortree])
      @warn "Style $style is unknown. Defaulted to :fulltree."
      style = :fulltree
    end
    curved ∈ (:none, :minor, :both) ||
        error("curved must be :none, :minor, or :both; got " * repr(curved))

    R"plot"(node_x[leaves], node_y[leaves], type="n",
         xlim=[xmin,xmax], ylim=[ymin,ymax],
         axes=false, xlab="", ylab="")
    R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)
    for (ie,e,le) in zip(1:nedges, net.edge, edgewidth_vec)
        if !e.hybrid || !e.ismajor || curved != :both
            R"segments"(edge_xB[ie], edge_yB[ie], edge_xE[ie], edge_yE[ie],
                col=eCol[ie], lwd=le)
        else
            draw_quadraticbezier(edge_xB[ie], edge_yB[ie], edge_xE[ie], edge_yE[ie],
                missing, missing, arrowlen, eCol[ie], le, "solid")
                # fixit: replace missings by correct control point coordinates
        end
    end
    for (ie,ce,le) in zip(1:nminor, hybmincol_vec, hybridedgewidth_vec)
        if curved==:none
            R"arrows"(hybridedge_xB[ie], hybridedge_yB[ie], hybridedge_xE[ie], hybridedge_yE[ie],
                length=arrowlen, angle=20, col=ce, lty=minorlinetype, lwd=le)
        else
            draw_quadraticbezier(hybridedge_xB[ie], hybridedge_yB[ie],
                hybridedge_xE[ie], hybridedge_yE[ie], hybridedge_xC[ie], hybridedge_yC[ie],
                arrowlen, ce, le, minorlinetype)
        end
    end
#=
    curved_major_parent_y2 = Dict{Int, Vector{Float64}}()
    hybrid_offset = Float64(bend) * (ymax - ymin) / max(net.numtaxa, 1)
    minor_edge_indices = [i for i in 1:length(net.edge) if !net.edge[i].ismajor]
    overlapping_minors = Set{Int}()
    if curved != :none
        for j in 1:length(hybridedge_xB)
            min_seg = _minor_hybrid_segment(j, hybridedge_xB, hybridedge_yB,
                                            hybridedge_xE, hybridedge_yE)
            child_j  = PhyloNetworks.getchild(net.edge[minor_edge_indices[j]])
            maj_part = _major_partner_segment(child_j, net, edge_xB, edge_xE, edge_yE, node_y)
            if maj_part !== nothing && _segs_overlap(min_seg, maj_part)
                push!(overlapping_minors, j)
            end
        end
    end

    if curved == :none
        R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
        R"arrows"(hybridedge_xB, hybridedge_yB, hybridedge_xE, hybridedge_yE,
                  length=arrowlen, angle=20, col=hybmincol_vec, lty=minorlinetype,
                  lwd=hybridedgewidth_vec)
    else
        if curved == :both
            straight_idx = [i for i in 1:length(net.edge)
                            if !net.edge[i].hybrid || !net.edge[i].ismajor]
            ew_straight = isa(edgewidth_vec, Number) ? edgewidth_vec : edgewidth_vec[straight_idx]
            R"segments"(edge_xB[straight_idx], edge_yB[straight_idx],
                        edge_xE[straight_idx], edge_yE[straight_idx],
                        col=eCol[straight_idx], lwd=ew_straight)
            for i in findall(e -> e.hybrid && e.ismajor, net.edge)
                seg = _major_hybrid_segment(i, net, edge_xB, edge_xE, edge_yE, node_y)
                _, _, _, y2_i = seg
                parent_ni_i = findfirst(x -> x === PhyloNetworks.getparent(net.edge[i]), net.node)
                child_node_i = PhyloNetworks.getchild(net.edge[i])
                minor_j = let jcount = 0, found = nothing
                    for ee in net.edge
                        if !ee.ismajor
                            jcount += 1
                            if PhyloNetworks.getchild(ee) === child_node_i
                                found = jcount
                                break
                            end
                        end
                    end
                    found
                end
                lwd_i = isa(edgewidth_vec, AbstractVector) ? edgewidth_vec[i] : edgewidth_vec
                if minor_j !== nothing && minor_j in overlapping_minors
                    # major overlaps minor: draw major straight so the minor arc is visible
                    x0s, y0s, x2s, y2s = seg
                    eff_al = Float64(arrowlen) > 0 ? Float64(arrowlen) : 0.1
                    R"arrows"(x0s, y0s, x2s, y2s, length=eff_al, angle=20,
                              col=eCol[i], lty="solid", lwd=lwd_i)
                else
                    # no overlap: curved major, update node-bar gap tracking
                    if !haskey(curved_major_parent_y2, parent_ni_i)
                        curved_major_parent_y2[parent_ni_i] = Float64[]
                    end
                    push!(curved_major_parent_y2[parent_ni_i], y2_i)
                    bow = if minor_j !== nothing
                        partner_seg = _minor_hybrid_segment(minor_j, hybridedge_xB, hybridedge_yB,
                                                            hybridedge_xE, hybridedge_yE)
                        _hybrid_bow_sign(seg; offset=hybrid_offset, partner=partner_seg)
                    else
                        _hybrid_bow_sign(seg; offset=hybrid_offset, partner=nothing)
                    end
                    draw_quadraticbezier(seg; bow_sign=bow, offset=hybrid_offset, bend=bend,
                                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                         col=eCol[i], lwd=lwd_i, arrowlen=arrowlen, linetype="solid",
                                         force_bow=false)
                end
            end
        else
            R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
        end
        for j in 1:length(hybridedge_xB)
            seg = _minor_hybrid_segment(j, hybridedge_xB, hybridedge_yB,
                                        hybridedge_xE, hybridedge_yE)
            is_overlap = j in overlapping_minors
            bow = if is_overlap
                -_bow_left_sign(seg..., hybrid_offset)
            else
                partner = curved == :both ?
                    _major_partner_segment(PhyloNetworks.getchild(net.edge[minor_edge_indices[j]]),
                                           net, edge_xB, edge_xE, edge_yE, node_y) : nothing
                _hybrid_bow_sign(seg; offset=hybrid_offset, partner=partner)
            end
            col_j = isa(hybmincol_vec, AbstractVector) ? hybmincol_vec[j] : hybmincol_vec
            lwd_j = isa(hybridedgewidth_vec, AbstractVector) ? hybridedgewidth_vec[j] :
                    hybridedgewidth_vec
            draw_quadraticbezier(seg; bow_sign=bow, offset=hybrid_offset, bend=bend,
                                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                                 col=col_j, lwd=lwd_j, arrowlen=arrowlen, linetype=minorlinetype,
                                 force_bow=is_overlap)
        end
    end
    if isempty(curved_major_parent_y2)
        R"segments"(node_x, node_yB, node_x, node_yE, col=defaultedgecolor)
    else
        normal_node_idx = [i for i in 1:length(net.node) if !haskey(curved_major_parent_y2, i)]
        R"segments"(node_x[normal_node_idx], node_yB[normal_node_idx],
                    node_x[normal_node_idx], node_yE[normal_node_idx],
                    col=defaultedgecolor)
        for (pni, y2_list) in curved_major_parent_y2
            y_parent = node_y[pni]
            exclusions = Tuple{Float64, Float64}[]
            for y2_i in y2_list
                push!(exclusions, (min(y_parent, y2_i), max(y_parent, y2_i)))
            end
            sort!(exclusions)
            merged = Tuple{Float64, Float64}[]
            for (lo, hi) in exclusions
                if isempty(merged) || lo > merged[end][2]
                    push!(merged, (lo, hi))
                else
                    merged[end] = (merged[end][1], max(merged[end][2], hi))
                end
            end
            tree_ys = sort!([node_y[findfirst(x -> x === PhyloNetworks.getchild(e), net.node)]
                             for e in net.edge
                             if PhyloNetworks.getparent(e) === net.node[pni] &&
                                !(e.hybrid && e.ismajor)])
            pos = node_yB[pni]
            for (lo, hi) in merged
                if pos < lo
                    R"segments"(node_x[pni], pos, node_x[pni], lo, col=defaultedgecolor)
                end
                draw_from = max(pos, lo)
                for ty in tree_ys
                    if draw_from < ty <= hi
                        R"segments"(node_x[pni], draw_from, node_x[pni], ty, col=defaultedgecolor)
                        draw_from = ty
                    end
                end
                pos = max(pos, hi)
            end
            if pos < node_yE[pni]
                R"segments"(node_x[pni], pos, node_x[pni], node_yE[pni],
                            col=defaultedgecolor)
            end
        end
    end
    =#
    if showtiplabel
      R"text"(node_x[leaves] .+ tipoffset, node_y[leaves],
              tiplabels(net), adj=0, font=3, cex=tipcex)
    end
    if shownodelabel
      R"text"(ndf[.!ndf[!,:lea],:x], ndf[.!ndf[!,:lea],:y],
              ndf[.!ndf[!,:lea],:name], font=3, cex=tipcex, adj=[.5,0])
    end
    if shownodenumber
      R"text"(ndf[!,:x], ndf[!,:y], ndf[!,:num], adj=1)
    end
    if labelnodes
      R"text"(ndf[!,:x], ndf[!,:y], ndf[!,:lab], adj=nodelabeladj,
              col=nodelabelcolor, cex=nodecex)
    end
    labeledges, edf = prepare_edgedataframe(net, edgelabel, style,
                        edge_xB, edge_xE, edge_yB, edge_yE,
                        hybridedge_xB, hybridedge_xE, hybridedge_yB, hybridedge_yE)
    if labeledges
      R"text"(edf[!,:x], edf[!,:y], edf[!,:lab], adj=edgelabeladj,
              col=edgelabelcolor, cex=edgecex)
    end
    if showedgelength
      R"text"(edf[!,:x], edf[!,:y], edf[!,:len], adj=[.5,1.])
    end
    if showgamma && net.numhybrids>0
      im = edf[!,:hyb] .& edf[!,:min]
      iM = edf[!,:hyb] .& .!edf[!,:min]
      R"text"(edf[im,:x], edf[im,:y], edf[im,:gam],
              adj=[.5,1], col=minorhybridedgecolor)
      R"text"(edf[iM,:x], edf[iM,:y], edf[iM,:gam],
              adj=[.5,1], col=majorhybridedgecolor)
    end
    if showedgenumber
      R"text"(edf[!,:x], edf[!,:y], edf[!,:num], adj=[.5,0], col=edgenumbercolor)
    end
    return (xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
      node_x=node_x, node_y=node_y,
      node_y_lo=node_yB, node_y_hi=node_yE,
      edge_x_lo=edge_xB, edge_x_hi=edge_xE,
      edge_y_lo=edge_yB, edge_y_hi=edge_yE,
      arrow_x_lo=hybridedge_xB, arrow_x_hi=hybridedge_xE,
      arrow_y_lo=hybridedge_yB, arrow_y_hi=hybridedge_yE,
      node_data=ndf, edge_data=edf)
end

function _major_hybrid_segment(i::Int, net::HybridNetwork,
                               edge_xB, edge_xE, edge_yE, node_y)
    parent_ni = findfirst(x -> x === PhyloNetworks.getparent(net.edge[i]), net.node)
    return (edge_xB[i], node_y[parent_ni], edge_xE[i], edge_yE[i])
end

_minor_hybrid_segment(j::Int, hxB, hyB, hxE, hyE) = (hxB[j], hyB[j], hxE[j], hyE[j])

function _minor_partner_segment(x2::Float64, y2::Float64, hxB, hyB, hxE, hyE)
    j = findfirst(k -> abs(hxE[k] - x2) < 1e-10 && abs(hyE[k] - y2) < 1e-10, 1:length(hxB))
    j === nothing ? nothing : _minor_hybrid_segment(j, hxB, hyB, hxE, hyE)
end

function _major_partner_segment(child_node, net::HybridNetwork,
                                edge_xB, edge_xE, edge_yE, node_y)
    i = findfirst(e -> e.hybrid && e.ismajor &&
                   PhyloNetworks.getchild(e) === child_node, net.edge)
    i === nothing ? nothing : _major_hybrid_segment(i, net, edge_xB, edge_xE, edge_yE, node_y)
end

function draw_quadraticbezier(
    x0,y0, x2,y2, x1,y1,
    arrowlen,
    ce, # color for the edge
    le, # line width for the edge
    linetype,
    nsegments=20,
)
    if ismissing(x1)
        R"arrows"(x0, y0, x2, y2, length=arrowlen, angle=20,
                  col=ce, lwd=le, lty=linetype)
    else
        ts = range(0, 1, step=1/nsegments)
        xs = [(1-t)^2 * x0 + 2*(1-t)*t * x1 + t^2 * x2 for t in ts]
        ys = [(1-t)^2 * y0 + 2*(1-t)*t * y1 + t^2 * y2 for t in ts]
        for i in 1:(nsegments-1)
            R"segments"(xs[i], ys[i], xs[i+1], ys[i+1], col=ce, lwd=le, lty=linetype)
        end
        R"arrows"(xs[nsegments], ys[nsegments], x2, y2,
            length=arrowlen, angle=20, col=ce, lwd=le, lty=linetype)
    end
end
