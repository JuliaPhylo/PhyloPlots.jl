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
- `edgewidth=1`: width of horizontal (not diagonal) edges. To vary them,
  use a dictionary to map the number of each edge to its desired width.
- `xlim`, `ylim`: array of 2 values, to determine the axes limits.

## tip annotations:

- `showtiplabel = true`: if true, taxon labels (names) are shown.
- `tipoffset = 0`: to offset tip labels.
- `tipcex = 1`: character expansion for tip and internal node names.

## nodes & edges annotations:

- `shownodelabel = false`: if true, internal nodes are labelled with their names.
  Useful for hybrid nodes, which do have tags like 'H1'.
- `shownodenumber = false`: if true, nodes are labelled with the number used internally.
- `showedgenumber = false`: if true, edges are labelled with the number used internally.
- `showedgelength = false`: if true, edges are labelled with their length (above).
- `showgamma = false`: if true, hybrid edges are labelled with their heritability (below).
- `edgelabel = DataFrame()`: dataframe with two columns: the first with edge numbers,
  the second with labels (like bootstrap values) to annotate edges. empty by default.
- `nodelabel = DataFrame()`: dataframe with two columns: the first with node numbers,
  the second with labels (like bootstrap values for hybrid relationships)
  to annotate nodes. empty by default.
- `nodecex = 1`: character expansion for labels in the `nodelabel` data frame.
- `edgecex = 1`: character expansion for labels in the `edgelabel` data frame.

## colors:

- `edgecolor = "black"`: color for tree edges.
- `majorhybridedgecolor = "deepskyblue4"`: color for major hybrid edges.
- `minorhybridedgecolor = "deepskyblue"`: color for minor hybrid edges.
- `edgenumbercolor = "grey"`: color for edge numbers.
- `edgelabelcolor = "black"`: color for labels in the `edgelabel` data frame.
- `nodelabelcolor = "black"`: color for labels in the `nodelabel` data frame.

Output the following named tuple, that can be used for downstream plot annotations
with RCall:

```
(xmin, xmax, ymin, ymax,
 node_x,    node_y,    node_y_lo, node_y_hi,
 edge_x_lo, edge_x_hi, edge_y_lo, edge_y_hi,
 node_data, edge_data)
```

1. `:xmin`: minimum x value of the plot
2. `:xmax`: maximum x value of the plot
3. `:ymin`: minimum y value of the plot
4. `:ymax`: maximum y value of the plot
5. `:node_x`: x values of the nodes in net.node in their respective order
6. `:node_y`: y values of the nodes
7. `:node_y_lo`: y value of the beginning of the vertical bar representing the clade at each node
8. `:node_y_hi`: y value of the end of the vertical bar
9. `:edge_x_lo`: x value of the beginning of the edges in `net.edge` in their respective order
10. `:edge_x_hi`: x value of the end of the edges
11. `:edge_y_lo`: y value of the beginning of the edges
12. `:edge_y_hi`: y value of the end of the edges
13. `:node_data`: node data frame: see section [Adding labels](@ref) for more
14. `:edge_data`: edge data frame

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
    edgecolor::AbstractString="black",
    majorhybridedgecolor::AbstractString="deepskyblue4",
    minorhybridedgecolor::AbstractString="deepskyblue",
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
    edgewidth = 1,
    edgenumbercolor = "grey", # don't limit the type because R accepts many types
    edgelabelcolor = "black", # and these colors are used as is
    nodelabelcolor = "black",
    preorder::Bool=true,
)
    if getroot(net).leaf
        @warn """The network is rooted at a leaf: the plot won't look good.
            Try rooting the network on the edge adjacent to that leaf, with
            rootonedge!(network_name, $(getroot(net).edge[1].number))"""
    end
    (edge_xB, edge_xE, edge_yB, edge_yE, node_x, node_y, node_yB, node_yE,
     hybridedge_xB, hybridedge_xE, hybridedge_yB, hybridedge_yE,
     xmin, xmax, ymin, ymax) = edgenode_coordinates(
        net, useedgelength, style==:majortree, preorder)
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
    eCol = fill(edgecolor, length(net.edge))
    eCol[ [ e.hybrid  for e in net.edge] ] .= majorhybridedgecolor
    eCol[ [!e.ismajor for e in net.edge] ] .= minorhybridedgecolor

    if isa(edgewidth, Number)
      edgewidth_vec = edgewidth
      hybridedgewidth_vec = edgewidth
    elseif isa(edgewidth, AbstractDict)
      ewtype = valtype(edgewidth)
      ewtype <: Number || error("edgewidth should be numerical")
      edgewidth_vec = Vector{ewtype}(undef,length(edge_xB))
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
    arrowstyle = style==:majortree ? "solid" : "longdash"

    if !(style in [:fulltree, :majortree])
      @warn "Style $style is unknown. Defaulted to :fulltree."
      style = :fulltree
    end

    R"""
    plot($(node_x[leaves]), $(node_y[leaves]), type='n',
         xlim=c($xmin,$xmax), ylim=c($ymin,$ymax),
         axes=FALSE, xlab='', ylab='')
    """
    R"segments"(edge_xB, edge_yB, edge_xE, edge_yE, col=eCol, lwd=edgewidth_vec)
    R"arrows"(hybridedge_xB, hybridedge_yB, hybridedge_xE, hybridedge_yE,
              length=arrowlen, angle=20, col=minorhybridedgecolor, lty=arrowstyle,
              lwd=hybridedgewidth_vec)
    R"segments"(node_x, node_yB, node_x, node_yE, col=edgecolor)
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
      R"text"(ndf[!,:x], ndf[!,:y], ndf[!,:lab], adj=1,
              col=nodelabelcolor, cex=nodecex)
    end
    labeledges, edf = prepare_edgedataframe(net, edgelabel, style,
                        edge_xB, edge_xE, edge_yB, edge_yE,
                        hybridedge_xB, hybridedge_xE, hybridedge_yB, hybridedge_yE)
    if labeledges
      R"text"(edf[!,:x], edf[!,:y], edf[!,:lab], adj=[.5,0],
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
      node_data=ndf, edge_data=edf)
end
