"""
    plot(net::HybridNetwork, method::Symbol)

Plot a network using R graphics.
`method` should be `:R` (actually, any symbol would do, for now!).

optional arguments, shared with the Gadfly-based plot function:
- `useEdgeLength = false` : if true, the tree edges and major hybrid edges are
  drawn proportionally to their length. Minor hybrid edges are not, however.
  Note that edge lengths in coalescent units may scale very poorly with time.
- `showTipLabel = true` : if true, taxon labels are shown. You may need to zoom out to see them.
- `showNodeNumber = false` : if true, nodes are labelled with the number used internally.
- `showEdgeLength = false` : if true, edges are labelled with their length (above)
- `showGamma = false` : if true, hybrid edges are labelled with their heritability (below)
- `edgeColor = "black"` : color for tree edges.
- `majorHybridEdgeColor = "deepskyblue4"` : color for major hybrid edges
- `minorHybridEdgeColor = "deepskyblue"` : color for minor hybrid edges
- `showEdgeNumber = false` : if true, edges are labelled with the number used internally.
- `showIntNodeLabel = false` : if true, internal nodes are labelled with their names.
  Useful for hybrid nodes, which do have tags like 'H1'.
- `edgeLabel = DataFrame()` : dataframe with two columns: the first with edge numbers, the second with labels
  (like bootstrap values) to annotate edges. empty by default.
- `nodeLabel = DataFrame()` : dataframe with two columns: the first with node numbers, the second with labels
  (like bootstrap values for hybrid relationships) to annotate nodes. empty by default.
- `style = :fulltree` : symbol indicating the style of the diagram
  * `:majortree` will simply draw minor edges onto the major tree.
  * `:fulltree` will draw minor edges as their own branches in the tree (like in icytree.org), 
    usefull for overlapping or confusing networks.
- `arrowlen` : the length of the arrow tips in the full tree style. if `style = :fulltree`, then 
  `arrowlen = 0.2`. otherwise, `arrowlen = 0`, which makes the arrows appear as segments.

optional arguments specific to this function:
- `xlim`, `ylim` : array of 2 values
- `tipOffset = 0.0` : to offset tip labels

plot() returns the following tuple: 
`(xmin, xmax, ymin, ymax, node_x, node_y, node_yB, node_yE, 
edge_xB, edge_xE, edge_yB, edge_yE, ndf, edf)`

1. `xmin` : the minimum x value of the plot
2. `xmax` : the maximum x value of the plot
3. `ymin` : the minimum y value of the plot
4. `ymax` : the maximum y value of the plot
5. `node_x` : the x values of the nodes in HybridNetWork.node in their respective order
6. `node_y` : the y values of the nodes in HybridNetWork.node in their respective order
7. `node_yB` : the y value of the beginning of the verticle bar coresponding to each node in HybridNetwork.node
8. `node_yE` : the y value of the end of the verticle bar coresponding to each node in HybridNetwork.node
9. `node_yE` : the y value of the end of the verticle bar coresponding to each node in HybridNetwork.node
10. `edge_xB` : the x value of the beginning of the line coresponding to each edge in HybridNetwork.edge
11. `edge_xE` : the x value of the end of the line coresponding to each edge in HybridNetwork.edge
12. `edge_yB` : the y value of the beginning of the line coresponding to each edge in HybridNetwork.edge
13. `edge_yE` : the y value of the end of the line coresponding to each edge in HybridNetwork.edge
14. `ndf` : the node data frame: see here for more
15. `edf` : the edge data frame: see here for more

Note that `plot` actually modifies some (minor) attributes of the network,
as it calls `directEdges!` and `preorder!`.

If hybrid edges cross tree and major edges, you may choose to rotate some tree
edges to eliminate crossing edges, using `rotate!`
(in [`PhyloNetworks`](http://crsl4.github.io/PhyloNetworks.jl/latest/lib/public/#PhyloNetworks.rotate!)).

**Alternative**: a tree or network can be exported with [`sexp`](@ref)
and then displayed with R's "plot" and all its options.
"""
function plot(net::HybridNetwork, ::Symbol; useEdgeLength=false::Bool,
    mainTree=false::Bool, showTipLabel=true::Bool, showNodeNumber=false::Bool,
    showEdgeLength=false::Bool, showGamma=false::Bool,
    edgeColor="black"::String,
    majorHybridEdgeColor="deepskyblue4"::String,
    minorHybridEdgeColor="deepskyblue"::String,
    showEdgeNumber=false::Bool, showIntNodeLabel=false::Bool,
    edgeLabel=DataFrame()::DataFrame, nodeLabel=DataFrame()::DataFrame,
    xlim=Float64[]::Array{Float64,1}, ylim=Float64[]::Array{Float64,1},
    tipOffset=0.0::Float64, tipcex=1.0::Float64,
    style=:fulltree::Symbol, arrowlen=(style==:majortree ? 0 : 0.1)::Real)

    (edge_xB, edge_xE, edge_yB, edge_yE, node_x, node_y, node_yB, node_yE,
     hybridedge_xB, hybridedge_xE, hybridedge_yB, hybridedge_yE,
     xmin, xmax, ymin, ymax) = getEdgeNodeCoordinates(net, useEdgeLength, style==:majortree)
    labelnodes, nodeLabel = checkNodeDataFrame(net, nodeLabel)
    ndf = prepareNodeDataFrame(net, nodeLabel, showNodeNumber,
            showIntNodeLabel, labelnodes, node_x, node_y)
    if (showTipLabel || showNodeNumber || showIntNodeLabel || labelnodes)
        expfac = 0.1  # force 10% more space to show tip/node/root name
        expfacy = 0.5 # additive expansion for y axis
        xmin -= (xmax-xmin)*expfac
        xmax += (xmax-xmin)*expfac
        ymin -= expfacy
        ymax += expfacy
    end
    xmax += tipOffset
    if length(xlim)==2
        xmin=xlim[1]; xmax=xlim[2]
    end
    if length(ylim)==2
        ymin=ylim[1]; ymax=ylim[2]
    end
    leaves = [n.leaf for n in net.node]
    eCol = fill(edgeColor, length(net.edge))
    eCol[ [ e.hybrid  for e in net.edge] ] .= majorHybridEdgeColor
    eCol[ [!e.isMajor for e in net.edge] ] .= minorHybridEdgeColor

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
    segments($edge_xB, $edge_yB, $edge_xE, $edge_yE, col=$eCol)
    arrows($hybridedge_xB, $hybridedge_yB, $hybridedge_xE, $hybridedge_yE, length=$arrowlen, angle = 20, col=$minorHybridEdgeColor, lty=$arrowstyle)
    segments($node_x, $node_yB, $node_x, $node_yE, col=$edgeColor,)
    """
    if showTipLabel
      R"text"(node_x[leaves] .+ tipOffset, node_y[leaves],
              tipLabels(net), adj=0, font=3, cex=tipcex)
    end
    if showIntNodeLabel
      R"text"(ndf[.!ndf[!,:lea],:x], ndf[.!ndf[!,:lea],:y],
              ndf[.!ndf[!,:lea],:name], font=3, cex=tipcex, adj=[.5,0])
    end
    if showNodeNumber
      R"text"(ndf[!,:x], ndf[!,:y], ndf[!,:num], adj=1)
    end
    if labelnodes
      R"text"(ndf[!,:x], ndf[!,:y], ndf[!,:lab], adj=1)
    end
    labeledges, edf = prepareEdgeDataFrame(net, edgeLabel, mainTree,
                        edge_xB, edge_xE, edge_yB, edge_yE,
                        hybridedge_xB, hybridedge_xE, hybridedge_yB, hybridedge_yE)
    if labeledges
      R"text"(edf[!,:x], edf[!,:y], edf[!,:lab], adj=[.5,0])
    end
    if showEdgeLength
      R"text"(edf[!,:x], edf[!,:y], edf[!,:len], adj=[.5,1.])
    end
    if (showGamma && net.numHybrids>0)
      im = edf[!,:hyb] .& edf[!,:min]
      iM = edf[!,:hyb] .& .!edf[!,:min]
      R"text"(edf[im,:x], edf[im,:y], edf[im,:gam],
              adj=[.5,1], col=minorHybridEdgeColor)
      R"text"(edf[iM,:x], edf[iM,:y], edf[iM,:gam],
              adj=[.5,1], col=majorHybridEdgeColor)
    end
    if showEdgeNumber
      R"text"(edf[!,:x], edf[!,:y], edf[!,:num], adj=[.5,0])
    end
    return (xmin, xmax, ymin, ymax, node_x, node_y, node_yB, node_yE,
      edge_xB, edge_xE, edge_yB, edge_yE, ndf, edf)
end
