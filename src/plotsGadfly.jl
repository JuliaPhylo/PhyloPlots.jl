"""
    plot(net::HybridNetwork; useEdgeLength=false, mainTree=false, showTipLabel=true,
         showNodeNumber=false, showEdgeLength=false, showGamma=false, edgeColor=colorant"black",
         majorHybridEdgeColor=colorant"deepskyblue4", minorHybridEdgeColor=colorant"deepskyblue",
         showEdgeNumber=false, showIntNodeLabel=false, edgeLabel=[], nodeLabel=[])

Plots a network, from left to right.

- useEdgeLength: if true, the tree edges and major hybrid edges are
  drawn proportionally to their length. Minor hybrid edges are not, however.
  Note that edge lengths in coalescent units may scale very poorly with time.
- mainTree: if true, the minor hybrid edges are ommitted.
- showTipLabel: if true, taxon labels are shown. You may need to zoom out to see them.
- showNodeNumber: if true, nodes are labelled with the number used internally.
- showEdgeLength: if true, edges are labelled with their length (above)
- showGamma: if true, hybrid edges are labelled with their heritability (below)
- edgeColor: color for tree edges. black by default.
- majorHybridEdgeColor: color for major hybrid edges
- minorHybridEdgeColor: color for minor hybrid edges
- showEdgeNumber: if true, edges are labelled with the number used internally.
- showIntNodeLabel: if true, internal nodes are labelled with their names.
  Useful for hybrid nodes, which do have tags like '#H1'.
- edgeLabel: dataframe with two columns: the first with edge numbers, the second with labels
  (like bootstrap values) to annotate edges. empty by default.
- nodeLabel: dataframe with two columns: the first with node numbers, the second with labels
  (like bootstrap values for hybrid relationships) to annotate nodes. empty by default.

Note that `plot` actually modifies some (minor) attributes of the network,
as it calls `directEdges!`, `preorder!` and `cladewiseorder!`.

If hybrid edges cross tree and major edges, you may choose to rotate some tree
edges to eliminate crossing edges, using `rotate!`.
"""
function Gadfly.plot(net::HybridNetwork; useEdgeLength=false::Bool,
        mainTree=false::Bool, showTipLabel=true::Bool, showNodeNumber=false::Bool,
        showEdgeLength=false::Bool, showGamma=false::Bool,
        edgeColor=colorant"black"::ColorTypes.Colorant,
        majorHybridEdgeColor=colorant"deepskyblue4"::ColorTypes.Colorant,
        minorHybridEdgeColor=colorant"deepskyblue"::ColorTypes.Colorant,
        showEdgeNumber=false::Bool, showIntNodeLabel=false::Bool,
        edgeLabel=DataFrame()::DataFrame, nodeLabel=DataFrame()::DataFrame)

    (edge_xB, edge_xE, edge_yB, edge_yE, node_x, node_y, node_yB, node_yE,
     xmin, xmax, ymin, ymax) = getEdgeNodeCoordinates(net, useEdgeLength)

    !net.node[net.root].leaf ||
        warn("the root is leaf $(net.node[net.root].name): the plot will look weird...")


    mylayers = Layer[] # gadfly layers
    # one layers for each edge
    for i=1:net.numEdges
        if (!mainTree || net.edge[i].isMajor)
            col = edgeColor
            if net.edge[i].hybrid
              if (net.edge[i].isMajor) col = majorHybridEdgeColor;
              else col = minorHybridEdgeColor; end
            end
            push!(mylayers,
              layer(x = [edge_xB[i],edge_xE[i]],
                    y = [edge_yB[i],edge_yE[i]], Geom.line,
                    Theme(default_color=col))[1])
        end
    end
    # one layer for each (vertical) clade
    for i=1:net.numNodes
        if (net.node[i].leaf) continue; end
        push!(mylayers,
              layer(y = [node_yB[i],node_yE[i]],
                    x = [node_x[i], node_x[i]], Geom.line)[1])
    end
    labelnodes, nodeLabel = checkNodeDataFrame(net, nodeLabel)
    if (showTipLabel || showNodeNumber || showIntNodeLabel || labelnodes)
      # white dot beyond tip labels/name and root node label to force enough zoom out
      expfac = 0.1
      expfacy = 0.5 # additive expansion for y axis
      push!(mylayers, layer(x=[xmin-(xmax-xmin)*expfac,xmax+(xmax-xmin)*expfac],
                            y=[ymin- expfacy, ymax+ expfacy],
               Geom.point, Theme(default_color=colorant"white"))[1])
      # data frame to place tip names and node annotations (labels)
      ndf = prepareNodeDataFrame(net, nodeLabel, showNodeNumber,
              showIntNodeLabel, labelnodes, node_x, node_y)
      if (showTipLabel)
        push!(mylayers, layer(ndf[ndf[:lea], [:x,:y,:name]], y="y", x="x", label="name",
            Geom.label(position=:right ;hide_overlaps=true))[1])
      end
      if (showIntNodeLabel)
        push!(mylayers, layer(ndf[.!ndf[:lea], [:x,:y,:name]], y="y", x="x", label="name",
            Geom.label(position=:above ;hide_overlaps=true))[1])
      end
      if (showNodeNumber)
        push!(mylayers, layer(ndf, y="y", x="x", label="num",
            Geom.label(position=:dynamic ;hide_overlaps=true))[1])
      end
      if labelnodes
        push!(mylayers, layer(ndf[:,[:x,:y,:lab]], y="y", x="x", label="lab",
            Geom.label(position=:left ;hide_overlaps=false))[1])
      end
    end
    # data frame for edge annotations.
    labeledges, edf = prepareEdgeDataFrame(net, edgeLabel, mainTree,
                        edge_xB, edge_xE, edge_yB, edge_yE)
        if labeledges
            push!(mylayers, layer(edf[:,[:x,:y,:lab]], y="y", x="x", label="lab",
                  Geom.label(position=:above ;hide_overlaps=false))[1])
        end
        if (showEdgeLength)
            push!(mylayers, layer(edf[:,[:x,:y,:len]], y="y", x="x", label="len",
                  Geom.label(position=:below ;hide_overlaps=false))[1])
        end
        if (showGamma && net.numHybrids>0)
            if !mainTree
            push!(mylayers, layer(edf[edf[:hyb] .& edf[:min], [:x,:y,:gam]], y="y", x="x",label="gam",
                  Geom.label(position=:below ;hide_overlaps=true),
                  Theme(point_label_color=minorHybridEdgeColor))[1])
            end
            push!(mylayers, layer(edf[edf[:hyb] .& .!edf[:min],[:x,:y,:gam]], y="y", x="x",label="gam",
                  Geom.label(position=:below ;hide_overlaps=true),
                  Theme(point_label_color=majorHybridEdgeColor))[1])
        end
        if (showEdgeNumber)
            push!(mylayers, layer(edf[:,[:x,:y,:num]], y="y", x="x", label="num",
                  Geom.label(position=:dynamic ;hide_overlaps=false))[1])
        end

    plot(mylayers, Guide.ylabel(nothing), Guide.xlabel(nothing), #("time"),
         Guide.xticks(ticks=:auto, label=false), # ticks=[xmin,xmin,xmax,xmax*1.1],
         Guide.yticks(ticks=:auto, label=false), # ticks=[ymin,ymax],
         Theme(default_color=edgeColor,grid_color=colorant"white",grid_line_width=0pt))
end
