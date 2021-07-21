```@setup better_edges
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Better edges

## Different hybrid edge styles

We can use the `style` option to visualize minor hybrid edges as simple lines,
unlike the [icytree](https://icytree.org/) style visualization. `style` is by default `:fulltree`,
but by switching it to `:majortree`, we can draw minor hybrid edges as diagonal lines.

```@example better_edges
R"svg"(figname("style_example.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, style=:majortree);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/style_example.svg)

## Using edge lengths

We can use `useEdgeLength=true` to draw a plot that uses the network's edge lengths to determine the lengths of the
lines. For this, we'll use a network that has branch lengths:

```@example better_edges
R"svg"(figname("edge_len_example.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2]) # hide
net = readTopology("(A:3.3,((B:1.5,#H1:0.5):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);")
df = DataFrame(Number=[-3, 3], Label=["N", "H1"]); # hide
plot(net, :R, useEdgeLength=true, ylim = [-1, 5.5], nodeLabel = df); # hide
R"text"([3], [0], ["useEdgeLength=true"]) # hide
plot(net, :R, useEdgeLength=false, ylim = [-1, 5.5], nodeLabel = df); # hide
R"text"([3], [0], ["useEdgeLength=false"]) # hide
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_len_example.svg)

!!! note
    I used a DataFrame to add labels to the plot. For more on this,
    see the [Adding labels](@ref) section.

If branch lengths represent time, D could represent a fossil, or a virus strain sequenced
a year before the others. Seeing this visually is the advantage of `useEdgeLengths=true`

This network happens to be time consistent, because the distance
along the time (x) axis from node `N` to the hybrid node `H1` is
the same both ways.

!!! note "Time consistency"
    A network is time-consistent if all the paths between 2 given nodes all
    have the same length.
    Time inconsistency can occur when branch lengths are not measured in
    calendar time, such as if branch lengths are in substitutions per site
    (some paths might evolve with more substitutions than others), or in
    number of generations (some lineages might have 1 generation per year,
    others more or fewer generations per year), or in coalescent units
    (number of generations / effective population size).

    A time-consistent network may be ultrametric (the distance
    between the root and the tips is the same across all tips),
    or not like the network above.

Time inconsistent networks like these ones below might cause confusion:

```@example better_edges
R"svg"(figname("edge_len_example2.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2]) # hide
net1 = readTopology("(A:3.3,((B:1.5,#H1:1.2):1.5,((C:1.8)#H1:1,D:1.1):.2):0.3);");
net2 = readTopology("(A:3.3,((B:1.5,#H1:0.2):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);");
plot(net1, :R, useEdgeLength=true); # hide
plot(net2, :R, useEdgeLength=true); # hide
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/edge_len_example2.svg)

It may be useful to consider using `style=:majortree` if it causes
too much confusion, since the `:majortree` style doesn't visually represent
minor edge lengths. Because of this, I used the `showEdgeLength=true` option to
see the information anyway.

```@example better_edges
R"svg"(figname("edge_len_example3.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2])
plot(net1, :R, useEdgeLength=true, style = :majortree, showEdgeLength=true, arrowlen=0.1);
plot(net2, :R, useEdgeLength=true, style = :majortree, showEdgeLength=true, arrowlen=0.1);
R"dev.off()" # hide
nothing # hide
```
![example4](../assets/figures/edge_len_example3.svg)

I also used the `arrowlen=0.1` option to show the arrow tips to show the direction of minor edges,
which are hidden by default when using the `style=:majortree` option.
