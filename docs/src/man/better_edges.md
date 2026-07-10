```@setup better_edges
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Better edges

## Different hybrid edge styles

We can use the `style` and `curved` options to visualize hybrid edges
in various ways.
We will see examples later when hybrid edges are drawn straight with
`curved=:none`, or only minor hybrid edges are curved with `curved=:minor`.
By default `curved=:both` so both minor and major edges are curved.

- The default `style = :majortree` has the advantage of drawing each
  minor edge as a single segment, but the disadvantage of being unable
  to draw it proportional to the edge length (because the segment
  connects its 2 nodes, whose placements are dictated by the other edges).
- The `:fulltree` style draws each minor edge as 2 segments: one straight
  whose length can represent the edge length, and another segment
  (diagonal straight or curved) connecting to the hybrid child.


```@example better_edges
using RCall # to add annotations to the R-based plot
net = readnewick("(A,((B,#H1),(C,(D)#H1)));") # hide
R"svg"(figname("style_example.svg"), width=6, height=3) # hide
R"layout"([1 2]) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
plot(net, style=:majortree); # default
R"mtext"("style = :majortree", side=1, line=-2);
plot(net, style=:fulltree);
R"mtext"("style = :fulltree", side=1, line=-2);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/style_example.svg)

## Using edge lengths

We can use `useedgelength=true` to draw a plot that uses the
network's edge lengths to determine the lengths of the lines.
For this, we'll use a network that has branch lengths:

```@example better_edges
net = readnewick("(A:3.3,((B:1.5,#H1:0.5):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);")
df = DataFrame(number=[-3], label=["N"]); # hide
R"svg"(figname("edge_len_example.svg"), width=6, height=6) # hide
R"layout"([1 3; 2 4]) # hide
R"par"(mar=[0,0,0,0], oma=[0,0,.3,0]) # hide
plot(net, useedgelength=false, nodelabel=df, nodelabeladj=[1.2,-.2]); # hide
R"mtext"("useedgelength = false (default)", side=3, line=-1.5); # hide
R"mtext"("style = :majortree (default)", side=2, line=-1.5, las=0); # hide
plot(net, useedgelength=false, style=:fulltree, nodelabel=df, nodelabeladj=[1.2,-.2]); # hide
R"mtext"("style = :fulltree", side=2, line=-1.5, las=0); # hide
plot(net, useedgelength=true, curved=:none, showedgelength=true, nodelabel=df, nodelabeladj=[1.2,-.2]); # hide
R"mtext"("useedgelength = true, curved = :none", side=3, line=-1.5); # hide
plot(net, useedgelength=true, curved=:none, style=:fulltree, showedgelength=true, nodelabel=df, nodelabeladj=[1.2,-.2]); # hide
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_len_example.svg)

!!! note "node N"
    We used a DataFrame (not shown) to add the label "N".
    For more on this, see the section on [Adding labels](@ref).

If branch lengths represent time, D could represent a fossil, or a virus strain sequenced
a year before the others. Seeing this visually is the advantage of `useedgelength=true`.

This network happens to be time consistent, because the distance
along the time (x) axis from node `N` to the hybrid node is
the same both ways: the "upper" path has length 0.2 + 1.8 = 2,
which is the same along the "lower" path, 1.5 + 0.5 = 2.
We used option `showedgelength=true` to annotate the edges with their length.

!!! note "time consistency"
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

Time-inconsistent networks like these ones below might cause confusion.
Below we use the `:fulltree` style for the minor hybrid edges to have a
straight segment showing (proportional to) their length.

```@example better_edges
net1 = readnewick("(A:3.3,((B:1.5,#H1:1.2):1.5,((C:1.8)#H1:1,D:1.1):.2):0.3);");
net2 = readnewick("(A:3.3,((B:1.5,#H1:0.2):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);");
R"svg"(figname("edge_len_example2.svg"), width=6, height=3) # hide
R"layout"([1 2]) # hide
R"par"(mar=[0,0,0,0], cex=0.8) # hide
plot(net1, style=:fulltree, curved=:minor,
     useedgelength=true, showedgelength=true);
R"mtext"("net1", side=3, line=-2); # hide
plot(net2, style=:fulltree, curved=:minor,
     useedgelength=true, showedgelength=true);
R"mtext"("net2", side=3, line=-2); # hide
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/edge_len_example2.svg)

The default `style=:majortree` simplifies the visualization, as it
does not visually represent minor edge lengths.
Because of this, the option `showedgelength=true` to annotate each edge
with its length gives us the information anyway.

```@example better_edges
R"svg"(figname("edge_len_example3.svg"), width=6, height=3) # hide
R"layout"([1 2]) # hide
R"par"(mar=[0,0,0,0], cex=0.8) # hide
plot(net1, useedgelength=true, showedgelength=true);
R"mtext"("net1", side=3, line=-2); # hide
plot(net2, useedgelength=true, showedgelength=true);
R"mtext"("net2", side=3, line=-2); # hide
R"dev.off()" # hide
nothing # hide
```
![example4](../assets/figures/edge_len_example3.svg)

## Varying edge widths

We can vary edge widths to show population sizes for example.
First we need to map each edge number to the desired width for that edge.
We do this with a dictionary.

```@repl better_edges
R"svg"(figname("edge_len_example5.svg"), width=6, height=3) # hide
using RCall # to send any command to R, to modify the plot
R"par"(mar=[.1,0,0,0]); R"layout"([1 2]);
plot(net1, showedgenumber=true);
R"mtext"("edge numbers, used\nas keys in edgewidth", side=1, line=-1);
# below: population sizes on the log scale
log_populationsize = Dict(e.number => log10(1_000) for e in net1.edge);
log_populationsize[9] = log10(100_000); # larger populations on edge 9
log_populationsize[1] = log10(100_000); #                and on edge 1
log_populationsize
plot(net1, edgewidth=log_populationsize);
R"dev.off()"; # hide
nothing # hide
```
![example5](../assets/figures/edge_len_example5.svg)

## Customization

Check out the list of [`plot`](@ref) options.

In the example below,
we first highlight in orange the edges on the 2 paths from the root to C.
Then we change the type of the minor edge (to hide it).

```@repl better_edges
ecols = Dict(i => "black" for i in 1:9); # make all black
for i in [9,8,6,5, 4,3] # except for edges ancestral to C
  ecols[i] = "orangered"
end
ecols
```
```@example better_edges
R"svg"(figname("edge_len_example6.svg"), width=6, height=3) # hide
R"par"(mar=[.1,0,0,0]); R"layout"([1 2]); # hide
plot(net1, edgecolor=ecols, defaultedgecolor="grey80",
     minorlinetype="solid");
R"mtext"("default curved = :both", side=3, line=-2); # hide
plot(net1, majorhybridedgecolor="red", defaultedgecolor="grey80",
     minorlinetype="blank", curved=:minor); # make minor edges (arrows) of type 'blank'
R"mtext"("curved = :minor", side=3, line=-2); # hide
R"mtext"("minor hybrid edge is\nhidden: 'blank' type", side=1, line=-1); # hide
R"dev.off()"; # hide
nothing # hide
```
![example6](../assets/figures/edge_len_example6.svg)
