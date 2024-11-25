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
net = readnewick("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, style=:majortree);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/style_example.svg)

## Using edge lengths

We can use `useedgelength=true` to draw a plot that uses the network's edge lengths to determine the lengths of the
lines. For this, we'll use a network that has branch lengths:

```@example better_edges
R"svg"(figname("edge_len_example.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2]) # hide
net = readnewick("(A:3.3,((B:1.5,#H1:0.5):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);")
df = DataFrame(number=[-3,3], label=["N","H1"]); # hide
plot(net, useedgelength=true, ylim = [-1, 5.5], nodelabel = df); # hide
R"text"([3], [0], ["useedgelength=true"]) # hide
plot(net, useedgelength=false, ylim = [-1, 5.5], nodelabel = df); # hide
R"text"([3], [0], ["useedgelength=false"]) # hide
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_len_example.svg)

!!! note
    I used a DataFrame (not shown) to add the label "N" to the plot.
    For more on this, see the [Adding labels](@ref) section.

If branch lengths represent time, D could represent a fossil, or a virus strain sequenced
a year before the others. Seeing this visually is the advantage of `useedgelength=true`.

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

Time-inconsistent networks like these ones below might cause confusion:

```@example better_edges
R"svg"(figname("edge_len_example2.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2]) # hide
net1 = readnewick("(A:3.3,((B:1.5,#H1:1.2):1.5,((C:1.8)#H1:1,D:1.1):.2):0.3);");
net2 = readnewick("(A:3.3,((B:1.5,#H1:0.2):1.5,((C:1)#H1:1.8,D:1.1):.2):0.3);");
plot(net1, useedgelength=true); # hide
plot(net2, useedgelength=true); # hide
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/edge_len_example2.svg)

It may be useful to consider using `style=:majortree` if it causes
too much confusion, since the `:majortree` style doesn't visually represent
minor edge lengths. Because of this, I used the `showedgelength=true` option to
see the information anyway.

```@example better_edges
R"svg"(figname("edge_len_example3.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
R"layout"([1 2])
plot(net1, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1);
plot(net2, useedgelength=true, style = :majortree, showedgelength=true, arrowlen=0.1);
R"dev.off()" # hide
nothing # hide
```
![example4](../assets/figures/edge_len_example3.svg)

I also used the `arrowlen=0.1` option to show the arrow tips to show the direction of minor edges,
which are hidden by default when using the `style=:majortree` option.

## Varying edge widths

We can vary edge widths to show population sizes for example.
First we need to map each edge number to the desired width for that edge.
We do this with a dictionary.

```@repl better_edges
R"svg"(figname("edge_len_example5.svg"), width=6, height=3) # hide
using RCall # to send any command to R, to modify the plot
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]);
plot(net1, showedgenumber=true);
R"mtext"("edge numbers, used\nas keys in edgewidth", side=1, line=-1);
log_populationsize = Dict(e.number => log10(1_000) for e in net1.edge); # pop size on log scale
log_populationsize[9] = log10(100_000); # larger populations on edges 9 and 1
log_populationsize[1] = log10(100_000);
log_populationsize
plot(net1, edgewidth=log_populationsize);
R"dev.off()"; # hide
nothing # hide
```
![example5](../assets/figures/edge_len_example5.svg)
