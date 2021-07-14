```@setup better_edges
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Better edges

## Different hybrid edge styles

We can use the `style` parameter to visualize minor hybrid edges as simple lines, unlike the icytree.org style visualization. `style` is by default `:fulltree`, but by switching it 
to `:majortree`, we can draw minor hybir edges as simple lines.

```@example better_edges
R"svg"(figname("style_example.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, style=:majortree);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/style_example.svg)

## Using edge lengths

We can use the `useEdgeLength=true` parameter to draw a plot that uses the network's edge lengths to determine the lengths of the lines. For this, we'll use a network that can be found [here](https://github.com/nkarimi/Adansonia_HybSeq/blob/master/trait-evolution/BestH1_372g_calibrated.tre).

```@example better_edges
R"svg"(figname("edge_len_example.svg"), width=6, height=6) # hide
net = readTopology("((Smi165:1.6261423761885154,Pcr070:1.6261423761885154):0.0345640579033647,(((#H18:0.23725347915651637::0.12440503333556951,(Adi001:0.2497235156848997,(Adi003:0.22937089525518586,Adi002:0.22937089525518586):0.02035262042971384):0.36375851550073335):1.3644833483534315e-9,((((Asu001:0.3743115977039037,(Aga001:0.24655560654629088,Aga002:0.24655560654629088):0.1277559911576128):0.0019169543252129813)#H18:0.10531506732808882::0.8755949666644305,(Aza135:0.35341877971983404,((Aza037:0.2766129548639923,(Ama018:0.2198029990484657,Ama006:0.2198029990484657):0.056809955815526614):0.030764173677767633,(Ape009:0.2219108012530084,Ape001:0.2219108012530084):0.08546632728875153):0.04604165117807412):0.12812483963737142):0.016691432891123632,(Aru001:0.2294802905252357,Aru127:0.2294802905252357):0.2687547617230934):0.11524698030178727):0.3865179674498836,Age001:1.0):0.6607064340918802);") # hide

plot(net, :R, useEdgeLength=true);
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_len_example.svg)