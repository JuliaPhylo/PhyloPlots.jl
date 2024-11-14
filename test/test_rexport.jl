# fixit: also test @rput net2

@testset "rexport and sexp tests" begin

useape = false; # these extra tests would require the installation of ape by the CI platform,
# but ape itself requires a number of R packages.
# would make the tests taken longer. overkill.
# use it locally to run more checks when desired.

# on a tree with some edge lengths missing
s = "(A,(B:1.0,(C:1.0,D:1.0):1.0):1.0);";
tree1 = PhyloPlots.rexport(readnewick(s)); # for testing rexport
net2 = readnewick(s); # for testing sexp on clean network
if useape
## check for correct unrooted topology and for identical lengths
## of (unrooted) internal edges, but requires t
    R"library(ape)"
    R"tree2 = read.tree(text = $s)";
    @test convert(Bool, R"dist.topo($tree1, tree2, method='score') == 0")
    @test convert(Bool, R"dist.topo($net2, tree2, method='score') == 0")
end
R"""
tree1r = list(Nnode=3L, edge=matrix(c(5,5,6,6,7,7, 1,6,2,7,3,4), 6,2),
             tip.label=c("A","B","C","D"), edge.length=c(NA,1,1,1,1,1))
storage.mode(tree1r$edge) <- "integer"
class(tree1r) = "phylo"
tree2r = tree1r[c(2,1,4,3)]
class(tree2r) = "phylo"
""";
@test convert(Bool, R"isTRUE(all.equal($tree1,tree1r))")
@test convert(Bool, R"isTRUE(all.equal($net2, tree2r))")
@rput net2
@test convert(Bool, R"isTRUE(all.equal(net2, tree2r))")

# network, h=1, missing gamma values
s = "(((A:4.0,(B:1.0)#H1:1.1):0.5,(C:0.6,#H1:1.0):1.0):3.0,D:5.0);";
phy1 = PhyloPlots.rexport(readnewick(s));
net2 = readnewick(s);
if useape # needs ape version > 4.1 with read.evonet (not in 4.1)
    R"phyr = read.evonet(text = $s)"
    @test convert(Bool, R"dist.topo($phy1, phyr, method='score') == 0")
    @test convert(Bool, R"dist.topo($net2, phyr, method='score') == 0")
end
R"""
phy1r = list(Nnode=5L, edge=matrix(c(5,5,6,6,7,8,8,9, 6,4,8,7,3,1,9,2), 8,2),
             tip.label=c("A","B","C","D"), edge.length=c(3,5,.5,1,.6,4,1.1,1),
             reticulation=matrix(c(7L,9L),1,2), reticulation.length=1)
storage.mode(phy1r$edge) <- "integer"
class(phy1r) = c("evonet","phylo")
phy2r = phy1r[c(2,6,1,4,5,3)]
class(phy2r) = c("evonet","phylo")
""";
@test convert(Bool, R"isTRUE(all.equal($phy1, phy1r))")
@test convert(Bool, R"isTRUE(all.equal($net2, phy2r))")
@rput net2
@test convert(Bool, R"isTRUE(all.equal( net2, phy2r))")

# network, h=1, maintree=true; minor hybrid edge length missing
s = "(((A:4.0,(B)#H1:1.1::0.9):0.5,(C:0.6,#H1):1.0):3.0,D:5.0);";
stree = "(((A:4.0,B):0.5,C:1.6):3.0,D:5.0);"; # main tree
phy1 = PhyloPlots.rexport(readnewick(s); maintree=true);
net2 = readnewick(s);
tree2 = majorTree(readnewick(s));
if useape
    R"tree2r = read.tree(text = $stree)";
    @test convert(Bool, R"dist.topo($phy1, tree2r, method='score') == 0")
    @test convert(Bool, R"dist.topo($tree2, tree2r, method='score') == 0")
end
R"""
tree1r = list(Nnode=3L, edge=matrix(c(5,5,6,6,7,7, 6,4,7,3,1,2), 6,2),
             tip.label=c("A","B","C","D"), edge.length=c(3,5,.5,1.6,4,NA))
storage.mode(tree1r$edge) <- "integer"
class(tree1r) = "phylo"
phy2r = list(edge=matrix(c(5,5,6,6,7,8,8,9, 6,4,8,7,3,1,9,2), 8,2),
             Nnode=5L, edge.length=c(3,5,.5,1,.6,4,1.1,NA),
             reticulation=matrix(c(7L,9L),1,2),
             reticulation.gamma=0.1,
             tip.label=c("A","B","C","D"))
storage.mode(phy2r$edge) <- "integer"
class(phy2r) = c("evonet","phylo")
""";
@test convert(Bool, R"isTRUE(all.equal($phy1, tree1r))")
@test convert(Bool, R"isTRUE(all.equal($net2, phy2r))")
@rput net2
@test convert(Bool, R"isTRUE(all.equal( net2, phy2r))")

# on a network, h=1, with useedgelength=false
s = "(((A:1.0,(B:.5)#H1:.2::0.9):.8,(C:1.5,#H1:.01::0.1):2):.6,D:5.0);";
phy1 = PhyloPlots.rexport(readnewick(s), useedgelength=false);
if useape # needs ape version > 4.1 with read.evonet (not in 4.1), also dist.topo for networks
    R"phyr = read.evonet(text = $s)"
    @test convert(Bool, R"dist.topo($phy1, phyr) == 0")
end
R"""
phy1r = list(Nnode=5L, edge=matrix(c(5,5,6,6,7,8,8,9, 6,4,8,7,3,1,9,2), 8,2),
             tip.label=c("A","B","C","D"), reticulation=matrix(c(7L,9L),1,2),
             reticulation.gamma=0.1)
storage.mode(phy1r$edge) <- "integer"
class(phy1r) = c("evonet","phylo")
""";
@rput phy1
@test convert(Bool, R"isTRUE(all.equal(phy1, phy1r))")

# on a network with h=2 hybridizations
s = "(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,
      (((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);";
phy1 = PhyloPlots.rexport(readnewick(s));
net2 = readnewick(s);
R"""
phy1r = list(Nnode=14L, edge=matrix(c(13,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20,21,21,22,23,24,25,25,26,26,
             15,14,12,10,11,18,16,17,9,24,8,1,19,20,21,23,2,22,3,4,25,26,7,5,6), 25,2),
             tip.label=c("Ag","Ak","E","M","Az","Ag2","As","Ap","Ar","P","20","165"),
             edge.length=c(5.943,1.863,NA,NA,NA,.026,.723,.187,NA,0,NA,NA,2.214,2.490,.078,NA,NA,
                           .023,.08,0,1.697,2.11,2.027,.002,.023),
             reticulation=matrix(c(19,22, 24,23),2,2),
             reticulation.gamma=c(0.056,0.004),
             reticulation.length=c(7.159,0.0))
storage.mode(phy1r$edge) <- "integer"
storage.mode(phy1r$reticulation) <- "integer"
class(phy1r) = c("evonet","phylo")
phy2r = phy1r[c(2,7,1,4,5,6,3)]
class(phy2r) = c("evonet","phylo")
""";
@test convert(Bool, R"identical($phy1, phy1r)")
@test convert(Bool, R"identical($net2, phy2r)")

end
