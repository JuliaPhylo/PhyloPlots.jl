using RCall # before using PhyloNetworks and PhyloPlots
using PhyloPlots
using Base.Test

using DataFrames
using PhyloNetworks # to read topologies
using Colors # is REQUIREd by Gadfly

@testset "PhyloPlots Tests" begin
  # each file should be its own testset, to run them all even if one has a failure
  include("test_phylonetworkPlots.jl")
  include("test_plotsGadfly.jl") # uses Colors
  include("test_plotsRCall.jl")  # uses DataFrames to test annotations
end
