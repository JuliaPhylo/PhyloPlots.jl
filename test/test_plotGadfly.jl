@testset "Gadfly-based plot Test" begin
  # fist on a tree, then on a network with h=2 and various options"
  # keep ";" below to avoid creating a new file
  tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
  net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");
  @test_nowarn plot(tre);
  @test_nowarn plot(net);
  @test_warn "At least one non-missing edge length: plotting any missing length as 1.0" plot(net, useEdgeLength=true);
  @test_nowarn plot(net, mainTree=true);
  @test_nowarn plot(net, showTipLabel=false);
  @test_nowarn plot(net, showNodeNumber=true);
  @test_nowarn plot(net, showEdgeLength=false, showEdgeNumber=true);
  @test_nowarn plot(net, showGamma=true);
  @test_nowarn plot(net, edgeColor=colorant"olive", # uses Colors
          minorHybridEdgeColor=colorant"tan",
          majorHybridEdgeColor=colorant"skyblue");
end
