@testset "RCall-based plot Test" begin
  # testing for absence of errors, not for correctness

  #net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);")
  net = readTopology("(((Ag,(#H1:7.159::0.056,((Ak,(E:0.08,#H2:0.0::0.004):0.023):0.078,(M:0.0)#H2:::0.996):2.49):2.214):0.026,(((((Az:0.002,Ag2:0.023):2.11,As:2.027):1.697)#H1:0.0::0.944,Ap):0.187,Ar):0.723):5.943,(P,20):1.863,165);");
  @test_logs plot(net,:RCall);
  @test_logs (:warn, "At least one non-missing edge length: plotting any missing length as 1.0") plot(net,:RCall, useEdgeLength=true);
  @test_logs plot(net,:RCall, showTipLabel=false);
  @test_logs plot(net,:RCall, showNodeNumber=true, showIntNodeLabel=true);
  @test_logs plot(net,:RCall, tipOffset=1, showGamma=true);
  @test_logs plot(net,:RCall, showEdgeLength=true, showEdgeNumber=true);
  @test_logs plot(net,:RCall, edgeColor="tomato4",minorHybridEdgeColor="skyblue",
          majorHybridEdgeColor="tan");
  dat = DataFrame(node=[-5,-10,-1],bs=["90","95","100"],edge=[11,22,26]);
  @test_logs (:warn, "Some node numbers in the nodeLabel data frame are not found in the network:\n -1") plot(net,:RCall, nodeLabel=dat);
  @test_logs plot(net,:RCall, edgeLabel=dat[!,[:edge,:bs]]);

  # plot based on RCall and ape:
  tre = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
  # fixit: plot(tre, :ape)
end
