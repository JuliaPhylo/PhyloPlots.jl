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

  # example with weird plot
  cui3str = "(Xmayae,((Xhellerii,(((Xclemenciae_F2,Xmonticolus):1.458,(((((Xmontezumae,(Xnezahuacoyotl)#H26:0.247::0.804):0.375,((Xbirchmanni_GARC,Xmalinche_CHIC2):0.997,Xcortezi):0.455):0.63,(#H26:0.0::0.196,((Xcontinens,Xpygmaeus):1.932,(Xnigrensis,Xmultilineatus):1.401):0.042):2.439):2.0)#H7:0.787::0.835,(Xmaculatus,(Xandersi,(Xmilleri,((Xxiphidium,#H7:9.563::0.165):1.409,(Xevelynae,(Xvariatus,(Xcouchianus,(Xgordoni,Xmeyeri):0.263):3.532):0.642):0.411):0.295):0.468):0.654):1.022):0.788):1.917)#H27:0.149::0.572):0.668,Xalvarezi):0.257,(Xsignum,#H27:1.381::0.428):4.669);"
  net3  = readTopology(cui3str);
  deleteleaf!(net3,"Xhellerii"); deleteleaf!(net3,"Xsignum");
  deleteleaf!(net3,"Xmayae", simplify=false);
  @test_nowarn plot(net3); # looks weird: k=2 cycle at the root. 3 root edges:
  # one to a leaf, 1 major & 1 minor hybrid edge to the same child.

end
