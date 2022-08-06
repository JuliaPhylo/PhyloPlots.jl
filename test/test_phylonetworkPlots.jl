@testset "Test setup for plotting PhyloNetworks objects" begin

  net = readTopology("(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);")
  # in this order:
  # edge_xB, edge_xE, edge_yB, edge_yE,
  # node_x, node_y, node_yB, node_yE,
  # minoredge_xB, minoredge_xE, minoredge_yB, minoredge_yE,
  # xmin, xmax, ymin, ymax
  @test PhyloPlots.edgenode_coordinates(net, true, false) == (
    [1.0, 2.5, 2.5, 1.5, 2.5, 3.0, 2.5, 1.5, 1.0],
    [3.5, 3.5, 3.0, 2.5, 3.5, 3.5, 3.0, 2.5, 1.5],
    [1.0, 2.0, 3.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5],
    [1.0, 2.0, 3.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5],
    [3.5, 3.5, 2.5, 3.5, 3.5, 3.0, 2.5, 1.5, 1.0],
    [1.0, 2.0, 2.5, 4.0, 5.0, 5.0, 4.5, 3.5, 2.25],
    [1.0, 2.0, 2.0, 4.0, 5.0, 5.0, 4.0, 2.5, 1.0],
    [1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 5.0, 4.5, 3.5],
    [3.0], [3.0], [3.0], [5.0], 1.0, 3.5, 1.0, 5.0)
  @test PhyloPlots.edgenode_coordinates(net, true, true) == (
    [1.0, 2.5, 2.5, 1.5, 2.5, 3.0, 2.5, 1.5, 1.0],
    [3.5, 3.5, 2.5, 2.5, 3.5, 3.5, 3.0, 2.5, 1.5],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75],
    [3.5, 3.5, 2.5, 3.5, 3.5, 3.0, 2.5, 1.5, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75, 1.875],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.5, 2.75],
    [2.5], [3.0], [2.0], [4.0], 1.0, 3.5, 1.0, 4)
  @test PhyloPlots.edgenode_coordinates(net, false, true) == (
    [1.0, 3.0, 3.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0],
    [5.0, 5.0, 3.0, 3.0, 5.0, 5.0, 4.0, 3.0, 2.0],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75],
    [1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75],
    [5.0, 5.0, 3.0, 5.0, 5.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.5, 2.75, 1.875],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0],
    [1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.5, 2.75],
    [3.0], [4.0], [2.0], [4.0], 1.0, 5.0, 1.0, 4)
  dat = DataFrame(node=[-5,-3,-4,5,100],bs=["90","95","99","mytip","bogus"],edge=[8,9,4,6,200]);
  @test_logs (:warn, "Some node numbers in the nodelabel data frame are not found in the network:\n 100") PhyloPlots.check_nodedataframe(net, dat);
  @test_logs (:warn, "nodelabel should have 2+ columns, the first one giving the node numbers (Integer)") PhyloPlots.check_nodedataframe(net, dat[!,2:3])
  dat = DataFrame(node=[-5,-3,missing,5],
                  bs=["90","95","99","mytip"],edge=[8,9,4,6]);
  @test PhyloPlots.check_nodedataframe(net, dat) == (true,
    DataFrame(node=[-5,-3,5],bs=["90","95","mytip"],edge=[8,9,6]))
  dat = DataFrame(node=[-5,-3,-4,5],bs=["90","95","99","mytips"]);
  @test PhyloPlots.prepare_nodedataframe(net,dat,true,true,true,collect(1.:9),collect(10.:18)) ==
    DataFrame(name=["A","B","","C","D","H1","","",""],
    num=["1","2","-4","4","5","3","-5","-3","-2"],
    lab=["","","99","","mytips","","90","95",""],
    lea=[true,true,false,true,true,false,false,false,false],
    x=collect(1.:9), y=collect(10.:18))
  dat = DataFrame(edge=[8,9,4,6,200],bs=["90","95","99","mytips","bogus"]);
  @test_logs (:warn, "Some edge numbers in the edgelabel data frame are not found in the network:\n 200") PhyloPlots.prepare_edgedataframe(
    net,dat,true,collect(1.:9),collect(10.:18),collect(19.:27),collect(28.:36),[6.5],[8.5],[24.5],[26.5]);
  dat = DataFrame(edge=[8,9,4,6],bs=[missing,"95","99","mytips"]);
  @test PhyloPlots.prepare_edgedataframe(net,dat,false,collect(1.:9),collect(11.:19),
    collect(21.:29),collect(31.:39),[7.],[9.],[27.],[29.]) == (true, DataFrame(
    len=["2.5","1","0.5","1","1","0.5","0.5","1","0.5"],
    gam=["1","1","0.1","1","1","1","0.9","1","1"],
    num=["1","2","3","4","5","6","7","8","9"],
    lab=["","","","99","","mytips","","","95"], # if not missing: "90" second to last (row 8)
    hyb=[false,false,true,false,false,false,true,false,false],
    min=[false,false,true,false,false,false,false,false,false],
    x=collect(6.:14),y=collect(26.:34)))

  # example with level-2 network, non-tree child:
  # one hybrid node ends up as a leaf in the major tree.
  # no major child edge to follow to set coordinates
  net = readTopology("((((B)#H1:::0.2)#H2,((D,C,#H2:::0.8)S1,(#H1,A)S2)S3)S4);")
  @test_logs plot(net, shownodenumber=true, showgamma=true);
  @test PhyloPlots.edgenode_coordinates(net, false, false) == (
    [5.0, 4.0, 1.0, 3.0, 3.0, 3.0, 2.0, 4.0, 4.0, 2.0, 1.0],
    [6.0, 5.0, 4.0, 6.0, 6.0, 4.0, 3.0, 5.0, 6.0, 4.0, 2.0],
    [5.0, 4.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.25],
    [5.0, 4.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.25],
    [6.0, 5.0, 4.0, 6.0, 6.0, 3.0, 6.0, 4.0, 2.0, 1.0],
    [5.0, 5.0, 4.0, 2.0, 3.0, 3.0, 6.0, 5.5, 4.25, 2.625],
    [5.0, 5.0, 4.0, 2.0, 3.0, 2.0, 6.0, 5.0, 3.0, 1.0],
    [5.0, 5.0, 4.0, 2.0, 3.0, 4.0, 6.0, 6.0, 5.5, 4.25],
    [5.0, 4.0], [5.0, 4.0], [4.0, 1.0], [5.0, 4.0],
    1.0, 6.0, 1.0, 6.0)
  net = readTopology("((((B)#H1:::0.2)#H2,((D,C,#H2)S1,(#H1,A)S2)S3)S4);")
  @test_logs plot(net, shownodenumber=true, showgamma=true);
  @test PhyloPlots.edgenode_coordinates(net, false, false) == (
    [5.0, 4.0, 1.0, 3.0, 3.0, 3.0, 2.0, 4.0, 4.0, 2.0, 1.0],
    [6.0, 5.0, 4.0, 6.0, 6.0, 4.0, 3.0, 5.0, 6.0, 4.0, 2.0],
    [5.0, 1.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.25],
    [5.0, 1.0, 1.0, 2.0, 3.0, 4.0, 3.0, 5.0, 6.0, 5.5, 4.25],
    [6.0, 5.0, 4.0, 6.0, 6.0, 3.0, 6.0, 4.0, 2.0, 1.0],
    [5.0, 5.0, 1.0, 2.0, 3.0, 3.0, 6.0, 5.5, 4.25, 2.625],
    [5.0, 5.0, 1.0, 2.0, 3.0, 2.0, 6.0, 5.0, 3.0, 1.0],
    [5.0, 5.0, 1.0, 2.0, 3.0, 4.0, 6.0, 6.0, 5.5, 4.25],
    [5.0, 4.0], [5.0, 4.0], [1.0, 4.0], [5.0, 1.0],
    1.0, 6.0, 1.0, 6.0)
end
