@testset "Test setup for plotting PhyloNetworks objects" begin

  net = readTopology("(A:2.5,((B:1,#H1:0.5::0.1):1,(C:1,(D:0.5)#H1:0.5::0.9):1):0.5);")
  @test PhyloPlots.getEdgeNodeCoordinates(net, true) == (
      [1., 2.5,2.5,1.5,2.5,3., 2.5,1.5, 1.],
      [3.5,3.5,3., 2.5,3.5,3.5,3., 2.5, 1.5],
      [4., 3., 3., 3., 2., 1., 1., 1.5, 2.25],
      [4., 3., 1., 3., 2., 1., 1., 1.5, 2.25],
      [3.5,3.5,2.5,3.5,3.5,3., 2.5,1.5, 1.],
      [4., 3., 3., 2., 1., 1., 1.5,2.25,3.125],
      [0., 0., 3., 0., 0., 1., 1., 1.5, 2.25],
      [0., 0., 3., 0., 0., 1., 2., 3.,  4.],
      1., 3.5, 1., 4.)
  # in this order:
  # edge_xB, edge_xE, edge_yB, edge_yE,
  # node_x, node_y, node_yB, node_yE,
  # xmin, xmax, ymin, ymax
  @test PhyloPlots.getEdgeNodeCoordinates(net, false) == (
      [1.,3.,3.,2.,3.,4.,3.,2., 1.],
      [5.,5.,4.,3.,5.,5.,4.,3., 2.],
      [4.,3.,3.,3.,2.,1.,1.,1.5,2.25],
      [4.,3.,1.,3.,2.,1.,1.,1.5,2.25],
      [5.,5.,3.,5.,5.,4.,3., 2.,  1.],
      [4.,3.,3.,2.,1.,1.,1.5,2.25,3.125],
      [0.,0.,3.,0.,0.,1.,1., 1.5, 2.25],
      [0.,0.,3.,0.,0.,1.,2., 3.,  4.],
      1., 5.,1., 4.)
  dat = DataFrame(node=[-5,-3,-4,5,100],bs=["90","95","99","mytip","bogus"],edge=[8,9,4,6,200]);
  @test_warn "Some node numbers in the nodeLabel data frame are not found in the network:\n 100" PhyloPlots.checkNodeDataFrame(net, dat);
  @test_warn "nodeLabel should have 2+ columns, the first one giving the node numbers (Integer)" PhyloPlots.checkNodeDataFrame(net, dat[2:3])
  dat = DataFrame(node=[-5,-3,missing,5],
                  bs=["90","95","99","mytip"],edge=[8,9,4,6]);
  @test PhyloPlots.checkNodeDataFrame(net, dat) == (true,
    DataFrame(node=[-5,-3,5],bs=["90","95","mytip"],edge=[8,9,6]))
  dat = DataFrame(node=[-5,-3,-4,5],bs=["90","95","99","mytips"]);
  @test PhyloPlots.prepareNodeDataFrame(net,dat,true,true,true,collect(1.:9),collect(10.:18)) ==
    DataFrame(name=["A","B","","C","D","#H1","","",""],
    num=["1","2","-4","4","5","3","-5","-3","-2"],
    lab=["","","99","","mytips","","90","95",""],
    lea=[true,true,false,true,true,false,false,false,false],
    x=collect(1.:9), y=collect(10.:18))
  dat = DataFrame(edge=[8,9,4,6,200],bs=["90","95","99","mytips","bogus"]);
  @test_warn "Some edge numbers in the edgeLabel data frame are not found in the network:\n 200" PhyloPlots.prepareEdgeDataFrame(
    net,dat,true,collect(1.:9),collect(10.:18),collect(19.:27),collect(28.:36));
  dat = DataFrame(edge=[8,9,4,6],bs=[missing,"95","99","mytips"]);
  @test PhyloPlots.prepareEdgeDataFrame(net,dat,false,collect(1.:9),collect(11.:19),
    collect(21.:29),collect(31.:39)) == (true, DataFrame(
    len=["2.5","1","0.5","1","1","0.5","0.5","1","0.5"],
    gam=["1","1","0.1","1","1","1","0.9","1","1"],
    num=["1","2","3","4","5","6","7","8","9"],
    lab=["","","","99","","mytips","","","95"], # if not missing: "90" second to last (row 8)
    hyb=[false,false,true,false,false,false,true,false,false],
    min=[false,false,true,false,false,false,false,false,false],
    x=collect(6.:14),y=collect(26.:34)))

end
