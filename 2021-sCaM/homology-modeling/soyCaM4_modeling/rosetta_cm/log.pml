load 3clnA.pdb
load S_0001.pdb
load S_0002.pdb
load S_0004.pdb
load S_0003.pdb
load S_0004.pdb
load S_0005.pdb
load S_0006.pdb
load S_0007.pdb
load S_0008.pdb
load S_0009.pdb
load S_0010.pdb
show cartoon, 
hide lines, 
cmd.align("polymer and name ca and (S_0001)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0001_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0001)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0001_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0004)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0004_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0003)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0003_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0005)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0005_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0006)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0006_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0007)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0007_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0008)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0008_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0009)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0009_to_3clnA",reset=1)
cmd.align("polymer and name ca and (S_0010)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0010_to_3clnA",reset=1)
cmd.center("(S_0004`1882)",state=-1)
cmd.center("(S_0003`1797)",state=-1)
cmd.align("polymer and name ca and (S_0002)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0002_to_3clnA",reset=1)
cmd.disable('all')
cmd.enable('3clnA',1)
color red, S_00*
cmd.enable('S_0001',1)
cmd.enable('S_0002',1)
cmd.enable('S_0003',1)
cmd.disable('S_0003')
cmd.enable('S_0003',1)
cmd.disable('S_0003')
cmd.enable('S_0003',1)
cmd.disable('S_0003')
cmd.enable('S_0003',1)
cmd.align("polymer and name ca and (S_0003)","polymer and name ca and (3clnA)",quiet=0,object="aln_S_0003_to_3clnA",reset=1)
cmd.disable('S_0003')
cmd.enable('S_0003',1)
cmd.enable('S_0005',1)
cmd.disable('S_0003')
cmd.enable('S_0004',1)
cmd.enable('S_0006',1)
cmd.disable('S_0006')
cmd.enable('S_0006',1)
cmd.disable('S_0006')
cmd.enable('S_0006',1)
cmd.disable('S_0006')
cmd.enable('S_0007',1)
cmd.disable('S_0007')
cmd.enable('S_0007',1)
cmd.disable('S_0007')
cmd.enable('S_0007',1)
cmd.disable('S_0007')
cmd.enable('S_0008',1)
cmd.disable('S_0008')
cmd.enable('S_0008',1)
cmd.enable('S_0009',1)
cmd.enable('S_0010',1)
cmd.disable('S_0010')
color cyan, 3clnA
cmd.disable('aln_S_0003_to_3clnA')
bg_color white
cmd.center("(3clnA`1091)",state=-1)
cmd.center("(S_0005`1171)",state=-1)
cmd.center("(S_0002`1197)",state=-1)
cmd.center("(S_0009`1235)",state=-1)
set cartoon_transparency, 0.1, S_0001
set cartoon_transparency, 0.2, S_0002
set cartoon_transparency, 0.3, S_0003
set cartoon_transparency, 0.4, S_0004
set cartoon_transparency, 0.4, S_0005
set cartoon_transparency, 0.4, S_0006
set cartoon_transparency, 0.4, S_0007
set cartoon_transparency, 0.4, S_0008
set cartoon_transparency, 0.4, S_0009
set cartoon_transparency, 0.4, S_0010
cmd.disable('S_0004')
cmd.disable('S_0002')
cmd.disable('S_0001')
cmd.enable('S_0003',1)
cmd.disable('S_0005')
cmd.enable('S_0006',1)
cmd.enable('S_0007',1)
cmd.disable('S_0008')
cmd.disable('S_0009')
cmd.enable('S_0010',1)
cmd.align("polymer and name ca and (S_0006)","polymer and name ca and (S_0003)",quiet=0,object="aln_S_0006_to_S_0003",reset=1)
cmd.align("polymer and name ca and (S_0007)","polymer and name ca and (S_0003)",quiet=0,object="aln_S_0007_to_S_0003",reset=1)
cmd.align("polymer and name ca and (S_0010)","polymer and name ca and (S_0003)",quiet=0,object="aln_S_0010_to_S_0003",reset=1)
cmd.disable('aln_S_0006_to_S_0003')
cmd.disable('aln_S_0007_to_S_0003')
cmd.enable('aln_S_0007_to_S_0003',1)
cmd.disable('aln_S_0007_to_S_0003')
cmd.disable('aln_S_0010_to_S_0003')
cmd.disable('3clnA')
cmd.disable('S_0006')
cmd.disable('S_0007')
cmd.disable('S_0010')
cmd.enable('3clnA',1)
cmd.disable('3clnA')
cmd.enable('S_0006',1)
cmd.disable('S_0006')
cmd.enable('S_0006',1)
cmd.enable('S_0009',1)
cmd.disable('S_0009')
cmd.enable('S_0008',1)
cmd.disable('S_0008')
cmd.enable('S_0008',1)
cmd.disable('S_0008')
cmd.enable('S_0007',1)
cmd.enable('S_0010',1)
cmd.center("(S_0003`139)",state=-1)
cmd.disable('S_0006')
cmd.enable('S_0006',1)
