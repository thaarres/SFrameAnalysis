import sys
import os, commands
import shutil
from ROOT import *
import ROOT

indir = "/scratch/thaarres/VTopTagSF_MiniTuple/reweighted/HaddedOutput/"

ttsamples = ["powheg","herwig","madgraph"]
ttsamples = [""]

for sample in ttsamples:
  
  infiles = ['WJetsToLNu.root','ST.root','TT.root','VV.root']

          
  inTreeName = 'tree'
  targetFile = "pseudodata_weighted"
  addcmd = "hadd -f %s.root " %(targetFile)
  rmcmd = "rm "

  n = 0 
  for inFileName in infiles:
    tmpname = indir+"tmp_" + inFileName
    path = indir + inFileName
    print "copying file " ,tmpname
    shutil.copy2("%s" %path, "%s"%tmpname)
    filetmp = TFile.Open( "%s" %tmpname, "update" )
    myTree = filetmp.Get( inTreeName )
    print "BEFORE " ,myTree.GetWeight()
    weight = 1.
    if tmpname.find("TT") != -1:
		weight = 0.705
		print "Setting ttbar weight to " , weight
		myTree.SetWeight(weight)
    print "AFTER " , myTree.GetWeight()
    myTree.AutoSave()

    addcmd+= ' %s' %(tmpname)
    rmcmd += ' %s' %(tmpname)
    filetmp.Close()
    del myTree
  
  print addcmd
  os.system(addcmd)
  print "Removing temporary files: ..."
  print rmcmd
  os.system(rmcmd)
