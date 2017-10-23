#! /usr/bin/env python

import os, multiprocessing, math, sys
from array import array
from ROOT import TFile, TH1, TF1, TLorentzVector



import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--input', action='store', type='string', dest='origin', default='/scratch/thaarres/VTopTagSF_MiniTuple/')
parser.add_option('-o', '--output', action='store', type='string', dest='target', default='/scratch/thaarres/VTopTagSF_MiniTuple/reweighted')
parser.add_option('-f', '--filter', action='store', type='string', dest='filter', default='')
parser.add_option('-s', '--single', action='store_true', dest='single', default=True)
parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)

(options, args) = parser.parse_args()

origin      = options.origin
target      = options.target
filterset   = options.filter
singlecore  = options.single
verboseout  = options.verbose

LUMI    = 35867.
ANALYSIS = "VVanalysis"

if not os.path.exists(origin):
    print 'Origin directory', origin, 'does not exist, aborting...'
    exit()
if not os.path.exists(target):
    print 'Target directory', target,'does not exist, creating it...'
    cmd = "mkdir %s" %target
    os.system(cmd)


##############################

    
def processFile(sample_name, verbose=False):
    sample = sample_name.replace(ANALYSIS+".", "").replace(".root", "")
    isMC = not 'Single' in sample
    
    # Unweighted input
    ref_file_name = origin + '/' + sample_name
    if not os.path.exists(ref_file_name): 
        print '  WARNING: file', ref_file_name, 'does not exist, continuing'
        return True
    
    
    # Open old file
    ref_file = TFile(ref_file_name, 'READ')
    ref_hist = ref_file.Get('genEvents')
    totalEntries = ref_hist.GetBinContent(2)
    print "Using N generated events = " , totalEntries
    tree = ref_file.Get('tree')
    nev = tree.GetEntriesFast()
    tree.GetEntry(0)
    XS = tree.xsec
    print "Using cross section = " , XS
    
    # Weighted output
    new_file_name = target + '/' + sample + '.root'


    new_file = TFile(new_file_name, 'RECREATE')
    new_file.cd()

    Leq = LUMI*XS/totalEntries if totalEntries > 0 else 0.

    print sample, ": weight =", (Leq if isMC else "Data")
  
    # Variables declaration
    stitchWeight = array('f', [1.0])
    eventWeightLumi = array('f', [1.0])  # global event weight with lumi

    # Looping over file content
    for key in ref_file.GetListOfKeys():
        obj = key.ReadObj()
        
        if obj.IsA().InheritsFrom('TTree'):
            nev = obj.GetEntriesFast()
            new_file.cd()
            new_tree = obj.CopyTree("")
            # New branches
            stitchWeightBranch = new_tree.Branch('stitchWeight', stitchWeight, 'stitchWeight/F')
            eventWeightLumiBranch = new_tree.Branch('eventWeightLumi', eventWeightLumi, 'eventWeightLumi/F')

            # looping over events
            for event in range(0, obj.GetEntries()):
                if verbose and (event%10000==0 or event==nev-1): print ' = TTree:', obj.GetName(), 'events:', nev, '\t', int(100*float(event+1)/float(nev)), '%\r',
                #print '.',#*int(20*float(event)/float(nev)),#printProgressBar(event, nev)
                obj.GetEntry(event)
                # Initialize
                eventWeightLumi[0] = stitchWeight[0] = 1.
#                if obj.bTagWeight < 0.5 or obj.bTagWeight > 1.5:
#                    print obj.bTagWeight

                # Weights
                if isMC:
                    # MC stitching
                    if sample=='DYJetsToLL' or sample=='WJetsToLNu':
                        if obj.LheHT > 100.: stitchWeight[0] = 0.
                    eventWeightLumi[0] = LUMI*XS/totalEntries
                # Fill the branches
                stitchWeightBranch.Fill()
                eventWeightLumiBranch.Fill()

            new_file.cd()
            new_tree.Write()
            if verbose: print ' '


    new_file.Close()



jobs = []
for d in os.listdir(origin):
    if not '.root' in d: continue
    if len(filterset)>0 and not filterset in d: continue
    # if singlecore:
    print " -", d
    processFile(d, verboseout)
    # else:
    #     p = multiprocessing.Process(target=processFile, args=(d,verboseout,))
    #     jobs.append(p)
    #     p.start()
    #exit()
    
#print '\nDone.'

