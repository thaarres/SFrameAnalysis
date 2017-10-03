#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re

import thread
import subprocess
import mmap
import itertools

def waitForBatchJobs(runningJobs, listOfJobs, userName):
  print "waiting for %d job(s) in the queue" %(len(runningJobs))
  while not len(runningJobs)==0:
    time.sleep(30.)
    queryString="qstat -u %s | grep %s | awk {\'print $10\'}" %(userName, userName)
    lock=thread.allocate_lock()
    lock.acquire()
    compileProcess=subprocess.Popen(queryString, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    compileProcess.wait()
    lock.release()
    jobList = compileProcess.stdout.read()
    jobList = jobList.split("\n")
    for j in runningJobs:
      jobId=j[0]
      if (jobId not in jobList):
        for l in listOfJobs:
          if j[1]==(l[2]+l[3]):
            l[4]="."
        runningJobs.remove(j)
        print "waiting for %d job(s) in the queue" %(len(runningJobs))
		

def createJobs(i,f, outfolder, outname,channel='el',isData='false'):
	template=open("config/submitJobs.xml", 'r').read()
	template=template.replace('OUTPUT', ('<Cycle Name="VVanalysis" OutputDirectory="%s/" PostFix="" TargetLumi="1.0">')%outfolder)
	template=template.replace('INPUTHEADER', ('<InputData Lumi="0.0" NEventsMax="-1" NEventsSkip="0" Type="%s" Version="%s">')%(outname,i) )

        
	template=template.replace('INFILE', ('<In FileName="dcap://t3se01.psi.ch:22125/%s" Lumi="1.0" />')%f)
	template=template.replace('CHANNEL', ('<Item Name="Channel" Value="%s" />')%channel)
	template=template.replace('ISDATA', ('<Item Name="IsData" Value="%s"/>')%isData)
	try: os.stat("xmls") 
	except: os.mkdir("xmls")
	xml = "xmls/"+outname+""+str(i)+".xml"
	outFile=open(xml, 'w')
	outFile.write(template)
	outFile.close()
	
	cmd = 'sframe_main '+xml +'\n'
	print cmd
	jobs.write(cmd)
	return outfolder+'/VVanalysis.'+outname+'.'+str(i)+'.root'

def submitJobs(jobList, outfolder):
    print 'Reading joblist'
    jobListName = jobList
    print jobList
#    subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobListName)
    subCmd = 'qsub -o %s/logs/ %s %s' %(outfolder,"psibatch_runner.sh",jobListName)
    print 'Going to submit jobs with', subCmd
    os.system(subCmd)

    return 1


def split_seq(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))


def mapcount(filename):
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines


if __name__ == "__main__":
  
    print "Pass 2 arguments. First: Wjets, TT, ST, VV, SingleMu, SingleEl, Second: el, mu"
    channel = sys.argv[2]
    patterns = []
    sample = sys.argv[1] 
    
    if sample.find("TT")!=-1: patterns= ["TT_TuneCUETP8M2T4_13TeV-powheg-pythia8"]
    elif sample.find("Wjets")!=-1: patterns = ["WJetsToLNu_HT-100To200","WJetsToLNu_HT-200To400","WJetsToLNu_HT-400To600","WJetsToLNu_HT-600To800","WJetsToLNu_HT-800To1200","WJetsToLNu_HT-1200To2500","WJetsToLNu_HT-2500ToInf"]
    elif sample.find("VV")!=-1: patterns = ["WW_TuneCUETP8M1","WZ_TuneCUETP8M1","ZZ_TuneCUETP8M1"]
    elif sample.find("ST")!=-1: patterns = ["ST_s-channel_4f_leptonDecays","ST_t-channel_antitop_4f_inclusiveDecays","ST_t-channel_top_4f_inclusiveDecays","ST_tW_antitop_5f_inclusiveDecays","ST_tW_top_5f_inclusiveDecays"]
    elif sample.find("SingleMu")!=-1: patterns= ["SingleMuon"]
    elif sample.find("SingleEl")!=-1: patterns= ["SingleElectron"]

    else:
    	print "Please pass either: Wjets, TT, ST, VV, SingleMu, SingleEl"
    	sys.exit()
	

    
    outfolder = "Output"

    try: os.stat(outfolder) 
    except: os.mkdir(outfolder)
	
    try: os.stat(outfolder+'/logs') 
    except: os.mkdir(outfolder+'/logs')
   
    location = '/pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Summer16/Ntuple_80_20170206/'
    isData = 'false'
    if 'SingleMu' in sample or 'SingleEl' in sample:
      print 'Running data, I will pick muon channel for SingleMu dataset etc.'
      location = '/pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Moriond17/'
      if 'SingleMu' in sample:
        channel = 'mu'
        isData = 'true'
      if 'SingleEl' in sample:
        channel = 'el'
        isData = 'true'



    jobLists = []
    print patterns
    for pattern in patterns:
      filelist = glob.glob(location+'/'+pattern+'*/*/*/*/*.root')
      jobList = 'joblist_'+pattern+'_'+channel+'.txt'
      jobs = open(jobList, 'w')
      outs = []
      for i,f in enumerate(filelist):
        print f
        outname = f.split("/")[10]
        print outname
        fout = createJobs(i,f, outfolder, outname,channel,isData)
        outs.append(fout)
      print outs
	
      jobLists.append(jobList)
      jobs.close()

    print jobLists
    print "Have ", len(jobLists), " joblists to submit with "
    for jobList in jobLists:
      print "Name: ", jobList, " nJobs: ", mapcount(jobList)

    submit = raw_input("Do you also want to submit the jobs to the batch system? [y/n] ")
    if submit == 'y' or submit=='Y':

      for jobList in jobLists:
        lock=thread.allocate_lock()
        lock.acquire()
        print 'Reading joblist'
        jobListName = jobList
        print jobList
        subCmd = 'qsub -t 1-%s -o %s/logs/ psibatch_runner.sh %s' %(i,outfolder,jobListName)
        print 'Going to submit jobs with', subCmd
	    	
        subProcess=subprocess.Popen( subCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        subProcessStatus=subProcess.poll()
        if subProcessStatus==1:
          print "job '%s' failed" %(subCmd, subProcessStatus)
          lock.release()
          sys.exit()
        else:
          lock.release()
          print "job finished"
          userName=os.environ['USER']
          # mergeCmd='hadd -f %s.root %s.root && rm -rf %s.root'  %(mergeFileBaseName,  fileToMerge,  fileToMerge)
          #         print "mergeCmd is %s " %mergeCmd
            
            
      else:
        print "Not submitting jobs"


    
