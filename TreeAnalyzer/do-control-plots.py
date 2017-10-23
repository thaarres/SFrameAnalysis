import ROOT
from tdrstyle import *
setTDRStyle()
from  CMS_lumi import *
from time import sleep

ROOT.gROOT.SetBatch(True)

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.04*W_ref

CMS_lumi.lumi_13TeV = "35.9 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

lumi = "1"
dir = "/scratch/thaarres/VTopTagSF_MiniTuple/reweighted//HaddedOutput/"
cutL="1"
# cutL = "(jetAK8_softDrop_mass>30&&jetAK8_softDrop_mass<140)"
cutL = "(fabs(dr_ak8Lep)>1.5708&&fabs(dphi_ak8Et)>2.&&fabs(dphi_ak8Wlep)>2.&&Wlep_pt_2>200)"
vars = ["jetAK8_softDrop_mass", "jetAK8_pt", "jetAK8_tau21", "MET","jetAK8_tau32", "jetAK8_highestSubJetCSV","jetAK8_csv","lep_pt","Wlep_pt","Wlep_pt_2","dr_ak8Lep","fabs(dphi_ak8Et)","fabs(dphi_ak8Wlep)"]
bkgs = ["ST.root","VV.root","WJetsToLNu.root","QCD_reduced.root","TT.root"]
data = "SingleMuon.root"
fillcolor = [432,600,632,617,417,418]

def getPavetext():
  addInfo = ROOT.TPaveText(0.73010112,0.2566292,0.8202143,0.5523546,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  return addInfo
  	
def drawTH1(id,tree,var,cuts,bins,min,max,fillcolor,titlex = "",units = "",drawStyle = "HIST"):
	h = ROOT.TH1D("tmpTH1","",bins,min,max)
	h.Sumw2()
	# h.SetLineColor(1)
# 	h.SetLineWidth(2)
# 	h.SetFillStyle(1)
	h.SetFillColor(fillcolor)
	if units=="":
	    h.GetXaxis().SetTitle(titlex)
	else:
	    h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
	corrString='1'
	tree.Draw(var+">>tmpTH1","("+cuts+")*"+lumi+"*(eventWeightLumi)*("+corrString+")","goff")
	return h


def doCP(postfix=""):
	for var in vars:
		unit = "GeV"
		minx = 0.
		maxx = 200.
		binsx = 40
		
		if var.find("dr_ak8Lep")!=-1:
			minx = 0.
			maxx = 6.
		elif var.find("jetAK8_pt")!=-1:
			minx = 200.
			maxx = 600.
		elif var.find("lep_pt")!=-1:
			minx = 50.
			maxx = 400.
		elif var.find("tau")!=-1 or var.find("csv")!=-1 or var.find("CSV")!=-1:
			minx = 0.
			maxx = 1.
			unit = ""
			binsx = 20
		elif var.find("MET")!=-1:
			minx = 0.
			maxx = 500.
		elif var.find("dphi")!=-1:
			minx = 0.0
			maxx = 3.2
			unit = ""
		name = var
		canvas = ROOT.TCanvas(name+"_c1",name+"_c1",50,50,W,H)
		canvas.SetFillColor(0)
		canvas.SetBorderMode(0)
		canvas.SetFrameFillStyle(0)
		canvas.SetFrameBorderMode(0)
		canvas.SetLeftMargin( L/W )
		canvas.SetRightMargin( R/W )
		canvas.SetTopMargin( T/H )
		canvas.SetBottomMargin( B/H )
		canvas.SetTickx(0)
		canvas.SetTicky(0)
		
		ROOT.gStyle.SetOptStat(0)
		ROOT.gStyle.SetOptTitle(0)
		
		legend = ROOT.TLegend(0.62,0.7,0.92,0.9,"","brNDC")
		legend.SetBorderSize(0)
		legend.SetLineColor(1)
		legend.SetLineStyle(1)
		legend.SetLineWidth(1)
		legend.SetFillColor(0)
		legend.SetFillStyle(0)
		legend.SetTextFont(42)
		
		dataf = ROOT.TFile(dir+data, 'READ')
		treeD = dataf.Get('tree')
		datahist = drawTH1("data",treeD,var,cutL,binsx,minx,maxx,1,var.replace("_", " "),unit)
		datahist.SetName("data")
		files=[]
		hists=[]
		stack = ROOT.THStack("stack","")
		ttint			= 0
		totalMinoInt 	= 0
		for i,bg in enumerate(bkgs):
			print "opening file " ,dir+bg
			ref_file = ROOT.TFile(dir+bg, 'READ')
			files.append(ref_file)
		for i,file in enumerate(files):	
			tree = file.Get('tree')
			hist = drawTH1(str(i),tree,var,cutL,binsx,minx,maxx,fillcolor[i],var.replace("_", " "),unit)
			hist.SetName(name+str(i))
			ROOT.SetOwnership(hist,True)
			print hist.GetName()
			legend.AddEntry(hist,bkgs[i].split(".")[0],"F")
			hist.SetFillColor(fillcolor[i])
			if file.GetName().find("TT")!=-1: 
				ttint = hist.Integral()
				hist.Scale(0.750762237782)
				# hist.Scale(0.714578847292)
			else: totalMinoInt += hist.Integral()
		    
			
			hist.Draw()
			# sleep(4)
			stack.Add(hist)
			hists.append(hist)
		
		scale = datahist.Integral()-totalMinoInt
		print "DATA/MC" ,scale/ttint
		canvas.cd()
		datahist.GetYaxis().SetRangeUser(0, datahist.GetMaximum()*1.6);
		if var.find("jetAK8_pt")!=-1:
			datahist.GetYaxis().SetRangeUser(0.1, datahist.GetMaximum()*1000);
			canvas.SetLogy()
		datahist.Draw("ME")
		stack.Draw("HISTsame")
		datahist.Draw("MEsame")
		legend.Draw("SAME")
		CMS_lumi(canvas, iPeriod, iPos)
		canvas.Update()
		
		canvas.SaveAs(var.replace("fabs(","").replace(")","")+postfix+".png")
		del files
		del hists
		del dataf
		del treeD
		# sleep(200)
		# ["jetAK8_softDrop_mass","jetAK8_softDrop_mass_unCorr","jetAK8_gen_softDrop_mass","jetAK8_tau1",  "jetAK8_tau2",  "jetAK8_tau3",  "jetAK8_tau21","jetAK8_tau32",
		# "jetAK8_highestSubJetCSV", "jetAK8_pt","jetAK8_eta","jetAK8_gen_pt","jetAK8_csv","lep_pt","Wlep_pt","lep_eta","lep_phi"]

def doEff(cutstring,cuttext,postfix=""):
	iPeriod = 0
	unit = "GeV"
	minx = 30.
	maxx = 140.
	binsx = 14
	var = "jetAK8_softDrop_mass"
	canvas = ROOT.TCanvas("c1","c1",50,50,W,H)
	canvas.SetFillColor(0)
	canvas.SetBorderMode(0)
	canvas.SetFrameFillStyle(0)
	canvas.SetFrameBorderMode(0)
	canvas.SetLeftMargin( L/W )
	canvas.SetRightMargin( R/W )
	canvas.SetTopMargin( T/H )
	canvas.SetBottomMargin( B/H )
	canvas.SetTickx(0)
	canvas.SetTicky(0)
	
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)
	
	legend = ROOT.TLegend(0.62,0.7,0.92,0.9,"","brNDC")
	legend.SetBorderSize(0)
	legend.SetLineColor(1)
	legend.SetLineStyle(1)
	legend.SetLineWidth(1)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	
	files=[]
	hists=[]
	stack = ROOT.THStack("stack","")
	ttumINT 		= 0
	ttint			= 0
	totalMinoInt	= 0
	for i,bg in enumerate(bkgs):
		print "opening file " ,dir+bg
		ref_file = ROOT.TFile(dir+bg, 'READ')
		files.append(ref_file)
	for i,file in enumerate(files):	
		tree = file.Get('tree')
		
		if file.GetName().find("TT")!=-1: 
			cutU = cutstring+"&&(mergedVTruth==0)"
			histB = drawTH1(str(i),tree,var,cutU,binsx,minx,maxx,fillcolor[i+1],var.replace("_", " "),unit)
			histB.SetName(var+str(i)+"_unmerged")
			histB.SetFillColor(fillcolor[i+1])
			ttumINT = histB.Integral()
			stack.Add(histB)
			legend.AddEntry(histB,"Unmerged W (from tt)","F")
			cutS = cutstring+"&&(mergedVTruth==1)"
			histS = drawTH1(str(i),tree,var,cutS,binsx,minx,maxx,fillcolor[i],var.replace("_", " "),unit)
			histS.SetFillColor(fillcolor[i])
			histS.SetName(var+str(i)+"_merged")
			stack.Add(histS)
			ttINT = histS.Integral()
			legend.AddEntry(histS,"Real merged W (from tt)","F")
		else:	
			hist = drawTH1(str(i),tree,var,cutstring,binsx,minx,maxx,fillcolor[i],var.replace("_", " "),unit)
			hist.SetFillColor(fillcolor[i])
			legend.AddEntry(hist,bkgs[i].split(".")[0]+ " (all)","F")
			hist.SetName(var+str(i))
			stack.Add(hist)

		if file.GetName().find("TT")!=-1:
			print "TT integral already set"
			ttint = ttINT+ttumINT
		
		# 	hist.Scale(0.746992414682)
		# 	hist.Scale(0.89305440873)
		else: totalMinoInt += hist.Integral()
	    
		
		
		hists.append(hist)
	
	sb = ttint/totalMinoInt
	sb2 = ttint/(totalMinoInt+ttumINT)
	print "S/B = " ,sb
	canvas.cd()
	# stack.GetYaxis().SetRangeUser(0, stack.GetMaximum()*1.3);
	stack.Draw("HIST")
	
	legend.Draw("SAME")
	CMS_lumi(canvas, iPeriod, iPos)
	pt = getPavetext()
	pt.AddText(cuttext)
	pt.AddText("N_{sig} = %i"%ttint)
	pt.AddText("N_{bg} = %i"%totalMinoInt)
	pt.AddText("#frac{S}{B} = %.2f"%(sb2))
	pt.AddText("(w.o unmerged = %.2f)"%(sb))
	pt.Draw("same")	
	canvas.Update()
	# sleep(200)
	canvas.SaveAs("ttEff_"+cuttext.replace(" ","")+postfix+".png")
	
	
if __name__ == "__main__":
	doCP()

	purities = ["1","jetAK8_tau21>0.45","jetAK8_tau21<0.45"]
	label    = ["_all","_LP","_HP"]
	for i,purity in enumerate(purities):

		doEff(purity,"no cut",label[i])
		doEff("(fabs(dr_ak8Lep)>1.5708&&fabs(dphi_ak8Et)>2.&&fabs(dphi_ak8Wlep)>2&&%s)"%purity,"angular cuts",label[i])
		doEff("MET>150&&%s"%purity,"MET 150 GeV",label[i])
		doEff("MET>175&&%s"%purity,"MET 175 GeV",label[i])
		doEff("MET>200&&%s"%purity,"MET 200 GeV",label[i])