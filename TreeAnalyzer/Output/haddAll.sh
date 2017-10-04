#!/bin/sh

for ch in el mu
do
	test -e VVanalysis.WW_TuneCUETP8M1_13TeV-pythia8*_$ch*			&& echo WW.$ch.root 										VVanalysis.*WW_TuneCUETP8M1_13TeV-pythia8*_$ch*					&
	test -e VVanalysis.WZ_TuneCUETP8M1_13TeV-pythia8*_$ch*			&& echo WZ.$ch.root 										VVanalysis.*WZ_TuneCUETP8M1_13TeV-pythia8*_$ch*					&
	test -e VVanalysis.ZZ_TuneCUETP8M1_13TeV-pythia8*_$ch*			&& echo ZZ.$ch.root 										VVanalysis.*ZZ_TuneCUETP8M1_13TeV-pythia8*_$ch*					&
	test -e VVanalysis.TT_TuneCUETP8M2T4_13TeV-powheg-pythia8*_$ch* && echo TT.$ch.root 										VVanalysis.*TT_TuneCUETP8M2T4_13TeV-powheg-pythia8*_$ch*		&
	test -e VVanalysis.WJets_100To200*_$ch*    						&& echo WJets_100To200.$ch.root   							VVanalysis.*WJets_100To200*_$ch*    							&
	test -e VVanalysis.WJets_200To400*_$ch*    						&& echo WJets_200To400.$ch.root   							VVanalysis.*WJets_200To400*_$ch*    							&
	test -e VVanalysis.WJets_400To600*_$ch*    						&& echo WJets_400To600.$ch.root   							VVanalysis.*WJets_400To600*_$ch*    							&
	test -e VVanalysis.WJets_600To800*_$ch*    						&& echo WJets_600To800.$ch.root   							VVanalysis.*WJets_600To800*_$ch*    							&
	test -e VVanalysis.WJets_800To1200*_$ch*   						&& echo WJets_800To1200.$ch.root  							VVanalysis.*WJets_800To1200*_$ch*   							&
	test -e VVanalysis.WJets_1200To2500*_$ch*  						&& echo WJets_1200To2500.$ch.root 							VVanalysis.*WJets_1200To2500*_$ch*  							&
	test -e VVanalysis.WJets_2500ToInf*_$ch*   						&& echo WJets_2500ToInf.$ch.root  							VVanalysis.*WJets_2500ToInf*_$ch*   							&
	test -e VVanalysis.ST_s-channel_4f_leptonDecays*_$ch*           && echo ST_s-channel_4f_leptonDecays.$ch.root            	VVanalysis.*ST_s-channel_4f_leptonDecays*_$ch*           		&
	test -e VVanalysis.ST_t-channel_antitop_4f_inclusiveDecays*_$ch*&& echo ST_t-channel_antitop_4f_inclusiveDecays.$ch.root 	VVanalysis.*ST_t-channel_antitop_4f_inclusiveDecays*_$ch*		&
	test -e VVanalysis.ST_t-channel_top_4f_inclusiveDecays*_$ch*    && echo ST_t-channel_top_4f_inclusiveDecays.$ch.root     	VVanalysis.*ST_t-channel_top_4f_inclusiveDecays*_$ch*    		&
	test -e VVanalysis.ST_tW_antitop_5f_inclusiveDecays*_$ch*       && echo ST_tW_antitop_5f_inclusiveDecays.$ch.root        	VVanalysis.*ST_tW_antitop_5f_inclusiveDecays*_$ch*       		&
	test -e VVanalysis.ST_tW_top_5f_inclusiveDecays*_$ch*           && echo ST_tW_top_5f_inclusiveDecays.$ch.root            	VVanalysis.*ST_tW_top_5f_inclusiveDecays*_$ch*           		&
	
	
	
	# rm VVanalysis.*
done



