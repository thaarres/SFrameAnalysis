#include "include/LumiWeight.h"

#include <iostream>

//==============================================================================================
LumiWeight::LumiWeight( void ){

}

//==============================================================================================
LumiWeight::~LumiWeight( void ){

}

//==============================================================================================
double LumiWeight::getLumiWeight( TString sample ){


  if( sample.Contains( "WJetsToQQ_HT"                         ) ) return 95.14;
  if( sample.Contains( "ZJetsToQQ_HT600toInf"                 ) ) return 41.34;
  if( sample.Contains( "QCD_Pt_170to300_"                     ) ) return 117276.;
  if( sample.Contains( "QCD_Pt_300to470_"                     ) ) return 7823.;                   
  if( sample.Contains( "QCD_Pt_470to600_"                     ) ) return 648.2;                  
  if( sample.Contains( "QCD_Pt_600to800_"                     ) ) return 186.9;                 
  if( sample.Contains( "QCD_Pt_800to1000_"                    ) ) return 32.293;               
  if( sample.Contains( "QCD_Pt_1000to1400_"                   ) ) return 9.4183; 
  if( sample.Contains( "QCD_Pt_1400to1800_"                   ) ) return 0.84265; 
  if( sample.Contains( "QCD_Pt_1800to2400_"                   ) ) return 0.114943; 
  if( sample.Contains( "QCD_Pt_2400to3200_"                   ) ) return 0.006830; 
  if( sample.Contains( "QCD_Pt_3200toInf_"                    ) ) return 0.000165445; 
  if( sample.Contains( "QCD_HT100to200"                       ) ) return 2.785e07;
  if( sample.Contains( "QCD_HT200to300"                       ) ) return 1717000.;
  if( sample.Contains( "QCD_HT300to500"                       ) ) return 351300.0;
  if( sample.Contains( "QCD_HT500to700"                       ) ) return 31630.;
  if( sample.Contains( "QCD_HT700to1000"                      ) ) return 6802.;
  if( sample.Contains( "QCD_HT1000to1500"                     ) ) return 1206.0;
  if( sample.Contains( "QCD_HT1500to2000"                     ) ) return 120.4;
  if( sample.Contains( "QCD_HT2000toInf"                      ) ) return 25.25;
  if( sample.Contains( "QCD_Pt-15to7000" ) or sample.Contains( "QCD_Pt_15to7000" )) return  2.022100000e+09*60.5387252324;
  if( sample.Contains("BulkGrav") or sample.Contains("Qstar") or sample.Contains("Wprime") or sample.Contains("Zprime")) return 1.0;
  if( sample.Contains( "TT_TuneCUETP8M2T4"                    ) ) return  831.76      ;
  if( sample.Contains( "TT_TuneEE5C_13TeV-powheg-herwigpp"    ) ) return  831.76      ;
  if (sample.Contains("WJetsToLNu_HT-100To200"                ) ) return 1347*1.21    ;
  if (sample.Contains("WJetsToLNu_HT-200To400"                ) ) return 360*1.21     ;
  if (sample.Contains("WJetsToLNu_HT-400To600"                ) ) return 48.9*1.21    ;
  if (sample.Contains("WJetsToLNu_HT-600To800"                ) ) return 12.08*1.21   ;
  if (sample.Contains("WJetsToLNu_HT-800To1200"               ) ) return 5.26*1.21    ;
  if (sample.Contains("WJetsToLNu_HT-1200To2500"              ) ) return 1.33*1.21    ;
  if (sample.Contains("WJetsToLNu_HT-2500ToInf"               ) ) return 0.03089*1.21 ;
  if (sample.Contains("WW_TuneCUETP8M1"                       ) ) return 118.7        ;
  if (sample.Contains("WZ_TuneCUETP8M1"                       ) ) return 16.5         ;
  if (sample.Contains("ZZ_TuneCUETP8M1"                       ) ) return 47.13        ;
  if (sample.Contains("ST_s-channel_4f_leptonDecays"          ) ) return 47.13        ;
  if (sample.Contains("ST_t-channel_top_4f_leptonDecays"      ) ) return 136.02*0.322 ;
  if (sample.Contains("ST_t-channel_antitop_4f_leptonDecays"  ) ) return 80.95*0.322  ;
  if (sample.Contains("ST_tW_antitop_5f_inclusiveDecays"      ) ) return 35.6         ;
  if (sample.Contains("ST_tW_top_5f_inclusiveDecays_"         ) ) return 35.6         ;
  
  
  std::cout << "Cross section not defined for this sample!!" << std::endl;
  return 0;
  
}

















