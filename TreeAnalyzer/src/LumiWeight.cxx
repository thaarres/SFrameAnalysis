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
  
  xSec_ = 1.;
  if( sample.Contains( "WJetsToQQ_HT"          ) ) xSec_ = 95.14;
  if( sample.Contains( "ZJetsToQQ_HT600toInf"  ) ) xSec_ = 41.34;
  if( sample.Contains( "QCD_Pt_170to300_"      ) ) xSec_ = 117276.;
  if( sample.Contains( "QCD_Pt_300to470_"      ) ) xSec_ = 7823.;                   
  if( sample.Contains( "QCD_Pt_470to600_"      ) ) xSec_ = 648.2;                  
  if( sample.Contains( "QCD_Pt_600to800_"      ) ) xSec_ = 186.9;                 
  if( sample.Contains( "QCD_Pt_800to1000_"     ) ) xSec_ = 32.293;               
  if( sample.Contains( "QCD_Pt_1000to1400_"    ) ) xSec_ = 9.4183; 
  if( sample.Contains( "QCD_Pt_1400to1800_"    ) ) xSec_ = 0.84265; 
  if( sample.Contains( "QCD_Pt_1800to2400_"    ) ) xSec_ = 0.114943; 
  if( sample.Contains( "QCD_Pt_2400to3200_"    ) ) xSec_ = 0.006830; 
  if( sample.Contains( "QCD_Pt_3200toInf_"     ) ) xSec_ = 0.000165445; 
  if( sample.Contains( "QCD_HT100to200"        ) ) xSec_ = 2.785e07;
  if( sample.Contains( "QCD_HT200to300"        ) ) xSec_ = 1717000.;
  if( sample.Contains( "QCD_HT300to500"        ) ) xSec_ = 351300.0;
  if( sample.Contains( "QCD_HT500to700"        ) ) xSec_ = 31630.;
  if( sample.Contains( "QCD_HT700to1000"       ) ) xSec_ = 6802.;
  if( sample.Contains( "QCD_HT1000to1500"      ) ) xSec_ = 1206.0;
  if( sample.Contains( "QCD_HT1500to2000"      ) ) xSec_ = 120.4;
  if( sample.Contains( "QCD_HT2000toInf"       ) ) xSec_ = 25.25;
  if( sample.Contains( "QCD_Pt-15to7000" ) or sample.Contains( "QCD_Pt_15to7000" ))xSec_ =  2.022100000e+09;
  if( sample.Contains("BulkGrav") or sample.Contains("Qstar") or sample.Contains("Wprime") or sample.Contains("Zprime")) xSec_ = 1.0;
  else{
    std::cout << "Cross section not defined for this sample!!" << std::endl;
      return -999.;
  }
  return xSec_;

      
}
