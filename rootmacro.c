#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TSystem.h"


void rootmacro() { 
  //--------------Open the two reqd root files-----------------------------------------------------------------------------
  TFile *f_SM = new TFile("/afs/cern.ch/user/n/ngubbi/aQGC/aQGC_rivetAnalysis-SM-newv5-rivetroutine.root");//SM root file
  if (!f_SM->IsOpen())
    {
     cout<<"Error opening SM root file\n"<<endl;                                                                                                                                                 
     return;
    }
  TFile *f_aQGC = new TFile("/afs/cern.ch/user/n/ngubbi/aQGC/aQGC_rivetAnalysis-M5-bestlimits-newv5rivetroutine.root");//aQGC root file
  if (!f_aQGC->IsOpen())
    {
      cout<<"Error opening aQGC root file\n"<<endl;                                                                                                                                                      
      return;
    }
 
  string vars[34] ={"Ntotev","Njets","Nmu","Nel","Nlep","HT","MET","leadingelectron_pt","leadingmuon_pt", "jet1_pt","jet2_pt","leadingjet_pt","Zlep1_pt","Zlep2_pt","W_lep_pt","Z_m_l_l","WWZ_m","Z_pt_l_l","Whad_pt_q_q","Wlep_pt_l_nu","Z_dPhi_l_l","Whad_dPhi_q_q","Wlep_dPhi_l_MET","Z_dEta_l_l","Whad_dEta_q_q","Whad_m_q_q","Wlep_m_l_nu","Wlep_dEta_l_nu","Z_dR_l_l","Wlep_dR_l_nu","Whad_dR_q_q","singlejet_mass","Jet1_Mass","Jet2_Mass"};

    for(int i=0; i<34;++i){
    // Define the Canvas                                                                                                                                                                             
   TCanvas *c = new TCanvas("c", "canvas", 800, 600);
   gPad->SetRightMargin(0.05);
   gPad->SetTopMargin(0.05);
   gStyle->SetOptStat(0);
   //Define the variables as a String
   std::cout << "Retrieving " << vars[i] << std::endl;                                                                                                                                                                       
  // lets add the first histogram                                                                                                                                                                          
   TH1D* histo1 = (TH1D*)f_SM->Get( ("WWZ/"+vars[i]).c_str() );
   if(histo1==NULL){std::cout << "ERROR: null pointer for histo1 " << std::endl;}
  // lets add the second histogram                                                                                                                                                                         
   TH1D* histo2 = (TH1D*)f_aQGC->Get( ("WWZ/"+vars[i]).c_str() );
   if(histo2==NULL){std::cout << "ERROR: null pointer for histo2 " << std::endl;}


  //Processing histograms-  add overflow, underflow, Normalization and Rebining

   //Overflow and Underflow -Histo1
   int N = histo1->GetNbinsX();
   histo1->SetBinContent(N, histo1->GetBinContent(N)+histo1->GetBinContent(N+1));//Overflow
   histo1->SetBinContent(1, histo1->GetBinContent(1)+histo1->GetBinContent(0));//Underflow
   //Overflow and Underflow -Histo2                                                                                                                                                                        
   N = histo2->GetNbinsX();
   histo2->SetBinContent(N, histo2->GetBinContent(N)+histo2->GetBinContent(N+1));//Overflow                                                                                                                
   histo2->SetBinContent(1, histo2->GetBinContent(1)+histo2->GetBinContent(0));//Underflow 
   //Normalizing histograms to 1         
   //Want the cut flow histogram to be showing absolute number of events- so removing the feature to normalize the histograms to 1                                       
   if(vars[i]=="Ntotev" ){
     histo1->Scale();//In order to obtain the cut flow histo in terms of absolute number of events after respective cut!
     histo2->Scale();
   }
   else{         
     histo1->Scale(1./histo1->Integral());
     histo2->Scale(1./histo2->Integral());
   } //- plot distributions must be normalised to 1 (e.g. scaled to unity by doing h->Scale(1./h->Integral())) to compare shapes of the two samples                    
  
   //Rebinning                                                                                                                                                               
   int rebin=0;
   if(vars[i]=="HT"||vars[i]=="MET"||vars[i]=="leadingelectron_pt"||vars[i]=="leadingmuon_pt"||vars[i]=="leadingjet_pt"||vars[i]=="WWZ_m"||vars[i]=="jet1_pt"||vars[i]=="jet2_pt"||vars[i]=="Z_pt_l,l"||vars[i]=="Whad_pt_q,q"||vars[i]=="Whad_dPhi_q_q"||vars[i]=="Wlep_dPhi_l_nu"||vars[i]=="Z_dEta_l_l"||vars[i]=="Whad_dEta_q_q"||vars[i]=="Wlep_dEta_l_nu"||vars[i]=="Z_dR_l_l"||vars[i]=="Wlep_dR_l_nu"||vars[i]=="Whad_dR_q_q"){
     rebin=4;//}
     histo1->Rebin(rebin);
     histo2->Rebin(rebin);
    }

  //Cosmetics
  // Histo1 settings
  histo1->SetTitle("");                                                                                                                                                      
  histo1->SetFillStyle(0);
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2);
  //  histo1->SetMarkerColor();
  histo1->SetMarkerStyle(8);
  histo1->SetMarkerSize(0.5);
  //Setting Minimum and Maximum -Histo1
  histo1->SetMinimum(0);
  if(vars[i]=="Z_dEta_l_l"||vars[i]=="Whad_dEta_q_q"||vars[i]=="Wlep_dEta_l_nu"||vars[i]=="Z_dR_l_l"||vars[i]=="Wlep_dR_l_nu"||vars[i]=="Whad_dR_q_q"||vars[i]=="Z_dR_l_l"){
    histo1->SetMaximum(0.4);
  }
  if(vars[i]=="Whad_dPhi_q_q"||vars[i]=="Wlep_dPhi_l_nu"||vars[i]=="Z_dPhi_l_l"){
    histo1->SetMaximum(0.1);
  }

  //histo1 -Y axis                                                                                                                                                          
  //histo1->GetYaxis()->SetTitle("No. of events");
  //histo1->GetYaxis()->SetTitleSize(0.20);
  //histo1->GetYaxis()->SetTitleFont(43);
  histo1->GetYaxis()->SetTitleOffset(1.4);
  //histo1->GetYaxis()->SetLabelSize(0.03);

  //Histo1- Xaxis
  //histo1->GetXaxis()->SetTitle("Mass [GeV]");//I commented this since i dont want to forcefully add one label to all histos and the label must be picked as given in rivet routine while booking the histogram 
  //histo1->GetXaxis()->SetTitleSize(0.2);//Set 0 to remove the label for axis
  //histo1->GetXaxis()->SetTitleFont(43);
  histo1->GetXaxis()->SetTitleOffset(1.4);
  //histo1->GetXaxis()->SetLabelSize(0.03);

  // histo2 settings                                                                                                                                                                                   
  histo2->SetTitle("");                                                                                                                                                                                    
  histo2->SetFillStyle(0);
  histo2->SetLineColor(kBlue);
  histo2->SetLineWidth(2);
  // histo2->SetMarkerColor();
  histo2->SetMarkerStyle(8);
  histo2->SetMarkerSize(0.5);
  //Setting Minimum and Maximum -Histo2                                                                                                                                          
  histo2->SetMinimum(0);
  if(vars[i]=="Z_dEta_l_l"||vars[i]=="Whad_dEta_q_q"||vars[i]=="Wlep_dEta_l_nu"||vars[i]=="Z_dR_l_l"||vars[i]=="Wlep_dR_l_nu"||vars[i]=="Whad_dR_q_q"||vars[i]=="Z_dR_l_l"){
    histo2->SetMaximum(0.4);
  }
  if(vars[i]=="Whad_dPhi_q_q"||vars[i]=="Wlep_dPhi_l_nu"||vars[i]=="Z_dPhi_l_l"){
    histo2->SetMaximum(0.1);
  }

  //histo2 -Y axis                                                                                                                                                           
  //histo2->GetYaxis()->SetTitle("No. of events");                                                                                                                         
  //histo2->GetYaxis()->SetTitleSize(0.20);
  //histo2->GetYaxis()->SetTitleFont(43);
  histo2->GetYaxis()->SetTitleOffset(1.4);
  //histo2->GetYaxis()->SetLabelSize(0.03);

  //Histo2- Xaxis                                                                                                                                                                                          
  //histo1->GetXaxis()->SetTitle("Mass [GeV]");                                                                                                                                                            
  //histo2->GetXaxis()->SetTitleSize(0.1);//Set 0 to remove the label for axis                                                                                                                               
  //histo2->GetXaxis()->SetTitleFont(43);
  histo2->GetXaxis()->SetTitleOffset(1.4);
  //histo2->GetXaxis()->SetLabelSize(0.03);

  //Drawing both histos- including "HITSOE" in the way carlo wanted
  histo1->Draw("HISTOE");
  histo2->Draw("HISTOEsame");

  //Legend
  TLegend *leg = new TLegend(0.85,0.85,0.95,0.95);
  leg->SetBorderSize(0);  //0 means no border
  leg->SetNColumns(1);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(histo1,"SM","l");
  leg->AddEntry(histo2,"M5","l");
  leg->Draw("Same");
  c->Update();
  c->Draw();

  //c->SaveAs( ( vars[i]+".pdf".c_str() ) );                                                                                                                                                                
  c->SaveAs( (vars[i]+".pdf").c_str() );

  } //end of for loop
}//end of rootmacro()


