// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

//to compile this rivet file : 1. setupATLAS 2. asetup 20.20.8.3,here 3.rivet-buildplugin RivetWWZ-Carlo.so WWZ-Carlo.cc

namespace Rivet {

  double lep_pt_min = 10.0;
  double lep_maxabseta = 2.5;

  double jet_pt_min = 25.0;
  double jet_maxabseta = 5.0;

  double MZ_PDG = 91.1876;
  double MW_PDG = 80.385;
  double Pi = 3.14159265;
 
  //structure to be used for reconstruction of a X-> a, b decay
    struct reco_boson
  {
    FourMomentum p4;
    int flav;
    unsigned int a_idx;
    unsigned int b_idx;
    double mass_diff;
    reco_boson(FourMomentum p, int f, unsigned int a, unsigned int b, double m){p4=p, flav=f, a_idx=a, b_idx=b, mass_diff = m;}
  };


  class WWZ : public Analysis {
  public:
    int jet_counter=0,lep_counter=0,tau_neutrino_counter=0,neutrino_event_counter=0;
    long counter=0, sel_counter=0; 

    /// Constructor
    WWZ() : Analysis("WWZ") {    }

    /// @name Analysis methods
    /// Book histograms and initialise projections before the run
    void init() {
      
      Cut eta_full = (Cuts::abseta < 5.0) & (Cuts::pT >= 1.0*MeV);
      Cut lep_cuts = (Cuts::abseta < lep_maxabseta) & (Cuts::pT >= lep_pt_min*GeV);

      FinalState fs(eta_full);
  
      //photons for dressing
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // *** LEPTONS (MUONS+ELECTRONS) altogether ***
      IdentifiedFinalState lep_id(fs);
      lep_id.acceptIdPair(PID::ELECTRON);
      lep_id.acceptIdPair(PID::MUON);
      PromptFinalState leptons(lep_id);
      leptons.acceptTauDecays(true);
      declare(leptons, "leptons");
      //leptons meant for analysis
      DressedLeptons dressedleptons(photons, leptons, 0.1, lep_cuts, true, true);
      declare(dressedleptons, "dressedleptons");
      //leptons for Overlap Removal
      DressedLeptons ewdressedleptons(photons, leptons, 0.1, eta_full, true, true);
      declare(ewdressedleptons, "ewdressedleptons");

      // *** NEUTRINOS ***
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      addProjection(neutrinos, "neutrinos");

      // *** JET CLUSTERING ***
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedleptons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(false);
      declare(jets, "all_jets");

      // Book histograms   
      //Cut flow and Number of leptons
      _h_Ntotev = bookHisto1D("Ntotev",3,-0.5,2.5,"Cut Flow","Cut Flow","No. of events");
      _h_Njets  = bookHisto1D("Njets",12,0, 12, "Njets", "Number of jets", "(1/N) n. of events");
      _h_Nmu    = bookHisto1D("Nmu",4,-0.5,3.5, "Nmu", "Number of muons", "(1/N) n. of events");
      _h_Nel    = bookHisto1D("Nel",4,-0.5,3.5, "Nel", "Number of electrons", "(1/N) n. of events");
      _h_Nlep   = bookHisto1D("Nlep",6,-0.5,5.5, "Nel", "Number of leptons", "(1/N) n. of events");
      //HT and MET
      _h_HT      = bookHisto1D("HT",200, 0, 5000, "HT","H_{T} [GeV]","(1/N) n. of events");
      _h_MET     = bookHisto1D("MET",200, 0, 2000, "MET","E_{T}^{miss} [GeV]","(1/N) n. of events");
      //pT for Leading jet,mu,e
      _h_leadingelectron_pt = bookHisto1D("leadingelectron_pt",200, 0, 2000, "leadingelectron_pt","Leading electron p_{T} [GeV]","(1/N) n. of events");
      _h_leadingmuon_pt     = bookHisto1D("leadingmuon_pt",200, 0, 2000, "leadingmuon_pt","Leading muon p_{T}  [GeV]","(1/N) n. of events");
      _h_leadingjet_pt      = bookHisto1D("leadingjet_pt",200, 0, 2000, "leadingjet_pt","Leading jet p_{T}[GeV]","(1/N) n. of events");
      //pT of the individual jets and individual leptons from Z and other remaining lepton from W
      _h_jet1_pt   = bookHisto1D("jet1_pt",150, 0, 2000, "jet1_pt","jet1 p_{T}[GeV]","(1/N) n. of events");
      _h_jet2_pt   = bookHisto1D("jet2_pt",150, 0, 2000, "jet2_pt","jet2 p_{T}[GeV]","(1/N) n. of events");
      _h_Z_lep2_pt = bookHisto1D("Zlep2_pt",50, 0,2000, "Z_lep2_pt","Z_lep2 p_{T}[GeV]","(1/N) n. of events");
      _h_Z_lep1_pt = bookHisto1D("Zlep1_pt",50, 0,2000, "Z_lep1_pt","Z_lep1 p_{T}[GeV]","(1/N) n. of events");
      _h_W_lep_pt  = bookHisto1D("W_lep_pt",50, 0,2000, "W_lep_pt","W_lep p_{T}[GeV]","(1/N) n. of events");//pT of the other remaining lepton -to avoid confusion in the name with _h_Wlep_pt which the pT fo the W boson
      //Reconstructed Invariant masses
      _h_Z_m       = bookHisto1D("Z_m_l_l" ,20, 0, 140, "Z_m_l_l","Reconstructed Z(l,l) mass [GeV]","(1/N) n. of events");
      _h_Whad_m    = bookHisto1D("Whad_m_q_q" ,20, 0, 140, "Whad_m_q_q","Reconstructed W(q,q) mass [GeV]","(1/N) n. of events");
      _h_Wlep_m    = bookHisto1D("Wlep_m_l_nu",20, 0, 140, "Wlep_m_l_nu","Reconstructed W(l,nu) mass [GeV]","(1/N) n. of events");
      _h_InvMass_WWZ = bookHisto1D("WWZ_m",200, 0, 5000, "WWZ_m","Reconstructed WWZ mass [GeV]","(1/N) n. of events");
      _h_singlejet_mass = bookHisto1D("singlejet_mass" ,20, 0, 140, "singlejet_mass","SingleJet mass [GeV]","(1/N) n. of events");
      _h_Jet1_Mass = bookHisto1D("Jet1_Mass" ,20, 0, 100, "Jet1_Mass","Jet1 Mass [GeV]","(1/N) n. of events");
      _h_Jet2_Mass = bookHisto1D("Jet2_Mass" ,20, 0, 100, "Jet2_Mass","Jet2 Mass [GeV]","(1/N) n. of events");
      
      //pT of the Bosons
      _h_Z_pt      = bookHisto1D("Z_pt_l_l"    ,20, 0, 2000, "Z_pt_l_l"   ,"Reconstructed Z(l,l) p_{T} [GeV]","(1/N) n. of events");
      _h_Whad_pt   = bookHisto1D("Whad_pt_q_q" ,20, 0, 2000, "Whad_pt_q_q","Reconstructed W(q,q) p_{T} [GeV]","(1/N) n. of events");
      _h_Wlep_pt   = bookHisto1D("Wlep_pt_l_nu",20, 0, 2000, "Wlep_pt_l_nu","Reconstructed W(l,nu) p_{T} [GeV]","(1/N) n. of events");
      //Delta Phi
      _h_Z_dPhi    = bookHisto1D("Z_dPhi_l_l"     ,20, 0, 3.14,"Z_dPhi_l_l"     ,"#Delta#Phi(l,l) " ,"(1/N) n. of events");
      _h_Whad_dPhi = bookHisto1D("Whad_dPhi_q_q"  ,90, 0, 3.14,"Whad_dPhi_q_q"  ,"#Delta#Phi(j,j)"  ,"(1/N) n. of events");
      _h_Wlep_dPhi = bookHisto1D("Wlep_dPhi_l_MET",20, 0, 3.14,"Wlep_dPhi_l_MET","#Delta#Phi(l,MET)","(1/N) n. of events");
      //Delta Eta
      _h_Z_dEta    = bookHisto1D("Z_dEta_l_l"    ,45, 0,5,"Z_dEta_l_l"    ,"#Delta#eta(l,l)","(1/N) n. of events");
      _h_Whad_dEta = bookHisto1D("Whad_dEta_q_q" ,45, 0,5,"Whad_dEta_q_q" ,"#Delta#eta(j,j)","(1/N) n. of events");
      _h_Wlep_dEta = bookHisto1D("Wlep_dEta_l_nu",45, 0,5,"Wlep_dEta_l_nu","#Delta#eta(l,#nu)","(1/N) n. of events");
      //Delta R
      _h_Z_dR      = bookHisto1D("Z_dR_l_l"    ,45,0,5,  "Z_dR_l_l"    , "#DeltaR(l,l)"   ,"(1/N) n. of events" );
      _h_Wlep_dR   = bookHisto1D("Wlep_dR_l_nu",45,0,5,  "Wlep_dR_l_nu", "#DeltaR(l,#nu)","(1/N) n. of events" );
      _h_Whad_dR   = bookHisto1D("Whad_dR_q_q" ,45,0,5,  "Whad_dR_q_q" , "#DeltaR(j,j)","(1/N) n. of events" );    
      //2dim Histogram
      _h_HTvsTot_m    = bookHisto2D("HTvsTot_m"   ,20,0,2000,20,0,2000, "HTvsInvMass"  , "H_{T} [GeV]", "WWZ reconstructed mass [GeV]" );
    }
    void analyze(const Event& event) {
      ++counter;
      _h_Ntotev->fill(0);
      const double weight = event.weight();

      // *** CONTAINERS RETRIEVAL ***
      // Get the selected electrons and muons
      std::vector<DressedLepton> leps = sortByPt(applyProjection<DressedLeptons>(event, "dressedleptons").dressedLeptons());
      //Get the jets applying the baseline cut on pT and eta
      const Jets& jets = applyProjection<FastJets>(event, "all_jets").jetsByPt(Cuts::pT > jet_pt_min*GeV && Cuts::abseta < jet_maxabseta);
      //Get neutrinos and compute MET
      const Particles& nu = applyProjection<PromptFinalState>(event, "neutrinos").particlesByPt();
      FourMomentum p_miss;
      foreach(const Particle& p, nu) p_miss+=p.momentum();

      // DEFINITIONS OF VARIABLES ***                                                                                                                                                                            // GLOBAL event variables                                                                                                                                                                 
      double HT = 0., Nmu = 0., Nel = 0, InvMass_WWZ = 0.;

      //*** EVENT SELECTION ***
      //Jet selection -atleast 2jets
      if(jets.size() < 2){
	//To investigate a peak in Dijet Invariant mass histogram around 40GeV as to what happens to a single jet
	//cout<<"less than 2jets;"<<"singlejet mass= "<<singlejet_mass<<endl;
	if(!jets.empty()){_h_singlejet_mass->fill(jets[0].momentum().mass(),weight);}
	vetoEvent;
      }
      ++sel_counter;
      _h_Ntotev->fill(1);

      //Atleast 3lepton selection                                                                                                                                                                          
      if(leps.size() <3){vetoEvent;}
      ++lep_counter;
      _h_Ntotev->fill(2);

      //Calculate HT first and then place the cut
      // HT = scalar sum of the jets' and leptons' pt                                                                                                                                                      
      /* for(unsigned int l=0; l<leps.size(); ++l){
        HT+=leps.at(l).pt();
        if(leps.at(l).abspid() == 13) {Nmu+=1;}
        if(leps.at(l).abspid() == 11) {Nel+=1;}
      }
      for(unsigned int j=0; j<jets.size(); ++j){HT+=jets.at(j).pt();}
      // cout<<"HT="<<HT<<endl;
      // _h_HT->fill(HT, weight);
      //Placing the HT cut
      if(HT>3000*GeV){vetoEvent;}
      // cout<<"HT after cut="<<HT<<endl;
      _h_HT->fill(HT, weight); 
      */
      //*** EVENT RECONSTRUCTION ***
      // 1. best Z->ll candidate, using invariant mass criterion
      // 2. best W->qq candidate, using invariant mass criterion
      // 3. best W->lv candidate, using remaining lepton

      // Z->ll (l = e, mu)
      std::vector<reco_boson> Zll;
      if(!Zll.empty()){Zll.clear();}
      //build all candidates
      for(unsigned int l1=0; l1<leps.size(); ++l1){
        for(unsigned int l2=l1+1; l2<leps.size(); ++l2){
          if(leps.at(l1).abspid() == leps.at(l2).abspid() && leps.at(l1).charge() != leps.at(l2).charge()){
            FourMomentum candidate = leps.at(l1).momentum() + leps.at(l2).momentum();
            double mass_diff = fabs(candidate.mass() - MZ_PDG);
            Zll.push_back(reco_boson(candidate, leps.at(l1).abspid(), l1, l2, mass_diff));
          }  
        }
      }
      //if no candidates veto the event
      if(Zll.empty()){vetoEvent;}	
      //now sort the candidates according to the smallest invariant mass difference
      std::sort(Zll.begin(), Zll.end(), mass_diff_sorting);
      //Best candidate is the first!
      FourMomentum Leptonic_Z = Zll.at(0).p4;
      unsigned int Z_lep1_idx = Zll.at(0).a_idx;
      unsigned int Z_lep2_idx = Zll.at(0).b_idx;
      unsigned int W_lep_idx = 99;
      //pT of the lep1 and lep2 matched to the leptonic Z and pT of the lepton matched to the leptonic W
      FourMomentum Z_lep1 = leps[Z_lep1_idx];
      FourMomentum Z_lep2 = leps[Z_lep2_idx];
      _h_Z_lep1_pt->fill(Z_lep1.pt(),weight);
      _h_Z_lep2_pt->fill(Z_lep2.pt(),weight);
      // W->lv
      FourMomentum Leptonic_W;
      // only one candidate as there is only one lepton left...
      for(unsigned int l1=0; l1<leps.size(); ++l1){
        if(l1 == Z_lep1_idx || l1 == Z_lep2_idx){continue;}
        //remaining lepton + hardest neutrino (hardest in pT!)
        Leptonic_W = leps.at(l1).momentum() + nu[0].momentum();
        W_lep_idx = l1;
        break;
      }
      FourMomentum W_lep = leps[W_lep_idx];
      _h_W_lep_pt->fill(W_lep.pt(),weight);

      // W->qq 
      std::vector<reco_boson> Wqq;
      for(unsigned int j1=0; j1<jets.size(); ++j1){
        for(unsigned int j2=0; j2<jets.size(); ++j2){
          FourMomentum candidate = jets.at(j1).momentum() + jets.at(j2).momentum();
          double mass_diff = fabs(candidate.mass() - MW_PDG);
          Wqq.push_back(reco_boson(candidate, 0, j1, j2, mass_diff));
        }
      }  
      std::sort(Wqq.begin(), Wqq.end(), mass_diff_sorting);
      //Best candidate is the first!
      FourMomentum Hadronic_W = Wqq.at(0).p4;
      unsigned int W_j1_idx = Wqq.at(0).a_idx;
      unsigned int W_j2_idx = Wqq.at(0).b_idx;
      // jet1 and jet2 matched to the hadronic W
      FourMomentum jet1 = jets[W_j1_idx];
      FourMomentum jet2 = jets[W_j2_idx];
      //pT of the jet1 and jet2 matched to the hadronic W   
      double jet1_pt = jet1.pt();
      double jet2_pt = jet2.pt();
      _h_jet1_pt->fill(jet1_pt,weight);
      _h_jet2_pt->fill(jet2_pt,weight);
 
     //want to plot the invariant mass of jets[W_j1_idx] and jets[W_j2_idx] whenever the W invariant mass is within 10 and 40 GeV                                                                   
      double Jet1_Mass = jet1.mass();                      
      double W1 = fabs(Jet1_Mass - MW_PDG);

      double Jet2_Mass = jet2.mass();
      double W2 = fabs(Jet2_Mass- MW_PDG);
      
      if(Hadronic_W.mass()>10*GeV && Hadronic_W.mass()<40*GeV  ){
	if(W1<W2){
	  // cout << "Jet1_Mass selected= "<<Jet1_Mass<<endl;                                                                                                                                            
	  _h_Jet1_Mass->fill(Jet1_Mass,weight);//closer to the smaller peak
	}
	else{
	  // cout << "Jet2_Mass selected= "<<Jet2_Mass<<endl;
	  _h_Jet2_Mass->fill(Jet2_Mass,weight);//farther away from smaller peak
	}
      }
      
      // HT = scalar sum of the jets' and leptons' pt
      for(unsigned int l=0; l<leps.size(); ++l){
        HT+=leps.at(l).pt();
        if(leps.at(l).abspid() == 13) {Nmu+=1;}
        if(leps.at(l).abspid() == 11) {Nel+=1;}
      }
      for(unsigned int j=0; j<jets.size(); ++j){HT+=jets.at(j).pt();}
      cout<<"HT="<<HT<<endl;
      _h_HT->fill(HT, weight);
      
      //Muons, Electrons and Jets multiplicity
      _h_Nmu->fill(Nmu, weight);
      _h_Nel->fill(Nel, weight);
      _h_Nlep->fill(Nmu+Nel, weight);
      _h_Njets->fill(jets.size(),weight);

      //MET  
      _h_MET->fill(p_miss.pt(), weight);

      //pT of leading muon. 
      double leadingmuon_pt=0.;
      for(unsigned int l=0; l<leps.size(); ++l){
        if(leps.at(l).abspid() == 13){leadingmuon_pt = leps.at(l).pt(); break;}
      }
      if(Nmu>0){_h_leadingmuon_pt->fill(leadingmuon_pt, weight);}
 
     //pT of leading electron
      double leadingelectron_pt=0.;
      for(unsigned int l=0; l<leps.size(); ++l){
        if(leps.at(l).abspid() == 11){leadingelectron_pt = leps.at(l).pt(); break;}
      }
      if(Nel>0){_h_leadingelectron_pt->fill(leadingelectron_pt, weight);}

      //pT of leading jet
      _h_leadingjet_pt->fill(jets.at(0).pt(),weight);

      //RECONSTRUCTION variables

      //invariant masses
      _h_Z_m   ->fill(Leptonic_Z.mass(),weight);
      _h_Whad_m->fill(Hadronic_W.mass(),weight);
      _h_Wlep_m->fill(Leptonic_W.mass(),weight);

      //transverse momenta
      _h_Z_pt   ->fill(Leptonic_Z.pt(),weight);
      _h_Whad_pt->fill(Hadronic_W.pt(),weight);
      _h_Wlep_pt->fill(Leptonic_W.pt(),weight);

      //Filling 2dim histos
      _h_HTvsTot_m->fill(HT,(Leptonic_Z+Hadronic_W+Leptonic_W).mass(),weight);
 
     //Invariant Mass of WWZ system ---by summing .mass() for all the three respective four vector and then filling into histo                                                  
      InvMass_WWZ = (Leptonic_Z+Hadronic_W+Leptonic_W).mass();
      // cout<<"InvMass_WWZ ="<<InvMass_WWZ<<endl;
      _h_InvMass_WWZ -> fill(InvMass_WWZ, weight);

      //delta phi, for the leptonic W we use the MET instead of the neutrino
      double Z_dPhi = (Pi-fabs(fabs(leps[Z_lep1_idx].phi()-leps[Z_lep2_idx].phi())-Pi));
      double Whad_dPhi = (Pi-fabs(fabs(jets[W_j1_idx].phi()-jets[W_j2_idx].phi())-Pi));
      double Wlep_dPhi = (Pi-fabs(fabs(leps[W_lep_idx].phi()-p_miss.phi())-Pi));
      
      //Filling Histos for delta of the three bosons
      _h_Z_dPhi   ->fill(Z_dPhi,weight);
      _h_Whad_dPhi->fill(Whad_dPhi,weight);
      _h_Wlep_dPhi->fill(Wlep_dPhi,weight);

      //delta Eta, for the leptonic W the leading neutrino is used                                                                                                     
      double Z_dEta = fabs(leps[Z_lep1_idx].eta()-leps[Z_lep2_idx].eta());
      double Whad_dEta = fabs(jets[W_j1_idx].eta()-jets[W_j2_idx].eta());
      double Wlep_dEta = fabs(leps[W_lep_idx].eta()-nu[0].eta());

      //Filling Histos for deltaEta of the three bosons                                                                                                                        
      _h_Z_dEta   ->fill(Z_dEta,weight);
      _h_Whad_dEta->fill(Whad_dEta,weight);
      _h_Wlep_dEta->fill(Wlep_dEta,weight);
      
      //Delta R
      double Z_dR    = sqrt(( Z_dPhi*Z_dPhi) + (Z_dEta*Z_dEta) );
      double Whad_dR = sqrt(( Whad_dPhi*Whad_dPhi) + (Whad_dEta*Whad_dEta) );
      double Wlep_dR = sqrt(( Wlep_dPhi*Wlep_dPhi) + (Wlep_dEta*Wlep_dEta) );

      //Filling Histos for deltaR of the three bosons                                                                                                                          
      _h_Z_dR   ->fill(Z_dR,weight);
      _h_Whad_dR->fill(Whad_dR,weight);
      _h_Wlep_dR->fill(Wlep_dR,weight);

      //EFFICIENCY OF THE RECONSTRUCTION
      //i.e. comparison with truth information
      //bool correct_Zlep = false;
      //bool correct_Zlep = (leps.at(Z_lep1_idx).constituentLepton().hasAncestor(23) && leps.at(Z_lep2_idx).constituentLepton().hasAncestor(23));
      //work in progress... 


    }// end of analyze()

    void finalize() {
      //    scale(_h_2l2j, crossSection()/sumOfWeights()/femtobarn); // norm to cross section
    }


  private: 

    //function to sort a vector of reco_bosons according to the smallest inv_mass_diff
    static bool mass_diff_sorting(reco_boson & one, reco_boson & two) {return one.mass_diff < two.mass_diff;}

    Histo1DPtr _h_Ntotev;
    Histo1DPtr _h_Njets;
    Histo1DPtr _h_Nmu;  
    Histo1DPtr _h_Nel;  
    Histo1DPtr _h_Nlep;
    Histo1DPtr _h_MET;
    Histo1DPtr _h_HT;
    Histo1DPtr _h_leadingelectron_pt;
    Histo1DPtr _h_leadingmuon_pt;
    Histo1DPtr _h_leadingjet_pt;
    
    Histo1DPtr _h_Z_lep1_pt;
    Histo1DPtr _h_Z_lep2_pt;
    Histo1DPtr _h_W_lep_pt;
    
    Histo1DPtr _h_Z_m   ;
    Histo1DPtr _h_Whad_m;
    Histo1DPtr _h_Wlep_m;
    Histo1DPtr _h_InvMass_WWZ;
    Histo1DPtr _h_singlejet_mass;
    Histo1DPtr _h_Jet1_Mass;
    Histo1DPtr _h_Jet2_Mass;

    Histo1DPtr _h_Z_pt   ; 
    Histo1DPtr _h_Whad_pt; 
    Histo1DPtr _h_Wlep_pt; 
    Histo1DPtr _h_jet1_pt;
    Histo1DPtr _h_jet2_pt;
    Histo2DPtr _h_HTvsTot_m;


    Histo1DPtr _h_Z_dPhi   ; 
    Histo1DPtr _h_Whad_dPhi; 
    Histo1DPtr _h_Wlep_dPhi;

    Histo1DPtr _h_Z_dEta   ;
    Histo1DPtr _h_Whad_dEta;
    Histo1DPtr _h_Wlep_dEta;

    Histo1DPtr _h_Z_dR ;
    Histo1DPtr _h_Whad_dR;
    Histo1DPtr _h_Wlep_dR;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WWZ);
}

