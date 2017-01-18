//================================================================================================
//
// Demo kinematic fitter for dilepton ttbar events and plot distributions
//
//________________________________________________________________________________________________


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TH1D.h>                     // 1D histogram class
#include <TH2D.h>                     // 2D histogram class
#include <TF2.h>   
#include <TLorentzVector.h>           // 4-vector class
#include <TVector2.h>                 // 2-vector class
#include <TRandom.h>
#include <vector>                     // STL vector class
#include <iostream>                   // standard I/O
#include <iomanip>                    // functions to format standard I/O
#include <fstream>                    // functions for file I/O
#include <string>                     // C++ string class
#include <cmath>                      // C++ math library
#include <cassert>


#include "CPlot.hh"                   // helper class for plots
#include "KStyle.hh"                  // style settings for drawing
#include "CSample.hh"                 // helper class to manage samples


#include "../../Utils/interface/Solvettbar.hh"
#include "../../Utils/interface/JetResolutions.hh"
#include "../../Utils/interface/ttbarCandidate.hh"
#include "../../Utils/interface/TopCandidate.hh"
#include "../../Utils/interface/KinematicFitter.hh"

#endif


using namespace std;

//=== FUNCTION DECLARATIONS ======================================================================================

double deltaPhi(const double phi1, const double phi2) {
  double result = phi1 - phi2;
  if     (result >  TMath::Pi()) { result = result - 2*TMath::Pi(); }
  else if(result < -TMath::Pi()) { result = result + 2*TMath::Pi(); }
  return result;
}

TH1D* makeRatioHist(TH1D* hData, TH1D* hMC, const string name, const bool doBlind);

// make "standard" plot
void makePlot(TCanvas *c, const string outname, const string xlabel, const string ylabel, const bool doBlind,
              const vector<TH1D*>& histv, const vector<CSample*>& samplev, TH1D* hExp, TH1D* hRatio,
              const int ich, double lumi, const bool doLogy=false, const double legdx=0, const double legdy=0,
              const double ymin=-1, const double ymax=-1);

void makeHTML(const string outputDir);

//=== MAIN MACRO =================================================================================================

void demoSolveTTbar()
{
  bool verbose = false;

  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  string outputDir;
  outputDir = "demo";
  //
  // Constants
  //
  const double Z_MASS          = 91.1876;    
  const double WMASS           = 80.358;
  const double TOPMASS         = 172.5;
  //
  // Cuts
  //
  const double MET_CUT         = 80; 
  const double LEADLEPPT_CUT   = 30;
  const double TRAILLEPPT_CUT  = 10;
  const float  DPHIDILEPMET_CUT = 1.2;
  const int    NJETS_CUT        = 2;
  const int    NBJETS_CUT       = 1;
  const double DILEPMASS_CUT    = 20;
  //
  // Settings
  //
  const double SOLVE_WP =1.1;//SOLVE_WP = 170.0; //

  // Create output directory
  gSystem->mkdir(outputDir.c_str(), true);
  CPlot::sOutDir = outputDir;
  
  std::ofstream output("output.txt");
  //
  // samples
  //
  vector<CSample*> samplev;

  unsigned int isTT2L = samplev.size();
  samplev.push_back(new CSample("t#bar{t}(2l)", kGreen+1, kGreen+2));
  samplev.back()->fnamev.push_back("/tthome/ssevova/Analysis/11/CMSSW_8_0_20/src/DMSAna/ttDM/baconbits/Spring16_TTTo2L2Nu_powheg_dilepbits.root");
  unsigned int isPS100 = samplev.size(); 
  samplev.push_back(new CSample("PS M_{#phi}=100 M_{#chi}=1", kOrange+7, kOrange+8));
  samplev.back()->fnamev.push_back("../baconbits/Spring16_TTbarDMJets_pseudoscalar_Mchi-1_Mphi-100_dilepbits.root");
  unsigned int isS300 = samplev.size(); 
  samplev.push_back(new CSample("S M_{#phi}=300 M_{#chi}=1", kBlue+1, kBlue+2));
  samplev.back()->fnamev.push_back("../baconbits/Spring16_TTbarDMJets_scalar_Mchi-1_Mphi-300_dilepbits.root");
  
  // integrated lumi to scale MC
  const double LUMI = 36.459;
  // histograms for various corrections
  const string cmssw_base = getenv("CMSSW_BASE");
  const string puWeightFilename = cmssw_base + ("/src/DMSAna/Utils/data/PUWeights_2016.root");
  const int kMT2BINS = 11;
  const double mt2binning[] = {0,15,30,45,60,75,90,105,120,150,180,300};
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  //
  // Declare histograms
  //
  vector<TH1D*> hNSol[3];
  vector<TH1D*> hNPass[3];
  vector<TH1D*> hSolWeight[3];
  vector<TH1D*> hTotErr[3];
  vector<TH1D*> hResErr[3];
  vector<TH1D*> hB1PtErr[3];
  vector<TH1D*> hB2PtErr[3];
  vector<TH1D*> hMETphiErr[3];
  vector<TH1D*> hMETxErr[3];
  vector<TH1D*> hMETyErr[3];
  vector<TH1D*> hMETshift[3];
  vector<TH1D*> hMETfitted[3];
  vector<TH1D*> hMTop1Reco[3], hMTop2Reco[3];
  vector<TH1D*> hMTop1Gen[3],  hMTop2Gen[3];
  vector<TH1D*> hNu1Pt[3], hNu2Pt[3];
  vector<TH1D*> hNu1GenPt[3], hNu2GenPt[3];
  vector<TH1D*> hMT2ll[3];
  vector<TH1D*> hMtt[3];
  vector<TH1D*> hMlb[3];
  vector<TH1D*> hMET[3];
  vector<TH1D*> hMET_mt2[3];
  vector<TH1D*> hdPhiMETB[3];
  vector<TH1D*> hLep1Pt[3], hLep2Pt[3], hB1Pt[3], hB2Pt[3];
  TF2* fLandau;
  vector<TH1D*> hTopPtReco[3], hTopPtGen[3];
  vector <TH2D*> hMttvWeight[3];
  vector <TH2D*> hNuE2D[3];
  vector <TH2D*> hNu1vWeight[3];
  vector <TH2D*> hNu2vWeight[3];
  
  char hname[100];
  
  fLandau = new TF2("landau2D", "[0]*TMath::Landau(x,[1],[2],0)*TMath::Landau(y,[3],[4],0)",0,500,0,500); 
  fLandau -> SetParameters(30.7137,56.2880,23.0744,59.1015,24.9145);

  for(int ich=0; ich<3; ich++){
    for(unsigned int isam=0; isam<samplev.size(); isam++){
      
      sprintf(hname,"hNu1vWeight_%i_%i",ich,isam); hNu1vWeight[ich].push_back(new TH2D(hname,"",50,0,500,50,0,1)); hNu1vWeight[ich][isam]->Sumw2();
      sprintf(hname,"hNu2vWeight_%i_%i",ich,isam); hNu2vWeight[ich].push_back(new TH2D(hname,"",50,0,500,50,0,1)); hNu2vWeight[ich][isam]->Sumw2();
      
      sprintf(hname,"hNuE2D_%i_%i",ich,isam); hNuE2D[ich].push_back(new TH2D(hname,"",50,0,500,50,0,500)); hNuE2D[ich][isam]->Sumw2();
      sprintf(hname,"hNSol_%i_%i",ich,isam); hNSol[ich].push_back(new TH1D(hname,"",5,0,5)); hNSol[ich][isam]->Sumw2();
      sprintf(hname,"hNPass_%i_%i",ich,isam); hNPass[ich].push_back(new TH1D(hname,"",2,-0.1,1.1)); hNPass[ich][isam]->Sumw2();
      sprintf(hname,"hNu1Pt_%i_%i", ich, isam);      hNu1Pt[ich].push_back(new TH1D(hname,"",25,0,100));   hNu1Pt[ich][isam]->Sumw2();
      sprintf(hname,"hNu1GenPt_%i_%i", ich, isam);   hNu1GenPt[ich].push_back(new TH1D(hname,"",25,0,100));   hNu1GenPt[ich][isam]->Sumw2();
      sprintf(hname,"hNu2Pt_%i_%i", ich, isam);      hNu2Pt[ich].push_back(new TH1D(hname,"",25,0,100));   hNu2Pt[ich][isam]->Sumw2();
      sprintf(hname,"hNu2GenPt_%i_%i", ich, isam);   hNu2GenPt[ich].push_back(new TH1D(hname,"",25,0,100));   hNu2GenPt[ich][isam]->Sumw2();
      sprintf(hname,"hSolWeight_%i_%i",ich,isam);  hSolWeight[ich].push_back(new TH1D(hname,"",20,0,1)); hSolWeight[ich][isam]->Sumw2();
      sprintf(hname,"hTotErr_%i_%i",ich,isam);  hTotErr[ich].push_back(new TH1D(hname,"",50,0,5)); hTotErr[ich][isam]->Sumw2();
      sprintf(hname,"hResErr_%i_%i",ich,isam);  hResErr[ich].push_back(new TH1D(hname,"",40,0,200)); hResErr[ich][isam]->Sumw2();
      sprintf(hname,"hB1PtErr_%i_%i",ich,isam);  hB1PtErr[ich].push_back(new TH1D(hname,"",50,-5,5)); hB1PtErr[ich][isam]->Sumw2();
      sprintf(hname,"hB2PtErr_%i_%i",ich,isam);  hB2PtErr[ich].push_back(new TH1D(hname,"",50,-5,5)); hB2PtErr[ich][isam]->Sumw2();
      sprintf(hname,"hMETphiErr_%i_%i",ich,isam);  hMETphiErr[ich].push_back(new TH1D(hname,"",25,0,5)); hMETphiErr[ich][isam]->Sumw2();
      sprintf(hname,"hMETxErr_%i_%i",ich,isam);  hMETxErr[ich].push_back(new TH1D(hname,"",50,-5,5)); hMETxErr[ich][isam]->Sumw2();
      sprintf(hname,"hMETyErr_%i_%i",ich,isam);  hMETyErr[ich].push_back(new TH1D(hname,"",50,-5,5)); hMETyErr[ich][isam]->Sumw2();
      sprintf(hname,"hMETshift_%i_%i",ich,isam); hMETshift[ich].push_back(new TH1D(hname,"",25,0,1)); hMETshift[ich][isam]->Sumw2();
      sprintf(hname,"hMETfitted_%i_%i",ich,isam); hMETfitted[ich].push_back(new TH1D(hname,"",25,50,550)); hMETfitted[ich][isam]->Sumw2();
      sprintf(hname,"hMTop1Reco_%i_%i",ich,isam);  hMTop1Reco[ich].push_back(new TH1D(hname,"",20,150,190)); hMTop1Reco[ich][isam]->Sumw2();
      sprintf(hname,"hMTop1Gen_%i_%i",ich,isam);   hMTop1Gen[ich].push_back(new TH1D(hname,"",20,150,190));  hMTop1Gen[ich][isam]->Sumw2();
      sprintf(hname,"hMTop2Reco_%i_%i",ich,isam);  hMTop2Reco[ich].push_back(new TH1D(hname,"",20,150,190)); hMTop2Reco[ich][isam]->Sumw2();
      sprintf(hname,"hMTop2Gen_%i_%i",ich,isam);   hMTop2Gen[ich].push_back(new TH1D(hname,"",20,150,190));  hMTop2Gen[ich][isam]->Sumw2();
      sprintf(hname,"hMttvWeight_%i_%i",ich,isam); hMttvWeight[ich].push_back(new TH2D(hname,"",10,300,500,20,0,1)); hMttvWeight[ich][isam]->Sumw2();
      sprintf(hname,"hMtt_%i_%i",ich,isam);        hMtt[ich].push_back(new TH1D(hname,"",25,300,800)); hMtt[ich][isam]->Sumw2();
      sprintf(hname,"hMlb_%i_%i",ich,isam);        hMlb[ich].push_back(new TH1D(hname,"",40,0,200)); hMlb[ich][isam]->Sumw2(); 
      sprintf(hname,"hTopPtReco_%i_%i",ich,isam);  hTopPtReco[ich].push_back(new TH1D(hname,"",25,0,500)); hTopPtReco[ich][isam]->Sumw2();
      sprintf(hname,"hTopPtGen_%i_%i",ich,isam);   hTopPtGen[ich].push_back(new TH1D(hname,"",25,0,500));  hTopPtGen[ich][isam]->Sumw2();
      sprintf(hname,"hMT2ll_%i_%i",ich,isam);      hMT2ll[ich].push_back(new TH1D(hname,"",kMT2BINS, mt2binning)); hMT2ll[ich][isam]->Sumw2();
      
      sprintf(hname,"hMET_%i_%i",ich,isam);        hMET[ich].push_back(new TH1D(hname,"",25,50,550)); hMET[ich][isam]->Sumw2();
      sprintf(hname,"hMET_mt2_%i_%i",ich,isam);    hMET_mt2[ich].push_back(new TH1D(hname,"",25,50,550)); hMET_mt2[ich][isam]->Sumw2();
      sprintf(hname,"hdPhiMETB_%i_%i",ich,isam);   hdPhiMETB[ich].push_back(new TH1D(hname,"",50,0,3)); hdPhiMETB[ich][isam]->Sumw2();
      //per solution num
      sprintf(hname,"hLep1Pt_%i_%i",ich,isam);      hLep1Pt[ich].push_back(new TH1D(hname,"",25,0,150)); hLep1Pt[ich][isam]->Sumw2();
      sprintf(hname,"hLep2Pt_%i_%i",ich,isam);      hLep2Pt[ich].push_back(new TH1D(hname,"",25,0,150)); hLep2Pt[ich][isam]->Sumw2();
      sprintf(hname,"hB1Pt_%i_%i",ich,isam);      hB1Pt[ich].push_back(new TH1D(hname,"",25,0,150)); hB1Pt[ich][isam]->Sumw2();
      sprintf(hname,"hB2Pt_%i_%i",ich,isam);      hB2Pt[ich].push_back(new TH1D(hname,"",25,0,150)); hB2Pt[ich][isam]->Sumw2();

    }
  }
  //
  // variables to read in bacon bits
  //
  unsigned int runNum, lumiSec, evtNum;          // event ID
  unsigned int metfilter;                        // MET filter bits
  unsigned int npv;                              // number of PV
  unsigned int njets, njetsC;                    // jet multiplicity
  unsigned int nbjetsL, nbjetsM, nbjetsT;        // b-tag multiplicity
  unsigned int nEle, nMu, nTau, nPho;            // loose lepton multiplicity
  float scale1fb;                                // cross-section scale factor per 1/fb
  float npu;                                     // mean expected PU
  float pfmet, pfmetphi;                         // PF MET
  float mt2ll;
  int lep1Id, lep2Id;                            // lepton PDG ID
  TLorentzVector *lep1=0, *lep2=0, *dilep=0;     // lepton 4-vector

  //generator level info
  float           genHT;
  int             genId11, genId12, genId13;
  int             genId21, genId22, genId23;
  TLorentzVector *genpar11=0, *genpar12=0, *genpar13=0;
  TLorentzVector *genpar21=0, *genpar22=0, *genpar23=0;
  
  int             genTop1Id, genTop2Id;
  TLorentzVector *genparTop1=0, *genparTop2=0;

  float           rhojet;
  float           jet1csv,  jet2csv,  jet3csv,  jet4csv;
  TLorentzVector *jet1=0,  *jet2=0,  *jet3=0,  *jet4=0;
  int lep1hlt, lep2hlt;

  TFile *infile=0;
  TTree *intree=0;

  double npassWP[3]={0,0,0};    double nfailWP[3]={0,0,0};
  double npassMT2WP[3]={0,0,0}; double nfailMT2WP[3]={0,0,0};
  double numsol[3]={0,0,0};     double nosol[3]={0,0,0};
  double nevts [3]={0,0,0};

  for(unsigned int isam=0; isam<samplev.size(); isam++) {
    CSample *sample = samplev[isam];
    cout << "Sample: " << sample->label << endl;

    unsigned int NEVT=0;
    if(isam==0) NEVT = 10000;
    else        NEVT = 5000;

    for(unsigned int ifile=0; ifile<sample->fnamev.size(); ifile++) {
      string infilename = sample->fnamev[ifile];
      cout << " ==> Processing " << infilename << "... " << endl;
      infile = new TFile(infilename.c_str()); assert(infile);
      intree = (TTree*)infile->Get("Events"); assert(intree);

      intree->SetBranchAddress("runNum",       &runNum);
      intree->SetBranchAddress("lumiSec",      &lumiSec);
      intree->SetBranchAddress("evtNum",       &evtNum);
      intree->SetBranchAddress("scale1fb",     &scale1fb);
      intree->SetBranchAddress("metfilter",    &metfilter);
      intree->SetBranchAddress("npv",          &npv);
      intree->SetBranchAddress("njets",        &njets);
      intree->SetBranchAddress("njetsC",       &njetsC);
      intree->SetBranchAddress("nbjetsL",      &nbjetsL);
      intree->SetBranchAddress("nbjetsM",      &nbjetsM);
      intree->SetBranchAddress("nbjetsT",      &nbjetsT);
      intree->SetBranchAddress("nEle",         &nEle);
      intree->SetBranchAddress("nMu",          &nMu);
      intree->SetBranchAddress("nTau",         &nTau);
      intree->SetBranchAddress("npu",          &npu);
      intree->SetBranchAddress("pfmet",        &pfmet);
      intree->SetBranchAddress("pfmetphi",     &pfmetphi);
      intree->SetBranchAddress("mt2ll",        &mt2ll);
      intree->SetBranchAddress("jet1csv",      &jet1csv);
      intree->SetBranchAddress("jet2csv",      &jet2csv);
      intree->SetBranchAddress("jet3csv",      &jet3csv);
      intree->SetBranchAddress("jet4csv",      &jet4csv);
      intree->SetBranchAddress("jet1",         &jet1);
      intree->SetBranchAddress("jet2",         &jet2);
      intree->SetBranchAddress("jet3",         &jet3);
      intree->SetBranchAddress("jet4",         &jet4);
      intree->SetBranchAddress("rhojet",       &rhojet);
      intree->SetBranchAddress("lep1Id",    &lep1Id);
      intree->SetBranchAddress("lep1",      &lep1);
      intree->SetBranchAddress("lep2Id",    &lep2Id);
      intree->SetBranchAddress("lep2",      &lep2);
      intree->SetBranchAddress("dilep",     &dilep);      
      intree->SetBranchAddress("lep1hlt",   &lep1hlt);
      intree->SetBranchAddress("lep2hlt",   &lep2hlt);
      intree->SetBranchAddress("genHT",        &genHT);
      intree->SetBranchAddress("genId11",      &genId11);
      intree->SetBranchAddress("genId12",      &genId12);
      intree->SetBranchAddress("genId13",      &genId13);
      intree->SetBranchAddress("genId21",      &genId21);
      intree->SetBranchAddress("genId22",      &genId22);
      intree->SetBranchAddress("genId23",      &genId23);
      intree->SetBranchAddress("genpar11",     &genpar11);
      intree->SetBranchAddress("genpar12",     &genpar12);
      intree->SetBranchAddress("genpar13",     &genpar13);
      intree->SetBranchAddress("genpar21",     &genpar21);
      intree->SetBranchAddress("genpar22",     &genpar22);
      intree->SetBranchAddress("genpar23",     &genpar23);
      /*      
      intree->SetBranchAddress("genparTop1",   &genparTop1);
      intree->SetBranchAddress("genparTop2",   &genparTop2);
      intree->SetBranchAddress("genTop1Id",    &genTop1Id);
      intree->SetBranchAddress("genTop2Id",    &genTop2Id);
      */

      TRandom * mysmear = new TRandom();	

      //instantiate solver object pointer
      Solvettbar * solver = new Solvettbar(160,180,1); 
      // instantiate solution storing object
      FitSolution sol;
      

      int nsol2=0, nsol4=0, nsol0=0;

      std::cout << "start loop ... " << std::endl;
      // for(unsigned int ientry=0; ientry<intree->GetEntries()/1000; ientry++) {
      for(unsigned int ientry=0; ientry<NEVT; ientry++) {
	sol.reset();
        intree->GetEntry(ientry);
	//	if(ientry%1000 == 0)
	//  std::cout << "entry: " << ientry << std::endl;

	TLorentzVector vMET;
        vMET.SetPtEtaPhiM(pfmet,0,pfmetphi,0);  
	TLorentzVector vDilep;
	vDilep = (*lep1) + (*lep2);

	float DphiDilepMET  = fabs(vDilep.DeltaPhi(vMET));	
	int ich=-1;
	if     (abs(lep1Id)==11 && abs(lep2Id)==11) { ich=0; }  // ee
	else if(abs(lep1Id)==13 && abs(lep2Id)==13) { ich=2; }  // mm
	else                                        { ich=1; }  // em

	if(lep1->Pt() < LEADLEPPT_CUT)      continue;
	if(lep2->Pt() < TRAILLEPPT_CUT)     continue; 
	if(vMET.Pt()  < MET_CUT)            continue;
	if(njets      < NJETS_CUT)          continue;
	if(nbjetsM    < NBJETS_CUT)         continue; 
	// if(DphiDilepMET < DPHIDILEPMET_CUT) continue; 
	if(vDilep.M()   < DILEPMASS_CUT)    continue;
	if(ich==0 || ich==2){
	  if(fabs(Z_MASS - vDilep.M())<15)  continue; 
	}
	if((nEle+nMu) != 2)                 continue;
	cout << "================================ new event =========================" << runNum<<" "<<evtNum<<endl;

	double weight = 1;
	// weight *= LUMI*scale1fb; 

	if(mt2ll>100){
	  npassMT2WP[isam]+=weight;
	  hMET_mt2[ich][isam]->Fill(pfmet);
	} else {
          nfailMT2WP[isam]+=weight;
	}

	vector<double> vjetcsv;
	vector<TLorentzVector> vjet;
	vjetcsv.push_back(jet1csv); vjetcsv.push_back(jet2csv); vjetcsv.push_back(jet3csv); vjetcsv.push_back(jet4csv);
	vjet.push_back(*jet1); 	    vjet.push_back(*jet2); 	vjet.push_back(*jet3);	    vjet.push_back(*jet4);

	if(verbose){
	  cout<<"SORTING by CSV"<<endl;
	  cout<<"-------UNSORTED-----"<<endl;
	  for(unsigned int i=0; i<vjet.size(); i++)
	    cout<<i<<" CSV "<<vjetcsv[i]<<endl;
	}

	double tempcsv=-999;
	TLorentzVector tempJet;

	for(unsigned int i=0; i<vjet.size(); i++){
	  for(unsigned int j=i+1; j<vjet.size(); j++){
	    if(vjetcsv[i] > vjetcsv[j]){
	      tempcsv = vjetcsv[i];     tempJet = vjet[i];
	      vjetcsv[i] = vjetcsv[j];  vjet[i] = vjet[j];
	      vjetcsv[j] = tempcsv;     vjet[j] = tempJet;
	    }
	  }
	}

	if(verbose){
	  cout<<"-------SORTED-----"<<endl;
	  for(unsigned int i=0; i<vjet.size(); i++)
	    cout<<i<<" CSV "<<vjetcsv[i]<<endl;
	}

	for(unsigned int i=0; i<vjet.size(); i++){
	  for(unsigned int j=i+1; j<vjet.size(); j++){
	    assert(vjetcsv[i]<=vjetcsv[j]);
	  }//for(unsigned int j=i+1; j<vjet.size(); j++)
	}//for(unsigned int i=0; i<vjet.size(); i++)


	if(verbose) cout << "done bjet" << endl;

	TLorentzVector vb1, vb2;
	vb1.SetPtEtaPhiM(vjet[vjet.size()-1].Pt(), vjet[vjet.size()-1].Eta(), vjet[vjet.size()-1].Phi(), vjet[vjet.size()-1].M());
	vb2.SetPtEtaPhiM(vjet[vjet.size()-2].Pt(), vjet[vjet.size()-2].Eta(), vjet[vjet.size()-2].Phi(), vjet[vjet.size()-2].M());
		
	TLorentzVector vnull; vnull.SetPtEtaPhiM(0,0,0,0);
	ttbarCandidate::ttbarCandidateParticle l1,n1,b1,l2,n2,b2;
	n1 = ttbarCandidate::ttbarCandidateParticle(vnull,0);
	n2 = ttbarCandidate::ttbarCandidateParticle(vnull,0);
	
	if(verbose) cout << "gen id 11: "<<genId11<< "\tgen id 21: "<<genId21<<endl;

	// save original 4-vecs before smearing
	TLorentzVector olep1, olep2, oB1, oB2, oMET;
	olep1.SetPtEtaPhiM(lep1->Pt(), lep1->Eta(), lep1->Phi(), lep1->M());
	olep2.SetPtEtaPhiM(lep2->Pt(), lep2->Eta(), lep2->Phi(), lep2->M());
	oB1.SetPtEtaPhiM(vb1.Pt(), vb1.Eta(), vb1.Phi(), vb1.M());
	oB2.SetPtEtaPhiM(vb2.Pt(), vb2.Eta(), vb2.Phi(), vb2.M());
	oMET.SetPtEtaPhiM(vMET.Pt(), vMET.Eta(), vMET.Phi(), vMET.M());


	// smear according to jet energy resolutions 

	JetResolutions myjer;
	TMatrixD covMat(2,2);
	//MT according to http://cms.cern.ch/iCMS/analysisadmin/cadilines?line=SMP-12-011&tp=an&id=862&ancode=SMP-12-011 recoil is defined as vectorial negative sum of everything once charged lepton(s) is subtracted for W/Z events. So for ttbar event I need to subtract all the visible objects => MET+visible objects
	// myjer.covMatrix(covMat,(oMET-olep1-olep2-oB1-oB2).Pt(),(oMET-olep1-olep2-oB1-oB2).Phi(),2);
	//
	myjer.covMatrix(covMat,(oMET+olep1+olep2+oB1+oB2).Pt(),(oMET+olep1+olep2+oB1+oB2).Phi(),2); //to be plotted
	if(verbose) covMat.Print();
	//	cout.flush();

	double b1ptsigma  = myjer.getPtSigmaSmear (oB1,rhojet); 
	double b2ptsigma  = myjer.getPtSigmaSmear (oB2,rhojet);
	double b1phisigma = myjer.getPhiSigmaSmear(oB1,rhojet);
	double b2phisigma = myjer.getPhiSigmaSmear(oB2,rhojet);
	
	float pxmiss_, pymiss_;
	int foundtrial=0;

	//RECO candidates
	l1 = ttbarCandidate::ttbarCandidateParticle(*lep1,0);
	l2 = ttbarCandidate::ttbarCandidateParticle(*lep2,0);

	TLorentzVector top1vec =  *genpar11 + *genpar12 + *genpar13;//*genparTop1;//
	TLorentzVector top2vec =  *genpar21 + *genpar22 + *genpar23;//*genparTop2;//

	hNu1GenPt[ich][isam]  ->Fill(genpar22->Pt());
	hNu2GenPt[ich][isam]  ->Fill(genpar12->Pt());
	hMTop1Gen[ich][isam]  ->Fill(top1vec.M());
	hMTop2Gen[ich][isam]  ->Fill(top2vec.M());
	hTopPtGen[ich][isam]  ->Fill(top1vec.Pt());

	TLorentzVector  iU, oU;
	std::cout << "oMET.Pt() : " << oMET.Pt() << "\toMET.Phi(): " << oMET.Phi() << std::endl;

	int ngoodtrials=0;
	double avg_recotop_pt=0;
	double avg_mtt = 0;
	double avg_Mlep1b1 = 0;
	double avg_Mlep2b2 =0;
	double avg_dphimetb =0;
	double avg_nu1_e =0;
	double avg_nu2_e =0;
	double avg_nu1_pt =0;
	double avg_nu2_pt =0;
	double avg_nu1_phi=0;
	double avg_nu2_phi=0;
	double avg_nu1_eta=0;
	double avg_nu2_eta=0;
	double avg_sol_w =0;
	double avg_toterr =0;
	double avg_reserr =0;
	double avg_top1_recoM =0;
	double avg_top2_recoM =0;
	double avg_b1pt_err=0;
	double avg_b2pt_err=0;
	double avg_metx_err=0;
	double avg_mety_err=0;
	double avg_metphi_err=0;
	
	for( int trial=0; trial<1000; trial++){
	  sol.reset();

	  unsigned int k=4;
	  vb1.SetPtEtaPhiM(mysmear->Gaus(oB1.Pt(), k*b1ptsigma*oB1.Pt()),oB1.Eta(),mysmear->Gaus(oB1.Phi(),k*b1phisigma*oB1.Phi()),oB1.M());
	  vb2.SetPtEtaPhiM(mysmear->Gaus(oB2.Pt(), k*b2ptsigma*oB2.Pt()),oB2.Eta(),mysmear->Gaus(oB2.Phi(),k*b2phisigma*oB2.Phi()),oB2.M());

	  double err_b1_pt  = (vb1.Pt()  - oB1.Pt()) /(k*b1ptsigma*oB1.Pt());
	  double err_b1_phi = (vb1.Phi() - oB1.Phi())/(k*b1ptsigma*oB1.Phi());
	  double err_b2_pt  = (vb2.Pt()  - oB2.Pt()) /(k*b2ptsigma*oB2.Pt());
	  double err_b2_phi = (vb2.Phi() - oB2.Phi())/(k*b2ptsigma*oB2.Phi());

	  // KH smear MET
	  // iU = oMET-olep1-olep2-vb1-vb2; //MT: replaced by line below (+ instead of -)
	  iU = oMET+olep1+olep2+vb1+vb2; 
	  double sX,sY;
	  myjer.getSmearedMET(&iU, &covMat, &oU, sX, sY);

	  if(verbose)   std::cout << "sX: " << sX << "\tsY: " << sY << std::endl;

	  // vMET = oU+olep1+olep2+vb1+vb2; //MT: replaced by line below (- instead of +)
	  vMET = oU-olep1-olep2-vb1-vb2;

	  double dMETx = vMET.Px() - oMET.Px();
	  double dMETy = vMET.Py() - oMET.Py();
	  
	  double err_metx = dMETx/(2*sX);
	  double err_mety = dMETy/(2*sY);
	  
	  double dMETPhi = vMET.Phi() - oMET.Phi();
	  double err_metphi = fabs(dMETPhi/oMET.Phi());

	  double dPhiMETB = fabs(vMET.DeltaPhi(vb1));

	  if(verbose)    cout << "vMET.Pt(): " << vMET.Pt() <<  "\tvMET.Phi(): " << vMET.Phi() << endl;

	  double toterr = // sqrt( err_b1_pt*err_b1_pt + err_b1_phi*err_b1_phi + err_b2_pt*err_b2_pt + err_b2_phi*err_b2_phi
	    sqrt( err_metx*err_metx + err_mety*err_mety );
	    // TMath::Exp(-err_metx*err_metx)*TMath::Exp(-err_mety*err_mety);
	  // cout << toterr << endl;
	  double res_meterr = sqrt(dMETx*dMETx + dMETy*dMETy);
	  
	  for(int i=0; i<2; i++){
	    if(i==1 && sol.numSol!=0) continue;
	    
	    if(i==0){
	      b1 = ttbarCandidate::ttbarCandidateParticle(vb1,0);
	      b2 = ttbarCandidate::ttbarCandidateParticle(vb2,0);
	    } else {
	      b1 = ttbarCandidate::ttbarCandidateParticle(vb2,0);
	      b2 = ttbarCandidate::ttbarCandidateParticle(vb1,0);
	    }

	    //must be passed in this order
	    // lep1, neu1 (null), b1; lep2, neu2 (mull), b2
	    ttbarCandidate ttbar(l1,n1,b1,l2,n2,b2);                                                                            
	    
	    pxmiss_ = ttbar.lep1vec.Px() + ttbar.B1vec.Px() + ttbar.lep2vec.Px() + ttbar.B2vec.Px() + vMET.Px();
	    pymiss_ = ttbar.lep1vec.Py() + ttbar.B1vec.Py() + ttbar.lep2vec.Py() + ttbar.B2vec.Py() + vMET.Py();
	    
	    solver->SetConstraints(pxmiss_,pymiss_);
	    solver->SetGenTopMass(TOPMASS, TOPMASS);//top1vec, top2vec);
	    solver->SetWMass(WMASS, WMASS);
	    solver->SolveNu(sol,ttbar);
	    
	    if(sol.numSol>0){
	      if(verbose){
		cout<<"SORTING SOL by WEIGHT from low to high"<<endl;
		cout<<"-------UNSORTED-----"<<endl;
		for(int i=0; i<sol.numSol; i++)
		  cout<<i<<" weight "<<sol.weight[i]<<endl;
	      }
      
	      // sort solutions according to weight
	      float temp_solw=-999; TLorentzVector temp_vNu1, temp_vNu2; float temp_topm1, temp_topm2;
	      for(int i=0; i < sol.numSol; i++){
		for(int j=i+1; j < sol.numSol; j++){
		  if(sol.weight[i] > sol.weight[j]){
		    temp_solw = sol.weight[i];     temp_vNu1 = sol.nu1[i];  temp_vNu2 = sol.nu2[i];  temp_topm1 = sol.mass1[i];   temp_topm2 = sol.mass2[i];
		    sol.weight[i] = sol.weight[j]; sol.nu1[i] = sol.nu1[j]; sol.nu2[i] = sol.nu2[j]; sol.mass1[i] = sol.mass1[j]; sol.mass2[i] = sol.mass2[j];
		    sol.weight[j] = temp_solw;     sol.nu1[j] = temp_vNu1;  sol.nu2[j] = temp_vNu2;  sol.mass1[j] = temp_topm1;   sol.mass2[j] = temp_topm2;
		  }
		}
	      }
	      
	      if(verbose){
		cout<<"-------SORTED-----"<<endl;
		for(int i=0; i<sol.numSol; i++)
		  cout<<i<<" weight "<<sol.weight[i]<<" "<<(int)(100000*sol.weight[i])<<endl;
	      }
	      
	      //comparing weights using 5 decimal points (assuming weight<=1)
	      //sometimes it fails for some rounding problem???
	      // for(unsigned int i=0; i<sol.numSol; i++){
	      // 	for(unsigned int j=i+1; j<sol.numSol; j++){
	      // 	  if((int)(1000000*sol.weight[i])<=int(1000000*(sol.weight[j])));
	      // 	  else{
	      // 	    cout<<"sol.weight["<<i<<"] "<<(int)(1000000*sol.weight[i])<<" more than sol.weight["<<j<<"] "<<(int)(1000000*sol.weight[j])<<endl;
	      // 	    exit(1);
	      // 	  }
	      // 	  //		  assert((int)(100000*sol.weight[i])<=int(100000*(sol.weight[j]))); 
	      // 	}//for(unsigned int i=0; i<sol.numSol; i++)
	      // }//for(unsigned int j=i+1; j<sol.numSol; j++)
	      
	      if(verbose) std::cout << "sol.mass1: " << sol.mass1[sol.numSol-1] << "\tsol.mass2: " << sol.mass2[sol.numSol-1] << std::endl;

	      double M_lep1_B1 = (olep1 + oB1).M();
	      double M_lep2_B2 = (olep2 + oB2).M();

	      TLorentzVector recotop1, recotop2; 
	      recotop1 = olep1 + sol.nu1[sol.numSol-1] + oB1;
	      recotop2 = olep2 + sol.nu2[sol.numSol-1] + oB2;
	      avg_recotop_pt += recotop1.Pt();
	      avg_mtt        += (recotop1 + recotop2).M();//sol.weight[sol.numSol-1];	  
	      avg_Mlep1b1    += M_lep1_B1;
	      avg_Mlep2b2    += M_lep2_B2;
	      avg_dphimetb   += dPhiMETB;
	      avg_nu1_e      += sol.nu1[sol.numSol-1].E();
	      avg_nu2_e      += sol.nu2[sol.numSol-1].E();
	      avg_nu1_pt     += sol.nu1[sol.numSol-1].Pt();
	      avg_nu2_pt     += sol.nu2[sol.numSol-1].Pt();
	      avg_nu1_phi    += sol.nu1[sol.numSol-1].Phi();
	      avg_nu2_phi    += sol.nu2[sol.numSol-1].Phi();
	      avg_nu1_eta    += sol.nu1[sol.numSol-1].Eta();
	      avg_nu2_eta    += sol.nu2[sol.numSol-1].Eta();
	      avg_sol_w      += sol.weight[sol.numSol-1];
	      avg_toterr     += toterr;
	      avg_reserr     += res_meterr;
	      avg_b1pt_err   += err_b1_pt;
	      avg_b2pt_err   += err_b2_pt;
	      avg_metx_err   += err_metx;//sX;
	      avg_mety_err   += err_mety;//sY;
	      avg_metphi_err += err_metphi;
	      avg_top1_recoM += sol.mass1[sol.numSol-1];
	      avg_top2_recoM += sol.mass2[sol.numSol-1];
	    }  
	  }
	  if(sol.numSol!=0){
	    foundtrial = trial;
	    ngoodtrials++;
	  }
	}

	if ( ngoodtrials > 0 ) { 
	  avg_recotop_pt /= ngoodtrials;
	  avg_mtt        /= ngoodtrials;	  
	  avg_Mlep1b1    /= ngoodtrials;
	  avg_Mlep2b2    /= ngoodtrials;
	  avg_dphimetb   /= ngoodtrials;
	  avg_nu1_e      /= ngoodtrials;
	  avg_nu2_e      /= ngoodtrials;
	  avg_nu1_pt     /= ngoodtrials;
	  avg_nu2_pt     /= ngoodtrials; 
	  avg_nu1_phi    /= ngoodtrials;
	  avg_nu2_phi    /= ngoodtrials;
	  avg_nu1_eta    /= ngoodtrials;
	  avg_nu2_eta    /= ngoodtrials;
	  avg_sol_w      /= ngoodtrials;
	  avg_toterr     /= ngoodtrials;
	  avg_reserr     /= ngoodtrials;
	  avg_top1_recoM /= ngoodtrials;
	  avg_top2_recoM /= ngoodtrials;
	  avg_b1pt_err   /= ngoodtrials;
	  avg_b2pt_err   /= ngoodtrials;
	  avg_metx_err   /= ngoodtrials;
	  avg_mety_err   /= ngoodtrials;
	  avg_metphi_err /= ngoodtrials;

	  hMT2ll[ich][isam]      ->Fill(mt2ll);
	  hMlb[ich][isam]        ->Fill(avg_Mlep1b1);
	  hMlb[ich][isam]        ->Fill(avg_Mlep2b2);
	  hMtt[ich][isam]        ->Fill(avg_mtt);
	  hNu1vWeight[ich][isam] ->Fill(avg_nu1_e,avg_sol_w);
	  hNu2vWeight[ich][isam] ->Fill(avg_nu2_e,avg_sol_w);
	  hNuE2D[ich][isam]      ->Fill(avg_nu1_e,avg_nu2_e);
	  hNu1Pt[ich][isam]      ->Fill(avg_nu1_pt);
	  hNu2Pt[ich][isam]      ->Fill(avg_nu2_pt);
	  hSolWeight[ich][isam]  ->Fill(avg_sol_w);
	  hTotErr[ich][isam]     ->Fill(avg_toterr);
	  hResErr[ich][isam]     ->Fill(avg_reserr);
	  hB1PtErr[ich][isam]    ->Fill(avg_b1pt_err);
	  hB2PtErr[ich][isam]    ->Fill(avg_b2pt_err);
	  hMETxErr[ich][isam]    ->Fill(avg_metx_err);
	  hMETyErr[ich][isam]    ->Fill(avg_mety_err);
	  hMETphiErr[ich][isam]  ->Fill(avg_metphi_err);
	  hMTop1Reco[ich][isam]  ->Fill(avg_top1_recoM);
	  hMTop2Reco[ich][isam]  ->Fill(avg_top2_recoM);
	  hTopPtReco[ich][isam]  ->Fill(avg_recotop_pt);
	  hdPhiMETB[ich][isam]   ->Fill(avg_dphimetb);


	  TLorentzVector avg_nu1v, avg_nu2v;
	  avg_nu1v.SetPtEtaPhiE(avg_nu1_pt, avg_nu1_eta, avg_nu1_phi, avg_nu1_e);
	  avg_nu2v.SetPtEtaPhiE(avg_nu2_pt, avg_nu2_eta, avg_nu2_phi, avg_nu2_e);

	  TLorentzVector fMET;
	  fMET = avg_nu1v + avg_nu2v;

	  std::cout << "fMET.Pt() : " << fMET.Pt() << "\tfMET.Phi(): " << fMET.Phi() << std::endl;
	  
	  double metshift = fabs(oMET.Pt() - fMET.Pt())/oMET.Pt();
	  hMETshift [ich][isam]->Fill(metshift); 
	  hMETfitted[ich][isam]->Fill(fMET.Pt());
	  hNPass    [ich][isam]->Fill(1);
	  numsol[isam]         += weight;
	  
	  if(avg_toterr > SOLVE_WP){//if(avg_reserr > SOLVE_WP){//
	    npassWP[isam]      += weight; 	  
	    hMET   [ich][isam]->Fill(pfmet);
	  } else {
	    nfailWP[isam]      +=weight;
	  }
	} else {//no solution
	  hMET     [ich][isam]->Fill(pfmet);
	  hNPass   [ich][isam]->Fill(0);
	  nosol[isam]         +=weight;
	  npassWP[isam]       +=weight;
	}

	hNSol      [ich][isam]->Fill(sol.numSol);
	/*
	if(sol.numSol==0){ nsol0++; }//isol=0; }
	if(sol.numSol==2){ nsol2++; }//isol=1; }
	if(sol.numSol==4){ nsol4++; }//isol=2; }
	*/
	nevts[isam]++;
	string issol;
	issol = ngoodtrials > 0 ? "yes" : "no";
	cout << "evt           : " << nevts[isam] << endl;
	cout << "solution?     : " << issol << endl;
	cout << "ngoodtrials   : " << ngoodtrials << endl;
	cout << "top1vec.M(): " << TOPMASS << "\ttop2vec.M(): " << TOPMASS << endl;
	if(issol.compare("yes")==0) cout << "reco mtop1 : " << avg_top1_recoM << "\treco mtop2 : " << avg_top2_recoM << endl;
      }
      delete infile;
      infile=0;
      intree=0;      
    }
  }

  output << "========== Output Stats ==========" << endl;
  output << "MT2(ll) > 100 S/sqrt(B): " << npassMT2WP[isPS100]/sqrt(npassMT2WP[isTT2L]) << endl;
  output << "MT2(ll) < 100 S/sqrt(B): " << nfailMT2WP[isPS100]/sqrt(nfailMT2WP[isTT2L]) << endl;
  output << "---> evts " << samplev[isPS100]->label << ": " << nevts[isPS100] << endl;
  output << "---> sig acc: " << npassMT2WP[isPS100]/nevts[isPS100] << endl; 
  output << "---> bkg acc: " << npassMT2WP[isTT2L]/nevts[isTT2L] << endl; 
  output << "----------------------------------" << endl;
  output << "tot err > " << SOLVE_WP << " S/sqrt(B): " << npassWP[isPS100]/sqrt(npassWP[isTT2L]) << endl;
  output << "tot err < " << SOLVE_WP << " S/sqrt(B): " << nfailWP[isPS100]/sqrt(nfailWP[isTT2L]) << endl;
  output << "---> evts " << samplev[isPS100]->label << ": " << nevts[isPS100] << endl;
  output << "---> sig acc: " << npassWP[isPS100]/nevts[isPS100] << endl;
  output << "---> bkg acc: " << npassWP[isTT2L]/nevts[isTT2L] << endl; 
  output << "----------------------------------" << endl;
  for(unsigned int isam=0; isam<samplev.size(); isam++){
    output << "*** sample  : " << samplev[isam]->label << " *** " << endl;
    output << "yes solution: " << numsol[isam] << endl;
    output << "no solution : " << nosol[isam] << endl; 
    output << "efficiency  : " << numsol[isam]/(numsol[isam]+nosol[isam]) << endl;
  }
  output.close();
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  gStyle->SetPalette(1);
  TCanvas *c = MakeCanvas("c","c",800,600);
  c->SetTickx(1);
  c->SetTicky(1);

  char pname[100];
  char xlabel[100];
  char ylabel[100];
  string suffix;

  for(unsigned int ich=0; ich<3; ich++){
    for(unsigned int isam=0; isam<samplev.size(); isam++){
      hMlb    [ich][isam]->Scale(1./hMlb[ich][isam]->Integral());
      hMtt    [ich][isam]->Scale(1./hMtt[ich][isam]->Integral());
      hTotErr [ich][isam]->Scale(1./hTotErr[ich][isam]->Integral());
      hResErr [ich][isam]->Scale(1./hResErr[ich][isam]->Integral());
      hB1PtErr[ich][isam]->Scale(1./hB1PtErr[ich][isam]->Integral());
      hB2PtErr[ich][isam]->Scale(1./hB2PtErr[ich][isam]->Integral());
      hMETxErr[ich][isam]->Scale(1./hMETxErr[ich][isam]->Integral());
      hMETyErr[ich][isam]->Scale(1./hMETyErr[ich][isam]->Integral());
      hMETphiErr[ich][isam]->Scale(1./hMETphiErr[ich][isam]->Integral());
      hMT2ll  [ich][isam]->Scale(1./hMT2ll  [ich][isam]->Integral());
      hMETshift[ich][isam]->Scale(1./hMETshift[ich][isam]->Integral());
      hMET     [ich][isam]->Scale(1./hMET[ich][isam]->Integral());
      hMET_mt2 [ich][isam]->Scale(1./hMET_mt2[ich][isam]->Integral());
      hMETfitted[ich][isam]->Scale(1./hMETfitted[ich][isam]->Integral());
      hdPhiMETB[ich][isam] ->Scale(1./hdPhiMETB[ich][isam]->Integral());
    }

    if(ich==0) { suffix = "ee"; }
    if(ich==1) { suffix = "em"; }
    if(ich==2) { suffix = "mm"; }

    ////// 2D ///////
    sprintf(pname,"mtt_weight_%s",suffix.c_str());
    CPlot plotmttw(pname,"","M_{t#bar{t}} [GeV]","weight");
    plotmttw.AddHist2D(hMttvWeight[ich][isTT2L], "COLZ");
    plotmttw.Draw(c,true,"png");

    sprintf(pname,"nu1E_weight_%s",suffix.c_str());
    CPlot plotnu1w(pname,"","#nu_{1} E [GeV]","weight");
    plotnu1w.AddHist2D(hNu1vWeight[ich][isTT2L], "COLZ");
    plotnu1w.Draw(c,true,"png");

    sprintf(pname,"nu2E_weight_%s",suffix.c_str());
    CPlot plotnu2w(pname,"","#nu_{2} E [GeV]","weight");
    plotnu2w.AddHist2D(hNu2vWeight[ich][isTT2L], "COLZ");
    plotnu2w.Draw(c,true,"png");

    sprintf(pname, "nu1_nu2_E_%s",suffix.c_str());
    CPlot plotNuE(pname,"","#bar{#nu} E", "#nu E");
    plotNuE.AddHist2D(hNuE2D[ich][isTT2L],"COLZ");
    plotNuE.Draw(c,true,"png");
    /////////////////////////////////
 
    sprintf(pname,"lep1_pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hLep1Pt[ich][0]->GetBinWidth(1)));
    CPlot plotlep1pt(pname,"","lep1 p_{T} [GeV]",ylabel);
    plotlep1pt.AddHist1D(hLep1Pt[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotlep1pt.Draw(c,true,"png");

    sprintf(pname,"lep2_pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hLep2Pt[ich][0]->GetBinWidth(1)));
    CPlot plotlep2pt(pname,"","lep2 p_{T} [GeV]",ylabel);
    plotlep2pt.AddHist1D(hLep2Pt[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotlep2pt.Draw(c,true,"png");

    sprintf(pname,"B1_pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hB1Pt[ich][0]->GetBinWidth(1)));
    CPlot plotB1pt(pname,"","B1 p_{T} [GeV]",ylabel);
    plotB1pt.AddHist1D(hB1Pt[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotB1pt.Draw(c,true,"png");

    sprintf(pname,"B2_pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hB2Pt[ich][0]->GetBinWidth(1)));
    CPlot plotB2pt(pname,"","B2 p_{T} [GeV]",ylabel);
    plotB2pt.AddHist1D(hB2Pt[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotB2pt.Draw(c,true,"png");


    sprintf(pname,"GenTop1_m_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hMTop1Gen[ich][0]->GetBinWidth(1)));
    CPlot plotGenTop1M(pname,"","GenTop1 Mass [GeV]",ylabel);
    plotGenTop1M.AddHist1D(hMTop1Gen[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotGenTop1M.Draw(c,true,"png");


    sprintf(pname,"GenTop2_m_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i GeV",int(hMTop2Gen[ich][0]->GetBinWidth(1)));
    CPlot plotGenTop2M(pname,"","GenTop2 Mass [GeV]",ylabel);
    plotGenTop2M.AddHist1D(hMTop2Gen[ich][isTT2L], samplev[isTT2L]->label,"histE", samplev[isTT2L]->linecolor);
    plotGenTop2M.Draw(c,true,"png");
      
    sprintf(pname,"num_solutions_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hNSol[ich][0]->GetBinWidth(1)));
    CPlot plotNSol(pname,"","N Solutions", ylabel);
    plotNSol.AddHist1D(hNSol[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotNSol.Draw(c,true,"png");

    sprintf(pname,"num_pass_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hNPass[ich][0]->GetBinWidth(1)));
    CPlot plotNPass(pname,"","N Pass", ylabel);
    plotNPass.AddHist1D(hNPass[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotNPass.Draw(c,true,"png");

    sprintf(pname,"weight_solutions_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hSolWeight[ich][0]->GetBinWidth(1)));
    CPlot plotSolW(pname,"","weight", ylabel);
    plotSolW.AddHist1D(hSolWeight[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotSolW.Draw(c,true,"png");

    sprintf(ylabel,"a.u.");
    sprintf(pname,"dphimetb_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hdPhiMETB[ich][0]->GetBinWidth(1)));
    CPlot plotdPhiMETB(pname,"","#Delta#phi(MET,B)", ylabel);
    plotdPhiMETB.AddHist1D(hdPhiMETB[ich][isTT2L] ,samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotdPhiMETB.AddHist1D(hdPhiMETB[ich][isPS100],samplev[isPS100]->label,"histE", samplev[isPS100]->linecolor,1,0,3);
    plotdPhiMETB.AddHist1D(hdPhiMETB[ich][isS300] ,samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotdPhiMETB.Draw(c,true,"png");

    sprintf(pname,"fitted_met_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETfitted[ich][0]->GetBinWidth(1)));
    CPlot plotMetfitted(pname,"","MET", ylabel);
    plotMetfitted.AddHist1D(hMETfitted[ich][isTT2L] ,samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMetfitted.AddHist1D(hMETfitted[ich][isPS100],samplev[isPS100]->label,"histE", samplev[isPS100]->linecolor,1,0,3);
    plotMetfitted.AddHist1D(hMETfitted[ich][isS300] ,samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMetfitted.Draw(c,true,"png");

    sprintf(pname,"met_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMET[ich][0]->GetBinWidth(1)));
    CPlot plotMet(pname,"","MET", ylabel);
    plotMet.AddHist1D(hMET[ich][isTT2L] ,samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMet.AddHist1D(hMET[ich][isPS100],samplev[isPS100]->label,"histE", samplev[isPS100]->linecolor,1,0,3);
    plotMet.AddHist1D(hMET[ich][isS300] ,samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMet.Draw(c,true,"png");

    sprintf(pname,"met_mt2cat_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMET[ich][0]->GetBinWidth(1)));
    CPlot plotMetMt2(pname,"","MET", ylabel);
    plotMetMt2.AddHist1D(hMET_mt2[ich][isTT2L] ,samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMetMt2.AddHist1D(hMET_mt2[ich][isPS100],samplev[isPS100]->label,"histE", samplev[isPS100]->linecolor,1,0,3);
    plotMetMt2.AddHist1D(hMET_mt2[ich][isS300] ,samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMetMt2.Draw(c,true,"png");

    sprintf(pname,"metshift_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETshift[ich][0]->GetBinWidth(1)));
    CPlot plotMetShift(pname,"","MET error", ylabel);
    plotMetShift.AddHist1D(hMETshift[ich][isTT2L] ,samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMetShift.AddHist1D(hMETshift[ich][isPS100],samplev[isPS100]->label,"histE", samplev[isPS100]->linecolor,1,0,3);
    plotMetShift.AddHist1D(hMETshift[ich][isS300] ,samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMetShift.Draw(c,true,"png");

    sprintf(pname,"mt2ll_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hTotErr[ich][0]->GetBinWidth(1)));
    CPlot plotMT2ll(pname,"","M_{T2}(ll)", ylabel);
    plotMT2ll.AddHist1D(hMT2ll[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMT2ll.AddHist1D(hMT2ll[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMT2ll.AddHist1D(hMT2ll[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMT2ll.Draw(c,true,"png");


    sprintf(pname,"toterr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hTotErr[ich][0]->GetBinWidth(1)));
    CPlot plotTotErr(pname,"","total error", ylabel);
    plotTotErr.AddHist1D(hTotErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotTotErr.AddHist1D(hTotErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotTotErr.AddHist1D(hTotErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotTotErr.Draw(c,true,"png");

    sprintf(pname,"reserr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hResErr[ich][0]->GetBinWidth(1)));
    CPlot plotResErr(pname,"","res error", ylabel);
    plotResErr.AddHist1D(hResErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotResErr.AddHist1D(hResErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotResErr.AddHist1D(hResErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotResErr.Draw(c,true,"png");


    sprintf(pname,"toterrlog_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hTotErr[ich][0]->GetBinWidth(1)));
    CPlot plotTotErrLog(pname,"","total error",ylabel);
    plotTotErrLog.AddHist1D(hTotErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotTotErrLog.AddHist1D(hTotErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotTotErrLog.AddHist1D(hTotErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotTotErrLog.SetLogy();
    plotTotErrLog.Draw(c,true,"png");

    sprintf(pname,"b1pterr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hB1PtErr[ich][0]->GetBinWidth(1)));
    CPlot plotB1PtErr(pname,"","b1 p_{T} error", ylabel);
    plotB1PtErr.AddHist1D(hB1PtErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotB1PtErr.AddHist1D(hB1PtErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotB1PtErr.AddHist1D(hB1PtErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotB1PtErr.Draw(c,true,"png");

    sprintf(pname,"b1pterrlog_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hB1PtErr[ich][0]->GetBinWidth(1)));
    CPlot plotB1PtErrLog(pname,"","b1 p_{T} error", ylabel);
    plotB1PtErrLog.AddHist1D(hB1PtErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotB1PtErrLog.AddHist1D(hB1PtErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotB1PtErrLog.AddHist1D(hB1PtErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotB1PtErrLog.SetLogy();
    plotB1PtErrLog.Draw(c,true,"png");

    sprintf(pname,"b2pterr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hB2PtErr[ich][0]->GetBinWidth(1)));
    CPlot plotB2PtErr(pname,"","b2 p_{T} error", ylabel);
    plotB2PtErr.AddHist1D(hB2PtErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotB2PtErr.AddHist1D(hB2PtErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotB2PtErr.AddHist1D(hB2PtErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotB2PtErr.Draw(c,true,"png");

    sprintf(pname,"b2pterrlog_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hB2PtErr[ich][0]->GetBinWidth(1)));
    CPlot plotB2PtErrLog(pname,"","b2 p_{T} error", ylabel);
    plotB2PtErrLog.AddHist1D(hB2PtErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotB2PtErrLog.AddHist1D(hB2PtErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotB2PtErrLog.AddHist1D(hB2PtErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotB2PtErrLog.SetLogy();
    plotB2PtErrLog.Draw(c,true,"png");

    sprintf(pname,"metphierr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETphiErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETphiErr(pname,"","E^{miss}_{T} #phi error", ylabel);
    plotMETphiErr.AddHist1D(hMETphiErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETphiErr.AddHist1D(hMETphiErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETphiErr.AddHist1D(hMETphiErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETphiErr.Draw(c,true,"png");

    sprintf(pname,"metphierrlog_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETphiErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETphiErrlog(pname,"","E^{miss}_{T} #phi error", ylabel);
    plotMETphiErrlog.AddHist1D(hMETphiErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETphiErrlog.AddHist1D(hMETphiErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETphiErrlog.AddHist1D(hMETphiErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETphiErrlog.SetLogy();
    plotMETphiErrlog.Draw(c,true,"png");

    sprintf(pname,"metxerr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETxErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETxErr(pname,"","E^{miss}_{x} error", ylabel);
    plotMETxErr.AddHist1D(hMETxErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETxErr.AddHist1D(hMETxErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETxErr.AddHist1D(hMETxErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETxErr.Draw(c,true,"png");

    sprintf(pname,"metxerrlog_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETxErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETxErrLog(pname,"","E^{miss}_{x} error", ylabel);
    plotMETxErrLog.AddHist1D(hMETxErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETxErrLog.AddHist1D(hMETxErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETxErrLog.AddHist1D(hMETxErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETxErrLog.SetLogy();
    plotMETxErrLog.Draw(c,true,"png");

    sprintf(pname,"metyerr_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETyErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETyErr(pname,"","E^{miss}_{y} error", ylabel);
    plotMETyErr.AddHist1D(hMETyErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETyErr.AddHist1D(hMETyErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETyErr.AddHist1D(hMETyErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETyErr.Draw(c,true,"png");

    sprintf(pname,"metyerrlogz_%s",suffix.c_str());
    // sprintf(ylabel, "Events /%i", int(hMETyErr[ich][0]->GetBinWidth(1)));
    CPlot plotMETyErrLog(pname,"","E^{miss}_{y} error", ylabel);
    plotMETyErrLog.AddHist1D(hMETyErr[ich][isTT2L],samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMETyErrLog.AddHist1D(hMETyErr[ich][isPS100],samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMETyErrLog.AddHist1D(hMETyErr[ich][isS300],samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMETyErrLog.SetLogy();
    plotMETyErrLog.Draw(c,true,"png");

    sprintf(pname,"mtt_%s",suffix.c_str());
    // sprintf(ylabel, "Events / %i", int(hMtt[ich][0]->GetBinWidth(1)));
    CPlot plotMtt(pname,"","M_{t#bar{t}} [GeV]",ylabel);
    plotMtt.AddHist1D(hMtt[ich][isTT2L], samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMtt.AddHist1D(hMtt[ich][isPS100], samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMtt.AddHist1D(hMtt[ich][isS300], samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMtt.Draw(c,true,"png");
    
    sprintf(pname,"mlb_%s",suffix.c_str());
    // sprintf(ylabel, "Events / %i", int(hMlb[ich][0]->GetBinWidth(1)));
    CPlot plotMlb(pname,"","M_{(l,b)} [GeV]",ylabel);
    plotMlb.AddHist1D(hMlb[ich][isTT2L], samplev[isTT2L]->label, "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMlb.AddHist1D(hMlb[ich][isPS100], samplev[isPS100]->label, "histE", samplev[isPS100]->linecolor,1,0,3);
    plotMlb.AddHist1D(hMlb[ich][isS300], samplev[isS300]->label, "histE", samplev[isS300]->linecolor,1,0,3);
    plotMlb.Draw(c,true,"png");

    sprintf(pname,"toppt_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hTopPtReco[ich][0]->GetBinWidth(1)));
    CPlot plotTopPt(pname,"","top p_{T} [GeV]", ylabel);
    // plotTopPt.AddHist1D(hTopPtReco[ich][isPS100], "S 1 100 signal smear","histE", samplev[isPS100]->linecolor,1,0,3);
    plotTopPt.AddHist1D(hTopPtReco[ich][isTT2L],"t#bar{t} reco", "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotTopPt.AddHist1D(hTopPtGen[ich][isTT2L],"t#bar{t} gen", "histE", kMagenta+2,1,0,3);
    plotTopPt.Draw(c,true,"png");
    
    sprintf(pname,"mtop1_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hMTop1Reco[ich][0]->GetBinWidth(1)));
    CPlot plotMTop1(pname,"","M_{#bar{t}} [GeV]", ylabel);
    // plotMTop1.AddHist1D(hMTop1Reco[ich][isPS100], "S 1 100 signal smear","histE", samplev[isPS100]->linecolor,1,0,3);
    plotMTop1.AddHist1D(hMTop1Reco[ich][isTT2L],"t#bar{t} reco", "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMTop1.AddHist1D(hMTop1Gen[ich][isTT2L],"t#bar{t} gen", "histE", kMagenta+2,1,0,3);
    plotMTop1.Draw(c,true,"png");

    sprintf(pname,"mtop2_%s",suffix.c_str());
    sprintf(ylabel, "Events /%i", int(hMTop2Reco[ich][0]->GetBinWidth(1)));
    CPlot plotMTop2(pname,"","M_{t} [GeV]", ylabel);
    // plotMTop2.AddHist1D(hMTop2Reco[ich][isPS100], "S 1 100 signal smear","histE", samplev[isPS100]->linecolor,1,0,3);
    plotMTop2.AddHist1D(hMTop2Reco[ich][isTT2L],"t#bar{t} reco", "histE", samplev[isTT2L]->linecolor,1,0,3);
    plotMTop2.AddHist1D(hMTop2Gen[ich][isTT2L],"t#bar{t} gen", "histE", kMagenta+2,1,0,3);
    plotMTop2.Draw(c,true,"png");

    sprintf(pname,"nu1pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i", int(hNu1Pt[ich][0]->GetBinWidth(1)));
    CPlot plotNu1Pt(pname, "", "#bar{#nu} p_{T} [GeV]",ylabel);
    plotNu1Pt.AddHist1D(hNu1Pt[ich][isTT2L],"ttbar smear","histE",samplev[isTT2L]->linecolor,1,0,3);
    plotNu1Pt.AddHist1D(hNu1GenPt[ich][isTT2L],"ttbar gen","histE",kPink+7,1,0,3);
    plotNu1Pt.Draw(c,true,"png");

    sprintf(pname,"nu2pt_%s",suffix.c_str());
    sprintf(ylabel,"Events / %i", int(hNu2Pt[ich][0]->GetBinWidth(1)));
    CPlot plotNu2Pt(pname, "", "#nu p_{T} [GeV]",ylabel);
    plotNu2Pt.AddHist1D(hNu2Pt[ich][isTT2L],"ttbar smear","histE",samplev[isTT2L]->linecolor,1,0,3);
    plotNu2Pt.AddHist1D(hNu2GenPt[ich][isTT2L],"ttbar gen","histE",kPink+7,1,0,3);
    plotNu2Pt.Draw(c,true,"png");

    
  }
  
  for(unsigned int isol=0; isol<3; isol++){
    for(unsigned int isam=0; isam<samplev.size(); isam++){
      hLep1Pt[isol][isam]->Scale(1./hLep1Pt[isol][isam]->Integral());
      hLep2Pt[isol][isam]->Scale(1./hLep2Pt[isol][isam]->Integral());
      hB1Pt[isol][isam]->Scale(1./hB1Pt[isol][isam]->Integral());
      hB2Pt[isol][isam]->Scale(1./hB2Pt[isol][isam]->Integral());
    }
  }
  
  sprintf(ylabel,"Events / %i GeV",int(hLep1Pt[0][0]->GetBinWidth(1)));
  CPlot plotLep1Pt("lep1_pt","","Lep1 p_{T} [GeV]",ylabel);
  plotLep1Pt.AddHist1D(hLep1Pt[0][isTT2L],"0 sol","histE",kPink+7);
  plotLep1Pt.AddHist1D(hLep1Pt[1][isTT2L],"2 sol","histE",kBlue+1);
  plotLep1Pt.AddHist1D(hLep1Pt[2][isTT2L],"4 sol","histE",kGreen+2);
  plotLep1Pt.Draw(c,true,"png");

  sprintf(ylabel,"Events / %i GeV",int(hLep2Pt[0][0]->GetBinWidth(1)));
  CPlot plotLep2Pt("lep2_pt","","Lep2 p_{T} [GeV]",ylabel);
  plotLep2Pt.AddHist1D(hLep2Pt[0][isTT2L],"0 sol","histE",kPink+7);
  plotLep2Pt.AddHist1D(hLep2Pt[1][isTT2L],"2 sol","histE",kBlue+1);
  plotLep2Pt.AddHist1D(hLep2Pt[2][isTT2L],"4 sol","histE",kGreen+2);
  plotLep2Pt.Draw(c,true,"png");

  sprintf(ylabel,"Events / %i GeV",int(hB1Pt[0][0]->GetBinWidth(1)));
  CPlot plotB1Pt("b1_pt","","B1 p_{T} [GeV]",ylabel);
  plotB1Pt.AddHist1D(hB1Pt[0][isTT2L],"0 sol","histE",kPink+7);
  plotB1Pt.AddHist1D(hB1Pt[1][isTT2L],"2 sol","histE",kBlue+1);
  plotB1Pt.AddHist1D(hB1Pt[2][isTT2L],"4 sol","histE",kGreen+2);
  plotB1Pt.Draw(c,true,"png");

  sprintf(ylabel,"Events / %i GeV",int(hB2Pt[0][0]->GetBinWidth(1)));
  CPlot plotB2Pt("b2_pt","","B2 p_{T} [GeV]",ylabel);
  plotB2Pt.AddHist1D(hB2Pt[0][isTT2L],"0 sol","histE",kPink+7);
  plotB2Pt.AddHist1D(hB2Pt[1][isTT2L],"2 sol","histE",kBlue+1);
  plotB2Pt.AddHist1D(hB2Pt[2][isTT2L],"4 sol","histE",kGreen+2);
  plotB2Pt.Draw(c,true,"png");
  
  
  
}
