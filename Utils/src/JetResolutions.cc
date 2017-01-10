#include "../interface/JetResolutions.hh"
//----------------------------------------------------------------------------------
double JetResolutions::getPtSigmaSmear(TLorentzVector & vJet, const double JetRho){
  double sigma_pt;

  double JetEta = vJet.Eta();
  std::vector<double> rholow; std::vector<double> rhohigh;
  rholow.push_back(0.0); rhohigh.push_back(6.73);
  rholow.push_back(6.73); rhohigh.push_back(12.51);
  rholow.push_back(12.51); rhohigh.push_back(18.29);
  rholow.push_back(18.29); rhohigh.push_back(24.08);
  rholow.push_back(24.08); rhohigh.push_back(44.31);

  std::vector<double> etalow; std::vector<double> etahigh;
  etahigh.push_back(4.7); etalow.push_back(3.2);
  etahigh.push_back(3.2); etalow.push_back(3.0);
  etahigh.push_back(3.0); etalow.push_back(2.8);
  etahigh.push_back(2.8); etalow.push_back(2.5);
  etahigh.push_back(2.5); etalow.push_back(2.3);
  etahigh.push_back(2.3); etalow.push_back(2.1);
  etahigh.push_back(2.1); etalow.push_back(1.9);
  etahigh.push_back(1.9); etalow.push_back(1.7);
  etahigh.push_back(1.7); etalow.push_back(1.3);
  etahigh.push_back(1.3); etalow.push_back(1.1);
  etahigh.push_back(1.1); etalow.push_back(0.8);
  etahigh.push_back(0.8); etalow.push_back(0.5);
  etahigh.push_back(0.5); etalow.push_back(0.0);

  int eta_bin = -1; int rho_bin=-1;
  for(unsigned int i=0; i<etahigh.size(); i++){
    if(fabs(JetEta) < etahigh[i] && fabs(JetEta) > etalow[i]) eta_bin = i;
    else eta_bin = 0;
  }
  for(unsigned int i=0; i<rhohigh.size(); i++){
    if(JetRho > rholow[i] && JetRho < rhohigh[i]) rho_bin = i;
    else rho_bin = 0;
  }

  fPtSmear->SetParameter(0,fPtPar0[eta_bin][rho_bin]);
  fPtSmear->SetParameter(1,fPtPar1[eta_bin][rho_bin]);
  fPtSmear->SetParameter(2,fPtPar2[eta_bin][rho_bin]);
  fPtSmear->SetParameter(3,fPtPar3[eta_bin][rho_bin]);

  fPhiSmear->SetParameter(0,fPhiPar0[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(1,fPhiPar1[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(2,fPhiPar2[eta_bin][rho_bin]);

  sigma_pt = fPtSmear->Eval(vJet.Pt());
  return sigma_pt;

}
//----------------------------------------------------------------------------------
double JetResolutions::getPhiSigmaSmear(TLorentzVector & vJet, const double JetRho){
  double sigma_phi;
  double JetEta = vJet.Eta();
  std::vector<double> rholow; std::vector<double> rhohigh;
  rholow.push_back(0.0); rhohigh.push_back(6.73);
  rholow.push_back(6.73); rhohigh.push_back(12.51);
  rholow.push_back(12.51); rhohigh.push_back(18.29);
  rholow.push_back(18.29); rhohigh.push_back(24.08);
  rholow.push_back(24.08); rhohigh.push_back(44.31);

  std::vector<double> etalow; std::vector<double> etahigh;
  etahigh.push_back(4.7); etalow.push_back(3.2);
  etahigh.push_back(3.2); etalow.push_back(3.0);
  etahigh.push_back(3.0); etalow.push_back(2.8);
  etahigh.push_back(2.8); etalow.push_back(2.5);
  etahigh.push_back(2.5); etalow.push_back(2.3);
  etahigh.push_back(2.3); etalow.push_back(2.1);
  etahigh.push_back(2.1); etalow.push_back(1.9);
  etahigh.push_back(1.9); etalow.push_back(1.7);
  etahigh.push_back(1.7); etalow.push_back(1.3);
  etahigh.push_back(1.3); etalow.push_back(1.1);
  etahigh.push_back(1.1); etalow.push_back(0.8);
  etahigh.push_back(0.8); etalow.push_back(0.5);
  etahigh.push_back(0.5); etalow.push_back(0.0);

  int eta_bin = -1; int rho_bin=-1;
  for(unsigned int i=0; i<etahigh.size(); i++){
    if(fabs(JetEta) < etahigh[i] && fabs(JetEta) > etalow[i]) eta_bin = i;
    else eta_bin = 0;
  }
  for(unsigned int i=0; i<rhohigh.size(); i++){
    if(JetRho > rholow[i] && JetRho < rhohigh[i]) rho_bin = i;
    else rho_bin = 0;
  }

  fPtSmear->SetParameter(0,fPtPar0[eta_bin][rho_bin]);
  fPtSmear->SetParameter(1,fPtPar1[eta_bin][rho_bin]);
  fPtSmear->SetParameter(2,fPtPar2[eta_bin][rho_bin]);
  fPtSmear->SetParameter(3,fPtPar3[eta_bin][rho_bin]);

  fPhiSmear->SetParameter(0,fPhiPar0[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(1,fPhiPar1[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(2,fPhiPar2[eta_bin][rho_bin]);

  sigma_phi = fPhiSmear->Eval(vJet.Pt());
  return sigma_phi;
}
//----------------------------------------------------------------------------------
void JetResolutions::setSmear(TLorentzVector & vJet, const double JetRho, TLorentzVector vSmear) {

  double JetEta = vJet.Eta();
  std::vector<double> rholow; std::vector<double> rhohigh;
  rholow.push_back(0.0); rhohigh.push_back(6.73);
  rholow.push_back(6.73); rhohigh.push_back(12.51);
  rholow.push_back(12.51); rhohigh.push_back(18.29);
  rholow.push_back(18.29); rhohigh.push_back(24.08);
  rholow.push_back(24.08); rhohigh.push_back(44.31);

  std::vector<double> etalow; std::vector<double> etahigh;
  etahigh.push_back(4.7); etalow.push_back(3.2);
  etahigh.push_back(3.2); etalow.push_back(3.0);
  etahigh.push_back(3.0); etalow.push_back(2.8);
  etahigh.push_back(2.8); etalow.push_back(2.5);
  etahigh.push_back(2.5); etalow.push_back(2.3);
  etahigh.push_back(2.3); etalow.push_back(2.1);
  etahigh.push_back(2.1); etalow.push_back(1.9);
  etahigh.push_back(1.9); etalow.push_back(1.7);
  etahigh.push_back(1.7); etalow.push_back(1.3);
  etahigh.push_back(1.3); etalow.push_back(1.1);
  etahigh.push_back(1.1); etalow.push_back(0.8);
  etahigh.push_back(0.8); etalow.push_back(0.5);
  etahigh.push_back(0.5); etalow.push_back(0.0);

  int eta_bin = -1; int rho_bin=-1;
  for(unsigned int i=0; i<etahigh.size(); i++){
    if(fabs(JetEta) < etahigh[i] && fabs(JetEta) > etalow[i]) eta_bin = i;
    else eta_bin = 0;
  }
  for(unsigned int i=0; i<rhohigh.size(); i++){
    if(JetRho > rholow[i] && JetRho < rhohigh[i]) rho_bin = i;
    else rho_bin = 0;
  }

  fPtSmear->SetParameter(0,fPtPar0[eta_bin][rho_bin]);
  fPtSmear->SetParameter(1,fPtPar1[eta_bin][rho_bin]);
  fPtSmear->SetParameter(2,fPtPar2[eta_bin][rho_bin]);
  fPtSmear->SetParameter(3,fPtPar3[eta_bin][rho_bin]);

  fPhiSmear->SetParameter(0,fPhiPar0[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(1,fPhiPar1[eta_bin][rho_bin]);
  fPhiSmear->SetParameter(2,fPhiPar2[eta_bin][rho_bin]);

  double unc_pt, unc_phi;
  unc_pt = fPtSmear->Eval(vJet.Pt())*vJet.Pt();
  unc_phi = fPhiSmear->Eval(vJet.Pt());
  vSmear.SetPtEtaPhiM(unc_pt, 0, unc_phi, 0);
  
  // fEtaSmear->SetParameter(0,fEtaPar0[eta_bin][rho_bin]);
  // fEtaSmear->SetParameter(1,fEtaPar1[eta_bin][rho_bin]);
  // fEtaSmear->SetParameter(2,fEtaPar2[eta_bin][rho_bin]);
  // fEtaSmear->SetParameter(3,fEtaPar3[eta_bin][rho_bin]);
  // fEtaSmear->SetParameter(4,fEtaPar4[eta_bin][rho_bin]);
};
//----------------------------------------------------------------------------------
void JetResolutions::getUncertainties( TLorentzVector & v, const double rho, double & upx, double & upy, double & um ) //double & upz
{
  double pt  = v.Pt();
  // double eta = v.Eta();
  double phi = v.Phi();
  /*
  std::vector<double> rholow; std::vector<double> rhohigh;
  rholow.push_back(0.0); rhohigh.push_back(6.73);
  rholow.push_back(6.73); rhohigh.push_back(12.51);
  rholow.push_back(12.51); rhohigh.push_back(18.29);
  rholow.push_back(18.29); rhohigh.push_back(24.08);
  rholow.push_back(24.08); rhohigh.push_back(44.31);

  std::vector<double> etalow; std::vector<double> etahigh;
  etahigh.push_back(4.7); etalow.push_back(3.2);
  etahigh.push_back(3.2); etalow.push_back(3.0);
  etahigh.push_back(3.0); etalow.push_back(2.8);
  etahigh.push_back(2.8); etalow.push_back(2.5);
  etahigh.push_back(2.5); etalow.push_back(2.3);
  etahigh.push_back(2.3); etalow.push_back(2.1);
  etahigh.push_back(2.1); etalow.push_back(1.9);
  etahigh.push_back(1.9); etalow.push_back(1.7);
  etahigh.push_back(1.7); etalow.push_back(1.3);
  etahigh.push_back(1.3); etalow.push_back(1.1);
  etahigh.push_back(1.1); etalow.push_back(0.8);
  etahigh.push_back(0.8); etalow.push_back(0.5);
  etahigh.push_back(0.5); etalow.push_back(0.0);
  int eta_bin = -1; int rho_bin=-1;
  for(unsigned int i=0; i<etahigh.size(); i++){
    if(fabs(eta) < etahigh[i] && fabs(eta) > etalow[i]) eta_bin = i;
    else eta_bin = 0;
  }
  for(unsigned int i=0; i<rhohigh.size(); i++){
    if(rho > rholow[i] && rho < rhohigh[i]) rho_bin = i;
    else rho_bin = 0;
  }
  */
  TLorentzVector vsmear;
  setSmear(v,rho,vsmear);

  /*
    int eta_bin = 9-int(fabs(eta*2));
    if( eta_bin < 0 ) eta_bin = 0.;
    setSmear(eta_bin);
  */
  double uncert_pt  = fPtSmear->Eval(pt)*pt;
  // double uncert_eta = fEtaSmear->Eval(pt);
  double uncert_phi = fPhiSmear->Eval(pt);
  //  double uncert_m = uncert_pt;
  
  upx 
    = sqrt(TMath::Cos(phi)*uncert_pt*TMath::Cos(phi)*uncert_pt 
	   + TMath::Sin(phi)*pt*uncert_phi*TMath::Sin(phi)*pt*uncert_phi);
  upy 
    = sqrt(TMath::Sin(phi)*uncert_pt *TMath::Sin(phi)*uncert_pt 
	   + TMath::Cos(phi)*pt*uncert_phi*TMath::Cos(phi)*pt*uncert_phi);
  // upz 
  //   = sqrt(TMath::SinH(eta)*uncert_pt*TMath::SinH(eta)*uncert_pt 
  // 	   + TMath::CosH(eta)*pt*uncert_eta*TMath::CosH(eta)*pt*uncert_eta);
  um  = uncert_pt;
  //um  = sqrt(uncert_pt*uncert_pt + upz*upz);
  
  float inflation=1.0;
  upx *= inflation;
  upy *= inflation;
  // upz *= inflation;
  um *= inflation;
};

void JetResolutions::covMatrix(TMatrixD &iCov, double iGenPt, double iGenPhi, int njet) {
  std::vector<TF1*> fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit;
  readRecoil(fF1U1Fit,fF1U1RMSSMFit,fF1U1RMS1Fit,fF1U1RMS2Fit,fF1U2Fit,fF1U2RMSSMFit,fF1U2RMS1Fit,fF1U2RMS2Fit);
  // const int fJet;
  // fJet = njet; if(njet > 2) fJet = 2;  
  // if(fJet >= int(fF1U1Fit.size())) fJet = 0; 
  double lSigma1 = fF1U1RMSSMFit[0]->Eval(iGenPt);
  double lSigma2 = fF1U2RMSSMFit[0]->Eval(iGenPt);
  double lCovU1  = lSigma1*lSigma1;
  double lCovU2  = lSigma2*lSigma2;
  double lUPhi   = iGenPhi;
  TMatrixD lRotate(2,2);
  lRotate(0,0) = cos(lUPhi);   lRotate(0,1) = -sin(lUPhi);
  lRotate(1,1) = cos(lUPhi);   lRotate(1,0) =  sin(lUPhi);
  TMatrixD lCov(2,2);
  lCov(0,0) = lCovU1;
  lCov(1,1) = lCovU2;
  TMatrixD lIRotate = lRotate.Invert();
  // TMatrixD lM = lRotate*lCov*lIRotate;
  TMatrixD lM = lCov; 
  (iCov)(0,0)   =  lM(0,0);
  (iCov)(1,0)   =  lM(1,0);
  (iCov)(0,1)   =  lM(0,1);
  (iCov)(1,1)   =  lM(1,1);
}

void JetResolutions::getSmearedMET(const TLorentzVector* iU, const TMatrixD* CovMat, TLorentzVector* oU,
				   double & sigmaX, double & sigmaY){
  /*
  float x = (*CovMat)(0,0);
  float y = (*CovMat)(1,0);
  float z = (*CovMat)(2,0);
  float E = TMath::Sqrt(x*x + y*y + z*z);
  TLorentzVector* sMET = new TLorentzVector(x, y, z, E);
  return sMET;
  */
  double lUPhi = iU->Phi();

  TVector2 tmpU;
  tmpU.SetX(iU->X());
  tmpU.SetY(iU->Y());
  tmpU.Rotate(-lUPhi);
  
  double metX = myRand->Gaus(tmpU.X(), sqrt((*CovMat)(0,0)));
  double metY = myRand->Gaus(tmpU.Y(), sqrt((*CovMat)(1,1)));
  
  tmpU.SetX(metX);
  tmpU.SetY(metY);
  tmpU.Rotate(lUPhi);
  
  oU->SetX(metX);
  oU->SetY(metY);  

}

void JetResolutions::readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,
                                 std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit) {
  //  if(!getenv("CMSSW_BASE")) {
  //    printf("error! RecoilCorrector called without input files. Define CMSSW_BASE or add by hand.\n");
  //    assert(0);
  //  }
  TFile *lFile  = new TFile("/tthome/ssevova/Analysis/11/CMSSW_8_0_20/src/DMSAna/Utils/data/recoilfit_zmmMC_pf_v2.root");
  int lNJet = 0;
  std::string iPrefix = "PF";
  std::stringstream lSS; lSS << iPrefix << "u1Mean_" << lNJet;
  while(lFile->FindObjectAny(lSS.str().c_str()) != 0) { lSS.str("");
    lSS << iPrefix << "u1Mean_"    << lNJet; iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lSS << iPrefix << "u1MeanRMS_" << lNJet; iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
    lSS << iPrefix << "u1RMS1_"    << lNJet; iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
    lSS << iPrefix << "u1RMS2_"    << lNJet; iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str(""); 
    lSS << iPrefix << "u2Mean_"    << lNJet; iU2Fit    .push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lSS << iPrefix << "u2MeanRMS_" << lNJet; iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lSS << iPrefix << "u2RMS1_"    << lNJet; iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lSS << iPrefix << "u2RMS2_"    << lNJet; iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lSS << iPrefix << "u2RMS2_"    << lNJet; iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(lSS.str().c_str())); lSS.str("");
    lNJet++; lSS << iPrefix << "u1Mean_" << lNJet;
  }
  lFile->Close();
}

