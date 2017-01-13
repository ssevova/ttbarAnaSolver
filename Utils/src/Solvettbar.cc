/* 

Analytical solution of ttbar dilepton equations 

Theory: https://arxiv.org/pdf/hep-ph/0603011v3.pdf Author: Lars Sonnenschein

Get roots of quartic equation as solutions to p_nu_x //neutrino x momentum

Constraints
-----------
W+/- and t tbar masses
Neutrino energy/mass relation

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TF2.h>
#include <TF1.h>

#include "../interface/JetResolutions.hh"
#include "../interface/FitSolution.hh"
#include "../interface/Solvettbar.hh"
#include "../interface/ttbarCandidate.hh"


Solvettbar::Solvettbar():
  topmass_begin(0), 
  topmass_end(0),
  topmass_step(0),
  mw1(80.4),
  mw2(80.4),
  mb(4.8),
  pxmiss(0),
  pymiss(0)
{}

Solvettbar::Solvettbar(const float b, const float e, const float s, const double mW, const double mB):
  topmass_begin(b),
  topmass_end(e), 
  topmass_step(s),
  mw1(mW),
  mw2(mW),
  mb(mB),
  pxmiss(0),
  pymiss(0)
{
  EventShape_ = new TF2("landau2D","[0]*TMath::Landau(x,[1],[2],0)*TMath::Landau(y,[3],[4],0)",0,500,0,500);
  EventShape_ ->SetParameters(30.7137,56.2880,23.0744,59.1015,24.9145);//nupars[0], nupars[1], nupars[2], nupars[3], nupars[4]);
}

//
//destructor
//
Solvettbar::~Solvettbar(){
  delete EventShape_;
}

FitSolution Solvettbar::SolveNu( FitSolution & fitsol, ttbarCandidate & ttbarcand){

 // ttbarcand.reset();

  TLorentzVector l1vec = ttbarcand.lep1vec;
  TLorentzVector l2vec = ttbarcand.lep2vec;
  TLorentzVector b1vec = ttbarcand.B1vec;
  TLorentzVector b2vec = ttbarcand.B2vec;
  
  // double discriminant = -999;

  TLorentzVector n1vec, n2vec;

  double q_coeff[5], q_sol[4]; double q_disc[1];

  FindCoeff(l1vec, l2vec, b1vec, b2vec, genTop1Mass, genTop2Mass, q_coeff);

  //get the number of solutions taking into account the coefficients
  int Nsol = quartic(q_coeff, q_sol, q_disc);

  fitsol.numSol   = Nsol;  
  //loop over all solutions
  for(int isol = 0; isol < Nsol; isol++){
    //reconstruct ttbar and get neutrinos
    TopRec(l1vec, l2vec, b1vec, b2vec, q_sol[isol]);

    n1vec.SetPtEtaPhiM(neu1vec.Pt(), neu1vec.Eta(), neu1vec.Phi(), neu1vec.M());
    n2vec.SetPtEtaPhiM(neu2vec.Pt(), neu2vec.Eta(), neu2vec.Phi(), neu2vec.M());
    
    double weight = WeightSolfromMC();//WeightSolfromShape();//   
    fitsol.weight.push_back(weight);
    fitsol.cudisc.push_back(q_disc[0]);
    fitsol.nu1.push_back(n1vec); 
    fitsol.nu2.push_back(n2vec);  
        
    double mtt = (tt_top1vec + tt_top2vec).M();
    fitsol.mass1.push_back(tt_top1vec.M());
    fitsol.mass2.push_back(tt_top2vec.M());
    fitsol.mtt.push_back(mtt);
    
  }

  return fitsol;
}

//
// Set the x and y momentum constraints: Emiss_x, Emiss_y
//
void Solvettbar::SetConstraints(const double xx, const double yy)
{
  pxmiss = xx;
  pymiss = yy;
  // std::cout << "pxmiss constraint: " << pxmiss << std::endl;
  // std::cout << "pymiss constraint: " << pymiss << std::endl;
}

void Solvettbar::SetGenTopMass(const double genMtop1, const double genMtop2)
{
  genTop1Mass = genMtop1; 
  genTop2Mass = genMtop2;
}

void Solvettbar::SetWMass(const double w1mass, const double w2mass){
  mw1 = w1mass;
  mw2 = w2mass;
}

void Solvettbar::FindCoeff(const TLorentzVector lep1vec,
			   const TLorentzVector lep2vec,
			   const TLorentzVector B1vec, 
			   const TLorentzVector B2vec, 
			   const double mtop1, 
			   const double mtop2,
			   double * coeff)
{
  double E;
  double c1, c2, c3;
  double d1, d2, d3;
  double k11, k21, k31, k41;
  double l11, l21, l31, l41, l51, l61;
  double k1, k2, k3, k4, k5, k6;
  double l1, l2, l3, l4, l5, l6;
  double k15, k25, k35, k45;
  
  C = - lep1vec.Px() - B1vec.Px() - lep2vec.Px() - B2vec.Px() + pxmiss;
  D = - lep1vec.Py() - B1vec.Py() - lep2vec.Py() - B2vec.Py() + pymiss;

  //// top leg ////
  E = ((sqr(mtop1)-sqr(mw1)-sqr(mb))/(2*B1vec.E()) - sqr(mw1)/(2*lep1vec.E()) - lep1vec.E() 
       + lep1vec.Px()*B1vec.Px()/B1vec.E() 
       + lep1vec.Py()*B1vec.Py()/B1vec.E() 
       + lep1vec.Pz()*B1vec.Pz()/B1vec.E());

  a1 = lep1vec.Px()/lep1vec.E() - B1vec.Px()/B1vec.E();
  a2 = lep1vec.Py()/lep1vec.E() - B1vec.Py()/B1vec.E();
  a3 = lep1vec.Pz()/lep1vec.E() - B1vec.Pz()/B1vec.E();

  G = E - a1*C - a2*D;

  c1 = sqr(lep1vec.Px()) - sqr(lep1vec.E());
  c2 = sqr(lep1vec.Py()) - sqr(lep1vec.E());
  c3 = sqr(lep1vec.Pz()) - sqr(lep1vec.E());

  b1 = lep2vec.Px()/lep2vec.E() - B2vec.Px()/B2vec.E();
  b2 = lep2vec.Py()/lep2vec.E() - B2vec.Py()/B2vec.E();
  b3 = lep2vec.Pz()/lep2vec.E() - B2vec.Pz()/B2vec.E();
  
  d1 = sqr(lep2vec.Px()) - sqr(lep2vec.E());
  d2 = sqr(lep2vec.Py()) - sqr(lep2vec.E());
  d3 = sqr(lep2vec.Pz()) - sqr(lep2vec.E());


    
  k11 = (1/sqr(lep1vec.E())*(pow(mw1,4)/4 + sqr(C)*c1 + sqr(D)*c2 + c3*sqr(G)/sqr(a3) 
			     + sqr(mw1)*(lep1vec.Px()*C + lep1vec.Py()*D + lep1vec.Pz()*G/a3)
			     + 2*lep1vec.Px()*lep1vec.Py()*C*D 
			     + 2*lep1vec.Px()*lep1vec.Pz()*C*G/a3
			     + 2*lep1vec.Py()*lep1vec.Pz()*D*G/a3)); 
  
  k21 = (1/sqr(lep1vec.E())*(- 2*C*a3*b3*c1 + 2*c3*b3*a1*G/a3
			     - sqr(mw1)*a3*b3*lep1vec.Px()
			     + sqr(mw1)*a1*b3*lep1vec.Pz()
			     - 2*lep1vec.Px()*lep1vec.Py()*D*a3*b3
			     + 2*lep1vec.Px()*lep1vec.Pz()*C*a1*b3
			     - 2*lep1vec.Px()*lep1vec.Pz()*b3*G 
			     + 2*lep1vec.Py()*lep1vec.Pz()*D*a1*b3));

  k31 = (1/sqr(lep1vec.E())*(- 2*D*a3*b3*c2 + 2*c3*b3*a2*G/a3
			     - sqr(mw1)*a3*b3*lep1vec.Py()
			     + sqr(mw1)*a2*b3*lep1vec.Pz()
			     - 2*lep1vec.Px()*lep1vec.Py()*C*a3*b3
			     + 2*lep1vec.Px()*lep1vec.Pz()*C*a2*b3
			     - 2*lep1vec.Py()*lep1vec.Pz()*b3*G
			     + 2*lep1vec.Py()*lep1vec.Pz()*D*a2*b3));
  
  k41 =  (1/sqr(lep1vec.E())*(2*c3*a1*a2*sqr(b3)
			     + 2*lep1vec.Px()*lep1vec.Py()*sqr(a3)*sqr(b3)
			     - 2*lep1vec.Px()*lep1vec.Pz()*a2*a3*sqr(b3)
			     - 2*lep1vec.Py()*lep1vec.Pz()*a1*a3*sqr(b3)));

  k51 = (1/sqr(lep1vec.E())*(c1*sqr(a3)*sqr(b3)
			     + c3*sqr(a1)*sqr(b3)
			     - 2*lep1vec.Px()*lep1vec.Pz()*a1*a3*sqr(b3)));

  k61 = (1/sqr(lep1vec.E())*(c2*sqr(a3)*sqr(b3)
			     + c3*sqr(a2)*sqr(b3)
			     - 2*lep1vec.Py()*lep1vec.Pz()*a2*a3*sqr(b3)));

  k1 = k11*k61; 
  k2 = k61*k21/k51;
  k3 = k31; 
  k4 = k41/k51;
  k5 = k61/k51;
  k6 = 1;

  //// antitop leg ////
  F = ((sqr(mtop2)-sqr(mw2)-sqr(mb))/(2*B2vec.E()) - sqr(mw2)/(2*lep2vec.E()) - lep2vec.E() 
       + lep2vec.Px()*B2vec.Px()/B2vec.E()
       + lep2vec.Py()*B2vec.Py()/B2vec.E()
       + lep2vec.Pz()*B2vec.Pz()/B2vec.E());

  l11 = (1/sqr(lep2vec.E())*(pow(mw2,4)/4
			     + d3*sqr(F)/sqr(b3)
			     + sqr(mw2)*lep2vec.Pz()*F/b3));
  
  l21 = (1/sqr(lep2vec.E())*(- 2*d3*F*a3*b1/b3
			     + sqr(mw2)*(lep2vec.Px()*a3*b3 - lep2vec.Pz()*b1*a3)
			     + 2*lep2vec.Px()*lep2vec.Pz()*F*a3));

  l31 = (1/sqr(lep2vec.E())*(- 2*d3*F*a3*b2/b3
			     + sqr(mw2)*(lep2vec.Py()*a3*b3 - lep2vec.Pz()*b2*a3)
			     + 2*lep2vec.Py()*lep2vec.Pz()*F*a3));

  l41 = (1/sqr(lep2vec.E())*(2*d3*b1*b2*sqr(a3)
			     + 2*lep2vec.Px()*lep2vec.Py()*sqr(a3)*sqr(b3)
			     - 2*lep2vec.Px()*lep2vec.Pz()*b2*b3*sqr(a3)
			     - 2*lep2vec.Py()*lep2vec.Pz()*b1*b3*sqr(a3)));

  l51 =(1/sqr(lep2vec.E())*(d1*sqr(a3)*sqr(b3) 
			    + d3*sqr(b1)*sqr(a3)
			    - 2*lep2vec.Px()*lep2vec.Pz()*b1*b3*sqr(a3)));
  
  l61 =(1/sqr(lep2vec.E())*(d2*sqr(a3)*sqr(b3) 
			    + d3*sqr(b2)*sqr(a3)
			    - 2*lep2vec.Py()*lep2vec.Pz()*b2*b3*sqr(a3)));
  
  l1 = l11*k61;
  l2 = l21*k61/k51;
  l3 = l31;
  l4 = l41/k51;
  l5 = l51*k61/(sqr(k51));
  l6 = l61/k61;
  
  k15 = k1*l5-l1*k5;
  k25 = k2*l5-l2*k5;
  k35 = k3*l5-l3*k5;
  k45 = k4*l5-l4*k5;
  
  k16 = k1*l6-l1*k6;
  k26 = k2*l6-l2*k6;
  k36 = k3*l6-l3*k6;
  k46 = k4*l6-l4*k6;
  k56 = k5*l6-l5*k6;

  coeff[0] = k15*sqr(k36)-k35*k36*k16-k56*sqr(k16);
  coeff[1] = 2*k15*k36*k46+k25*sqr(k36)+k35*(-k46*k16-k36*k26)-k45*k36*k16-2*k56*k26*k16;
  coeff[2] = k15*sqr(k46)+2*k25*k36*k46+k35*(-k46*k26-k36*k56)-k56*(sqr(k26)+2*k56*k16)-k45*(k46*k16+k36*k26);
  coeff[3] = k25*sqr(k46)-k35*k46*k56-k45*(k46*k26+k36*k56)-2*sqr(k56)*k26;
  coeff[4] = -k45*k46*k56-pow(k56,3);
  
  // normalization of coefficients
  int moc=(int(log10(fabs(coeff[0])))+int(log10(fabs(coeff[4]))))/2;
  coeff[0]=coeff[0]/TMath::Power(10,moc);
  coeff[1]=coeff[1]/TMath::Power(10,moc);
  coeff[2]=coeff[2]/TMath::Power(10,moc);
  coeff[3]=coeff[3]/TMath::Power(10,moc);
  coeff[4]=coeff[4]/TMath::Power(10,moc);
  
}

void Solvettbar::TopRec(const TLorentzVector lep1vec,
			const TLorentzVector lep2vec,
			const TLorentzVector B1vec,
			const TLorentzVector B2vec,
			const double sol)
{
  TVector3 tt_boost;
  TLorentzVector vtemp;
  
  double px1, py1, pz1;
  double px2, py2, pz2;
  
  px1 = sol*(a3*b3/k51);
  py1 = -(a3*b3/k61)*(k56*pow(sol,2) + k26*sol + k16)/(k36 + k46*sol);
  pz1 = -1/b3*(b1*px1 + b2*py1 - F);

  px2 = C - px1;
  py2 = D - py1;
  pz2 = 1/a3*(a1*px1 + a2*py1 + G);

  neu1vec.SetXYZM(px1, py1, pz1, 0.0);
  neu2vec.SetXYZM(px2, py2, pz2, 0.0);

  top1vec = B2vec + lep2vec + neu1vec;
  top2vec = B1vec + lep1vec + neu2vec;
  
  vtemp = (top1vec + top2vec);
  tt_boost = -vtemp.BoostVector();
  tt_top1vec = top1vec;
  tt_top2vec = top2vec;
  tt_top1vec.Boost(tt_boost);
  tt_top2vec.Boost(tt_boost);
}

		       
int Solvettbar::quartic(double *coeff, double *roots, double *disc) const 
{
  double w, n0, n1, n2;
  double cu[4];
  double p0, p1, h, t, z;
  double *px;

  // std::cout << "******* start by solving quartic problem *******" << std::endl;
  
  // cubic solution
  if(coeff[4] == 0.0){
    // std::cout << "look for cubic soln" << std::endl;
    return cubic(coeff, roots, disc);
  }
  // quartic solution
  // normalize all coefficients to coeff[4] 
  w = coeff[3]/(4*coeff[4]);
  // offset
  n2 = (-6*sqr(w) + coeff[2]/coeff[4]);
  // coeffs of shifted polynomial
  n1 = ((8*sqr(w) - 2*coeff[2]/coeff[4])*w + coeff[1]/coeff[4]);
  n0 = (((-3*sqr(w) + coeff[2]/coeff[4])*w - coeff[1]/coeff[4])*w + coeff[0]/coeff[4]);

  // solve with cubic  
  cu[3] = 1.0;
  cu[2] = n2; 
  cu[1] = -4*n0;
  cu[0] = sqr(n1) - 4*n0*n2;
  
  cubic(cu, roots, disc);

  z = roots[0];
 
  int nreal = 0;
  px = roots;
 
  t = sqrt(0.25*sqr(z) - n0);
  // std::cout << "t : " << t << std::endl;
  for(int i=-1; i<=1; i+=2){
    p0 = -0.5*z + i*t;

    // coeffs of quadratic factor
    p1 = (t!=0.0) ? -i*0.5*n1/t : i*sqrt(-z - n2);

    h  = 0.25*sqr(p1) - p0;
    /*
    std::cout << "p0 = -0.5*"<< z <<"+"<< i <<"*" << t << std::endl;
    std::cout << "p0 : " << p0 << std::endl;
    std::cout << "p1 = " << -i <<"*"<<"0.5*"<<n1<<"/"<<t<<std::endl;
    std::cout << "p1 : " << p1 << std::endl;
    std::cout << "h = 0.25*" <<z<< "+" << i <<"*"<<t<<std::endl; 
    std::cout << "h : " << h << std::endl;
    */
    if(h>=0.0){
      h = sqrt(h);
      nreal += 2;
      *px++ = -0.5*p1 - h - w;
      *px++ = -0.5*p1 + h - w;
    }
  }

  return nreal; 
}

int Solvettbar::cubic(const double *coeffs, double *roots, double *disc) const
{
  unsigned int nreal;
  double w, p, q, h, phi;
  double dis;
  // std::cout <<"cubic coefficients: " << coeffs[0] << ", " << coeffs[1] << ", " << coeffs[2] << ", " << coeffs[3] << std::endl;

  if (coeffs[3]!=0.0) {

    /* cubic problem? */
    w = coeffs[2]/(3*coeffs[3]);
    p = sqr(coeffs[1]/(3*coeffs[3])-sqr(w))*(coeffs[1]/(3*coeffs[3])-sqr(w));
    q = -0.5*(2*sqr(w)*w-(coeffs[1]*w-coeffs[0])/coeffs[3]);
    
    dis = sqr(q)+p;
    
    disc[0] = dis;
    /* discriminant */
    // std::cout << "cubic discriminant: " << dis << std::endl;
    if (dis<0.0) {
      /* 3 real solutions */
      h = q/sqrt(-p);

      if (h>1.0) h = 1.0;
      /* confine the argument of */
      if (h<-1.0) h = -1.0;
      /* acos to [-1;+1] */
      phi = acos(h);
      p = 2*TMath::Power(-p, 1.0/6.0);

      for(unsigned int i=0; i<3; i++) {
	roots[i] = p*cos((phi+2*i*TMath::Pi())/3.0) - w;
      }
     
      /* sort results */
      
      if (roots[1]<roots[0]) SWAP(roots[0], roots[1]);
      if (roots[2]<roots[1]) SWAP(roots[1], roots[2]);
      if (roots[1]<roots[0]) SWAP(roots[0], roots[1]);
      
      // std::cout << "root 0: " << roots[0] << ", root 1: " << roots[1] << ", root 2: " << roots[2] << std::endl;
      nreal = 3;
    }
    else {
      /* only one real solution */
      dis = sqrt(dis);
      h = TMath::Power(fabs(q+dis), 1.0/3.0);
      p = TMath::Power(fabs(q-dis), 1.0/3.0);
      roots[0] = ((q+dis>0.0)? h : -h) + ((q-dis>0.0)? p : -p) -  w;
      // std::cout << "root 0: " << roots[0] << std::endl;
      nreal = 1;
    }
    
    // Perform one step of a Newton iteration in order to minimize round-off errors
    for(unsigned int i=0; i<nreal; i++) {
      h = coeffs[1] + roots[i] * (2 * coeffs[2] + 3 * roots[i] * coeffs[3]);
      if (h != 0.0){
	roots[i] -= (coeffs[0] + roots[i] * (coeffs[1] + roots[i] * (coeffs[2] + roots[i] * coeffs[3])))/h;
      }
    }
    
  }

  else if (coeffs[2]!=0.0) {
    /* quadratic problem? */
    p = 0.5*coeffs[1]/coeffs[2];
    dis = sqr(p) - coeffs[0]/coeffs[2];
    if (dis>=0.0) {
      /* two real solutions */
      dis = sqrt(dis);
      roots[0] = -p - dis;
      roots[1] = -p + dis;
      nreal = 2;
    }
    else
      /* no real solution */
      nreal = 0;
  }

  else if (coeffs[1]!=0.0) {
    /* linear problem? */
    roots[0] = -coeffs[0]/coeffs[1];
    nreal = 1;
  }

  else
    /* no equation */
    nreal = 0;
  
  // std::cout << "cubic nreal : " << nreal << std::endl;
  return nreal;
}

void Solvettbar::SWAP(double& realone, double&realtwo) const
{
  if(realtwo < realone){
    double temp = realtwo;
    realtwo = realone;
    realone = temp;
  }
}

double Solvettbar::WeightSolfromMC() const
{
  double weight = 1;

  weight = 1-fabs(1-top1vec.M()/genTop1Mass)*fabs(1-top2vec.M()/genTop2Mass);
  // weight = (((top1vec.M() > genTop1.M()) ? genTop1.M()/top1vec.M() : top1vec.M()/genTop1.M())
  // 	    *((top2vec.M() > genTop2.M())? genTop2.M()/top2vec.M() : top2vec.M()/genTop2.M()));
  return weight;
}

double Solvettbar::WeightSolfromShape() const
{
  return EventShape_ -> Eval(neu1vec.E(),neu2vec.E());
}
