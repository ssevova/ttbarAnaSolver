#ifndef SOLVETTBAR
#define SOLVETTBAR

#include "TLorentzVector.h"
#include "TMath.h"
#include "TF2.h"
#include "TF1.h"

#include "JetResolutions.hh"
#include "ttbarCandidate.hh"
#include "FitSolution.hh"


class Solvettbar {

public:
  //default constructor
  Solvettbar();
  Solvettbar(const float, const float, const float, const double=80.4, const double=4.8);
  //destructor 
  ~Solvettbar();

  void SetConstraints( const double xx=0, const double yy=0 );
  // void SmearJet( TLorentzVector & jetvec, double & jetrho, JetResolutions & myres, double uncert_pt, double uncert_phi );
  // void SmearMet( TLorentzVector & metvec, JetResolutions & myres );
  // FitSolution Solve( FitSolution & fitsol, ttbarCandidate & ttbarcand);
  FitSolution SolveNu( FitSolution & fitsol, ttbarCandidate & ttbarcand);
  void SetGenTops(const TLorentzVector & gentop1, const TLorentzVector & gentop2);
  void SetWMass(const double w1mass, const double w2mass);
public:
  double dis;

private:
  
  void FindCoeff(const TLorentzVector lep1vec,
		 const TLorentzVector lep2vec,
		 const TLorentzVector B1vec,
		 const TLorentzVector B2vec,
		 const double mtop1,
		 const double mtop2,                                                                                                              
		 double* coeff);
  
  void TopRec(const TLorentzVector lep1vec,
	      const TLorentzVector lep2vec,
	      const TLorentzVector B1vec,
	      const TLorentzVector B2vec,
	      const double sol);
  
  double WeightSolfromMC() const;

  double WeightSolfromShape() const;

  int quartic(double* q_coeff, double* q_sol, double* disc) const;

  int cubic(const double* coeffs, double* roots, double* disc) const;
  
  double sqr(const double x) const { return (x*x); }

  void SWAP(double& realone, double& realtwo) const;

private:

  const float topmass_begin;
  const float topmass_end;
  const float topmass_step;

  TLorentzVector genTop1, genTop2;
  double mw1, mw2;
  const double mb;



  double pxmiss, pymiss;
  double C;
  double D;
  double F;
  double G;
  double k51;
  double k61;
  double k16;
  double k26;
  double k36;
  double k46;
  double k56;
  double a1, a2, a3;
  double b1, b2, b3;
  
  TLorentzVector neu1vec, neu2vec, top1vec, top2vec, tt_top1vec, tt_top2vec;

  TF2* EventShape_;  

};

#endif
