#ifndef FIT_SOLUTION
#define FIT_SOLUTION

class FitSolution { 
public:

  FitSolution():
    numSol(-1.)
  {}    
  ~FitSolution(){}
  std::vector<TLorentzVector> nu1; std::vector<TLorentzVector> nu2;
  int numSol;
  std::vector<double> mass1;
  std::vector<double> mass2;
  std::vector<double> weight;
  std::vector<double> cudisc;
  std::vector<double> mtt;

  // KH
  void reset() { 
    mass1.clear();
    mass2.clear();
    weight.clear();
    cudisc.clear();
    mtt.clear();
    nu1.clear();
    nu2.clear();
    numSol=0;
  }

};

#endif
