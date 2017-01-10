#ifndef JET_RESOLUTIONS
#define JET_RESOLUTIONS

#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMatrixD.h"

class JetResolutions { 
public:
  
  JetResolutions() {
    
    // KH
    myRand = new TRandom();
    myRand->SetSeed(0xBEEFBABE);

    // 2016 JER taken from https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Spring16_25nsV6_MC/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt    
    
    //Smear in Pt: fPtPar0[eta][rho] 
    
    fPtSmear = new TF1("JetPt","sqrt([0]*fabs([0])/(x*x)+[1]*[1]*pow(x,[3])+[2]*[2])");
    
    // std::vector<double> fPtPar0[13];// = new double[26][5];
    // std::vector<double> fPtPar1[13];// = new double[26][5];
    // std::vector<double> fPtPar2[13];// = new double[26][5];
    // std::vector<double> fPtPar3[13];// = new double[26][5];

    //13 Eta Bins, 5 Rho bins
    //-4.7 -3.2
    fPtPar0[0].push_back(2.557); 
    fPtPar0[0].push_back(3.111); 
    fPtPar0[0].push_back(3.506); 
    fPtPar0[0].push_back(3.976); 
    fPtPar0[0].push_back(-21.36);
    fPtPar1[0].push_back(0.1523);
    fPtPar1[0].push_back(0.1684);
    fPtPar1[0].push_back(0.1863);
    fPtPar1[0].push_back(0.2705);
    fPtPar1[0].push_back(21.45);
    fPtPar2[0].push_back(0.0588);
    fPtPar2[0].push_back(0.000138);
    fPtPar2[0].push_back(2.267e-05);
    fPtPar2[0].push_back(0.09184);
    fPtPar2[0].push_back(0.104);
    fPtPar3[0].push_back(-0.2175);
    fPtPar3[0].push_back(-0.1774);
    fPtPar3[0].push_back(-0.2249);
    fPtPar3[0].push_back(-0.6338);
    fPtPar3[0].push_back(-1.986);

    //-3.2 -3
    fPtPar0[1].push_back(1.402);
    fPtPar0[1].push_back(-0.02076); 
    fPtPar0[1].push_back(6.215);
    fPtPar0[1].push_back(6.949);
    fPtPar0[1].push_back(7.753);
    fPtPar1[1].push_back(3.892);
    fPtPar1[1].push_back(4.758);
    fPtPar1[1].push_back(0.154);
    fPtPar1[1].push_back(-1.793e-05);
    fPtPar1[1].push_back(-5.195e-07);
    fPtPar2[1].push_back(0.1354);
    fPtPar2[1].push_back(0.136);
    fPtPar2[1].push_back(0.0007419);
    fPtPar2[1].push_back(0.1397);
    fPtPar2[1].push_back(0.1407);
    fPtPar3[1].push_back(-1.908);
    fPtPar3[1].push_back(-1.898);
    fPtPar3[1].push_back(-0.04238);
    fPtPar3[1].push_back(-1.466);
    fPtPar3[1].push_back(-1.47);

    //-3 -2.8
    fPtPar0[2].push_back(4.573);
    fPtPar0[2].push_back(5.63);  
    fPtPar0[2].push_back(6.06);  
    fPtPar0[2].push_back(7.013);
    fPtPar0[2].push_back(-102.9);
    fPtPar1[2].push_back(0.1926);
    fPtPar1[2].push_back(0.2078);
    fPtPar1[2].push_back(0.3396);
    fPtPar1[2].push_back(0.2743);
    fPtPar1[2].push_back(102.4);
    fPtPar2[2].push_back(3.189e-06);
    fPtPar2[2].push_back(4.508e-06);
    fPtPar2[2].push_back(0.05938);
    fPtPar2[2].push_back(2.302e-06);
    fPtPar2[2].push_back(0.08302);
    fPtPar3[2].push_back(-0.2268);
    fPtPar3[2].push_back(-0.2741);
    fPtPar3[2].push_back(-0.5531);
    fPtPar3[2].push_back(-0.3832);
    fPtPar3[2].push_back(-1.995);

    //-2.8 -2.5
    fPtPar0[3].push_back(4.585);
    fPtPar0[3].push_back(5.318);
    fPtPar0[3].push_back(6.201);
    fPtPar0[3].push_back(7.217);
    fPtPar0[3].push_back(7.962);
    fPtPar1[3].push_back(0.2118);
    fPtPar1[3].push_back(0.2424);
    fPtPar1[3].push_back(0.2732);
    fPtPar1[3].push_back(0.2285);
    fPtPar1[3].push_back(0.3192);
    fPtPar2[3].push_back(6.684e-07);
    fPtPar2[3].push_back(2.644e-06);
    fPtPar2[3].push_back(2.313e-07);
    fPtPar2[3].push_back(2.705e-06);
    fPtPar2[3].push_back(5.103e-07);
    fPtPar3[3].push_back(-0.3778);
    fPtPar3[3].push_back(-0.4489);
    fPtPar3[3].push_back(-0.4922);
    fPtPar3[3].push_back(-0.4242);
    fPtPar3[3].push_back(-0.546);

    //-2.5 -2.3
    fPtPar0[4].push_back(4.205);
    fPtPar0[4].push_back(4.862);
    fPtPar0[4].push_back(3.554);
    fPtPar0[4].push_back(4.061);
    fPtPar0[4].push_back(7.087);
    fPtPar1[4].push_back(0.3201);
    fPtPar1[4].push_back(0.3633);
    fPtPar1[4].push_back(1.426);
    fPtPar1[4].push_back(1.696);
    fPtPar1[4].push_back(0.4703);
    fPtPar2[4].push_back(3.283e-06);
    fPtPar2[4].push_back(2.844e-07);
    fPtPar2[4].push_back(0.04178);
    fPtPar2[4].push_back(0.04539);
    fPtPar2[4].push_back(0.02352);
    fPtPar3[4].push_back(-0.586);
    fPtPar3[4].push_back(-0.6215);
    fPtPar3[4].push_back(-1.191);
    fPtPar3[4].push_back(-1.251);
    fPtPar3[4].push_back(-0.7083);
    //-2.3 -2.1
    fPtPar0[5].push_back(3.484);
    fPtPar0[5].push_back(3.938);
    fPtPar0[5].push_back(4.801);
    fPtPar0[5].push_back(5.449);
    fPtPar0[5].push_back(6.281);
    fPtPar1[5].push_back(0.4017);
    fPtPar1[5].push_back(0.5564);
    fPtPar1[5].push_back(0.4242);
    fPtPar1[5].push_back(0.4428);
    fPtPar1[5].push_back(0.4152);
    fPtPar2[5].push_back(0.02593);
    fPtPar2[5].push_back(0.03333);
    fPtPar2[5].push_back(0.01964);
    fPtPar2[5].push_back(0.01493);
    fPtPar2[5].push_back(-2.964e-06);
    fPtPar3[5].push_back(-0.6917);
    fPtPar3[5].push_back(-0.823);
    fPtPar3[5].push_back(-0.6759);
    fPtPar3[5].push_back(-0.6716);
    fPtPar3[5].push_back(-0.6261);
    //-2.1 -1.9 
    fPtPar0[6].push_back(3.405);
    fPtPar0[6].push_back(3.857);
    fPtPar0[6].push_back(4.363);
    fPtPar0[6].push_back(4.46);
    fPtPar0[6].push_back(5.394);
    fPtPar1[6].push_back(0.4581);
    fPtPar1[6].push_back(0.491);
    fPtPar1[6].push_back(0.5315);
    fPtPar1[6].push_back(0.6932);
    fPtPar1[6].push_back(0.5988);
    fPtPar2[6].push_back(0.03125);
    fPtPar2[6].push_back(0.03292);
    fPtPar2[6].push_back(0.03569);
    fPtPar2[6].push_back(0.03589);
    fPtPar2[6].push_back(0.03651);
    fPtPar3[6].push_back(-0.746);
    fPtPar3[6].push_back(-0.7656);
    fPtPar3[6].push_back(-0.7876);
    fPtPar3[6].push_back(-0.872);
    fPtPar3[6].push_back(-0.8103);
    //-1.9 -1.7
    fPtPar0[7].push_back(2.663);
    fPtPar0[7].push_back(3.108);
    fPtPar0[7].push_back(3.444);
    fPtPar0[7].push_back(2.351);
    fPtPar0[7].push_back(4.126);
    fPtPar1[7].push_back(0.7106);
    fPtPar1[7].push_back(0.8232);
    fPtPar1[7].push_back(0.8742);
    fPtPar1[7].push_back(1.452);
    fPtPar1[7].push_back(1.004);
    fPtPar2[7].push_back(0.03652);
    fPtPar2[7].push_back(0.04012);
    fPtPar2[7].push_back(0.0407);
    fPtPar2[7].push_back(0.04581);
    fPtPar2[7].push_back(0.04237);
    fPtPar3[7].push_back(-0.9051);
    fPtPar3[7].push_back(-0.9596);
    fPtPar3[7].push_back(-0.9711);
    fPtPar3[7].push_back(-1.161);
    fPtPar3[7].push_back(-0.9977);
    //-1.7 -1.3
    fPtPar0[8].push_back(-1.259);
    fPtPar0[8].push_back(-0.4585);
    fPtPar0[8].push_back(2.314);
    fPtPar0[8].push_back(3.078);
    fPtPar0[8].push_back(3.526);
    fPtPar1[8].push_back(1.39);
    fPtPar1[8].push_back(1.512);
    fPtPar1[8].push_back(1.441);
    fPtPar1[8].push_back(1.414);
    fPtPar1[8].push_back(1.522);
    fPtPar2[8].push_back(0.05809);
    fPtPar2[8].push_back(0.05827);
    fPtPar2[8].push_back(0.05851);
    fPtPar2[8].push_back(0.05949);
    fPtPar2[8].push_back(0.06019);
    fPtPar3[8].push_back(-1.155);
    fPtPar3[8].push_back(-1.177);
    fPtPar3[8].push_back(-1.153);
    fPtPar3[8].push_back(-1.144);
    fPtPar3[8].push_back(-1.163);
    //-1.3 -1.1
    fPtPar0[9].push_back(1.783);
    fPtPar0[9].push_back(2.566);
    fPtPar0[9].push_back(2.988);
    fPtPar0[9].push_back(3.449);
    fPtPar0[9].push_back(3.802);
    fPtPar1[9].push_back(0.7418);
    fPtPar1[9].push_back(0.7314);
    fPtPar1[9].push_back(0.7753);
    fPtPar1[9].push_back(0.7868);
    fPtPar1[9].push_back(0.8917);
    fPtPar2[9].push_back(0.05023);
    fPtPar2[9].push_back(0.04992);
    fPtPar2[9].push_back(0.04999);
    fPtPar2[9].push_back(0.04921);
    fPtPar2[9].push_back(0.05064);
    fPtPar3[9].push_back(-0.8928);
    fPtPar3[9].push_back(-0.8785);
    fPtPar3[9].push_back(-0.8917);
    fPtPar3[9].push_back(-0.8871);
    fPtPar3[9].push_back(-0.9323);
    //-1.1 -0.8
    fPtPar0[10].push_back(1.85);
    fPtPar0[10].push_back(2.383);
    fPtPar0[10].push_back(2.937);
    fPtPar0[10].push_back(3.259);
    fPtPar0[10].push_back(3.965);
    fPtPar1[10].push_back(0.4623);
    fPtPar1[10].push_back(0.4783);
    fPtPar1[10].push_back(0.5203);
    fPtPar1[10].push_back(0.5784);
    fPtPar1[10].push_back(0.5726);
    fPtPar2[10].push_back(0.0341);
    fPtPar2[10].push_back(0.03276);
    fPtPar2[10].push_back(0.03549);
    fPtPar2[10].push_back(0.03626);
    fPtPar2[10].push_back(0.03623);
    fPtPar3[10].push_back(-0.705);
    fPtPar3[10].push_back(-0.7081);
    fPtPar3[10].push_back(-0.7447);
    fPtPar3[10].push_back(-0.7783);
    fPtPar3[10].push_back(-0.7726);
    //-0.8 -0.5
    fPtPar0[11].push_back(1.598);
    fPtPar0[11].push_back(2.209);
    fPtPar0[11].push_back(2.857);
    fPtPar0[11].push_back(3.384);
    fPtPar0[11].push_back(3.661);
    fPtPar1[11].push_back(0.368);
    fPtPar1[11].push_back(0.4084);
    fPtPar1[11].push_back(0.41);
    fPtPar1[11].push_back(0.431);
    fPtPar1[11].push_back(0.5082);
    fPtPar2[11].push_back(0.01819);
    fPtPar2[11].push_back(0.02308);
    fPtPar2[11].push_back(0.02228);
    fPtPar2[11].push_back(0.02329);
    fPtPar2[11].push_back(0.02756);
    fPtPar3[11].push_back(-0.614);
    fPtPar3[11].push_back(-0.6567);
    fPtPar3[11].push_back(-0.6559);
    fPtPar3[11].push_back(-0.671);
    fPtPar3[11].push_back(-0.7339);
    //-0.5 0
    fPtPar0[12].push_back(1.308);
    fPtPar0[12].push_back(1.962);
    fPtPar0[12].push_back(2.584);
    fPtPar0[12].push_back(3.215);
    fPtPar0[12].push_back(3.728);
    fPtPar1[12].push_back(0.3597);
    fPtPar1[12].push_back(0.4033);
    fPtPar1[12].push_back(0.4265);
    fPtPar1[12].push_back(0.4443);
    fPtPar1[12].push_back(0.4898);
    fPtPar2[12].push_back(0.01513);
    fPtPar2[12].push_back(0.01816);
    fPtPar2[12].push_back(0.02032);
    fPtPar2[12].push_back(0.02137);
    fPtPar2[12].push_back(0.02347);
    fPtPar3[12].push_back(-0.6041);
    fPtPar3[12].push_back(-0.6447);
    fPtPar3[12].push_back(-0.6648);
    fPtPar3[12].push_back(-0.6788);
    fPtPar3[12].push_back(-0.7146);


    //Smear in Phi
    fPhiSmear = new TF1("JetPhi","[0]+[1]*exp(-x/[2])");
    // std::vector<double> fPhiPar0[13];// = new double[13][5];
    // std::vector<double> fPhiPar1[13];// = new double[13][5];
    // std::vector<double> fPhiPar2[13];// = new double[13][5];
    
    //-4.7 -3.2
    fPhiPar0[0].push_back(0.02332);
    fPhiPar0[0].push_back(0.0254);
    fPhiPar0[0].push_back(0.02592);
    fPhiPar0[0].push_back(0.02616);
    fPhiPar0[0].push_back(0.02611);
    fPhiPar1[0].push_back(0.08121);
    fPhiPar1[0].push_back(0.1465);
    fPhiPar1[0].push_back(0.1447);
    fPhiPar1[0].push_back(0.1744);
    fPhiPar1[0].push_back(0.1703);
    fPhiPar2[0].push_back(10.29);
    fPhiPar2[0].push_back(8.57);
    fPhiPar2[0].push_back(9.493);
    fPhiPar2[0].push_back(9.628);
    fPhiPar2[0].push_back(10.7);


    fPhiPar0[1].push_back(0.0256);
    fPhiPar0[1].push_back(0.02598);
    fPhiPar0[1].push_back(0.02745);
    fPhiPar0[1].push_back(0.02668);
    fPhiPar0[1].push_back(0.02814);
    fPhiPar1[1].push_back(0.1153);
    fPhiPar1[1].push_back(0.1221);
    fPhiPar1[1].push_back(0.1578);
    fPhiPar1[1].push_back(0.1339);
    fPhiPar1[1].push_back(0.1798);
    fPhiPar2[1].push_back(14.98);
    fPhiPar2[1].push_back(16.01);
    fPhiPar2[1].push_back(14.33);
    fPhiPar2[1].push_back(17.67);
    fPhiPar2[1].push_back(15.07);

    fPhiPar0[2].push_back(0.01952);
    fPhiPar0[2].push_back(0.02003);
    fPhiPar0[2].push_back(0.01947);
    fPhiPar0[2].push_back(0.02011);
    fPhiPar0[2].push_back(0.022);
    fPhiPar1[2].push_back(0.09488);
    fPhiPar1[2].push_back(0.1009);
    fPhiPar1[2].push_back(0.08778);
    fPhiPar1[2].push_back(0.09804);
    fPhiPar1[2].push_back(0.1132);
    fPhiPar2[2].push_back(22.05);
    fPhiPar2[2].push_back(21.61);
    fPhiPar2[2].push_back(27.12);
    fPhiPar2[2].push_back(27.46);
    fPhiPar2[2].push_back(26.41);

    fPhiPar0[3].push_back(0.0135);
    fPhiPar0[3].push_back(0.01329);
    fPhiPar0[3].push_back(0.01368);
    fPhiPar0[3].push_back(0.01396);
    fPhiPar0[3].push_back(0.01457);
    fPhiPar1[3].push_back(0.06913);
    fPhiPar1[3].push_back(0.06797);
    fPhiPar1[3].push_back(0.06985);
    fPhiPar1[3].push_back(0.07725);
    fPhiPar1[3].push_back(0.08012);
    fPhiPar2[3].push_back(29.44);
    fPhiPar2[3].push_back(31.01);
    fPhiPar2[3].push_back(32.46);
    fPhiPar2[3].push_back(33.6);
    fPhiPar2[3].push_back(35.82);

    fPhiPar0[4].push_back(0.009148);
    fPhiPar0[4].push_back(0.009558);
    fPhiPar0[4].push_back(0.01024);
    fPhiPar0[4].push_back(0.01146);
    fPhiPar0[4].push_back(0.01085);
    fPhiPar1[4].push_back(0.03953);
    fPhiPar1[4].push_back(0.04439);
    fPhiPar1[4].push_back(0.05175);
    fPhiPar1[4].push_back(0.06791);
    fPhiPar1[4].push_back(0.05936);
    fPhiPar2[4].push_back(44.86);
    fPhiPar2[4].push_back(40.92);
    fPhiPar2[4].push_back(36.51);
    fPhiPar2[4].push_back(28.96);
    fPhiPar2[4].push_back(40.49);

    fPhiPar0[5].push_back(0.007723);
    fPhiPar0[5].push_back(0.007738);
    fPhiPar0[5].push_back(0.008745);
    fPhiPar0[5].push_back(0.00929);
    fPhiPar0[5].push_back(0.009744);
    fPhiPar1[5].push_back(0.03159);
    fPhiPar1[5].push_back(0.03778);
    fPhiPar1[5].push_back(0.0451);
    fPhiPar1[5].push_back(0.05474);
    fPhiPar1[5].push_back(0.06141);
    fPhiPar2[5].push_back(58.13);
    fPhiPar2[5].push_back(53.51);
    fPhiPar2[5].push_back(43.43);
    fPhiPar2[5].push_back(36.36);
    fPhiPar2[5].push_back(35.04);

    fPhiPar0[6].push_back(0.007595);
    fPhiPar0[6].push_back(0.008105);
    fPhiPar0[6].push_back(0.008014);
    fPhiPar0[6].push_back(0.008295);
    fPhiPar0[6].push_back(0.008872);
    fPhiPar1[6].push_back(0.03543);
    fPhiPar1[6].push_back(0.03738);
    fPhiPar1[6].push_back(0.04165);
    fPhiPar1[6].push_back(0.04993);
    fPhiPar1[6].push_back(0.05531);
    fPhiPar2[6].push_back(52.48);
    fPhiPar2[6].push_back(50.34);
    fPhiPar2[6].push_back(48.25);
    fPhiPar2[6].push_back(42.22);
    fPhiPar2[6].push_back(40.96);
    
    fPhiPar0[7].push_back(0.007854);
    fPhiPar0[7].push_back(0.008684);
    fPhiPar0[7].push_back(0.008432);
    fPhiPar0[7].push_back(0.008178);
    fPhiPar0[7].push_back(0.008973);
    fPhiPar1[7].push_back(0.03468);
    fPhiPar1[7].push_back(0.04505);
    fPhiPar1[7].push_back(0.04922);
    fPhiPar1[7].push_back(0.0483);
    fPhiPar1[7].push_back(0.05905);
    fPhiPar2[7].push_back(58.09);
    fPhiPar2[7].push_back(43.6);
    fPhiPar2[7].push_back(43.67);
    fPhiPar2[7].push_back(49.66);
    fPhiPar2[7].push_back(39.83);

    fPhiPar0[8].push_back(0.008264);
    fPhiPar0[8].push_back(0.008579);
    fPhiPar0[8].push_back(0.009563);
    fPhiPar0[8].push_back(0.009542);
    fPhiPar0[8].push_back(0.01032);
    fPhiPar1[8].push_back(0.03564);
    fPhiPar1[8].push_back(0.04355);
    fPhiPar1[8].push_back(0.05593);
    fPhiPar1[8].push_back(0.06466);
    fPhiPar1[8].push_back(0.07494);
    fPhiPar2[8].push_back(61.14);
    fPhiPar2[8].push_back(51.01);
    fPhiPar2[8].push_back(39.32);
    fPhiPar2[8].push_back(36.94);
    fPhiPar2[8].push_back(33.85);

    fPhiPar0[9].push_back(0.007255);
    fPhiPar0[9].push_back(0.007838);
    fPhiPar0[9].push_back(0.008754);
    fPhiPar0[9].push_back(0.009131);
    fPhiPar0[9].push_back(0.009631);
    fPhiPar1[9].push_back(0.03034);
    fPhiPar1[9].push_back(0.03813);
    fPhiPar1[9].push_back(0.04878);
    fPhiPar1[9].push_back(0.05651);
    fPhiPar1[9].push_back(0.07024);
    fPhiPar2[9].push_back(67.08);
    fPhiPar2[9].push_back(52.11);
    fPhiPar2[9].push_back(40.49);
    fPhiPar2[9].push_back(36.84);
    fPhiPar2[9].push_back(32.77);

    fPhiPar0[10].push_back(0.006828);
    fPhiPar0[10].push_back(0.007306);
    fPhiPar0[10].push_back(0.007676);
    fPhiPar0[10].push_back(0.008184);
    fPhiPar0[10].push_back(0.00859);
    fPhiPar1[10].push_back(0.02736);
    fPhiPar1[10].push_back(0.03378);
    fPhiPar1[10].push_back(0.04057);
    fPhiPar1[10].push_back(0.05021);
    fPhiPar1[10].push_back(0.06082);
    fPhiPar2[10].push_back(70.3);
    fPhiPar2[10].push_back(56.4);
    fPhiPar2[10].push_back(48.08);
    fPhiPar2[10].push_back(40.48);
    fPhiPar2[10].push_back(35.84);

    fPhiPar0[11].push_back(0.006604);
    fPhiPar0[11].push_back(0.006738);
    fPhiPar0[11].push_back(0.007114);
    fPhiPar0[11].push_back(0.007778);
    fPhiPar0[11].push_back(0.007987);
    fPhiPar1[11].push_back(0.02718);
    fPhiPar1[11].push_back(0.03036);
    fPhiPar1[11].push_back(0.03698);
    fPhiPar1[11].push_back(0.04757);
    fPhiPar1[11].push_back(0.05785);
    fPhiPar2[11].push_back(67.03);
    fPhiPar2[11].push_back(63.87);
    fPhiPar2[11].push_back(53.1);
    fPhiPar2[11].push_back(42.44);
    fPhiPar2[11].push_back(38.29);

    fPhiPar0[12].push_back(0.006662);
    fPhiPar0[12].push_back(0.00665);
    fPhiPar0[12].push_back(0.007047);
    fPhiPar0[12].push_back(0.007271);
    fPhiPar0[12].push_back(0.007744);
    fPhiPar1[12].push_back(0.02627);
    fPhiPar1[12].push_back(0.02818);
    fPhiPar1[12].push_back(0.03473);
    fPhiPar1[12].push_back(0.04078);
    fPhiPar1[12].push_back(0.05205);
    fPhiPar2[12].push_back(66.67);
    fPhiPar2[12].push_back(68.9);
    fPhiPar2[12].push_back(55.96);
    fPhiPar2[12].push_back(50.09);
    fPhiPar2[12].push_back(41.09);
  };
  ~JetResolutions() {
    //KH
    delete myRand;

    delete fPtSmear;
    delete fPhiSmear;
    // delete[] fPtPar0; delete[] fPtPar1; delete[] fPtPar2; delete[] fPtPar3;
    // delete[] fPhiPar0; delete[] fPhiPar1; delete[] fPhiPar2;

  }
  double getPtSigmaSmear(TLorentzVector & vJet, const double JetRho);
  double getPhiSigmaSmear(TLorentzVector & vJet, const double JetRho);
  void setSmear(TLorentzVector & vJet, const double JetRho, TLorentzVector vSmear);
  void getUncertainties( TLorentzVector & v, const double rho, double & upx, double & upy, double & um );//double & upz, 
  void readRecoil(std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,        
		  std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit);
  void covMatrix(TMatrixD & iCov, double iGenPt, double iGenPhi, int njet);
  void getSmearedMET(const TLorentzVector * iMET, const TMatrixD* CovMat, TLorentzVector * oMET, double & sigmaX, double & sigmaY);
  // double ptres(double x);
  // double etares(double x);
  // double phires(double x);

  TRandom * myRand;

  TF1 * fPtSmear; 
  // TF1 * fEtaSmear; 
  TF1 * fPhiSmear; 

  std::vector<double> fPtPar0[13] ;
  std::vector<double> fPtPar1[13] ;
  std::vector<double> fPtPar2[13] ;
  std::vector<double> fPtPar3[13] ;
  /*
  std::vector<double*>fEtaPar0 ;
  std::vector<double*>fEtaPar1 ;
  std::vector<double*>fEtaPar2 ;
  std::vector<double*>fEtaPar3 ;
  std::vector<double*>fEtaPar4 ;
  */
  std::vector<double> fPhiPar0[13] ;
  std::vector<double> fPhiPar1[13] ;
  std::vector<double> fPhiPar2[13] ;

};


#endif
