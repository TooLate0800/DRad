#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"
#include "TFoamIntegrand.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <cstdlib>
#include <cstdio>
#include <time.h>

#include "Math/Functor.h"
#include "Math/Interpolator.h"
#include "Math/GSLIntegrator.h"

#include <iostream>



const double Pi = TMath::Pi();
const double deg = Pi / 180.; // Degree to radian conversion
const double m = 0.51099893e-3; // Mass of the electron/positron (in GeV)
const double m2 = TMath::Power(m, 2);
const double M = 1875.612928e-3; // Mass of the deuteron (in GeV)
const double M2 = TMath::Power(M, 2);
const double alpha = 1. / 137.036; // Fine-structure constant
const double e = TMath::Sqrt(4. * Pi *alpha);   // Electron charge magnitude
const double fm = 0.197327; // GeV^{-1} to fm conversion
const double mkb = 389.379404; // GeV^{-2} to mkbarn conversion

const double rd = 1.0 * 2.130 / fm; // GeV^{-1}
const double mepsilon = 0.936 * 0.0022246 / 0.197 / 0.197;

 
Double_t theta_max = 6.0;
Double_t theta_min = 0.6;

Double_t sqr(Double_t x){return x*x;};
Double_t Camel2(Int_t nDim, Double_t *Xarg){
     // 1-dimensional distribution for FOAM, normalized to one (within 1e-5)
     Double_t theta=(theta_min + Xarg[0] * (theta_max - theta_min))*deg;
     Double_t Ei_1 = 1.1;
     Double_t cost2 = cos(theta / 2.0);
     Double_t mott = TMath::Power((alpha * cost2 / (2. * Ei_1 * (1. - cost2 * cost2))),2);
     return mott*sin(theta); 
}// Camel2

void Mott(){
      gSystem->Load("libFoam");
      TH1D *hst_x = new TH1D("hst_x" , "theta plot", 50 , theta_min, theta_max);//TH1D(name of histogram, title, number of bin, lower limit, upper limit)
      Double_t *MCvect =new Double_t[1]; // 1-dim vector generated in the MC run
      TRandom3 *PseRan = new TRandom3(); // Create random number generator
      PseRan->SetSeed(4357); // Set seed
      TFoam *FoamX = new TFoam("FoamX"); // Create Simulator
      FoamX->SetkDim(1); // No. of dimensions, obligatory!
      FoamX->SetnCells(500); // No. of cells, can be omitted, default=2000
      FoamX->SetRhoInt(Camel2); // Set 1-dim distribution, included below
      FoamX->SetPseRan(PseRan); // Set random number generator
      FoamX->Initialize(); // Initialize simulator, takes a few seconds...
      // From now on FoamX is ready to generate events according to Camel2(x,y)
      for(Long_t loop=0; loop<100000; loop++){
          FoamX->MakeEvent(); // generate MC event
          FoamX->GetMCvect( MCvect); // get generated vector (x,y)
          Double_t theta=theta_min + MCvect[0] * (theta_max - theta_min);
          hst_x->Fill(theta); // fill 1-D histogram
      }// loop
      Double_t mcResult, mcError;
      FoamX->GetIntegMC( mcResult, mcError); // get MC integral, should be one
      std::cout << " mcResult= " << mcResult << " +- " << mcError <<std::endl;
      // now hst_x will be plotted visualizing generated distribution
      TCanvas *c1 = new TCanvas("cKanwa","Canvas for plotting",600,600);
      c1->cd();
      hst_x->Draw("");
}//Mott
