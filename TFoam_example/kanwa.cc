#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TFoam.h"
#include "TRandom2.h"
Double_t sqr(Double_t x){return x*x;};

Double_t Camel2(Int_t nDim, Double_t *Xarg){
     // 2-dimensional distribution for FOAM, normalized to one (within 1e-5)
     Double_t x=Xarg[0];
     Double_t y=Xarg[1];
     Double_t GamSq= sqr(0.100e0);
     Double_t Dist=exp(-(sqr(x-1./3) +sqr(y-1./3))/GamSq)/GamSq/TMath::Pi();
     Dist +=exp(-(sqr(x-2./3) +sqr(y-2./3))/GamSq)/GamSq/TMath::Pi();
     return Dist;
}// Camel2
 
void kanwa(){
      gSystem->Load("libFoam");
      TH2D *hst_xy = new TH2D("hst_xy" , "x-y plot", 50,0,1.0, 50,0,1.0);
      Double_t *MCvect =new Double_t[2]; // 2-dim vector generated in the MC run
      TRandom3 *PseRan = new TRandom3(); // Create random number generator
      PseRan->SetSeed(4357); // Set seed
      TFoam *FoamX = new TFoam("FoamX"); // Create Simulator
      FoamX->SetkDim(2); // No. of dimensions, obligatory!
      FoamX->SetnCells(500); // No. of cells, can be omitted, default=2000
      FoamX->SetRhoInt(Camel2); // Set 2-dim distribution, included below
      FoamX->SetPseRan(PseRan); // Set random number generator
      FoamX->Initialize(); // Initialize simulator, takes a few seconds...
      // From now on FoamX is ready to generate events according to Camel2(x,y)
      for(Long_t loop=0; loop<100000; loop++){
          FoamX->MakeEvent(); // generate MC event
          FoamX->GetMCvect( MCvect); // get generated vector (x,y)
          Double_t x=MCvect[0];
          Double_t y=MCvect[1];
          if(loop<10) std::cout<<"(x,y) = ( "<< x <<", "<< y <<" )"<<std::endl;
          hst_xy->Fill(x,y); // fill scattergram
      }// loop
      Double_t mcResult, mcError;
      FoamX->GetIntegMC( mcResult, mcError); // get MC integral, should be one
      std::cout << " mcResult= " << mcResult << " +- " << mcError <<std::endl;
      // now hst_xy will be plotted visualizing generated distribution
      TCanvas *cKanwa = new TCanvas("cKanwa","Canvas for plotting",600,600);
      cKanwa->cd();
      hst_xy->Draw("lego2");
}//kanwa

