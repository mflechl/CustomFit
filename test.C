#include "CustomFit.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include <stdlib.h>

void test(int cat=0){
  gStyle->SetOptStat(0);

  //  gROOT->ProcessLine(".L CustomFit.C+");
  CustomFit cf;

  int ret=cf.setInputHisto( "FF_corr_Wjets_MCsum_noGen.root" , "c_t" );
  if ( ret != 0 ) return;
  TH1D *h=cf.returnInputHisto();

  //qcd
  const int nbins=8;
  double a_bins[nbins]   = {20.,22.5,25.,27.5,30.,35.,40.,50.};

  //w+jets
  //const int nbins=9;
  //double a_bins[nbins] = {20.,22.5,25.,27.5,30.,35.,40.,50.,60.};


  cf.set_fitFunc( "landau(0)+pol0(2)" );
  //  cf.set_fitFunc( "expo(0)+pol0(2)" );
  //  cf.set_fitFunc( "pol1(0)*expo(2)+pol0(4)" );
  cf.set_fitFromBin( 1+cat*nbins );
  cf.set_fitToBin( cat*(nbins+1) );
  cf.set_fitMin( 20 );
  cf.set_fitMax( 200 );
  cf.set_histo_bins( 200 );

  std::vector<double> bins;
  //  for (int i=0; i<nbins; i++) bins.push_back(a_bins[i] );
  for (int i=0; i<(nbins-1); i++) bins.push_back( (a_bins[i+1]+a_bins[i])/2 ); //for now
  bins.push_back( a_bins[nbins-1]+  (a_bins[nbins-1]-a_bins[nbins-2])/2 );
  //  for (int i=0; i<bins.size(); i++) std::cout << bins.at(i) << " " << flush; std::cout << std::endl;
  cf.set_bin_centers( bins );

  cf.fitHisto();

  TGraph *g_fit_input=cf.returnFitInputGraph(); //the input to the fit: data points in the range given
  TF1 *f_fit=cf.returnFitForm();                //the fit result (function)
  TH1D *h_fit=cf.returnFitHisto();              //the fit result binned (histo)

  TCanvas *c2=new TCanvas();
  //  h->Draw("E");     //to set the axis
  //  h->GetYaxis()->SetRangeUser(0.4,2.0);
  //  g_fit->SetMarkerStyle(21);
  f_fit->SetTitle("");
  f_fit->Draw();
  f_fit->GetXaxis()->SetTitle("p_{T} [GeV]");
  f_fit->GetYaxis()->SetTitle("Fake factor");
  f_fit->SetMinimum(0);
  //  h_fit->SetLineWidth(3);
  h_fit->SetLineColor(602);
  h_fit->Draw("E same");
  g_fit_input->SetLineWidth(2);
  g_fit_input->Draw("P same");

  c2->SaveAs("fit.png");

}

#ifndef __CINT__
int main(int argc, char* argv[]) {
  int channel = argc > 1 ? atoi(argv[1]) : 0;
  test(channel);    
}
#endif
