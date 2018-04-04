#include "CustomFit.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TString.h"
#include <stdlib.h>

void test(int cat=0, TString fname="FF_corr_QCD_MCsum_noGen.root"){
  gStyle->SetOptStat(0);

  //  gROOT->ProcessLine(".L CustomFit.C+");
  CustomFit cf;

  //  TString fname="FF_corr_Wjets_MCsum_noGen.root";
  //  TString fname="FF_corr_QCD_MCsum_noGen.root";
  TString hname="c_t";
  TString wname="bins_weighted";

  TString proc="unknown";
  if ( fname.Contains("_QCD") ) proc="qcd";
  if ( fname.Contains("_Wj") ) proc="wjets";
  if ( fname.Contains("_TT") ) proc="tt";

  int ret=cf.setInputHisto( fname , hname );
  if ( ret != 0 ){
    std::cout << "Cannot open input file " << fname << " with histogram " << hname << std::endl;
    return;
  }
  TH1D *h=cf.returnInputHisto();

  //et, ICHEP
  //  if ( ( cat==1 && proc=="tt" )     ||     ( (cat==1 || cat==3) && proc=="qcd" ) ){
  if ( (cat==1 || cat==3) && proc=="qcd" ){
    cf.set_fitFunc( "landau(0)+pol0(2)" );
    //    cf.set_errFunc( "landau(0)+pol1(2)" ); // opt 2
    //    cf.set_err_scale( 1.0 ); //opt 2
    cf.set_err_scale( 3.0 ); //opt 1
    cf.set_err_cl( 0 );
  } else{
    cf.set_fitFunc( "landau(0)+pol1(2)" );
    cf.set_err_scale( 1.2 );
    if ( (cat==1 || cat==3) && proc=="wjets" ) cf.set_histMaxFrac( 0.65 );
  }

  //mt, ICHEP
  /*
  if ( cat==1 && ( proc=="qcd" || proc=="tt" ) ){
    cf.set_fitFunc( "landau(0)+pol0(2)" );
    //    cf.set_errFunc( "landau(0)+pol1(2)" ); // opt 2
    //    cf.set_err_scale( 1.0 ); //opt 2
    cf.set_err_scale( 3.0 ); //opt 1
    cf.set_err_cl( 0 );
  } else{
    cf.set_fitFunc( "landau(0)+pol1(2)" );
    cf.set_err_scale( 1.2 );
    if ( (cat==1 || cat==3) && proc=="wjets" ) cf.set_histMaxFrac( 0.75 );
  }
  */

  //  cf.set_fitFunc( "landau(0)+pol0(2)" );
  //  cf.set_err_scale( 2.0 );

  int nbins=h->GetNbinsX()/4;

  //  cf.set_fitFunc( "expo(0)+pol1(2)" );
  //  cf.set_fitFunc( "pol1(0)*expo(2)+pol0(4)" );
  cf.set_fitFromBin( 1+cat*nbins );
  //  cf.set_fitToBin( cat*(nbins+1) );
  cf.set_fitMin( 20 );
  cf.set_fitMax( 200 );
  cf.set_histo_bins( 200 );


  //Option 1: set bins explicitly

  /*
  //qcd
  //  const int nbins=9;
  //  double a_bins[nbins]   = {20.,22.5,25.,27.5,30.,35.,40.,50.};

  //w+jets
  //const int nbins=9;
  //double a_bins[nbins] = {20.,22.5,25.,27.5,30.,35.,40.,50.,60.};

  std::vector<double> bins;
  for (int i=0; i<(nbins-1); i++) bins.push_back( (a_bins[i+1]+a_bins[i])/2 ); //for now
  bins.push_back( a_bins[nbins-1]+  (a_bins[nbins-1]-a_bins[nbins-2])/2 );
  cf.set_bin_centers( bins );
  */

  //use histo to set bin centers
  cf.set_bin_centers( fname , wname, nbins ); //fitFromBin  needs to be set BEFORE

  cf.fitHisto();

  TGraph *g_fit_input=cf.returnFitInputGraph(); //the input to the fit: data points in the range given
  TF1 *f_fit=cf.returnFitForm();                //the fit result (function)
  TGraphAsymmErrors *g_fit=cf.returnFitGraph();              //the fit result binned (histo)

  TH1D *h_fit=cf.returnFitHisto();              //the fit result binned (histo)
  TH1D *h_fit_lo=cf.returnFitHisto(-1);         //the fit result binned (histo), lower uncertainty band
  TH1D *h_fit_hi=cf.returnFitHisto(+1);         //the fit result binned (histo), higher uncertainty band

  TString pfile=proc+"_cat_"; pfile+=cat; pfile+="_";

  //plot fit
  TCanvas *c2=new TCanvas();
  //  h->Draw("E");     //to set the axis
  //  h->GetYaxis()->SetRangeUser(0.4,2.0);
  //  g_fit->SetMarkerStyle(21);
  f_fit->SetTitle("");
  f_fit->Draw();
  f_fit->GetXaxis()->SetTitle("p_{T} [GeV]");
  f_fit->GetYaxis()->SetTitle("Fake factor");
  f_fit->SetMinimum(0);
  f_fit->SetMaximum(f_fit->GetMaximum()*1.6);
  //  //  h_fit->SetLineWidth(3);
  //  h_fit->SetLineColor(602);
  //  h_fit->Draw("E same");
  g_fit->SetLineColor(602);
  g_fit->Draw("Ez same");
  g_fit_input->SetLineWidth(2);
  g_fit_input->Draw("P same");

  h_fit_lo->SetLineColor(kGreen);
  h_fit_hi->SetLineColor(kGreen+3);
  h_fit_lo->SetLineWidth(2);
  h_fit_hi->SetLineWidth(2);
  h_fit_lo->Draw("same");
  h_fit_hi->Draw("same");

  c2->SaveAs(pfile+"fit.png");

  //plot error fits
  const int NTOYS=200;
  TF1 *f[NTOYS]; //NTOYS, set in CustomFit.C !!!
  for (int i=0; i<NTOYS; i++) f[i]=cf.returnErrFits(i);

  TCanvas *c_err=new TCanvas();
  f_fit->SetLineWidth(3);
  f_fit->SetLineColor(kBlack);
  f_fit->Draw();

  for (unsigned i=0; i<NTOYS; i++){
    f[i]->Draw("same");
  }
  f_fit->Draw("same"); //draw again on top
  g_fit_input->Draw("P same");

  c_err->SaveAs(pfile+"err.png");

}

#ifndef __CINT__
int main(int argc, char* argv[]) {
  int channel = argc > 1 ? atoi(argv[1]) : 0;
  TString file = argc > 2 ? argv[2] : "FF_corr_QCD_MCsum_noGen.root";

  if ( file=="all" ){
    for (int i=0; i<4; i++) test(i, "FF_corr_QCD_MCsum_noGen.root");
    for (int i=0; i<4; i++) test(i, "FF_corr_Wjets_MCsum_noGen.root");
    for (int i=0; i<4; i++) test(i, "FF_TT_J_only_SR.root");
  } else{
    if (channel<10) test(channel, file);
    else{
      for (int i=0; i<4; i++) test(i, file);
    }
  }

}
#endif
