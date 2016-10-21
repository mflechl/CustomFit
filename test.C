
void test(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L CustomFit.C+");
  CustomFit cf;

  int ret=cf.setInputHisto( "FF_corr_Wjets_MCsum_noGen.root" , "c_t" );
  TH1D *h=cf.returnInputHisto();
  if ( ret != 0 ) return;

  cf.set_fitFunc( "p2" );
  
  cf.set_fitFromBin( 1 );
  cf.set_fitToBin( 9 );

  //qcd
  //  const int nbins=8;
  //  double a_bins[nbins]   = {20.,22.5,25.,27.5,30.,35.,40.,50.};

  //w+jets
  const int nbins=9;
  double a_bins[nbins] = {20.,22.5,25.,27.5,30.,35.,40.,50.,60.};

  std::vector<double> bins;
  for (int i=0; i<nbins; i++) bins.push_back(a_bins[i] );
  cf.set_bins( bins );

  cf.getFitHisto();

  TH1D *h2=cf.returnFitHisto();

  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  //  h->GetYaxis()->SetRangeUser(0.4,2.0);
  h2->Draw("same LP");
  h->Draw("E same");

}
