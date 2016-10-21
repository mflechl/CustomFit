
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
  //  for (int i=0; i<nbins; i++) bins.push_back(a_bins[i] );
  for (int i=0; i<(nbins-1); i++) bins.push_back( (a_bins[i+1]+a_bins[i])/2 ); //for now
  bins.push_back( a_bins[nbins-1]+  (a_bins[nbins-1]-a_bins[nbins-2])/2 );
  for (int i=0; i<bins.size(); i++) std::cout << bins.at(i) << " " << flush; std::cout << std::endl;
  cf.set_bin_centers( bins );

  cf.fitHisto();

  TGraph *g_fit=cf.returnFitGraph();

  TCanvas *c2=new TCanvas();
  //  h->Draw("E");     //to set the axis
  //  h->GetYaxis()->SetRangeUser(0.4,2.0);
  //  g_fit->SetMarkerStyle(21);
  g_fit->Draw("AP");
  //  h->Draw("E same");

}
