#include "CustomFit.h"

ClassImp(CustomFit)

CustomFit::CustomFit(){
  this->fitFunc="p2";
  //  this->bins.push_back(0);
  this->fitFromBin=0;
  this->fitToBin=3;
}

void CustomFit::fitHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: CustomFit::fitHisto: input histo does not exist." << std::endl;

  this->g_fit=this->makeFitGraph( h_in );

  this->h_out = (TH1D*)this->h_in->Clone(this->h_in->GetName()+(TString)"_smoothed");
  this->h_out->SetDirectory(0);

  Int_t nbins = this->h_out->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++){
    //    double x = this->h_out->GetBinCenter(ib);
    double smooth_val = this->h_out->GetBinContent(ib);
    this->h_out->SetBinContent(ib,smooth_val);
  }
}

TGraphAsymmErrors* CustomFit::makeFitGraph(TH1D* h_in){

  TGraphAsymmErrors* g;
  int nbins=bin_centers.size();  

  double *x=new double[nbins];
  double *y=new double[nbins];
  double *ey_u=new double[nbins];
  double *ey_d=new double[nbins];

  for (int i=0; i<nbins; i++){
    x[i]=bin_centers.at(i);
    y[i]=h_in->GetBinContent( this->fitFromBin+i );
    ey_u[i]=h_in->GetBinErrorUp( this->fitFromBin+i );
    ey_d[i]=h_in->GetBinErrorLow( this->fitFromBin+i );
    std::cout << i << " " << x[i] << " " << y[i] << std::endl;
  }

  g=new TGraphAsymmErrors( bin_centers.size() , x , y , 0 , 0 , ey_d , ey_u );

  return g;
}

double CustomFit::std_dev( std::vector<double> v ){
  unsigned size=v.size();

  double mean=0;
  for (unsigned i=0; i<size; i++) mean+=v.at(i)/size;

  double std_dev=0;
  for (unsigned i=0; i<size; i++){
    double diff = v.at(i) - mean;
    std_dev+=diff*diff;
  }
  std_dev/=size;
  std_dev=sqrt(std_dev);
 
  return std_dev;
}

void CustomFit::createTestHisto(){

  TF1 *p1 = new TF1("p1","[0]+x*[1]+x*x*[2]");
  p1->SetParameters(5., 10., 0.);

  //  std::cout << p1->Eval(-5) << "\t" << p1->Eval(0) << "\t" << p1->Eval(5) << std::endl;

  this->h_in = new TH1D("h_in","",100,1,11);
  this->h_in->FillRandom("p1",1000);
  this->h_in->SetMarkerColor(kBlack);
  this->h_in->SetLineColor(kBlack);

}

