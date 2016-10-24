#include "CustomFit.h"

ClassImp(CustomFit)

static const int NTOYS = 100;

CustomFit::CustomFit(){
  this->fitFunc="p2";
  //  this->bins.push_back(0);
  this->fitFromBin=0;
  this->fitToBin=3;
  this->histo_bins=100;
}

void CustomFit::fitHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: CustomFit::fitHisto: input histo does not exist." << std::endl;

  this->g_fit_input=this->makeFitGraph( h_in );

  this->f_fit = new TF1("f1",this->fitFunc,fitMin,fitMax);
  //  g_fit_input->Fit(fitFunc,"","",fitMin,fitMax);
  g_fit_input->Fit(f_fit,"N");

  TF1 *f_toys[NTOYS];
  for (int i=0; i<NTOYS; i++){
    TGraphAsymmErrors *g_toys=this->fluctuateGraph();
    f_toys[i]=new TF1("f_toys",this->fitFunc,fitMin,fitMax);
    g_toys->Fit(f_toys[i],"NQ");
  }
  //  std::cout << f_toys[3]->Eval(27.2) << " " << f_toys[37]->Eval(27.2) << std::endl;

  //  this->h_fit = new TH1D("h_fit","",this->histo_bins,this->bin_centers.front(),this->bin_centers.back());
  this->h_fit = new TH1D("h_fit","",this->histo_bins,this->fitMin,this->fitMax);
  for (int i=1; i<=this->histo_bins; i++){
    this->h_fit->SetBinContent( i , f_fit->Eval( h_fit->GetBinCenter(i) ) );
    double err = 2*this->std_dev( f_toys , NTOYS , h_fit->GetBinCenter(i) );
    this->h_fit->SetBinError( i , err );
  }
  //  this->h_fit->SetDirectory(0);

}

TGraphAsymmErrors* CustomFit::fluctuateGraph(){
  TGraphAsymmErrors *g = (TGraphAsymmErrors*)this->g_fit_input->Clone((TString)this->g_fit_input->GetName()+"_toy");

  Int_t nbins = g->GetN();

  double x,y;
  double el,eh,err;
  for (int ib=0; ib<nbins; ib++){
    g->GetPoint(ib,x,y);
    el=g->GetErrorYlow(ib);
    eh=g->GetErrorYlow(ib);
    err=(el+eh)/2.;
    double newval = this->rand.Gaus(y,err);
    g->SetPoint(ib,x,newval);
  }
  return g;
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
    //    std::cout << i << " " << x[i] << " " << y[i] << std::endl;
  }

  g=new TGraphAsymmErrors( bin_centers.size() , x , y , 0 , 0 , ey_d , ey_u );

  return g;
}

double CustomFit::std_dev( const std::vector<double> v ){
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

double CustomFit::std_dev( TF1* f[] , const unsigned fsize ,  const double val ){
  std::vector<double> v;
  for (unsigned i=0; i<fsize; i++) v.push_back( f[i]->Eval(val)  );

  return this->std_dev( v );
}

