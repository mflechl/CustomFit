#include "CustomFit.h"

//ClassImp(CustomFit)

//static const int NTOYS = 100;

CustomFit::CustomFit(){
  this->fitFunc="";
  this->errFunc="";
  //  this->bins.push_back(0);
  this->fitFromBin=0;
  this->fitToBin=3;
  this->histo_bins=100;
  this->err_scale=2;
}

void CustomFit::fitHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: CustomFit::fitHisto: input histo does not exist." << std::endl;

  this->g_fit_input=this->makeFitGraph( h_in );

  this->f_fit = new TF1("f1",this->fitFunc,fitMin,fitMax);
  g_fit_input->Fit(f_fit,"N");

  for (int i=0; i<NTOYS; i++){
    int status;
    this->f_fit_err[i]=this->getFitErr(status);
  }
  this->h_fit = new TH1D("h_fit","",this->histo_bins,this->fitMin,this->fitMax);
  for (int i=1; i<=this->histo_bins; i++){
    double val=f_fit->Eval( h_fit->GetBinCenter(i) );
    this->h_fit->SetBinContent( i , val );
    //    double err = 2*this->std_dev( this->f_fit_err , NTOYS , h_fit->GetBinCenter(i) );
    //    double err = this->err_scale*this->std_dev( this->f_fit_err , NTOYS , h_fit->GetBinCenter(i) ); //standard method
    double err = this->err_scale*this->std_dev( this->f_fit_err , NTOYS , h_fit->GetBinCenter(i) , 1 ); //CL - 16% / 84% value
    if ( err>val ) err=val;
    this->h_fit->SetBinError( i , err );
  }

}

TF1* CustomFit::getFitErr(int &status){
  int fit_status=1;
  TF1* f;

  while (fit_status>0){
    TGraphAsymmErrors *g_toys=this->fluctuateGraph();
    f=new TF1("this->f_fit_err",this->errFunc,fitMin,fitMax);
    fit_status=g_toys->Fit(f,"NQ");
    //  if (fit_status>0) f=this->getFitErr(); //recursive
  }
  status=fit_status;
  return f;
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

//standard deviation
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

//instead of calculating std dev, return the 16% / 84% quantile
void CustomFit::std_dev( std::vector<double> v, double& err_lo, double& err_hi ){
  std::sort(v.begin(),v.end());
  int i_lo=(int)(0.16*v.size());
  int i_hi=(int)(0.84*v.size());
  float val_lo=v.at(i_lo);
  float val_hi=v.at(i_hi);
  //  err_lo=fabs(val_lo-val_hi)/(val_lo+val_hi); //rel symm error
  err_lo=fabs(val_lo-val_hi)/2;
  err_hi=err_lo;

  //  std::cout << "DEBUG : std_dev CL " << i_lo << " " << i_hi << " " << val_lo << " " << val_hi << " " << err_lo << " " << err_hi << std::endl;
}

//double CustomFit::std_dev( TF1* f[] , const unsigned fsize ,  const double val , const double binc ){
double CustomFit::std_dev( TF1* f[] , const unsigned fsize ,  const double val , const int cl){
  std::vector<double> v;
  for (unsigned i=0; i<fsize; i++){
    double y=f[i]->Eval(val);
    if ( !cl && y<0 ) y=0;
    //    if (y>2*binc) y=2*binc;
    v.push_back( y  );
  }

  double ret;
  if ( cl ) this->std_dev( v , ret , ret );
  else ret=this->std_dev( v );
  return ret;
}

