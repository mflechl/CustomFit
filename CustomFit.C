#include "CustomFit.h"

//ClassImp(CustomFit)

//static const int NTOYS = 100;

CustomFit::CustomFit(){
  this->fitFunc="";
  this->errFunc="";
  //  this->bins.push_back(0);
  this->fitFromBin=-1;
  //  this->fitToBin=-1;
  this->histo_bins=100;
  this->err_scale=2;
  this->err_cl=1;
  this->histMaxFrac=1;
  this->smoothFrac=-1;
  this->smoothMode="spline3"; //spline3 or simple
  this->smoothParam=2;
  this->smoothExp=0.7;
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
  //  this->h_fit = new TH1D("h_fit","",this->histo_bins,this->fitMin,this->fitMax);
  //  int nbins=this->histo_bins*this->histMaxFrac;  
  int nbins=this->histo_bins;  
  this->h_fit = new TH1D("h_fit","",nbins,this->fitMin,this->fitMax);
  this->h_fit_lo = new TH1D("h_fit_lo","",nbins,this->fitMin,this->fitMax);
  this->h_fit_hi = new TH1D("h_fit_hi","",nbins,this->fitMin,this->fitMax);

  double *x=new double[nbins]; double *y=new double[nbins];
  double *ey_lo=new double[nbins]; double *ey_hi=new double[nbins];

  //info used for spline
  const int np=4;
  double xx[np]={0}; xx[0]=-1;
  double yy[np]={0};
  TSpline3 *spline3;

  //create spline object if set
  if (this->smoothMode=="spline3" && this->smoothFrac>0){
    for (int i=1; i<=nbins; i++){
      double xv=h_fit->GetBinCenter(i);

      if( xv>this->fitMax*this->smoothFrac && xv <= this->fitMax*this->histMaxFrac ){
	if ( xx[0]<0 ){
	  xx[0]=h_fit->GetBinCenter(i-2); yy[0]=f_fit->Eval( xx[0] );
	  xx[1]=xv; yy[1]=f_fit->Eval( xv );
	}
	xx[np-2]=xv; yy[np-2]=f_fit->Eval( xv );
	xx[np-1]=h_fit->GetBinCenter(i+2); yy[np-1]=yy[np-2];
      }

    }
    spline3=new TSpline3("sp3",xx,yy,np);
  }

  //  for (int i=1; i<=this->histo_bins; i++){
  float smooth_term=0;
  for (int i=1; i<=nbins; i++){
    double xv=h_fit->GetBinCenter(i);
    double spline_val=-1;
    //to just continue after a certain value:
    if ( xv > this->fitMax*this->histMaxFrac ){
      double val=this->h_fit->GetBinContent(i-1);
      this->h_fit->SetBinContent( i , val );
      x[i-1]=xv;               y[i-1]=y[i-2];
      ey_lo[i-1]=ey_lo[i-2];   ey_hi[i-1]=ey_hi[i-2];
      this->h_fit->SetBinError( i , this->h_fit->GetBinError(i-1) );
      this->h_fit_lo->SetBinContent( i , val-ey_lo[i-1] );
      this->h_fit_hi->SetBinContent( i , val+ey_hi[i-1] );
      continue;
      //    } else if( this->smoothMode=="simple" && this->smoothFrac>0 && xv>this->fitMax*this->smoothFrac ){
    } else if( this->smoothFrac>0 && xv>this->fitMax*this->smoothFrac ){
      if( this->smoothMode=="simple" ){
	float a=this->fitMax*this->smoothFrac;
	float b=this->fitMax*this->histMaxFrac;

	const float s_s=f_fit->Eval( b ) - f_fit->Eval( a );
	float s_d=fabs( (xv-a)/(b-a) -0.5 ) -0.5;
	float sign_sd=(s_d > 0) - (s_d < 0);
	
	smooth_term=-s_s/this->smoothParam*pow(fabs(s_d),this->smoothExp)*sign_sd;
      } else if ( this->smoothMode=="spline3" ){
	spline_val=spline3->Eval( xv );
      }

    }
    double val=f_fit->Eval( xv ) + smooth_term;
    if ( spline_val>0 ){
      val=spline_val;
    }
    this->h_fit->SetBinContent( i , val );
    double err = this->err_scale*this->std_dev( this->f_fit_err , NTOYS , xv , ey_lo[i-1], ey_hi[i-1], this->err_cl );
    //    double err = this->err_scale*this->std_dev( this->f_fit_err , NTOYS , xv , ey_lo[i-1], ey_hi[i-1], 1 ); //CL - 16% / 84% value
    //double err = this->err_scale*this->std_dev( this->f_fit_err , NTOYS , xv , ey_lo[i-1], ey_hi[i-1]); //std dev
    ey_lo[i-1]*=this->err_scale;
    ey_hi[i-1]*=this->err_scale;
    x[i-1]=xv;
    y[i-1]=val;
    if ( err>val ) err=val;
    this->h_fit->SetBinError( i , err );
    this->h_fit_lo->SetBinContent( i , val-ey_lo[i-1] );
    this->h_fit_hi->SetBinContent( i , val+ey_hi[i-1] );


  } //end loop over x bins

  /*
  if (this->smoothMode=="spline3"){
    TCanvas *tc=new TCanvas();
    spline3->Draw();
    gPad->SaveAs("x.png");
  }
  */
  

  this->g_fit = new TGraphAsymmErrors( nbins , x , y , 0 , 0 , ey_lo , ey_hi );

}

TF1* CustomFit::getFitErr(int &status){
  int fit_status=11;
  TF1* f;
  //  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);

  while (fit_status>10){
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

  if (this->fitFromBin < 0){
    std::cout << "ERROR: fitFromBin has not been set!" << std::endl;
    return 0;
  }

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
double CustomFit::std_dev( std::vector<double> v, double& err_lo, double& err_hi , double central_value){
  std::sort(v.begin(),v.end());
  int i_lo=(int)(0.16*v.size());
  int i_hi=(int)(0.84*v.size());
  float val_lo=v.at(i_lo);
  float val_hi=v.at(i_hi);

  err_lo=central_value-val_lo;
  err_hi=val_hi-central_value;

 if (err_lo<err_hi/10.) err_lo=err_hi/10.;
 if (err_hi<err_lo/10.) err_hi=err_lo/10.;

  return fabs(val_lo-val_hi)/2; //symm. error
  //  std::cout << "DEBUG : std_dev CL " << i_lo << " " << i_hi << " " << val_lo << " " << val_hi << " " << err_lo << " " << err_hi << std::endl;
}

//double CustomFit::std_dev( TF1* f[] , const unsigned fsize ,  const double val , const double binc ){
double CustomFit::std_dev( TF1* f[] , const unsigned fsize ,  const double val , double& err_lo, double& err_hi, const int cl){
  std::vector<double> v;
  for (unsigned i=0; i<fsize; i++){
    double y=f[i]->Eval(val);
    //    if ( !cl && y<0 ) y=0;
    if ( y<0 ) y=0;
    //    if (y>2*binc) y=2*binc;
    v.push_back( y  );
  }

  double ret;
  if ( cl ) ret=this->std_dev( v , err_lo , err_hi , this->f_fit->Eval(val) );
  else{ ret=this->std_dev( v ); err_lo=ret; err_hi=ret; }

  return ret;
}

