#ifndef __CustomFit__
#define __CustomFit__

//#include <algorithm>
#include <stdlib.h>
#include<math.h>
#include<iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;

class CustomFit
{

  public:
  CustomFit();
  void set_fitFunc( TString m_fitFunc ){ 
    this->fitFunc = m_fitFunc; 
  }
  void set_bin_centers( std::vector<double> m_bin_centers ){
    this->bin_centers = m_bin_centers; 
  }
  void set_fitFromBin( int m_fitFromBin ){ 
    this->fitFromBin = m_fitFromBin;
  }
  void set_fitToBin( int m_fitToBin ){ 
    this->fitToBin = m_fitToBin;
  }
  void set_fitMin( float m_fitMin ){ 
    this->fitMin = m_fitMin;
  }
  void set_fitMax( float m_fitMax ){ 
    this->fitMax = m_fitMax;
  }
  void set_histo_bins( int m_histo_bins ){ 
    this->histo_bins = m_histo_bins;
  }

  void fitHisto();
  TGraphAsymmErrors* returnFitInputGraph(){ return this->g_fit_input; }
  TH1D* returnInputHisto(){ return this->h_in; }
  TH1D* returnFitHisto(){ return this->h_fit; }
  TF1* returnFitForm(){ return this->f_fit; }
  int setInputHisto(TString fname , TString hname ){ 
    TFile *f=new TFile( fname );
    this->h_in=(TH1D*) f->Get( hname );
    if ( this->h_in ) return 0;
    else{ cerr << "ERROR: " << hname << " not found in file " << fname << endl; return 1; }
  }
  int setInputHisto( TH1D* m_h ){ 
    this->h_in=(TH1D*) m_h->Clone();
    if ( this->h_in ) return 0;
    else{ cerr << "ERROR: Could not assign histogram with address " << m_h << endl; return 1; }
  }
  TGraphAsymmErrors* fluctuateGraph();

  private:
  TGraphAsymmErrors* makeFitGraph(TH1D* h_in);
  double std_dev( const std::vector<double> v );
  double std_dev( TF1* f[] , const unsigned fsize , double val );

  TH1D *h_in;
  TH1D *h_fit;
  TF1 *f_fit;
  TGraphAsymmErrors* g_fit_input;

  TString fitFunc;
  std::vector<double> bin_centers;
  int fitFromBin;
  int fitToBin;
  float fitMin;
  float fitMax;
  int histo_bins;

  TRandom3 rand;

  public:
  //    virtual ~CustomFit();

  public:
    ClassDef(CustomFit,0)

};

#endif

