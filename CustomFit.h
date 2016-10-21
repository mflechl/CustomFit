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
  void set_bins( std::vector<double> m_bins ){ 
    this->bins = m_bins; 
  }
  void set_fitFromBin( int m_fitFromBin ){ 
    this->fitFromBin = m_fitFromBin;
  }
  void set_fitToBin( int m_fitToBin ){ 
    this->fitToBin = m_fitToBin;
  }

  void getFitHisto();
  void createTestHisto();
  TH1D* returnFitHisto(){ return this->h_out; }
  TH1D* returnInputHisto(){ return this->h_in; }
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

  private:
  double std_dev( std::vector<double> v );

  TH1D *h_in;
  TH1D *h_out;

  TString fitFunc;
  std::vector<double> bins;
  int fitFromBin;
  int fitToBin;

  public:
  //    virtual ~CustomFit();

  public:
    ClassDef(CustomFit,0)

};

#endif

