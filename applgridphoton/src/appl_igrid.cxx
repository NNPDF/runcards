
//  appl_igrid.cxx        

//   grid class - all the functions needed to create and 
//   fill the grid from an NLO calculation program. 
//  
//  Copyright (C) 2007 Mark Sutton (sutt@hep.ucl.ac.uk)    

//   $Id: appl_igrid.cxx, v1.00 2007/10/16 17:01:39 sutt

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>


#include "appl_igrid.h"
#include "appl_grid/appl_grid.h"
#include "sparse_matrix3d_root_interface.hpp"

#include "Cache.h"
#include "Directory.h"
#include "hoppet_init.h"
#include "SparseMatrix3d.h"


#include "TFile.h"
#include "TObject.h"
#include "TObjString.h"
#include "TH3D.h"
#include "TVectorT.h"

#include "TFileString.h"

namespace
{

// fast -1^i function
int pow1(int i) { return ( 1&i ? -1 : 1 ); }

// fast factorial
double fac(int i) {
  int j;
  static int ntop = 4;
  static double f[34] = { 1, 1, 2, 6, 24 };
  if ( i<0 )  {  std::cerr << "igrid::fac() negative input"  << std::endl; return 0;  }
  if ( i>33 ) {  std::cerr << "igrid::fac() input too large" << std::endl; return 0;  }
  while ( ntop<i ) {
    j=ntop++;
    f[ntop]=f[j]*ntop;
  }
  return f[i];
}

// although not the best way to calculate interpolation coefficients,
// it may be the best for our use, where the "y" values of the nodes
// are not yet defined at the time of evaluation.
double fI(int i, int n, double u) {
  if ( n==0 && i==0 )     return 1.0;
  if ( std::fabs(u-i)<1e-8 ) return 1.0;
  //    if ( std::fabs(u-n)<u*1e-8 ) return 1.0;
  double product = pow1(n-i) / ( fac(i)*fac(n-i)*(u-i) );
  for( int z=0 ; z<=n ; z++ )  product *= (u-z);
  return product;
}

}

// splitting function code

void Splitting(const double& x, const double& Q, double* f);


// pdf reweighting
// bool   igrid::m_reweight   = false;
// bool   igrid::m_symmetrise = false;

// variable tranformation parameters
double appl::igrid::transvar = 5;

appl::igrid::exception::exception(const std::string& s)
{
    std::cerr << s << std::endl;
}

appl::igrid::exception::exception(std::ostream&)
{
    std::cerr << std::endl;
}

// standard constructor
appl::igrid::igrid(int NQ2, double Q2min, double Q2max, int Q2order, 
		   int Nx,  double xmin,  double xmax,  int xorder, 
		   std::string transform, int Nproc, bool disflag ):
  mfy(0),   mfx(0),
  m_Ny1(Nx),   m_Ny2( disflag ? 1 : Nx ),  m_yorder(xorder), 
  m_Ntau(NQ2), m_tauorder(Q2order), 
  m_Nproc(Nproc), 
  m_transform(transform), 
  //  m_transvarlocal(m_transvar),
  m_transvar(transvar),
  m_reweight(false),
  m_symmetrise(false), 
  m_optimised(false),
  m_weight(0),
  m_fg1(0),     m_fg2(0),  
  m_fsplit1(0), m_fsplit2(0),
  m_alphas(0),
  m_DISgrid(disflag)   
{
  //  std::cout << "igrid::igrid() transform=" << m_transform << std::endl;
  init_fmap();
  if ( m_fmap.find(m_transform)==m_fmap.end() ) throw exception("igrid::igrid() transform " + m_transform + " not found\n");
  mfx = m_fmap.find(m_transform)->second.mfx;
  mfy = m_fmap.find(m_transform)->second.mfy;
  


  double ymin1 = fy(xmax);
  double ymax1 = fy(xmin);

  m_y1min  = ymin1 ;
  m_y1max  = ymax1 ;
  m_y2min  = ymin1 ;
  m_y2max  = ymax1 ;

  if ( m_DISgrid ) m_y2min = m_y2max = 1;

  if ( m_Ny1>1 ) m_deltay1 = (m_y1max-m_y1min)/(m_Ny1-1);
  else           m_deltay1 = 0;
  
  if ( m_Ny2>1 ) m_deltay2 = (m_y2max-m_y2min)/(m_Ny2-1);
  else           m_deltay2 = 0;
  
  double taumin=ftau(Q2min);
  double taumax=ftau(Q2max);
  m_taumin   = taumin;
  m_taumax   = taumax;
  if ( m_Ntau>1 ) m_deltatau = (taumax-taumin)/(m_Ntau-1);
  else            m_deltatau = 0;

  if ( m_Ny1-1<m_yorder ) { 
    std::cerr << "igrid() not enough nodes for this interpolation order Ny1=" << m_Ny1 
	      << "\tyorder=" << m_yorder << std::endl;
 
    while ( m_Ny1-1<m_yorder ) m_yorder--;
  } 

  if ( !m_DISgrid ) { 
    if ( m_Ny2-1<m_yorder ) { 
      std::cerr << "igrid() not enough nodes for this interpolation order Ny2=" << m_Ny2 
		<< "\tyorder=" << m_yorder << std::endl;
      
      while ( m_Ny2-1<m_yorder ) m_yorder--;
    } 
  }
  
  if ( m_Ntau-1<m_tauorder ) { 
    std::cerr << "igrid() not enough nodes for this interpolation order Ntau=" << m_Ntau 
	      << "\ttauorder=" << m_tauorder << std::endl;
    
    while ( m_Ntau-1<m_tauorder ) m_tauorder--;
  } 
  
  m_weight = new SparseMatrix3d*[m_Nproc];
  construct();
}





// copy constructor
appl::igrid::igrid(const appl::igrid& g) : 
  mfy(0),  mfx(0),  
  m_Ny1(g.m_Ny1),     
  m_y1min(g.m_y1min),     m_y1max(g.m_y1max),     m_deltay1(g.m_deltay1),   
  m_Ny2(g.m_Ny2),     
  m_y2min(g.m_y2min),     m_y2max(g.m_y2max),     m_deltay2(g.m_deltay2),   
  m_yorder(g.m_yorder),   
  m_Ntau(g.m_Ntau), 
  m_taumin(g.m_taumin), m_taumax(g.m_taumax), m_deltatau(g.m_deltatau), m_tauorder(g.m_tauorder), 
  m_Nproc(g.m_Nproc),
  m_transform(g.m_transform), 
  //  m_transvarlocal(g.m_transvarlocal),
  m_transvar(g.m_transvar),
  m_reweight(g.m_reweight),
  m_symmetrise(g.m_symmetrise),
  m_optimised(g.m_optimised),
  m_weight(NULL),
  m_fg1(NULL),     m_fg2(NULL),
  m_fsplit1(NULL), m_fsplit2(NULL),
  m_alphas(NULL)   
{
  init_fmap();
  if ( m_fmap.find(m_transform)==m_fmap.end() ) throw exception("igrid::igrid() transform " + m_transform + " not found\n");
  mfx = m_fmap.find(m_transform)->second.mfx;
  mfy = m_fmap.find(m_transform)->second.mfy;


  m_weight = new SparseMatrix3d*[m_Nproc];
  for( int ip=0 ; ip<m_Nproc ; ip++ )   m_weight[ip] = new SparseMatrix3d(*g.m_weight[ip]);
  //  construct();
}


// read from a file 
appl::igrid::igrid(void* pf, const std::string& s) :
  mfy(0),  mfx(0),  
  m_Ny1(0),   m_y1min(0),   m_y1max(0),   m_deltay1(0),   
  m_Ny2(0),   m_y2min(0),   m_y2max(0),   m_deltay2(0),   
  m_yorder(0),   
  m_Ntau(0), m_taumin(0), m_taumax(0), m_deltatau(0), m_tauorder(0), 
  m_Nproc(0),
  m_transform(""), 
  //  m_transvarlocal(m_transvar),
  m_transvar(transvar),
  m_reweight(false),
  m_symmetrise(false),
  m_optimised(false),
  m_weight(NULL), 
  m_fg1(NULL),     m_fg2(NULL),
  m_fsplit1(NULL), m_fsplit2(NULL),    
  m_alphas(NULL) 
{ 
  //  std::cout << "igrid::igrid()" << std::endl;
  
  // using the title of a TH1D because I don't know 
  // how else to save a std::string in a root file 
  // TH1D* _transform = (TH1D*)f.Get((s+"/Transform").c_str());
  // m_transform = _transform->GetTitle();
  // delete _transform;
 
  // get the name of the transform pair
  TFile& f = *static_cast <TFile*> (pf);
  TFileString _tag = *(TFileString*)f.Get((s+"/Transform").c_str());
  m_transform = _tag[0];

  init_fmap();
  if ( m_fmap.find(m_transform)==m_fmap.end() ) throw exception("igrid::igrid() transform " + m_transform + " not found\n");
  mfx = m_fmap.find(m_transform)->second.mfx;
  mfy = m_fmap.find(m_transform)->second.mfy;


  // delete _transform;
  
  // retrieve setup parameters 
  TVectorT<double>* setup=(TVectorT<double>*)f.Get((s+"/Parameters").c_str());
  //  f.GetObject((s+"/Parameters").c_str(), setup);

  // NB: round integer variables to nearest integer 
  //     in case (unlikely) truncation error during 
  //     conversion to double when they were stored
  m_Ny1      = int((*setup)(0)+0.5);  
  m_y1min    = (*setup)(1);
  m_y1max    = (*setup)(2);

  m_Ny2      = int((*setup)(3)+0.5);  
  m_y2min    = (*setup)(4);
  m_y2max    = (*setup)(5);

  m_yorder   = int((*setup)(6)+0.5);

  m_Ntau     = int((*setup)(7)+0.5);
  m_taumin   = (*setup)(8);
  m_taumax   = (*setup)(9);
  m_tauorder = int((*setup)(10)+0.5);

  m_transvar = (*setup)(11);

  m_Nproc    = int((*setup)(12)+0.5);

  m_reweight   = ((*setup)(13)!=0 ? true : false );
  m_symmetrise = ((*setup)(14)!=0 ? true : false );
  m_optimised  = ((*setup)(15)!=0 ? true : false );
  m_DISgrid    =  ( setup->GetNoElements()>16 && (*setup)[16]!=0 ) ;


  delete setup;

  //  std::cout << "igrid::igrid() read setup" << std::endl;

  // create grids
  m_deltay1   = (m_y1max-m_y1min)/(m_Ny1-1);
  m_deltay2   = (m_y2max-m_y2min)/(m_Ny2-1);

  m_deltatau = (m_taumax-m_taumin)/(m_Ntau-1);
  
  //  int rawsize=0;
  //  int trimsize=0;

  m_weight = new SparseMatrix3d*[m_Nproc];

  for( int ip=0 ; ip<m_Nproc ; ip++ ) {
    char name[128];  sprintf(name,"/weight[%i]", ip );
    // get storage histogram
    TH3D* htmp = (TH3D*)f.Get((s+name).c_str()); 
 
    //    std::cout << "igrid::igrid() read " << name << std::endl;

    // create grid
    m_weight[ip]=sparse_matrix3d_from_th3d(*htmp);

    // save some space
    m_weight[ip]->trim();
    // delete storage histogram
    delete htmp;

    // DON'T trim on reading unless the user wants it!!
    // he can trim himself if need be!! 
    // rawsize += m_weight[ip]->size();
    // m_weight[ip]->trim(); // trim the grid and do some book keeping
    // trimsize += m_weight[ip]->size();
  }
}




// constructor common internals 
void appl::igrid::construct() 
{
  // Initialize histograms representing the weight grid
  for( int ip=0 ; ip<m_Nproc ; ip++ ) {
    m_weight[ip] = new SparseMatrix3d(m_Ntau, m_taumin,   m_taumax,    
    				      m_Ny1,  m_y1min,    m_y1max, 
    				      m_Ny2,  m_y2min,    m_y2max ); 
    
  }  
}


// destructor
appl::igrid::~igrid() {
  deletepdftable();
  deleteweights();
}



// clean up internal pdf table (could be moved to local variable)
void appl::igrid::deletepdftable() { 

  //  std::cout << "deleting pdf tables" << std::endl;

  if ( m_fg1 ) { 
    for ( int i=0 ; i<m_Ntau ; i++ ) {
      for ( int j=0 ; j<Ny1() ; j++ )  delete[] m_fg1[i][j];
      delete[] m_fg1[i];
    }
    delete[] m_fg1;
    m_fg1=NULL;
  }
  
  //  std::cout << "deleting splitting tables" << std::endl;

  if ( m_fsplit1 ) { 
    for ( int i=0 ; i<m_Ntau ; i++ ) {
      for ( int j=0 ; j<Ny1() ; j++ )  delete[] m_fsplit1[i][j];
      delete[] m_fsplit1[i];
    }
    delete[] m_fsplit1;
    m_fsplit1=NULL;
  }

  // if grid is symmetric, then don't need to deallocate 
  // x2 pdf values and splitting functions 
  if ( isSymmetric() ) { 
    m_fg2=NULL;
    m_fsplit2=NULL;
  }
  else { 
    if ( m_fg2 ) { 
      for ( int i=0 ; i<m_Ntau ; i++ ) {
	for ( int j=0 ; j<Ny2() ; j++ )  delete[] m_fg2[i][j];
	delete[] m_fg2[i];
      }
      delete[] m_fg2;
      m_fg2=NULL;
    }
    
    
    //  std::cout << "deleting splitting tables" << std::endl;
    if ( m_fsplit2 ) { 
      for ( int i=0 ; i<m_Ntau ; i++ ) {
	for ( int j=0 ; j<Ny2() ; j++ )  delete[] m_fsplit2[i][j];
	delete[] m_fsplit2[i];
      }
      delete[] m_fsplit2;
      m_fsplit2=NULL;
    }
  }

  //  std::cout << "delete alphas table" << std::endl;

  if ( m_alphas ) { 
    delete[] m_alphas;
    m_alphas=NULL;
  }
}





// delete internal grids
void appl::igrid::deleteweights() { 
  if ( m_weight ) { 
    for ( int ip=0 ; ip<m_Nproc ; ip++ ) if (m_weight[ip]) delete m_weight[ip];
    delete[] m_weight;
    // m_weight=NULL;
  }
}

int appl::igrid::fk1(double x) const {
  double y = fy(x);
  // make sure we are in the range covered by our binning
  if( y<y1min() || y>y1max() ) {
    if ( y<y1min() ) std::cerr <<"\tWarning: x1 out of range: x=" << x << "\t(y=" << y << ")\tBelow Delx=" << x-fx(y1min());
    else             std::cerr <<"\tWarning: x1 out of range: x=" << x << "\t(y=" << y << ")\tAbove Delx=" << x-fx(y1min());
    std::cerr << " ( " <<  fx(y1max())  << " - " <<  fx(y1min())  << " )"
    << "\ty=" << y << "\tDely=" << y-y1min() << " ( " << y1min() << " - " <<  y1max()  << " )" << std::endl;
    //     cerr << "\t" << m_weight[0]->yaxis() << "\n\t"
    //          << m_weight[0]->yaxis().transform(fx) << std::endl;
  }
  int k = (int)((y-y1min())/deltay1() - (m_yorder>>1)); // fast integer divide by 2
  if ( k<0 ) k=0;
  // shift interpolation end nodes to enforce range
  if(k+m_yorder>=Ny1())  k=Ny1()-1-m_yorder;
  return k;
}

int appl::igrid::fk2(double x) const {
  double y = fy(x);
  // make sure we are in the range covered by our binning
  if( y<y2min() || y>y2max() ) {
    if ( y<y2min() ) std::cerr <<"\tWarning: x2 out of range: x=" << x << "\t(y=" << y << ")\tBelow Delx=" << x-fx(y2min());
    else             std::cerr <<"\tWarning: x2 out of range: x=" << x << "\t(y=" << y << ")\tAbove Delx=" << x-fx(y2min());
    std::cerr << " ( " <<  fx(y2max())  << " - " <<  fx(y2min())  << " )"
     << "\ty=" << y << "\tdely=" << y-y2min() << " ( " << y2min() << " - " <<  y2max()  << " )" << std::endl;
    //     cerr << "\t" << m_weight[0]->yaxis() << "\n\t" << m_weight[0]->yaxis().transform(fx) << std::endl;
  }
  int k = (int)((y-y2min())/deltay2() - (m_yorder>>1)); // fast integer divide by 2
  if ( k<0 ) k=0;
  // shift interpolation end nodes to enforce range
  if(k+m_yorder>=Ny2())  k=Ny2()-1-m_yorder;
  return k;
}

int appl::igrid::fkappa(double Q2) const {
  double tau = ftau(Q2);
  // make sure we are in the range covered by our binning
  if( tau<taumin() || tau>taumax() ) {
    std::cerr << "\tWarning: Q2 out of range Q2=" << Q2
    << "\t ( " << fQ2(taumin()) << " - " << fQ2(taumax()) << " )" << std::endl;
    //      cerr << "\t" << m_weight[0]->xaxis() << "\n\t" << m_weight[0]->xaxis().transform(fQ2) << std::endl;
  }
  int kappa = (int)((tau-taumin())/deltatau() - (m_tauorder>>1)); // fast integer divide by 2
  // shift interpolation end nodes to enforce range
  if(kappa+m_tauorder>=Ntau()) kappa=Ntau()-1-m_tauorder;
  if(kappa<0) kappa=0;
  return kappa;
}


// write to file
void appl::igrid::write(const std::string& name) { 
  Directory d(name);
  d.push();

  // using a TH1D to store the transform pair tag since I don't know how 
  // to write a TString to a root file
  // TH1D* _transform = new TH1D("Transform", m_transform.c_str(), 1, 0, 1);
  //  _transform->Write();

  // write the name of the transform pair
  TFileString("Transform",m_transform).Write();


  TVectorT<double>* setup=new TVectorT<double>(20); // a few spare

  (*setup)(0)  = m_Ny1;
  (*setup)(1)  = m_y1min;
  (*setup)(2)  = m_y1max;

  (*setup)(3)  = m_Ny2;
  (*setup)(4)  = m_y2min;
  (*setup)(5)  = m_y2max;

  (*setup)(6)  = m_yorder;

  (*setup)(7)  = m_Ntau;
  (*setup)(8)  = m_taumin;
  (*setup)(9)  = m_taumax;
  (*setup)(10) = m_tauorder;

  (*setup)(11) = m_transvar;

  (*setup)(12) = m_Nproc;
 
  (*setup)(13) = ( m_reweight   ? 1 : 0 );
  (*setup)(14) = ( m_symmetrise ? 1 : 0 );
  (*setup)(15) = ( m_optimised  ? 1 : 0 );
  (*setup)(16) = ( m_DISgrid    ? 1 : 0 );

  setup->Write("Parameters");

  int igridsize     = 0;
  int igridtrimsize = 0;

  for ( int ip=0 ; ip<m_Nproc ; ip++ ) { 

    char hname[128];
    //    sprintf(hname,"%s[%d]", name.c_str(), ip);
    sprintf(hname,"weight[%d]", ip);

    //    int oldsize = m_weight[ip]->size();

    igridsize += m_weight[ip]->size();
 
    // trim it so that it's quicker to copy into the TH3D
    m_weight[ip]->trim();

    igridtrimsize += m_weight[ip]->size();

    TH3D* h=sparse_matrix3d_to_th3d(*(m_weight[ip]), hname);
    h->SetDirectory(0);
    h->Write();
    delete h; // is this dengerous??? will root try to delete it later? I guess not if we SetDirectory(0)
  }

#if 0
  //    std::cout << name << " trimmed" << std::endl;
  std::cout << name << "\tsize=" << igridsize << "\t-> " << igridtrimsize;
  if ( igridsize ) std::cout << "\t( " << igridtrimsize*100/igridsize << "% )";
  std::cout << std::endl;
#endif

  d.pop();
}





void appl::igrid::fill(const double x1, const double x2, const double Q2, const double* weight) 
{  

  // find preferred vertex for low end of interpolation range
  int k1=fk1(x1);
  int k2=fk2(x2);
  int k3=fkappa(Q2);
  
  double u_y1  = ( fy(x1)-gety1(k1) )/deltay1();
  double u_y2  = ( fy(x2)-gety2(k2) )/deltay2();
  double u_tau = ( ftau(Q2)-gettau(k3))/deltatau();
  
  
  double fillweight, fI_factor;
  
  // cache the interpolation coefficients so only 
  // have to calculate each one once 
  // NB: a (very naughty) hardcoded maximum interpolation order 
  //     of 16 for code simplicity.
  double _fI1[16];
  double _fI2[16];
  double _fI3[16];
  
  for( int i1=0 ; i1<=m_yorder   ; i1++ )  _fI1[i1] = fI(i1, m_yorder, u_y1);
  for( int i2=0 ; i2<=m_yorder   ; i2++ )  _fI2[i2] = fI(i2, m_yorder, u_y2);
  for( int i3=0 ; i3<=m_tauorder ; i3++ )  _fI3[i3] = fI(i3, m_tauorder, u_tau);
  
  double invwfun = 1/(weightfun(x1)*weightfun(x2));
	 
  for( int i3=0 ; i3<=m_tauorder ; i3++ ) { // "interpolation" loop in Q2

    for( int i1=0 ; i1<=m_yorder ; i1++ )   { // "interpolation" loop in x1      

      for( int i2=0 ; i2<=m_yorder ; i2++ ) { // "interpolation" loop in x2
	
	fI_factor = _fI1[i1] * _fI2[i2] * _fI3[i3];

	if ( m_reweight ) fI_factor *= invwfun;

	for( int ip=0 ; ip<m_Nproc ; ip++ ) {
	  
	  fillweight = weight[ip] * fI_factor;

	  // this method only works for grids where the x1, x2 axes are identical 
	  // ranges and binning 
	  //  if ( m_symmetrise ) { 
	  //	// symmetrise the grid
	  //    if ( abs(k1+i1)<abs(k2+i2) )  (*m_weight[ip])(k3+i3, k2+i2, k1+i1) += fillweight;   
	  //    else                          (*m_weight[ip])(k3+i3, k1+i1, k2+i2) += fillweight;  
	  //  } 
	  //  else {

	  //	  std::cout << "weight[" << ip << "]=" << weight[ip] << "\tfillweight=" << fillweight << std::endl;;

	  (*m_weight[ip])(k3+i3, k1+i1, k2+i2) += fillweight;	  
	  //  }	  
	  
	  //	  m_weight[ip]->print();	  

	  //	  m_weight[ip]->fill_fast(k3+i3, k1+i1, k2+i2) += fillweight/wfun;	  
	}
      }
    }
  }
}


void appl::igrid::fill_phasespace(const double x1, const double x2, const double Q2, const double* weight) { 

  int k1=fk1(x1);
  int k2=fk2(x2);
  int k3=fkappa(Q2);

  for ( int ip=0 ; ip<m_Nproc ; ip++ ) (*m_weight[ip])(k3, k1, k2) += weight[ip];

} 

// initialise the transform map - no longer shared between class members
void appl::igrid::init_fmap() {
  if ( m_fmap.size()==0 ) {
    add_transform( "f",  &igrid::_fx,  &igrid::_fy  );
    add_transform( "f0", &igrid::_fx0, &igrid::_fy0  );
    add_transform( "f1", &igrid::_fx1, &igrid::_fy1  );
    add_transform( "f2", &igrid::_fx2, &igrid::_fy2  );
    add_transform( "f3", &igrid::_fx3, &igrid::_fy3  );
    add_transform( "f4", &igrid::_fx4, &igrid::_fy4  );
  }
}

/// add a transform
void appl::igrid::add_transform(const std::string transform, transform_t __fx, transform_t __fy ) {
  if ( m_fmap.find(transform)!=m_fmap.end() ) {
    throw exception("igrid::add_fmap() transform "+transform+" already in std::map");
  }
  m_fmap[transform] = transform_vec( __fx, __fy );
}


double appl::igrid::weightfun(double x) { double n=(1-0.99*x); return std::sqrt(x*x*x)/(n*n*n); }

double appl::igrid::ftau(double Q2) { return std::log(std::log(Q2/0.0625)); }
double appl::igrid::fQ2(double tau) { return 0.0625*std::exp(std::exp(tau)); }

// define all these so that ymin=fy(xmin) rather than ymin=fy(xmax)
double appl::igrid::_fy(double x) const { return std::log(1/x-1); }
double appl::igrid::_fx(double y) const { return 1/(1+std::exp(y)); }

double appl::igrid::_fy0(double x) const { return -std::log(x); }
double appl::igrid::_fx0(double y) const { return  std::exp(-y); }

double appl::igrid::_fy1(double x) const { return std::sqrt(-std::log(x)); }
double appl::igrid::_fx1(double y) const { return std::exp(-y*y); }

double appl::igrid::_fy2(double x) const { return -std::log(x)+m_transvar*(1-x); }
double appl::igrid::_fx2(double y) const {
  // use Newton-Raphson: y = ln(1/x)
  // solve   y - yp - a(1 - exp(-yp)) = 0
  // deriv:  - 1 -a exp(-yp)

  if ( m_transvar==0 )  return std::exp(-y);

  const double eps  = 1e-12;  // our accuracy goal
  const int    imax = 100;    // for safety (avoid infinite loops)

  double yp = y;
  double x, delta, deriv;
  for ( int iter=imax ; iter-- ; ) {
    x = std::exp(-yp);
    delta = y - yp - m_transvar*(1-x);
    if ( std::fabs(delta)<eps ) return x; // we have found good solution
    deriv = -1 - m_transvar*x;
    yp  -= delta/deriv;
  }
  // exceeded maximum iterations
  std::cerr << "_fx2() iteration limit reached y=" << y << std::endl;
  std::cout << "_fx2() iteration limit reached y=" << y << std::endl;
  return std::exp(-yp);
}

double appl::igrid::_fy3(double x) const { return std::sqrt(-std::log10(x)); }
double appl::igrid::_fx3(double y) const { return std::pow(10,-y*y); }

// fastnlo dis transform
double appl::igrid::_fy4(double x) const { return -std::log10(x); }
double appl::igrid::_fx4(double y) const { return  std::pow(10,-y); }

double appl::igrid::gety1(int iy) const { return m_weight[0]->yaxis()[iy]; }

double appl::igrid::gety2(int iy) const { return m_weight[0]->zaxis()[iy]; }

double appl::igrid::gettau(int itau) const { return m_weight[0]->xaxis()[itau]; }

int    appl::igrid::Ny1()      const { return m_weight[0]->yaxis().N(); }
double appl::igrid::y1min()    const { return m_weight[0]->yaxis().min(); }
double appl::igrid::y1max()    const { return m_weight[0]->yaxis().max(); }
double appl::igrid::deltay1()  const { return m_weight[0]->yaxis().delta(); }

int    appl::igrid::Ny2()      const { return m_weight[0]->zaxis().N(); }
double appl::igrid::y2min()    const { return m_weight[0]->zaxis().min(); }
double appl::igrid::y2max()    const { return m_weight[0]->zaxis().max(); }
double appl::igrid::deltay2()  const { return m_weight[0]->zaxis().delta(); }

int    appl::igrid::Ntau()     const { return m_weight[0]->xaxis().N(); }
double appl::igrid::taumin()   const { return m_weight[0]->xaxis().min(); }
double appl::igrid::taumax()   const { return m_weight[0]->xaxis().max(); }
double appl::igrid::deltatau() const { return m_weight[0]->xaxis().delta(); }

void appl::igrid::setuppdf(double (*alphas)(const double&),
			   NodeCache* pdf0,
			   NodeCache* pdf1,
			   int _nloop,
			   double rscale_factor,
			   double fscale_factor,
			   double beam_scale ) 
{

  int nloop = std::fabs(_nloop);

  if ( pdf1==0 ) pdf1 = pdf0;

  bool initialise_hoppet = false;

#ifdef HAVE_HOPPET
  // check if we need to use the splitting function, and if so see if we 
  // need to initialise it again, and do so if required
  if ( fscale_factor!=1 ) {
    if ( pdf1 && pdf1!=pdf0 ) inititialise_hoppet = false;
  }
#endif

  void (*splitting)(const double& , const double&, double* ) = Splitting;
		     
  const int n_tau = Ntau();
  const int n_y1  = Ny1();
  const int n_y2  = Ny2();

  // pdf table
  m_fg1 = new double**[Ntau()];
  if ( isSymmetric() ) m_fg2 = m_fg1;
  else                 m_fg2 = new double**[Ntau()];

  const double invtwopi = 0.5/(M_PI);

  // splitting function table for nlo 
  // factorisation scale dependence
  if ( nloop==1 && fscale_factor!=1 ) { 
    m_fsplit1 = new double**[Ntau()];
    if ( isSymmetric() ) m_fsplit2 = m_fsplit1;
    else                 m_fsplit2 = new double**[Ntau()];
  }

  // alphas table
  m_alphas = new double[Ntau()];

  // allocate tables
  for ( int i=0 ; i<n_tau ; i++ ) {
    
    // pdf table
    m_fg1[i] = new double*[n_y1];
    for ( int j=0 ; j<n_y1 ; j++ ) m_fg1[i][j] = new double[14];
    if ( !isSymmetric() && !isDISgrid() ) { 
      m_fg2[i] = new double*[n_y2];
      for ( int j=0 ; j<n_y2 ; j++ ) m_fg2[i][j] = new double[14];   
    }    
    
    // splitting function table
    if ( nloop==1 && fscale_factor!=1 ) { 
      m_fsplit1[i] = new double*[n_y1];
      for ( int j=0 ; j<n_y1 ; j++ ) m_fsplit1[i][j] = new double[14];
      if ( !isSymmetric() && !isDISgrid() ) { 
	m_fsplit2[i] = new double*[n_y2];
	for ( int j=0 ; j<n_y2 ; j++ ) m_fsplit2[i][j] = new double[14];   
      }
    }
  }  
  
  bool scale_beams = false;
  if ( beam_scale!=1 ) scale_beams = true;

  if ( initialise_hoppet ) hoppet_init::assign( pdf0->pdf() );
  
  // set up pdf grid, splitting function grid 
  // and alpha_s grid
  //  for ( int itau=m_Ntau ; itau-- ; ) {
  for ( int itau=0 ; itau<m_Ntau ; itau++  ) {
    
    double tau = gettau(itau);
    double Q2  = fQ2(tau);
    double Q   = std::sqrt(Q2); 
    
    // alpha_s table 
    m_alphas[itau] = alphas(rscale_factor*Q)*invtwopi;

    //    std::cout << itau << "\ttau " << tau 
    //	      << "\tQ2 " << Q2 << "\tQ " << Q 
    //	      << "\talphas " << m_alphas[itau] << std::endl; 

#if 0
    int iymin1 = m_weight[0]->ymin();
    int iymax1 = m_weight[0]->ymax();
    
    if ( isSymmetric() ) { 
      int iymin2 = iymin1;
      int iymax2 = iymax1;
    }
    else {
      int iymin2 = m_weight[0]->zmin();
      int iymax2 = m_weight[0]->zmax();
    }  
#endif

    // grid not filled for iy<iymin || iy>iymax so no need to 
    // consider outside this range
    
    // y1 tables
    for ( int iy=n_y1 ; iy-- ;  ) { 
      
      double y = gety1(iy);
      double x = fx(y);
      double fun = 1;
      if ( m_reweight ) fun = weightfun(x);
      
      if ( scale_beams ) { 
	x *= beam_scale;
	if ( x>=1 ) { 
	  for ( int ip=0 ; ip<14 ; ip++ ) m_fg1[itau][iy][ip]=0; 
	  if ( nloop==1 && fscale_factor!=1 ) { 
	    for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit1[itau][iy][ip] = 0;
	  }
	  continue; 
	}    
      }

      pdf0->evaluate(x, fscale_factor*Q, m_fg1[itau][iy]);
      
      double invx = 1/x;
      for ( int ip=0 ; ip<14 ; ip++ ) m_fg1[itau][iy][ip] *= invx;
      if ( m_reweight ) for ( int ip=0 ; ip<14 ; ip++ ) m_fg1[itau][iy][ip] *= fun;
      
      // splitting function table
      if ( nloop==1 && fscale_factor!=1 ) { 
	splitting(x, fscale_factor*Q, m_fsplit1[itau][iy]);
	for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit1[itau][iy][ip] *= invx;
	if ( m_reweight ) for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit1[itau][iy][ip] *= fun;
      }
    }
  }

  if ( initialise_hoppet ) hoppet_init::assign( pdf1->pdf() );
  
  if ( ( !isSymmetric() && !isDISgrid() ) || ( pdf1 && pdf1!=pdf0 ) ) {
    
    for ( int itau=0 ; itau<m_Ntau ; itau++  ) {
    
      double tau = gettau(itau);
      double Q2  = fQ2(tau);
      double Q   = std::sqrt(Q2); 
      
      /// alpha_s table has already been filled  
      /// m_alphas[itau] = alphas(rscale_factor*Q)*invtwopi;
      
      // y2 tables
      //    for ( int iy=iymin2 ; iy<=iymax2 ; iy++ ) { 
      for ( int iy=n_y2 ; iy-- ;  ) { 
	
	double y = gety2(iy);
	double x = fx(y);
	double fun = 1;
	if (m_reweight) fun = weightfun(x);
	
	if ( scale_beams ) { 
	  x *= beam_scale;	
	  if ( x>=1 ) { 
	    for ( int ip=0 ; ip<14 ; ip++ ) m_fg2[itau][iy][ip]=0; 
	    if ( nloop==1 && fscale_factor!=1 ) { 
	      for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit2[itau][iy][ip] = 0;
	    }
	    continue; 
	  }
	}
	
	pdf1->evaluate(x, fscale_factor*Q, m_fg2[itau][iy]);
	
	double invx = 1/x;
	for ( int ip=0 ; ip<14 ; ip++ ) m_fg2[itau][iy][ip] *= invx;
	if ( m_reweight ) for ( int ip=0 ; ip<14 ; ip++ ) m_fg2[itau][iy][ip] *= fun;      
	
	// splitting functions
	if ( nloop==1 && fscale_factor!=1 ) { 
	  splitting(x, fscale_factor*Q, m_fsplit2[itau][iy]);
	  for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit2[itau][iy][ip] *= invx;
	  if ( m_reweight ) for ( int ip=0 ; ip<14 ; ip++ ) m_fsplit2[itau][iy][ip] *= fun;
	}
      }

    } // isSymmetric()

  } // loop over itau

}

/// this is the convolute routine for the amcatnlo convolution - essentially it 
/// is the same as for the standard calculation, but the amcatnlo calculation
/// stores weights for the NLO born contribution, and counterterms, so we need
/// more grids than the usual two, 
double appl::igrid::amc_convolute(NodeCache* pdf0,
				  NodeCache* pdf1,
				  appl_pdf*  genpdf,
				  double (*alphas)(const double& ), 
				  int     lo_order,  
				  int     nloop, 
				  double  rscale_factor,
				  double  fscale_factor,
				  double Escale) 
{ 

  //  m_transvar = m_transvarlocal;

  //char name[]="appl_grid:igrid::convolute(): ";
  //  static const double twopi = 2*M_PI;
  static const double eightpisquared = 8*M_PI*M_PI;
  // static const int nc = 3;
  //TC   const int nf = 6;
  // static const int nf = 5;
  //  static double beta0=(11.*nc-2.*nf)/(6.*twopi);
  //const bool debug=false;  

  double alphas_tmp = 0.;  
  double dsigma  = 0.; //, xsigma = 0.;
  double _alphas = 1.;
  //  double  alphaplus1 = 0.;
  // do the convolution  
  // if (debug) std::cout<<name<<" nloop= "<<nloop<<endl;
  //  std::cout << "\torder=" << lo_order << "\tnloop=" << nloop << std::endl;
  // is the grid empty
  int size=0;
   for ( int ip=0 ; ip<m_Nproc ; ip++ ) { 
    if ( !m_weight[ip]->trimmed() )  {
      //  std::cout << "igrid::convolute() naughty, naughty!" << std::endl;
      m_weight[ip]->trim();
    }
    size += m_weight[ip]->xmax() - m_weight[ip]->xmin() + 1;
  }



  // grid is empty
  if ( size==0 )  return 0;

  // 
  //  if ( m_fg1==NULL ) setuppdf(pdf);
  setuppdf( alphas, pdf0, pdf1, nloop, rscale_factor, fscale_factor, Escale);

  double* sig = new double[m_Nproc];  // weights from grid
  double* H   = new double[m_Nproc];  // generalised pdf  
  double* HA  = 0;  // generalised splitting functions
  double* HB  = 0;  // generalised splitting functions
  //  if ( nloop==1 && fscale_factor!=1 ) { 
  //    HA  = new double[m_Nproc];  // generalised splitting functions
  //    HB  = new double[m_Nproc];  // generalised splitting functions
  //  }

  // cross section for this igrid  

  // loop over the grid 
  // 
  for ( int itau=0 ; itau<Ntau() ; itau++  ) {
    _alphas  = 1;    
    alphas_tmp = m_alphas[itau]*eightpisquared;
    for ( int iorder=0 ; iorder<lo_order ; iorder++ ) _alphas *= alphas_tmp;
    //   alphaplus1 = _alphas*alphas_tmp;

    for ( int iy1=Ny1() ; iy1-- ;  ) {            
      for ( int iy2=Ny2() ; iy2-- ;  ) { 
 	// test if this element is actually filled
	// if ( !m_weight[0]->trimmed(itau,iy1,iy2) ) continue; 
	bool nonzero = false;
	// basic convolution order component for either the born level
	// or the convolution of the nlo grid with the pdf 
	for ( int ip=0 ; ip<m_Nproc ; ip++ ) {
	  if ( (sig[ip] = (*(const SparseMatrix3d*)m_weight[ip])(itau,iy1,iy2)) ) nonzero = true;
	}
	
	//	for ( int ip=0 ; ip<m_Nproc ; ip++ ) std::cout << "\t" << sig[ip]; 
	//	std::cout << std::endl;

	if ( nonzero ) { 	

	  // build the generalised pdfs from the actual pdfs
	  genpdf->evaluate( m_fg1[itau][iy1],  m_fg2[itau][iy2], H );
	
	  // do the convolution

          double xsigma=0.;
	  for ( int ip=0 ; ip<m_Nproc ; ip++ ) xsigma+=sig[ip]*H[ip];
	  dsigma += _alphas*xsigma;

#if 0	
	  /// JR scaling by R and F scale factors to grid::convoute() method 
	  // now do the convolution for the variation of factorisation and 
	  // renormalisation scales, proportional to the leading order weights

	  // renormalisation scale dependent bit
	  if ( rscale_factor!=1 ) { 
	    // nlo relative ln mu_R^2 term 
	    dsigma += _alphas*std::log(rscale_factor*rscale_factor)*xsigma;
	    
	  }
	  else if ( fscale_factor!=1 ) {
	    // factorisation scale dependent bit
	    // nlo relative ln mu_F^2 term 
	    dsigma += _alphas*std::log(fscale_factor*fscale_factor)*xsigma;
	    //if (debug) 
	    //cout <<name<<" fscale= " << fscale_factor << " dsigma= "<<dsigma << std::endl;
	  }
	  else { 
	    dsigma += _alphas*xsigma;
	  }
#endif	
	}  // nonzero
      }  // iy2
    }  // iy1
  }  // itau
  
  //if (debug)  std::cout << name<<"     convoluted dsigma=" << dsigma << std::endl; 
  
  delete[] sig;
  delete[] H;
  delete[] HA;
  delete[] HB;
  
  deletepdftable();
  
  //  std::cout << "dsigma " << dsigma << std::endl;

  // NB!!! the return value dsigma must be scaled by Escale*Escale which 
  // is done in grid::vconvolute. It would be better here, but is reduces 
  // the number of operations if done in grid. 
  return dsigma; 
}




bool appl::igrid::shrink( const std::vector<int>& keep ) {
 
  /// save the old grids
  int          Nproc = m_Nproc;
  SparseMatrix3d** w = m_weight;

  /// resize
  m_Nproc  = keep.size();
  m_weight = new SparseMatrix3d*[m_Nproc];

  /// copy across the grids we want to save
  for ( unsigned ip=0 ; ip<keep.size() ; ip++ ) {
 
    //    std::cout << "\tcopy " << keep[ip] << " -> " << ip << std::endl; 

    m_weight[ip] = w[keep[ip]]; /// move across the processes we want to keep
    w[keep[ip]] = 0; /// now flag as 0 
  }

  for ( int ip=0 ; ip<Nproc ; ip++ ) if ( w[ip]!=0 ) delete w[ip];

  delete[] w; 
  
  return true;
}




// FIXME: this has lots of legacy code and stuff commented which
// is no longer needed or even correct, so it should be tidied.
// the point of the algorithm here is to find the extent of the filled
// bines + 1 on either side and create a new grid with these limits
   
void appl::igrid::optimise(int NQ2, int Nx1, int Nx2) {     

  std::cout << "\tsize(untrimmed)=" << m_weight[0]->size();

  //  std::cout << "ymin=" << gety(0) << "\tymax=" << gety(m_Ny-1) 
  //       << "\txmin=" << fx(gety(m_Ny-1)) << "\txmax=" << fx(gety(0)) << std::endl;

  for ( int i=0 ; i<m_Nproc ; i++ )  m_weight[i]->trim();

  std::cout << "\tsize(trimmed)=" << m_weight[0]->size() << std::endl;


#if 1

  // overall igrid optimisation limits

  int _y1setmin = Ny1();
  int _y1setmax = -1;
  
  int _y2setmin = Ny2();
  int _y2setmax = -1;
  
  int _tausetmin = Ntau(); 
  int _tausetmax = -1; 
  
  //  double oldy1min = m_y1min;
  //  double oldy1max = m_y1max;
  
  //  double oldy2min = m_y2min;
  //  double oldy2max = m_y2max;

  // go through all the subprocess to get the limits
  // FIXME: this is actually redundant, at the moment, it assumes all the subprocesses 
  //        occupy the same phase space, so only the first need be tested, and if
  //        they don;t all occupy the same phase space, then they should each be
  //        optimised seperately

  // find limits for this grid
  // initialise to original grid limits
  int y1setmin = _y1setmin;
  int y1setmax = _y1setmax;
  
  int y2setmin = _y2setmin;
  int y2setmax = _y2setmax;
  
  int tausetmin = _tausetmin; 
  int tausetmax = _tausetmax; 
  
  for ( int ip=0 ; ip<m_Nproc ; ip++ ) { 
    
    // is it empty?
    if ( m_weight[ip]->xmax()-m_weight[ip]->xmin()+1 == 0 ) continue;

    //    m_weight[ip]->print();

    //    m_weight[ip]->yaxis().print(fx);

    //    std::cout << "initial ip=" << ip 
    //	       << "\tm_y2min=" << m_y2min << "\tm_y2max=" << m_y2max
    //         << "\tm_y1min=" << m_y1min << "\tm_y1max=" << m_y1max  
    //         << std::endl;  

    //    m_weight[ip]->print();

    // find actual limits

    // y1 optimisation
    int ymin1 = m_weight[ip]->ymin();
    int ymax1 = m_weight[ip]->ymax();
    
    if ( ymin1<y1setmin )                 y1setmin = ymin1;
    if ( ymin1<=ymax1 && y1setmax<ymax1 ) y1setmax = ymax1; 

    //    std::cout << "ip=" << ip << std::endl; // << "\tymin1=" << ymin1 << "\tymax1=" << ymax1 << std::endl;
    
    // y2 optimisation
    int ymin2 = m_weight[ip]->zmin();    
    int ymax2 = m_weight[ip]->zmax();

    if ( ymin2<y2setmin )                 y2setmin = ymin2;
    if ( ymin2<=ymax2 && y2setmax<ymax2 ) y2setmax = ymax2; 
    
    // tau optimisation
    int taumin = m_weight[ip]->xmin();    
    int taumax = m_weight[ip]->xmax();

    if ( taumin<tausetmin )                   tausetmin = taumin;
    if ( taumin<=taumax && taumax>tausetmax ) tausetmax = taumax;

  }
  
  // if grid is empty, do "nothing" ie create the grid with the same
  // limits as before but with the new required number of bins 
  if (  y1setmax==-1 || y2setmax==-1 || tausetmax==-1 ) { 
    m_Ny1  = Nx1;
    m_Ny2  = Nx2;
    m_Ntau = NQ2;
    m_deltay1  = (m_y1max-m_y1min)/(m_Ny1-1);
    m_deltay2  = (m_y2max-m_y2min)/(m_Ny2-1);
    m_deltatau = (m_taumax-m_taumin)/(m_Ntau-1);
  }
  else { 
 

    
    // y1 optimisation
    //    double oldy1min = m_y1min;
    //    double oldy1max = m_y1max;
   

    if ( isOptimised() ) { 
      // add a bit on each side
      if ( y1setmin>0 )        y1setmin--;
      if ( y1setmax<m_Ny1-1 )  y1setmax++;
    }
    else { 
      // not optimised yet so filled with phase space only so add 
      // the order to the max filled grid element position also 
      y1setmax += m_yorder+1;
      if ( y1setmin>0 )           y1setmin--;
      if ( y1setmax>=m_Ny1 )      y1setmax=m_Ny1-1; 
    }
    
    double _min = gety1(y1setmin); 
    double _max = gety1(y1setmax); 
    
    m_Ny1     = Nx1;
    m_y1min   = _min;
    m_y1max   = _max;
    m_deltay1 = (m_y1max-m_y1min)/(m_Ny1-1);
    
    
    
    // y2 optimisation
    //   double oldy2min = m_y2min;
    //   double oldy2max = m_y2max;
    
    if ( isOptimised() ) { 
      // add a bit on each side
      if ( y2setmin>0 )        y2setmin--;
      if ( y2setmax<m_Ny2-1 )  y2setmax++;
    }
    else { 
      // not optimised yet so add the order to the  
      // max filled grid element position also 
      y2setmax+=m_yorder+1;
      if ( y2setmin>0 )           y2setmin--;
      if ( y2setmax>=m_Ny2 )      y2setmax=m_Ny2-1; 
    }
    
    //  m_Ny2 =  y2setmax-y2setmin+1;
    _min = gety2(y2setmin); 
    _max = gety2(y2setmax); 
    
    m_Ny2     = Nx2;
    m_y2min   = _min;
    m_y2max   = _max;
    m_deltay2 = (m_y2max-m_y2min)/(m_Ny2-1);
    
    
    
    // tau optimisation
    //   double oldtaumin = m_taumin;
    //   double oldtaumax = m_taumax;
    
    // add a bit on each side
    if ( isOptimised() ) { 
      // add a bit on each side
      if ( tausetmin>0 )         tausetmin--;
      if ( tausetmax<m_Ntau-1 )  tausetmax++;
    }
    else { 
      // not optimised yet so add the order to the  
      // max filled grid element position also 
      tausetmax+=m_tauorder+1;
      if ( tausetmin>0 )            tausetmin--;
      if ( tausetmax>=m_Ntau )      tausetmax=m_Ntau-1; 
    }  
    
    _min   = gettau(tausetmin); 
    _max   = gettau(tausetmax);
    
    m_Ntau     = NQ2;
    m_taumin   = _min;
    m_taumax   = _max;
    m_deltatau = (m_taumax-m_taumin)/(m_Ntau-1);
    
    //    std::cout << "done ip=" << ip << std::endl; 
    
  }
  
#endif

  // now create the new subprocess grids with optimised limits
  for ( int ip=0 ; ip<m_Nproc ; ip++ ) { 
    delete m_weight[ip];
    m_weight[ip] = new SparseMatrix3d(m_Ntau, m_taumin,   m_taumax,    
				      m_Ny1,   m_y1min,   m_y1max, 
				      m_Ny2,   m_y2min,   m_y2max ); 
  }   
  
  m_optimised = true;
}

int appl::igrid::size() const {
  int _size = 0;
  for ( int i=0 ; i<m_Nproc ; i++ ) _size += m_weight[i]->size();
  return _size;
}

void appl::igrid::trim() {
  for ( int i=0 ; i<m_Nproc ; i++ )  m_weight[i]->trim();
}

void appl::igrid::untrim() {
  for ( int i=0 ; i<m_Nproc ; i++ ) m_weight[i]->untrim();
}

appl::igrid& appl::igrid::operator*=(const double& d) {
  for ( int ip=0 ; ip<m_Nproc ; ip++ ) if ( m_weight[ip] ) (*m_weight[ip]) *= d;
  return *this;
}

// should really check all the limits and *everything* is the same
appl::igrid& appl::igrid::operator+=(const igrid& g) {
  for ( int ip=0 ; ip<m_Nproc ; ip++ ) {
    if ( m_weight[ip] && g.m_weight[ip] ) {
  //if ( (*m_weight[ip]) == (*g.m_weight[ip]) ) (*m_weight[ip]) += (*g.m_weight[ip]);
  if ( m_weight[ip]->compare_axes( *g.m_weight[ip] ) ) (*m_weight[ip]) += (*g.m_weight[ip]);
  else {
    throw exception("igrid::operator+=() grids do not match");
  }
    }
  }
  return *this;
}
