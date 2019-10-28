// emacs: this is -*- c++ -*-
//
//   SparseMatrix3d.h        
//
//                   
// 
//   Copyright (C) 2007 M.Sutton (sutt@hep.ucl.ac.uk)    
//
//   $Id: SparseMatrix3d.h, v1.0   Wed Nov 14 14:23:49 GMT 2007 sutt


#ifndef __SPARSEMATRIX3D_H
#define __SPARSEMATRIX3D_H

#include <iostream>

#include "axis.h"
#include "tsparse3d.h"

#include "TH1D.h"
#include "TH3D.h"


class SparseMatrix3d : public tsparse3d<double> {

public:

  // constructors and destructor

  SparseMatrix3d( int Nx, double lx, double ux, 
		  int Ny, double ly, double uy, 
		  int Nz, double lz, double uz);

  SparseMatrix3d(const SparseMatrix3d& s); 

  ~SparseMatrix3d() { empty_fast(); } 
  

  // utilities for file access and storage

  TH3D* getTH3D(const std::string& s) const; 

  // axis accessors
  const axis<double>& xaxis() const { return m_xaxis; } 
  const axis<double>& yaxis() const { return m_yaxis; } 
  const axis<double>& zaxis() const { return m_zaxis; } 


  // trim to sparse structure 
  void trim() { empty_fast(); tsparse3d<double>::trim(); }
    
  // set up fast lookup table into the (untrimmed) 3d array.
  void setup_fast();

  // and clean up
  void empty_fast() { 
    if ( m_fastindex ) delete[] m_fastindex;
    m_fastindex = NULL;
  }

  // access using the fast (dangerous) methods
  double& fill_fast(int i, int j, int k)       { return *m_fastindex[(i*Ny()+j)*Nz()+k]; }
  double  fill_fast(int i, int j, int k) const { return *m_fastindex[(i*Ny()+j)*Nz()+k]; }

  void fill(double x, double y, double z, double w) { 

    int i=xaxis().bin(x);
    int j=yaxis().bin(y);
    int k=zaxis().bin(z);

    if ( i<0 || i>=Nx() || j<0 || j>=Ny() || k<0 || k>=Nz() ) return; 
   
    if ( m_fastindex ) fill_fast(i,j,k) += w; 
    else                 (*this)(i,j,k) += w;
  }
  
  /// check if the axes are all the same
  bool compare_axes(const SparseMatrix3d& s) const { 

#if 0

    std::cout << "SparseMatrix: " 
	      << "\tx " << m_xaxis 
	      << "\ty " << m_yaxis 
	      << "\tz " << m_zaxis << std::endl; 

    std::cout << "SparseMatrix: " 
	      << "\tx " << s.m_xaxis 
	      << "\ty " << s.m_yaxis 
	      << "\tz " << s.m_zaxis << std::endl; 

    std::cout << "SparseMatrix: " 
	      << "\tx " << ( m_xaxis == s.m_xaxis ) 
	      << "\ty " << ( m_yaxis == s.m_yaxis )
	      << "\tz " << ( m_zaxis == s.m_zaxis ) << std::endl; 


    std::cout << "SparseMatrix: " 
	      << "\tx " << ( m_xaxis.N() - s.m_xaxis.N() ) 
	      << "\tmin "   << m_xaxis.min()   << " " << ( m_xaxis.min() - s.m_xaxis.min() )
	      << "\tmax "   << m_xaxis.max()   << " " << ( m_xaxis.max() - s.m_xaxis.max() )
	      << "\tdelta " << m_xaxis.delta() << " " << ( m_xaxis.delta() - s.m_xaxis.delta() )
	      << "\tz " << std::endl; 

#endif

    return ( m_xaxis == s.m_xaxis &&  
	     m_yaxis == s.m_yaxis &&
	     m_zaxis == s.m_zaxis );
   
  }


  /// checks if the actual contents are the same
  bool operator==(const SparseMatrix3d& s) const;

  bool operator!=(const SparseMatrix3d& s) const { 
    return !( *this == s );
  }


private:

  axis<double> m_xaxis;
  axis<double> m_yaxis;
  axis<double> m_zaxis;
  
  double** m_fastindex;

};

#endif  // __SPARSEMATRIX3D_H 










