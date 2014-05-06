/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
// Author: Hadrien Courtecuisse

#ifndef SOFA_COMPONENT_LINEARSOLVER_SPARSEXXTSOLVER_INL
#define SOFA_COMPONENT_LINEARSOLVER_SPARSEXXTSOLVER_INL

#include <sofa/component/linearsolver/SparseXXTSolver.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ObjectFactory.h>
#include <iostream>
#include "sofa/helper/system/thread/CTime.h"
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/behavior/LinearSolver.h>
#include <math.h>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.inl>

#ifdef SOFA_HAVE_METIS
extern "C" {
#include <metis.h>
}
#endif

namespace sofa {

namespace component {

namespace linearsolver {

using namespace sofa::defaulttype;
using namespace sofa::core::behavior;
using namespace sofa::simulation;
using namespace sofa::core::objectmodel;
using sofa::helper::system::thread::CTime;
using sofa::helper::system::thread::ctime_t;
using std::cerr;
using std::endl;


template<class TMatrix, class TVector>
SparseXXTSolver<TMatrix,TVector>::SparseXXTSolver()
: f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
, f_dropTol( initData(&f_dropTol,(double) 0.0,"ic_dropTol","Drop tolerance use for incomplete factorization") )
{
}

template<class TMatrix, class TVector>
SparseXXTSolver<TMatrix,TVector>::~SparseXXTSolver() {}


template<class TMatrix, class TVector>
void SparseXXTSolver<TMatrix,TVector>::solve (Matrix& M, Vector& Z, Vector& R) {
    SparseXXTSolverInvertData * data = (SparseXXTSolverInvertData *) getMatrixInvertData(&M);

    //const int * colptr = &data->colptr[0];
    //const int * rowind = &data->rowind[0];
    const int * perm = &data->perm[0];
    const int * invperm = &data->invperm[0];
    //const Real * values = &data->values[0];
    const Real * D = &data->D[0];
    
    // permutation according to metis
    for (int i=0; i<data->n; i++) B[i] = R[perm[i]];    

    data->A.mul(Z,B);
    
    for (int j = 0 ; j < data->n ; j++) {
	Z[j] /= D[j] ;
    }
    
    data->A.mulTranspose(B,Z);

    for (int i=0; i<data->n; i++) Z[i] = B[invperm[i]];
}

template<class TMatrix, class TVector>
void SparseXXTSolver<TMatrix,TVector>::invert(Matrix& M) {
    SparseXXTSolverInvertData * data = (SparseXXTSolverInvertData *) getMatrixInvertData(&M);
	    
    //remplir A avec M
    data->n = M.colSize();// number of columns	  
    data->Mfiltered.clear();
    data->Mfiltered.resize(M.rowSize(),M.colSize());
    data->Mfiltered.copyNonZeros(M);
    data->Mfiltered.compress();	  

    int * M_colptr = (int *) &data->Mfiltered.getRowBegin()[0];
    int * M_rowind = (int *) &data->Mfiltered.getColsIndex()[0];
    Real * M_values = (Real *) &data->Mfiltered.getColsValue()[0];

    int * perm;
    int * invperm;
    int * colptr;
    int * rowind;
    Real * values;
    Real * D;
    
    data->perm.resize(data->n);     perm = &data->perm[0];
    data->invperm.resize(data->n);  invperm = &data->invperm[0];
    data->colptr.resize(data->n+1); colptr = &data->colptr[0];
    data->D.resize(data->n+1);      D = &data->D[0];
    
//     sofa::helper::AdvancedTimer::stepBegin("LDL_ordering");    
      LDL_ordering(data->n,M_colptr,M_rowind,perm,invperm);    
//     sofa::helper::AdvancedTimer::stepEnd("LDL_ordering");
//     sofa::helper::AdvancedTimer::stepBegin("LDL_symbolic");   
      LDL_symbolic(data->n,M_colptr,M_rowind,colptr,perm,invperm);    
//     sofa::helper::AdvancedTimer::stepEnd("LDL_symbolic");  
    data->nnz = colptr[data->n];
    
    data->rowind.resize(data->nnz); rowind = &data->rowind[0];
    data->values.resize(data->nnz); values = &data->values[0];	
    
//     sofa::helper::AdvancedTimer::stepBegin("LDL_numeric");
      LDL_numeric(data->n,M_colptr,M_rowind,M_values,colptr,rowind,values,D,perm,invperm);
//     sofa::helper::AdvancedTimer::stepEnd("LDL_numeric");
    
    //// Now compute the XXT factorisation
    
    data->A.clear();
    data->A.resize(data->n,data->n);
    
//     sofa::helper::AdvancedTimer::stepBegin("LDL_XXT");
    data->A.set(0,0,1.0);
    for (int j=1;j<data->n;j++) {
      const typename SparseMatrix<Real>::Line & lA = data->A[j-1];
      for (typename SparseMatrix<Real>::LElementConstIterator i1 = lA.begin(), i1end = lA.end(); i1 != i1end; ++i1) {
	  int a = i1->first;
	  Real T_ij = i1->second;
	  for (int i=colptr[j-1];i<colptr[j];i++) {
	    Real A_ij = - values[i];
	    data->A.add(rowind[i],a,A_ij * T_ij);
	  }
      }
      data->A.set(j,j,1.0);
    }
//     sofa::helper::AdvancedTimer::stepEnd("LDL_XXT");  
    
    
//     int nnz = 0;
//     const typename SparseMatrix<Real>::LineConstIterator jitend = data->A.end();    
//     for (typename SparseMatrix<Real>::LineConstIterator jit1 = data->A.begin(); jit1 != jitend; ++jit1) {
// 	for (typename SparseMatrix<Real>::LElementConstIterator i1 = jit1->second.begin(), i1end = jit1->second.end(); i1 != i1end; ++i1) {
// 	      nnz++;
// 	}
//     }
// 
//     double pldl = (data->nnz * 100.0) / (data->n * data->n * 0.5);
//     double pxxt = (nnz * 100.0) / (data->n * data->n * 0.5);
// 
//     printf("ldl %f xxt %f\n",pldl,pxxt);
//     
    B.resize(data->n);    
}


template<class TMatrix, class TVector>
void SparseXXTSolver<TMatrix,TVector>::LDL_ordering(int n,int * M_colptr,int * M_rowind,int * perm,int * invperm) {
    int  num_flag     = 0;
    int  options_flag = 0;
    
    xadj.resize(n+1);
    adj.resize(M_colptr[n]-n);
    
    int it = 0;
    for (int j=0; j<n; j++) {
      xadj[j] = M_colptr[j] - j;
      
      for (int ip = M_colptr[j]; ip < M_colptr[j+1]; ip++) {
	int i = M_rowind[ip];
	if (i != j) adj[it++] = i;
      }
    }
    xadj[n] = M_colptr[n] - n;

    METIS_NodeND(&n, &xadj[0],&adj[0], &num_flag, &options_flag, perm,invperm);
}


template<class TMatrix, class TVector>
void SparseXXTSolver<TMatrix,TVector>::LDL_symbolic (int n,int * M_colptr,int * M_rowind,int * colptr,int * perm,int * invperm) {
    Parent.clear();
    Lnz.clear();
    Flag.clear();
    Pattern.clear();
    Parent.resize(n);
    Lnz.resize(n);
    Flag.resize(n);
    Pattern.resize(n);
    
    for (int k = 0 ; k < n ; k++) {
	Parent [k] = -1 ;	    /* parent of k is not yet known */
	Flag [k] = k ;		    /* mark node k as visited */
	Lnz [k] = 0 ;		    /* count of nonzeros in column k of L */
	int kk = perm[k];  /* kth original, or permuted, column */
	for (int p = M_colptr[kk] ; p < M_colptr[kk+1] ; p++) {
	    /* A (i,k) is nonzero (original or permuted A) */
	    int i = invperm[M_rowind[p]];
	    if (i < k) {
		/* follow path from i to root of etree, stop at flagged node */
		for ( ; Flag [i] != k ; i = Parent [i]) {
		    /* find parent of i if not yet determined */
		    if (Parent [i] == -1) Parent [i] = k ;
		    Lnz [i]++ ;				/* L (k,i) is nonzero */
		    Flag [i] = k ;			/* mark i as visited */
		}
	    }
	}
    }

    colptr[0] = 0 ;
    for (int k = 0 ; k < n ; k++) colptr[k+1] = colptr[k] + Lnz[k] ;
}

template<class TMatrix, class TVector>
void SparseXXTSolver<TMatrix,TVector>::LDL_numeric(int n,int * M_colptr,int * M_rowind,Real * M_values,int * colptr,int * rowind,Real * values,Real * D,int * perm,int * invperm) {
    Real yi, l_ki ;
    int i, p, kk, len, top ;
    
    Y.resize(n);
    
    for (int k = 0 ; k < n ; k++) {
	Y [k] = 0.0 ;		    /* Y(0:k) is now all zero */
	top = n ;		    /* stack for pattern is empty */
	Flag [k] = k ;		    /* mark node k as visited */
	Lnz [k] = 0 ;		    /* count of nonzeros in column k of L */
	kk = perm[k];  /* kth original, or permuted, column */
	for (p = M_colptr[kk] ; p < M_colptr[kk+1] ; p++) {
	    i = invperm[M_rowind[p]];	/* get A(i,k) */
	    if (i <= k) {
		Y[i] += M_values[p] ;  /* scatter A(i,k) into Y (sum duplicates) */
		for (len = 0 ; Flag[i] != k ; i = Parent[i]) {
		    Pattern [len++] = i ;   /* L(k,i) is nonzero */
		    Flag [i] = k ;	    /* mark i as visited */
		}
		while (len > 0) Pattern[--top] = Pattern [--len] ;
	    }
	}
	/* compute numerical values kth row of L (a sparse triangular solve) */
	D[k] = Y [k] ;		    /* get D(k,k) and clear Y(k) */
	Y[k] = 0.0 ;
	for ( ; top < n ; top++) {
	    i = Pattern [top] ;	    /* Pattern [top:n-1] is pattern of L(:,k) */
	    yi = Y [i] ;	    /* get and clear Y(i) */
	    Y [i] = 0.0 ;
	    for (p = colptr[i] ; p < colptr[i] + Lnz [i] ; p++) {
		Y[rowind[p]] -= values[p] * yi ;
	    }
	    l_ki = yi / D[i] ;	    /* the nonzero entry L(k,i) */
	    D[k] -= l_ki * yi ;
	    rowind[p] = k ;	    /* store L(k,i) in column form of L */
	    values[p] = l_ki ;
	    Lnz[i]++ ;		    /* increment count of nonzeros in col i */
	}
	if (D[k] == 0.0) {
	  std::cerr << "SparseXXTSolver failure to factorize, D(k,k) is zero" << std::endl;
	  return;
	}
    }
}

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif
