/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004-2012 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/* 
mexAggSDPcones.cpp
[A1,c1] = mexAggSDPcones(A0,c0,K0s,Ks)
*/
#include <string> 
#include <iostream>
#include "mex.h"

using namespace std;

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{    
    /* [A,c,Ks] = mexAggregateSDPcone(A0,c0,K0s,Ks) */
    /* ************************************************** */
    
        /* Check for proper number of input and output arguments */    
    if (nrhs != 4) {
        mexErrMsgTxt("4 input arguments required.");
    } 
    if(nlhs != 2){
        mexErrMsgTxt("2 output arguments required");
    } 

    
    /* ---> */ 
    /* reading the sparse mex matrix A0 */
    const mxArray* A0_ptr = prhs[0];
    if (mxIsSparse(A0_ptr) == 0) {
        mexErrMsgTxt("A0 must be a sparse matrix.");
    }     
    // mwSize A0colSize = mxGetN(A0_ptr);
//    std::cout << "A0colSize = " << A0colSize << std::endl;    
    mwSize A0rowSize = mxGetM(A0_ptr);
//    std::cout << "A0rowSize = " << A0rowSize << std::endl; 
    mwIndex A0nnz = mxGetNzmax(A0_ptr);
//    std::cout << "A0nnz = " << A0nnz << std::endl;	
    mwIndex* A0jc = mxGetJc(A0_ptr);
//    std::cout << "A0jc = ";
//    for (mwSize j=0; j <= A0colSize; j++) 
//        std::cout << A0jc[j] << " "; 
//    std::cout << std::endl;
    mwIndex* A0ir = mxGetIr(A0_ptr);
//    std::cout << "A0ir = ";
//    for (mwSize i=0; i < A0nnz; i++) 
//        std::cout << A0ir[i] << " "; 
//    std::cout << std::endl;            
    double*  A0pr = mxGetPr(A0_ptr);
//    std::cout << "A0pr = ";
//    for (mwSize i=0; i < A0nnz; i++) 
//        std::cout << A0pr[i] << " "; 
//    std::cout << std::endl;
        
    /* reading the sparse mex matrix c0 --- a row vector */
    /* vector c0: Mex */ 
    const mxArray* c0_ptr = prhs[1];
    if (mxIsSparse(c0_ptr) == 0) {
        mexErrMsgTxt("c0 must be sparse.");
    }     

    mwSize c0colSize = mxGetN(c0_ptr); 
//    std::cout << "c0colSize = " << c0colSize << std::endl;
    if (c0colSize <= 1) {
        mexErrMsgTxt("c0 must be a row vector.");
    }     

    mwSize c0rowSize = mxGetM(c0_ptr); // = 1
//    std::cout << "c0rowSize = " << c0rowSize << std::endl;
    if (c0rowSize > 1) {
        mexErrMsgTxt("c0 must be a row vector.");
    }     
    mwIndex c0nnz = mxGetNzmax(c0_ptr);
//    std::cout << "c0nnz = " << c0nnz << std::endl;
	mwIndex* c0jc = mxGetJc(c0_ptr);
//    std::cout << "c0jc = ";
//    for (mwSize j=0; j <= c0colSize; j++) 
//        std::cout << c0jc[j] << " "; 
//    std::cout << std::endl;    
    mwIndex* c0ir = mxGetIr(c0_ptr);
//    std::cout << "c0ir = ";
//    for (mwSize i=0; i < c0nnz; i++) 
//        std::cout << c0ir[i] << " "; 
//    std::cout << std::endl;
    double*  c0pr = mxGetPr(c0_ptr);
//    std::cout << "c0pr = ";
//    for (mwSize i=0; i < c0nnz; i++) 
//        std::cout << c0pr[i] << " "; 
//    std::cout << std::endl;
    
    /* reading the dense mex vector K0s */ 
    const mxArray* K0s_ptr = prhs[2];
    if (mxIsSparse(K0s_ptr) == 1) {
        mexErrMsgTxt("K0s must be a dense column vector.");
    }
    mwSize  K0colSize = mxGetN(K0s_ptr);
    if (K0colSize > 1) {
        mexErrMsgTxt("K0s must be a column vector.");
    }     
    mwSize  K0sSize = mxGetM(K0s_ptr);
    double* K0spr = mxGetPr(K0s_ptr);
//    std::cout << "K0spr = ";
//    for (mwSize k0=0; k0 < K0sSize; k0++) 
//        std::cout << K0spr[k0] << " "; 
//    std::cout << std::endl;
    
    /* reading the dense mex vector K1s */ 
    const mxArray* K1s_ptr = prhs[3];
    if (mxIsSparse(K1s_ptr) == 1) {
        mexErrMsgTxt("K1s must be a dense column vector.");
    }
    mwSize  K1colSize = mxGetN(K1s_ptr);
    if (K1colSize > 1) {
        mexErrMsgTxt("K1s must be a column vector.");
    }     
    mwSize  K1sSize = mxGetM(K1s_ptr);
    double* K1spr = mxGetPr(K1s_ptr);
//    std::cout << "K1spr = ";
//    for (mwSize k1=0; k1 < K1sSize; k1++) 
//        std::cout << K1spr[k1] << " "; 
//    std::cout << std::endl;
    
    /* Constructing the sparse mex matrix A1 */
    
    /*   the column size of A1 be constructed 
       = the column size of c1 to be constructed */ 
    mwSize A1colSize = 0;
    for (mwSize k1 = 0; k1 < K1sSize; k1++)
        A1colSize += (mwSize)(K1spr[k1]*K1spr[k1]);
//    std::cout << "A1colSize = " << A1colSize << std::endl;
    /*   the row size of A1 be constructed */ 
    mwSize A1rowSize = A0rowSize;
    mwIndex A1nnz = A0nnz; 
    
    plhs[0] = mxCreateSparse(A1rowSize,A1colSize,A1nnz,mxREAL);
    mwIndex* A1jc = mxGetJc(plhs[0]);
    mwIndex* A1ir = mxGetIr(plhs[0]);
    double*  A1pr = mxGetPr(plhs[0]);
    
    mwSize c1colSize = A1colSize; 
    mwSize c1rowSize = c0rowSize;
    mwIndex c1nnz = c0nnz; 
    plhs[1] = mxCreateSparse(c1rowSize,c1colSize,c1nnz,mxREAL);
    mwIndex* c1jc = mxGetJc(plhs[1]);
    mwIndex* c1ir = mxGetIr(plhs[1]);
    double*  c1pr = mxGetPr(plhs[1]);
    
    
    int k1 = -1; 
    int topNDblks = 0; 
    int bottomNDblks = 0; 
    mwSize A0colPt = 0;
    mwSize A1colPt = 0; 
    A1jc[A1colPt] = 0; 
    c1jc[A1colPt] = 0; 
    for (mwSize k0 = 0; k0 < K0sSize; k0++) {
        bottomNDblks -= (int) K0spr[k0];
        if (bottomNDblks < 0) {
            k1++; 
            topNDblks = 0;
            bottomNDblks = (int) (K1spr[k1] - K0spr[k0]); 
        } 
        else {
            topNDblks += (int) K0spr[k0-1]; 
        }             
//        std::cout << " k0 = " << k0 << " topNDblks = " << topNDblks << " bottomNDblks =  " << bottomNDblks << std::endl;
        for (mwSize j0 = 0; (int) j0 < K0spr[k0]; j0++) {
            if (topNDblks > 0) {
                for (mwSize k1 = 0; (int) k1 < topNDblks; k1++) {
                    A1jc[A1colPt+1] = A1jc[A1colPt];
                    c1jc[A1colPt+1] = c1jc[A1colPt];
                    A1colPt++;
                }
            }
            for (mwSize i0 = 0; i0 < K0spr[k0]; i0++) {               
                A1jc[A1colPt+1] = A0jc[A0colPt+1];
                c1jc[A1colPt+1] = c0jc[A0colPt+1];
                A1colPt++;
                A0colPt++;
            }
            if (bottomNDblks > 0) {
                for (mwSize k1 = 0; (int) k1 < bottomNDblks; k1++) {
                    A1jc[A1colPt+1] = A1jc[A1colPt];
                    c1jc[A1colPt+1] = c1jc[A1colPt];
                    A1colPt++;
                }
            }
        }
    }
    
//    std::cout << "A1jc =  ";
//    for (mwSize k1 = 0; k1 <= A1colSize; k1++)
//        std::cout << A1jc[k1] << "  "; 
//    std::cout << std::endl;        
    /* A1ir = A0ir and A1pr = A0pr */
    for (mwSize p=0; p < A0nnz; p++) {
        A1ir[p] = A0ir[p];
        A1pr[p] = A0pr[p];
    }
    
    /* c1ir = c0ir and c1pr = c0pr */
    for (mwSize p=0; p < c0nnz; p++) {
        c1ir[p] = c0ir[p];
        c1pr[p] = c0pr[p];
    }

    return;
}
