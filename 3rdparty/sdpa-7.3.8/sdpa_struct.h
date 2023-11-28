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

// printing presicion of vectors and matrices
#define P_FORMAT ((char*)"%+8.3e")
#define NO_P_FORMAT "NOPRINT"


#ifndef __sdpa_struct_h__
#define __sdpa_struct_h__

#include "sdpa_include.h"
#include "sdpa_block.h"

#define DATA_CAPSULE 1
  // DATA_CAPSULE 0 : Three Arrays (row,column,sp_ele)
  // DATA_CAPSULE 1 : Capsuled data storage
  
namespace sdpa {

class Vector
{
public:
  int nDim;
  double* ele;

  Vector();
  Vector(int nDim, double value = 0.0);
  ~Vector();

  void initialize(int nDim, double value = 0.0);
  void initialize(double value);
  void terminate();

  void setZero();
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  void display(FILE* fpout,double scalar, char* printFormat = P_FORMAT);
  bool copyFrom(Vector& other);
};

class BlockVector
{
public:
  int  nBlock;
  int* blockStruct;

  Vector* ele;
  
  BlockVector();
  BlockVector(BlockStruct& bs, double value = 0.0);
  BlockVector(int nBlock, int* blockStruct, double value = 0.0);
  ~BlockVector();
  
  void initialize(BlockStruct& bs, double value = 0.0);
  void initialize(int nBlock, int* blockStruct, double value = 0.0);
  void initialize(double value);
  void terminate();

  void setZero();
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(BlockVector& other);
};

class SparseMatrix
{
public:
  int nRow, nCol;

  enum Type { SPARSE, DENSE};
  Type type;
  
  int NonZeroNumber;
  // for memory
  int NonZeroCount;
  // currentry stored
  int NonZeroEffect;
  // use for calculation of F1,F2,F3 

  // for Dense
  double* de_ele;

  // for Sparse ; 0:sparse 1:dense
  enum dsType {DSarrays, DScapsule};
  dsType DataStruct;

  // for Sparse Data1 // dsArrays
  int*    row_index;
  int*    column_index;
  double* sp_ele;

  // for Sparse Data2 // dsCapsule
  typedef struct{
    int vRow;
    int vCol;
    double vEle;
  } SparseElement __attribute__( (aligned (16)));

  SparseElement* DataS;

  SparseMatrix();
  SparseMatrix(int nRow,int nCol, Type type, int NonZeroNumber);
  ~SparseMatrix();

  #if DATA_CAPSULE 
  void initialize(int nRow,int nCol, Type type, int NonZeroNumber,
		  dsType DataStruct = DScapsule);
  #else
  void initialize(int nRow,int nCol, Type type, int NonZeroNumber,
		  dsType DataStruct = DSarrays);
  #endif
  void terminate();

  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(SparseMatrix& other);

  void changeToDense(bool forceChange = false);
  void setZero();
  void setIdentity(double scalar = 1.0);

  bool sortSparseIndex(int&i, int& j);
};

class DenseMatrix
{
public:
  int nRow, nCol;

  enum Type { DENSE, COMPLETION};
  Type type;
  
  double* de_ele;

  DenseMatrix();
  DenseMatrix(int nRow,int nCol, Type type);
  ~DenseMatrix();

  void initialize(int nRow,int nCol, Type type);
  void terminate();
  
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(DenseMatrix& other);
  bool copyFrom(SparseMatrix& other);

  void setZero();
  void setIdentity(double scalar = 1.0);
};

class SparseLinearSpace
{
public:
  int  SDP_sp_nBlock;
  int  SOCP_sp_nBlock;
  int  LP_sp_nBlock;

  int*  SDP_sp_index;
  int*  SOCP_sp_index;
  int*  LP_sp_index;

  SparseMatrix* SDP_sp_block;
  SparseMatrix* SOCP_sp_block;
  double* LP_sp_block;
  
  SparseLinearSpace();
  SparseLinearSpace(int SDP_nBlock, int* SDP_blockStruct, 
		    int* SDP_NonZeroNumber,
		    int SOCP_nBlock, int* SOCP_blockStruct,
		    int* SOCP_NonZeroNumber,
		    int LP_nBlock, bool* LP_NonZeroNumber);
  SparseLinearSpace(int SDP_sp_nBlock, 
                    int* SDP_sp_index,
                    int* SDP_sp_blockStruct, 
                    int* SDP_sp_NonZeroNumber,
                    int SOCP_sp_nBlock, 
                    int* SOCP_sp_index,
                    int* SOCP_sp_blockStruct,
                    int* SOCP_sp_NonZeroNumber,
                    int LP_sp_nBlock, 
                    int* LP_sp_index);
  ~SparseLinearSpace();

  // dense form of block index
  void initialize(int SDP_nBlock, int* SDP_blockStruct, 
		    int* SDP_NonZeroNumber,
		    int SOCP_nBlock, int* SOCP_blockStruct,
		    int* SOCP_NonZeroNumber,
		    int LP_nBlock, bool* LP_NonZeroNumber);
  // sparse form of block index      2008/02/27 kazuhide nakata
  void initialize(int SDP_sp_nBlock, 
                  int* SDP_sp_index,
                  int* SDP_sp_blockStruct, 
                  int* SDP_sp_NonZeroNumber,
                  int SOCP_sp_nBlock, 
                  int* SOCP_sp_index,
                  int* SOCP_sp_blockStruct,
                  int* SOCP_sp_NonZeroNumber,
                  int LP_sp_nBlock, 
                  int* LP_sp_index);
  void terminate();
  
  void changeToDense(bool forceChange=false);
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(SparseLinearSpace& other);
  
  void setElement_SDP(int block, int nCol, int nRow, double ele);
  void setElement_SOCP(int block, int nCol, int nRow, double ele);
  void setElement_LP(int block, double ele);

  void setZero();
  void setIdentity(double scalar = 1.0);
  // no check
  bool sortSparseIndex(int&l , int& i, int& j);
};

class DenseLinearSpace
{
 public:
  int  SDP_nBlock;
  int  SOCP_nBlock;
  int  LP_nBlock;

  DenseMatrix* SDP_block;
  DenseMatrix* SOCP_block;
  double* LP_block;

  DenseLinearSpace();
  DenseLinearSpace(BlockStruct& bs);
  ~DenseLinearSpace();
  void initialize(BlockStruct& bs);
  void terminate();

  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  void displaySolution(BlockStruct& bs, FILE* fpout = stdout,
		       char* printFormat = P_FORMAT);
  bool copyFrom(DenseLinearSpace& other);
  void setElement_SDP(int block, int nCol, int nRow, double ele);
  void setElement_SOCP(int block, int nCol, int nRow, double ele);
  void setElement_LP(int block, double ele);
  void setZero();
  void setIdentity(double scalar = 1.0);
};

} // end of namespace 'sdpa'

#endif // __sdpa_struct_h__
