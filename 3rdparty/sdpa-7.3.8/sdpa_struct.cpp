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

#include "sdpa_struct.h"
#include "sdpa_algebra.h"

namespace sdpa{

Vector::Vector()
{
  nDim = 0;
  ele  = NULL;
}

Vector::Vector(int nDim, double value)
{
  ele  = NULL;
  initialize(nDim,value);
}

Vector::~Vector()
{
  terminate();
}

void Vector::initialize(int nDim,double value)
{
  // rMessage("Vector initialize");
  if (nDim<=0) {
    rError("Vector:: nDim is nonpositive");
  }
  if (this->nDim!=nDim) {
    DeleteArray(ele);
  }
  this->nDim = nDim;
  if (ele == NULL) {
    NewArray(ele,double,nDim);
  }
  sdpa_dset(nDim,value,ele,IONE);
}

void Vector::initialize(double value)
{
  if (ele==NULL) {
    NewArray(ele,double,nDim);
  }
  sdpa_dset(nDim,value,ele,IONE);
}

void Vector::terminate()
{
  DeleteArray(ele);
}

void Vector::setZero()
{
  initialize(0.0);
}

void Vector::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    fprintf(fpout,printFormat,ele[j]);
    fprintf(fpout, ",");
  }
  if (nDim>0) {
    fprintf(fpout,printFormat,ele[nDim-1]);
    fprintf(fpout,"}\n");
  } else {
    fprintf(fpout,"  }\n");
  }
}

void Vector::display(FILE* fpout,double scalar, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    fprintf(fpout,printFormat,ele[j]*scalar);
    fprintf(fpout,",");
  }
  if (nDim>0) {
    fprintf(fpout,printFormat,ele[nDim-1]*scalar);
    fprintf(fpout,"}\n");
  } else {
    fprintf(fpout,"  }\n");
  }
}

bool Vector::copyFrom(Vector& other)
{
  if (this == &other) {
    return SDPA_SUCCESS;
  }
  if (other.nDim<=0) {
    rError("Vector:: nDim is nonpositive");
  }
  if (nDim != other.nDim) {
    DeleteArray(ele);
  }
  nDim = other.nDim;
  if (ele==NULL) {
    NewArray(ele,double,nDim);
  }
  dcopy_fc(&nDim,other.ele,&IONE,ele,&IONE);
  return SDPA_SUCCESS;
}

BlockVector::BlockVector()
{
  nBlock      = 0;
  blockStruct = NULL;
  ele         = NULL;
}

BlockVector::BlockVector(BlockStruct& bs, double value)
{
  initialize(bs.SDP_nBlock,bs.SDP_blockStruct,value);
}

BlockVector::BlockVector(int nBlock, int* blockStruct,
			   double value)
{
  initialize(nBlock,blockStruct,value);
}

BlockVector::~BlockVector()
{
  terminate();
}

void BlockVector::initialize(BlockStruct& bs, double value)
{
  initialize(bs.SDP_nBlock,bs.SDP_blockStruct,value);
}

void BlockVector::initialize(int nBlock, int* blockStruct,
			      double value)
{
  // rMessage("BlockVector initialize");
  if (nBlock<=0) {
    rError("BlockVector:: nBlock is nonpositive");
  }
  this->nBlock = nBlock;
  NewArray(this->blockStruct,int,nBlock);
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
  }

  NewArray(ele,Vector,nBlock);
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size<0) {
      size = -size;
    }
    ele[l].initialize(size,value);
  }
}

void BlockVector::initialize(double value)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].initialize(value);
    }
  }
}

void BlockVector::terminate()
{
  if (ele && blockStruct && nBlock>=0) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].terminate();
    }
    DeleteArray(ele);
    DeleteArray(blockStruct);
  }
}

void BlockVector::setZero()
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setZero();
    }
  }
}

void BlockVector::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout,"{ ");
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].display(fpout,printFormat);
    }
  }
  fprintf(fpout,"} \n");
}

bool BlockVector::copyFrom(BlockVector& other)
{
  if (this == &other) {
    return SDPA_SUCCESS;
  }
  
  if (other.nBlock<=0) {
    rError("BlockVector:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    DeleteArray(blockStruct);
    DeleteArray(ele);
  }
  if (blockStruct==NULL) {
    nBlock = other.nBlock;
    NewArray(blockStruct,int,nBlock);
    for (int l=0; l<nBlock; ++l) {
      blockStruct[l] = other.blockStruct[l];
    }
  }
  if (ele==NULL) {
    NewArray(ele,Vector,nBlock);
  }
  for (int l=0; l<nBlock; ++l) {
    ele[l].copyFrom(other.ele[l]);
  }
  return SDPA_SUCCESS;
}

SparseMatrix::SparseMatrix()
{
  nRow = 0;
  nCol = 0;
  type = SPARSE;

  NonZeroNumber = 0;
  
  de_ele = NULL;

  row_index     = NULL;
  column_index  = NULL;
  sp_ele        = NULL;

  DataStruct = DSarrays;
  DataS = NULL;

  NonZeroCount  = 0;
  NonZeroEffect = 0;
}

SparseMatrix::SparseMatrix(int nRow, int nCol,
			     SparseMatrix::Type type,
			     int NonZeroNumber)
{
  initialize(nRow, nCol, type, NonZeroNumber);
}

SparseMatrix::~SparseMatrix()
{
  terminate();
}

void SparseMatrix::initialize(int nRow, int nCol,
			      SparseMatrix::Type type,
			      int NonZeroNumber,
			      SparseMatrix::dsType DataStruct)
{
  // rMessage("SparseMatrix initialize");

  SparseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("SparseMatrix:: Dimensions are nonpositive");
  }
  this->nRow          = nRow;
  this->nCol          = nCol;
  this->type          = type;
  this->DataStruct    = DataStruct;

  int length;
  switch(type) {
  case SPARSE:
      this->NonZeroNumber  = NonZeroNumber;
      this->NonZeroCount   = 0;
      this->NonZeroEffect  = 0;
      if (NonZeroNumber > 0) {
	if (DataStruct == DSarrays) {
	  NewArray(row_index,int,NonZeroNumber);
	  NewArray(column_index,int,NonZeroNumber);
	  NewArray(sp_ele,double,NonZeroNumber);
	  if (row_index==NULL || column_index==NULL
	      || sp_ele==NULL) {
	    rError("SparseMatrix:: memory exhausted");
	  }
	}
	else {
	  NewArray(DataS, SparseElement, NonZeroNumber);
	  if (DataS == NULL) {
	    rError("SparseElement:: memory exhausted");
	  }
	}
      }
    break;
  case DENSE:
    this->NonZeroNumber = nRow*nCol;
    this->NonZeroCount  = nRow*nCol;
    this->NonZeroEffect = nRow*nCol;
    NewArray(de_ele,double,NonZeroNumber);
    if (de_ele==NULL) {
      rError("SparseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    // all elements are 0.
    break;
  }
}

void SparseMatrix::terminate()
{
  DeleteArray(de_ele);
  if (DataStruct == DSarrays) {
    DeleteArray(row_index);
    DeleteArray(column_index);
    DeleteArray(sp_ele);
  }
  else {
    DeleteArray(DataS);
  }
}

void SparseMatrix::display(FILE* fpout, char* printFormat)
{
  int i, j;
  double value;
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  switch(type) {
  case SPARSE:
    fprintf(fpout,"{");
    for (int index=0; index<NonZeroCount; ++index) {
      if (DataStruct == DSarrays) {
	i        = row_index[index];
	j        = column_index[index];
	value    = sp_ele[index];
      }
      else {
	i        = DataS[index].vRow;
	j        = DataS[index].vCol;
	value    = DataS[index].vEle;
      }
      fprintf(fpout,"val[%d,%d] = ", i,j);
      fprintf(fpout,printFormat,value);
      fprintf(fpout,"\n");
    }
    fprintf(fpout,"}\n");
    break;
  case DENSE:
    fprintf(fpout,"{\n");
    for (int i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
	fprintf(fpout,printFormat,de_ele[i+nCol*j]);
	fprintf(fpout, ",");
      }
      fprintf(fpout,printFormat,de_ele[i+nCol*(nCol-1)]);
      fprintf(fpout, " },\n");
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (int j=0; j<nCol-1; ++j) {
      fprintf(fpout,printFormat,de_ele[(nRow-1)+nCol*j]);
      fprintf(fpout, ",");
    }
    fprintf(fpout,printFormat,de_ele[(nRow-1)+nCol*(nCol-1)]);
    fprintf(fpout, " }");
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  }
}

bool SparseMatrix::copyFrom(SparseMatrix& other)
{
  if (type != other.type || nRow != other.nRow
      || nCol != other.nCol) {
    this->~SparseMatrix();
    initialize(other.nRow,other.nCol,other.type,
	       NonZeroNumber);
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length;
    switch(type) {
    case SPARSE:
      for (int index = 0; index<NonZeroCount;++index) {
	if (DataStruct == DSarrays) {
	  row_index[index]    = other.row_index[index];
	  column_index[index] = other.column_index[index];
	  sp_ele[index]       = other.sp_ele[index];
	}
	else {
	  DataS[index].vRow    = other.DataS[index].vRow;
	  DataS[index].vCol    = other.DataS[index].vCol;
	  DataS[index].vEle    = other.DataS[index].vEle;
	}
      }
      break;
    case DENSE:
      length = nRow*nCol;
      dcopy_fc(&length,other.de_ele,&IONE,de_ele,&IONE);
      break;
    }
  } else { // Sp_De_Di == other.Sp_De_Di
           // && nRow == other.nRow && nCol == other.nCol
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length;
    switch(type) {
    case SPARSE:
      if (NonZeroNumber!=other.NonZeroNumber) {
	if (DataStruct == DSarrays) {
	  DeleteArray(row_index);
	  DeleteArray(column_index);
	  DeleteArray(sp_ele);
	  NewArray(row_index   ,int   ,NonZeroNumber);
	  NewArray(column_index,int   ,NonZeroNumber);
	  NewArray(sp_ele      ,double,NonZeroNumber);
	}
	else {
	  NewArray(DataS, SparseElement,NonZeroNumber);
	}
      }
      for (int index = 0; index<NonZeroCount;++index) {
	if (DataStruct == DSarrays) {
	  row_index[index]    = other.row_index[index];
	  column_index[index] = other.column_index[index];
	  sp_ele[index]       = other.sp_ele[index];
	}
	else {
	  DataS[index].vRow    = other.DataS[index].vRow;
	  DataS[index].vCol    = other.DataS[index].vCol;
	  DataS[index].vEle    = other.DataS[index].vEle;
	}
      }
      break;
    case DENSE:
      length = nRow*nCol;
      dcopy_fc(&length,other.de_ele,&IONE,de_ele,&IONE);
      break;
    } // end of switch
  } // end of else
  return SDPA_SUCCESS;
}

void SparseMatrix::changeToDense(bool forceChange)
{
  if (type!=SPARSE) {
    return;
  }
  // if (false)
  // rMessage(" NonZeroCount " << NonZeroCount);
  // rMessage(" nRow*nCol*0.2 " << nRow*nCol*0.2);
  if (forceChange == false && NonZeroCount < (nRow*nCol) * 0.20) {
    // if the number of elements are less than 20 percent,
    // we don't change to Dense.
    return;
  }
  // rMessage("change");
  type = DENSE;
  de_ele = NULL;
  int length = nRow*nCol;
  NewArray(de_ele,double,length);
  sdpa_dset(length,DZERO,de_ele,IONE);
  // all elements are set 0.
  for (int index=0; index<NonZeroCount; ++index) {
#if DATA_CAPSULE
    int        i = DataS[index].vRow;
    int        j = DataS[index].vCol;
    double value = DataS[index].vEle;
#else
    int        i = row_index[index];
    int        j = column_index[index];
    double value = sp_ele[index];
#endif
    if (i==j) {
      de_ele[i+nCol*j] = value;
    } else {
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
  }
  NonZeroCount = NonZeroNumber = NonZeroEffect = length;
  if (DataStruct == DSarrays) {
    DeleteArray(row_index);
    DeleteArray(column_index);
    DeleteArray(sp_ele);
  }
  else {
    DeleteArray(DataS);
  }
}

void SparseMatrix::setZero()
{
  int length;
  switch(type) {
  case SPARSE:
    NonZeroCount  = 0;
    NonZeroEffect = 0;
    // No element is stored.
    break;
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    break;
  }
}

void SparseMatrix::setIdentity(double scalar)
{
  if (nRow != nCol) {
    rError("SparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step;
  switch(type) {
  case SPARSE:
    if (nCol > NonZeroNumber) {
      rError("SparseMatrix:: cannot store over NonZeroNumber");
      // the number of Diagonal elements equals nCol.
    }
    NonZeroCount  = nCol;
    NonZeroEffect = nCol;
    for (int index=0; index< NonZeroCount; ++index) {
      #if DATA_CAPSULE
      DataS[index].vRow   = index;
      DataS[index].vCol   = index;
      DataS[index].vEle   = scalar;
      #else
      row_index[index]    = index;
      column_index[index] = index;
      sp_ele[index]       = scalar;
      #endif
    }
    break;
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    step = nCol+1;
    sdpa_dset(nCol,scalar,de_ele,step);
    // only diagonal elements are set the value of scalar.
    break;
  }
}
    
bool SparseMatrix::sortSparseIndex(int& i, int& j)
{
  // if this matrix is not symmetric,
  // return the index(i,j) whose values are not symmetric.
  i = -1;
  j = -1;
  const double tolerance = 1.0e-8;
  switch(type) {
  case SPARSE:
    // Make matrix as Upper Triangluar
    for (int i1=0; i1<NonZeroCount; ++i1) {
#if DATA_CAPSULE
      int tmpi = DataS[i1].vRow;
      int tmpj = DataS[i1].vCol;
      if (tmpi>tmpj) {
	DataS[i1].vRow = tmpj;
	DataS[i1].vCol = tmpi;
      }
#else
      int tmpi = row_index[i1];
      int tmpj = column_index[i1];
      if (tmpi>tmpj) {
	row_index   [i1] = tmpj;
	column_index[i1] = tmpi;
      }
#endif
    }
    // simple sort
    for (int i1=0; i1<NonZeroCount; ++i1) {
      for (int i2=0; i2<i1; ++i2) {
	#if DATA_CAPSULE
	int index1 = DataS[i1].vRow + DataS[i1].vCol;
	int index2 = DataS[i2].vRow + DataS[i2].vCol;
	if (index1<index2) {
	  int         tmpi = DataS[i2].vRow;
	  int         tmpj = DataS[i2].vCol;
	  double      tmpv = DataS[i2].vEle;
	  DataS[i2].vRow   = DataS[i1].vRow;
	  DataS[i2].vCol   = DataS[i1].vCol;
	  DataS[i2].vEle   = DataS[i1].vEle;
	  DataS[i1].vRow   = tmpi;
	  DataS[i1].vCol   = tmpj;
	  DataS[i1].vEle   = tmpv;
	}
	#else
	int index1 = row_index[i1]+nCol*column_index[i1];
	int index2 = row_index[i2]+nCol*column_index[i2];
	if (index1<index2) {
	  int         tmpi = row_index   [i2];
	  int         tmpj = column_index[i2];
	  double      tmpv = sp_ele      [i2];
	  row_index   [i2] = row_index   [i1];
	  column_index[i2] = column_index[i1];
	  sp_ele      [i2] = sp_ele      [i1];
	  row_index   [i1] = tmpi;
	  column_index[i1] = tmpj;
	  sp_ele      [i1] = tmpv;
	}
	#endif
      }
    }
    // the process for the same index
    for (int i1=0; i1<NonZeroCount-1; ++i1) {
#if DATA_CAPSULE
      int index1 = DataS[i1].vRow + DataS[i1].vCol;
      int index2 = DataS[i1+1].vRow + DataS[i1+1].vCol;
      if (index1 == index2) {
	if (fabs(DataS[index1].vEle - DataS[index2].vEle) > tolerance) {
	  // Here must not be symmetric
	  if (i<0 || j<0) {
	    i = DataS[i1].vRow;
	    j = DataS[i1].vCol;
	  }
	}
	// remove redudunt
	for (int i2 = i1+1; i2<NonZeroCount-2;++i2) {
	  DataS[i2].vRow = DataS[i2+1].vRow;
	  DataS[i2].vCol = DataS[i2+1].vCol;
	  DataS[i2].vEle = DataS[i2+1].vEle;
	}
	NonZeroCount--;
	if (i==j) {
	  NonZeroEffect--;
	} else {
	  NonZeroEffect -= 2;
	}
      } // end of 'if (index1==index2)'
#else
      int index1 = row_index[i1  ]+nCol*column_index[i1  ];
      int index2 = row_index[i1+1]+nCol*column_index[i1+1];
      if (index1 == index2) {
	if (fabs(sp_ele[index1] - sp_ele[index2]) > tolerance) {
	  // Here must not be symmetric
	  if (i<0 || j<0) {
	    i = row_index   [i1];
	    j = column_index[i1];
	  }
	}
	// remove redudunt
	for (int i2 = i1+1; i2<NonZeroCount-2;++i2) {
	  row_index   [i2] = row_index   [i2+1];
	  column_index[i2] = column_index[i2+1];
	  sp_ele      [i2] = sp_ele      [i2+1];
	}
	NonZeroCount--;
	if (i==j) {
	  NonZeroEffect--;
	} else {
	  NonZeroEffect -= 2;
	}
      } // end of 'if (index1==index2)'
#endif
    }
    break;
  case DENSE:
    if (nRow!=nCol) {
      return SDPA_FAILURE;
    }
    for (j=1; j<nCol; ++j) {
      for (i=0; i<j; ++i) {
	if (fabs(de_ele[i+nCol*j]-de_ele[j+nCol*i]) > tolerance) {
	  return SDPA_FAILURE;
	}
      }
    }
    break;
  }
  return SDPA_SUCCESS;
}

DenseMatrix::DenseMatrix()
{
  nRow = 0;
  nCol = 0;
  type = DENSE;

  de_ele = NULL;
}

DenseMatrix::DenseMatrix(int nRow, int nCol,
			   DenseMatrix::Type type)
{
  initialize(nRow, nCol, type);
}

DenseMatrix::~DenseMatrix()
{
  terminate();
}

void DenseMatrix::initialize(int nRow, int nCol,
			     DenseMatrix::Type type)
{
  // rMessage("DenseMatrix::initialize");

  DenseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("DenseMatrix:: Dimensions are nonpositive");
  }
  int old_length = this->nRow*this->nCol;
  this->nRow  = nRow;
  this->nCol  = nCol;

  int length;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    if (de_ele && old_length!=length) {
      DeleteArray(de_ele);
    }
    if (de_ele==NULL) {
      NewArray(de_ele,double,length);
    }
    sdpa_dset(length,DZERO,de_ele,IONE);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

void DenseMatrix::terminate()
{
  DeleteArray(de_ele);
}

void DenseMatrix::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  switch(type) {
  case DENSE:
    fprintf(fpout,"{");
    for (int i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
	fprintf(fpout,printFormat,de_ele[i+nCol*j]);
	fprintf(fpout, ",");
      }
      fprintf(fpout,printFormat,de_ele[i+nCol*(nCol-1)]);
      fprintf(fpout, " },\n");
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (int j=0; j<nCol-1; ++j) {
      fprintf(fpout,printFormat,de_ele[(nRow-1)+nCol*j]);
      fprintf(fpout, ",");
    }
    fprintf(fpout,printFormat,de_ele[(nRow-1)+nCol*(nCol-1)]);
    fprintf(fpout, " }");
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

bool DenseMatrix::copyFrom(SparseMatrix& other)
{
  int length;
  switch(other.type) {
  case SparseMatrix::SPARSE:
    type = DENSE;
    DeleteArray(de_ele);
    nRow = other.nRow;
    nCol = other.nCol;
    NewArray(de_ele,double,nRow*nCol);
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    for (int index = 0; index<other.NonZeroCount; ++index) {
      #if DATA_CAPSULE
      int i = other.DataS[index].vRow;
      int j = other.DataS[index].vCol;
      double value = other.DataS[index].vEle;
      #else
      int i = other.row_index[index];
      int j = other.column_index[index];
      double value = other.sp_ele[index];
      #endif
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
    break;
  case SparseMatrix::DENSE:
    type = DENSE;
    if (other.nRow!=nRow || other.nCol!=nCol) {
      DeleteArray(de_ele);
    }
    nRow = other.nRow;
    nCol = other.nCol;
    NewArray(de_ele,double,nRow*nCol);
    length = nRow*nCol;
    dcopy_fc(&length,other.de_ele,&IONE,de_ele,&IONE);
    break;
  }
  return SDPA_SUCCESS;
}

bool DenseMatrix::copyFrom(DenseMatrix& other)
{
  if (this == &other) {
    return SDPA_SUCCESS;
  }
  int length;
  switch(other.type) {
  case DENSE:
    type = DENSE;
    if (other.nRow!=nRow || other.nCol!=nCol) {
      DeleteArray(de_ele);
    }
    nRow = other.nRow;
    nCol = other.nCol;
    if (de_ele==NULL) {
      NewArray(de_ele,double,nRow*nCol);
    }
    length = nRow*nCol;
    dcopy_fc(&length,other.de_ele,&IONE,de_ele,&IONE);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
  return SDPA_SUCCESS;
}


void DenseMatrix::setZero()
{
  int length;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}

void DenseMatrix::setIdentity(double scalar)
{
  if (nRow != nCol) {
    rError("SparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step;
  switch(type) {
  case DENSE:
    length = nRow*nCol;
    sdpa_dset(length,DZERO,de_ele,IONE);
    step = nCol+1;
    sdpa_dset(nCol,scalar,de_ele,step);
    break;
  case COMPLETION:
    rError("DenseMatrix:: no support for COMPLETION");
    break;
  }
}
    
SparseLinearSpace::SparseLinearSpace()
{
  SDP_sp_nBlock  = 0;
  SDP_sp_index   = NULL;
  SDP_sp_block   = NULL;
  SOCP_sp_nBlock = 0;
  SOCP_sp_index  = NULL;
  SOCP_sp_block  = NULL;
  LP_sp_nBlock   = 0;
  LP_sp_index    = NULL;
  LP_sp_block    = NULL;
}

SparseLinearSpace::SparseLinearSpace(int SDP_nBlock, 
				     int* SDP_blockStruct,
				     int* SDP_NonZeroNumber,
				     int SOCP_nBlock, 
				     int* SOCP_blockStruct,
				     int* SOCP_NonZeroNumber,
				     int LP_nBlock, 
				     bool* LP_NonZeroNumber)
{
  initialize(SDP_nBlock, SDP_blockStruct, SDP_NonZeroNumber,
	     SOCP_nBlock, SOCP_blockStruct, SOCP_NonZeroNumber,
	     LP_nBlock, LP_NonZeroNumber);
}

SparseLinearSpace::SparseLinearSpace(int SDP_sp_nBlock, 
                                     int* SDP_sp_index,
                                     int* SDP_sp_blockStruct, 
                                     int* SDP_sp_NonZeroNumber,
                                     int SOCP_sp_nBlock, 
                                     int* SOCP_sp_index,
                                     int* SOCP_sp_blockStruct,
                                     int* SOCP_sp_NonZeroNumber,
                                     int LP_sp_nBlock, 
                                     int* LP_sp_index)
{
  initialize(SDP_sp_nBlock, SDP_sp_index,
             SDP_sp_blockStruct, SDP_sp_NonZeroNumber,
	     SOCP_sp_nBlock, SOCP_sp_index,
             SOCP_sp_blockStruct, SOCP_sp_NonZeroNumber,
	     LP_sp_nBlock, LP_sp_index);
}

SparseLinearSpace::~SparseLinearSpace()
{
  terminate();
}

// dense form of block index
void SparseLinearSpace::initialize(int SDP_nBlock, 
				   int* SDP_blockStruct,
				   int* SDP_NonZeroNumber,
				   int SOCP_nBlock, 
				   int* SOCP_blockStruct,
				   int* SOCP_NonZeroNumber,
				   int LP_nBlock, 
				   bool* LP_NonZeroNumber)
{
  // rMessage("SparseLinearSpace::initialize");
  SDP_sp_nBlock = 0;
  SOCP_sp_nBlock = 0;
  LP_sp_nBlock = 0;
  int counter;

  // for SDP
  for (int l=0; l<SDP_nBlock; l++){
    if (SDP_NonZeroNumber[l] > 0){
      SDP_sp_nBlock++;
    }    
  }
  if (SDP_sp_nBlock > 0){
    NewArray(SDP_sp_index,int,SDP_sp_nBlock);
    NewArray(SDP_sp_block,SparseMatrix,SDP_sp_nBlock);
  }
  counter = 0;
  for (int l=0; l<SDP_nBlock; ++l) {
    if (SDP_NonZeroNumber[l] > 0){
      SDP_sp_index[counter] = l;
      int size = SDP_blockStruct[l];
      SDP_sp_block[counter].initialize(size,size,SparseMatrix::SPARSE,
				       SDP_NonZeroNumber[l]);
      counter++;
    }    
  }


  // for SOCP
#if 0
  for (int l=0; l<SOCP_nBlock; l++){
    if (SOCP_NonZeroNumber[l] > 0){
      SOCP_sp_nBlock++;
    }    
  }
  if (SOCP_sp_nBlock > 0){
    NewArray(SOCP_sp_index,int,SOCP_sp_nBLock);
    NewArray(SOCP_sp_block,SparseMatrix,SOCP_sp_nBLock);
  }
  counter = 0;
  for (int l=0; l<SOCP_nBlock; ++l) {
    if (SOCP_NonZeroNumber[l] > 0){
      SOCP_sp_index[counter] = l;
      int size = SOCP_blockStruct[l];
      SOCP_sp_block[counter].initialize(size,size,SparseMatrix::SPARSE,
					SOCP_NonZeroNumber[l]);
      counter++;
    }    
  }
#endif

  // for LP
  for (int l=0; l<LP_nBlock; l++){
    if (LP_NonZeroNumber[l] == true){
      LP_sp_nBlock++;
    }    
  }
  if (LP_sp_nBlock > 0){
    NewArray(LP_sp_index,int,LP_sp_nBlock);
    NewArray(LP_sp_block,double,LP_sp_nBlock);
  }
  counter = 0;
  for (int l=0; l<LP_nBlock; ++l) {
    if (LP_NonZeroNumber[l] == true){
      LP_sp_index[counter] = l;
      counter++;
    }    
  }
}

// sparse form of block index      2008/02/27 kazuhide nakata
void SparseLinearSpace::initialize(int SDP_sp_nBlock, 
                                   int* SDP_sp_index,
                                   int* SDP_sp_blockStruct, 
                                   int* SDP_sp_NonZeroNumber,
                                   int SOCP_sp_nBlock, 
                                   int* SOCP_sp_index,
                                   int* SOCP_sp_blockStruct,
                                   int* SOCP_sp_NonZeroNumber,
                                   int LP_sp_nBlock, 
                                   int* LP_sp_index)
{
  // rMessage("SparseLinearSpace::initialize");

  // for SDP
  this->SDP_sp_nBlock = SDP_sp_nBlock;
  if (SDP_sp_nBlock > 0){
    NewArray(this->SDP_sp_index,int,SDP_sp_nBlock);
    NewArray(this->SDP_sp_block,SparseMatrix,SDP_sp_nBlock);
  }
  for (int l=0; l<SDP_sp_nBlock; ++l) {
    this->SDP_sp_index[l] = SDP_sp_index[l];
    int size = SDP_sp_blockStruct[l];
    SDP_sp_block[l].initialize(size,size,SparseMatrix::SPARSE,
                               SDP_sp_NonZeroNumber[l]);
  }

  // for SOCP
#if 0
  this->SOCP_sp_nBlock = SOCP_sp_nBlock;
  if (SOCP_sp_nBlock > 0){
    NewArray(this->SOCP_sp_index,int,SOCP_sp_nBlock);
    NewArray(this->SOCP_sp_block,SparseMatrix,SOCP_sp_nBlock);
  }
  for (int l=0; l<SOCP_sp_nBlock; ++l) {
    this->SOCP_sp_index[l] = SOCP_sp_index[l];
    int size = SOCP_sp_blockStruct[l];
    SOCP_sp_block[l].initialize(size,size,SparseMatrix::SPARSE,
				SOCP_sp_NonZeroNumber[l]);
  }
#endif

  // for LP
  this->LP_sp_nBlock = LP_sp_nBlock;
  if (LP_sp_nBlock > 0){
    NewArray(this->LP_sp_index,int,LP_sp_nBlock);
    NewArray(this->LP_sp_block,double,LP_sp_nBlock);
  }
  for (int l=0; l<LP_sp_nBlock; ++l) {
    this->LP_sp_index[l] = LP_sp_index[l];
  }
}

void SparseLinearSpace::terminate()
{
  // for SDP
  if (SDP_sp_block && SDP_sp_index && SDP_sp_nBlock>=0) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].terminate();
    }
    DeleteArray(SDP_sp_block);
    DeleteArray(SDP_sp_index);
  }
  // for SOCP
#if 0
  if (SOCP_sp_block && SOCP_sp_index && SOCP_sp_nBlock>=0) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].terminate();
    }
    DeleteArray(SOCP_sp_block);
    DeleteArray(SOCP_sp_index);
  }
#endif 
  // for LP
  if (LP_sp_block && LP_sp_index && LP_sp_nBlock>=0) {
    DeleteArray(LP_sp_block);
    DeleteArray(LP_sp_index);
  }
}

void SparseLinearSpace::changeToDense(bool forceChange)
{
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].changeToDense(forceChange);
    }
  }
#if 0
  if (SOCP_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_sp_block[l].changeToDense(forceChange);
    }
  }
#endif

}

void SparseLinearSpace::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  // SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    fprintf(fpout,"SDP part{\n");
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      fprintf(fpout,"block %d\n",SDP_sp_index[l]);
      SDP_sp_block[l].display(fpout,printFormat);
    }
    fprintf(fpout,"} \n");
  }
  // for SOCP
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    fprintf(fpout,"SOCP part{\n");
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      fprintf(fpout,"block %d\n",SOCP_sp_index[l]);
      SOCP_sp_block[l].display(fpout,printFormat);
    }
    fprintf(fpout,"} \n");
  }
#endif
  // LP
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    fprintf(fpout,"LP part{\n");
    for (int l=0; l<LP_sp_nBlock; ++l) {
      fprintf(fpout,"index: %d, element ",LP_sp_index[l]);
      fprintf(fpout,printFormat,LP_sp_block[l]);
      fprintf(fpout,"\n");
    }
    fprintf(fpout,"} \n");
  }
}

bool SparseLinearSpace::copyFrom(SparseLinearSpace& other)
{
  bool total_judge;

  if (this == &other) {
    return SDPA_SUCCESS;
  }
  if (other.SDP_sp_nBlock+other.SOCP_sp_nBlock+LP_sp_nBlock < 0) {
    rError("SparseLinearSpace:: nBlock is negative");
  }

  // for SDP
  if (other.SDP_sp_nBlock < 0) {
    rError("SparseLinearSpace:: SDP_nBlock is negative");
  }
  if (SDP_sp_nBlock!=other.SDP_sp_nBlock) {
    DeleteArray(SDP_sp_index);
    DeleteArray(SDP_sp_block);
  }
  SDP_sp_nBlock = other.SDP_sp_nBlock;
  if ( SDP_sp_nBlock > 0 && SDP_sp_index==NULL ) {
    NewArray(SDP_sp_index,int,SDP_sp_nBlock);
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_index[l] = other.SDP_sp_index[l];
    }
  }
  if ( SDP_sp_nBlock > 0 && SDP_sp_block==NULL ) {
    NewArray(SDP_sp_block,SparseMatrix,SDP_sp_nBlock);
  }
  total_judge = SDPA_SUCCESS;
  for (int l=0; l<SDP_sp_nBlock; ++l) {
    total_judge = SDP_sp_block[l].copyFrom(other.SDP_sp_block[l]);
  }
  if (total_judge==SDPA_FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }


  // for SOCP
#if 0
  if (other.SOCP_sp_nBlock<0) {
    rError("SparseLinearSpace:: SOCP_nBlock is negative");
  }
  if (SOCP_sp_nBlock!=other.SOCP_sp_nBlock) {
    DeleteArray(SOCP_sp_index);
    DeleteArray(SOCP_sp_block);
  }
  SOCP_sp_nBlock = other.SOCP_sp_nBlock;
  if ( SOCP_sp_nBlock > 0 && SOCP_sp_index==NULL) {
    NewArray(SOCP_sp_index,int,SOCP_sp_nBlock);
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_index[l] = other.SOCP_sp_index[l];
    }
  }
  if ( SOCP_sp_nBlock > 0 && SOCP_sp_block==NULL) {
    NewArray(SOCP_sp_block,SparseMatrix,SOCP_sp_nBlock);
  }
  total_judge = SDPA_SUCCESS;
  for (int l=0; l<SOCP_sp_nBlock; ++l) {
    total_judge = SOCP_sp_block[l].copyFrom(other.SOCP_sp_block[l]);
  }
  if (total_judge==SDPA_FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }
#endif

  // for LP
  if (other.LP_sp_nBlock<0) {
    rError("SparseLinearSpace:: LP_nBlock is negative");
  }
  if (LP_sp_nBlock!=other.LP_sp_nBlock) {
    DeleteArray(LP_sp_index);
    DeleteArray(LP_sp_block);
  }
  LP_sp_nBlock = other.LP_sp_nBlock;
  if ( LP_sp_nBlock > 0 && LP_sp_index==NULL) {
    NewArray(LP_sp_index,int,LP_sp_nBlock);
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_index[l] = other.LP_sp_index[l];
    }
  }
  if ( LP_sp_nBlock > 0 && LP_sp_block==NULL) {
    NewArray(LP_sp_block,double,LP_sp_nBlock);
  }
  total_judge = SDPA_SUCCESS;
  for (int l=0; l<LP_sp_nBlock; ++l) {
    LP_sp_block[l] = other.LP_sp_block[l];
  }
  if (total_judge==SDPA_FAILURE) {
    rError("SparseLinearSpace:: copy miss");
  }


  return total_judge;
}

void SparseLinearSpace::setElement_SDP(int block,
				       int i, int j, double ele)
{
  int l;

  // seek block
  for (l=0; l<SDP_sp_nBlock; l++){
    if (SDP_sp_index[l] == block){
      break;
    }
  }
  if (l == SDP_sp_nBlock){
    rError("SparseLinearSpace::setElement no block");
  }

  // check range
  if (SDP_sp_block[l].NonZeroCount >= SDP_sp_block[l].NonZeroNumber){
    rError("SparseLinearSpace::setElement NonZeroCount >= NonZeroNumber");
  }
  if ((i >= SDP_sp_block[l].nRow) || (j >= SDP_sp_block[l].nCol)){
    rError("out of range in input data");
  }

  // set element
  int count = SDP_sp_block[l].NonZeroCount;
  #if DATA_CAPSULE
  SDP_sp_block[l].DataS[count].vRow   = i;
  SDP_sp_block[l].DataS[count].vCol   = j;
  SDP_sp_block[l].DataS[count].vEle   = ele;
  #else
  SDP_sp_block[l].row_index[count]    = i;
  SDP_sp_block[l].column_index[count] = j;
  SDP_sp_block[l].sp_ele[count]       = ele;
  #endif
  SDP_sp_block[l].NonZeroCount++;
  if (i==j){
    SDP_sp_block[l].NonZeroEffect++;
  } else {
    SDP_sp_block[l].NonZeroEffect += 2;
  }
  
}

void SparseLinearSpace::setElement_SOCP(int block, int i, int j,
					double ele)
{
  rError("DenseLinearSpace:: current version does not support SOCP");
}

void SparseLinearSpace::setElement_LP(int block, double ele)
{
  int l;

  for (l=0; l<LP_sp_nBlock; l++){
    if (LP_sp_index[l] == block){
      break;
    }
  }
  if (l == LP_sp_nBlock){
    rError("SparseLinearSpace::"
	   "setElement cannot find the appropriate block");
  }
  LP_sp_block[l] = ele;
}

void SparseLinearSpace::setZero()
{
  // for SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].setZero();
    }
  }
  // for SOCP
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].setZero();
    }
  }
#endif
  // for LP
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_block[l] = 0;
    }
  }
}

void SparseLinearSpace::setIdentity(double scalar)
{
  rError("SparseLinearSpace::setIdentity   no support");
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l=0; l<SDP_sp_nBlock; ++l) {
      SDP_sp_block[l].setIdentity(scalar);
    }
  }
#if 0
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l=0; l<SOCP_sp_nBlock; ++l) {
      SOCP_sp_block[l].setIdentity(scalar);
    }
  }
  if (LP_sp_nBlock>0 && LP_sp_index && LP_sp_block) {
    for (int l=0; l<LP_sp_nBlock; ++l) {
      LP_sp_block[l].setIdentity(scalar);
    }
  }
#endif
}

bool SparseLinearSpace::sortSparseIndex(int& l, int& i, int& j)
{
  bool total_judge = SDPA_SUCCESS;
  l = -1;
  int i_in,j_in; 
  // for SDP
  if (SDP_sp_nBlock>0 && SDP_sp_index && SDP_sp_block) {
    for (int l_in=0; l_in<SDP_sp_nBlock; ++l_in) {
      total_judge = SDP_sp_block[l_in].sortSparseIndex(i_in,j_in);
      if (total_judge==SDPA_FAILURE && l<0) {
	l = l_in;
	i = i_in;
	j = j_in;
      }
    }
  }
  // for SOCP
  l = -1;
  if (SOCP_sp_nBlock>0 && SOCP_sp_index && SOCP_sp_block) {
    for (int l_in=0; l_in<SOCP_sp_nBlock; ++l_in) {
      total_judge = SOCP_sp_block[l_in].sortSparseIndex(i_in,j_in);
      if (total_judge==SDPA_FAILURE && l<0) {
	l = l_in;
	i = i_in;
	j = j_in;
      }
    }
  }

  return total_judge;
}


DenseLinearSpace::DenseLinearSpace()
{
  SDP_nBlock  = 0;
  SDP_block   = NULL;
  SOCP_nBlock = 0;
  SOCP_block  = NULL;
  LP_nBlock   = 0;
  LP_block    = NULL;
}

DenseLinearSpace::DenseLinearSpace(BlockStruct& bs)
{
  initialize(bs);
}

DenseLinearSpace::~DenseLinearSpace()
{
  terminate();
}

void DenseLinearSpace::initialize(BlockStruct& bs)
{
  // First clean up, then initialzie
  terminate();
  
  this->SDP_nBlock  = bs.SDP_nBlock;
  this->SOCP_nBlock = bs.SOCP_nBlock;
  this->LP_nBlock   = bs.LP_nBlock;
  SDP_block  = NULL;
  SOCP_block = NULL;
  LP_block   = NULL;
  // rMessage("DenseLinearSpace::initialize");
  if (SDP_nBlock + SOCP_nBlock + LP_nBlock <= 0) {
    rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
  }

  // for SDP
  if (SDP_nBlock<0) {
    rError("DenseLinearSpace:: SDP_nBlock is negative");
  }
  if (SDP_nBlock > 0) {
    NewArray(SDP_block,DenseMatrix,SDP_nBlock);
  }
  for (int l=0; l<SDP_nBlock; ++l) {
    int size = bs.SDP_blockStruct[l];
    if (size>0) {
      SDP_block[l].initialize(size,size,DenseMatrix::DENSE);
    } else {
      rError("DenseLinearSpace:: SDP size is nonpositive");
    }
  }

  // for SOCP
  this->SOCP_nBlock = 0;
#if 0
  if (SOCP_nBlock<0) {
    rError("DenseLinearSpace:: SOCP_nBlock is negative");
  }
  if (SOCP_nBlock > 0) {
    NewArray(SOCP_block,DenseMatrix,SOCP_nBlock);
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    int size = SOCP_blockStruct[l];
    if (size>0) {
      SOCP_block[l].initialize(size,size,DenseMatrix::DENSE);
    } else {
      rError("DenseLinearSpace:: SOCP size is nonpositive");
    }
  }
#endif

  // for LP
  if (LP_nBlock<0) {
    rError("DenseLinearSpace:: LP_nBlock is negative");
  }
  if (LP_nBlock > 0) {
    NewArray(LP_block,double,LP_nBlock);
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_block[l] = 0.0;
  }
}

void DenseLinearSpace::terminate()
{
  // for SDP
  if (SDP_block && SDP_nBlock>0) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].terminate();
    }
    DeleteArray(SDP_block);
  }

  // SOCP
#if 0
  if (SOCP_block && SOCP_nBlock>0) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].terminate();
    }
    DeleteArray(SOCP_block);
  }
#endif

  // LP
  if (LP_block && LP_nBlock>0) {
    DeleteArray(LP_block);
  }
}

void DenseLinearSpace::display(FILE* fpout, char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  if (SDP_nBlock>0 && SDP_block) {
    fprintf(fpout,"SDP part{\n");
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }

#if 0
  if (SOCP_nBlock>0 && SOCP_block) {
    fprintf(fpout,"SOCP part{\n");
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].display(fpout);
    }
    fprintf(fpout,"} \n");
  }
#endif

  if (LP_nBlock>0 && LP_block) {
    fprintf(fpout,"LP part{\n");
    for (int l=0; l<LP_nBlock; ++l) {
      fprintf(fpout,printFormat,LP_block[l]);
      fprintf(fpout,", ");
    }
    fprintf(fpout,"} \n");
  }
}

void DenseLinearSpace::displaySolution(BlockStruct& bs,
				       FILE* fpout,
				       char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout,"{\n");
  for (int l=0; l<bs.nBlock; l++){
    if (bs.blockType[l] == BlockStruct::btSDP) {
      int l2 = bs.blockNumber[l];
      SDP_block[l2].display(fpout,printFormat);
    } else if (bs.blockType[l] == BlockStruct::btSOCP) {
      rError("io:: current version does not support SOCP");
      int l2 = bs.blockNumber[l];
      SOCP_block[l2].display(fpout,printFormat);
    } else if (bs.blockType[l] == BlockStruct::btLP) {
      fprintf(fpout,"{");
      int size  = bs.blockStruct[l];
      int start = bs.blockNumber[l];
      for (int l2=0; l2<size-1; ++l2) {
	fprintf(fpout,printFormat,LP_block[start+l2]);
	fprintf(fpout,",");
      }
      if (size > 0) {
	fprintf(fpout,printFormat,LP_block[start+size-1]);
	fprintf(fpout,"}\n");
      } else {
	fprintf(fpout,"  }\n");
      }
    } else {
      rError("io::displayDenseLinearSpaceLast not valid blockType");
    }
  }
  fprintf(fpout,"}\n");
}

bool DenseLinearSpace::copyFrom(DenseLinearSpace& other)
{
  if (this == &other) {
    return SDPA_SUCCESS;
  }

  if (other.SDP_nBlock+other.SOCP_nBlock+other.LP_nBlock<=0) {
    rError("DenseLinearSpace:: SDP + SOCP + LP Block is nonpositive");
  }
  bool total_judge = SDPA_SUCCESS;

  // for SDP
  if (other.SDP_nBlock<0) {
    rError("DenseLinearSpace:: SDP_nBlock is negative");
  }
  if (SDP_nBlock!=other.SDP_nBlock) {
    DeleteArray(SDP_block);
  }
  SDP_nBlock = other.SDP_nBlock;
  if (SDP_nBlock > 0 && SDP_block == NULL) {
    NewArray(SDP_block,DenseMatrix,SDP_nBlock);
  }
  for (int l=0; l<SDP_nBlock; ++l) {
    total_judge = SDP_block[l].copyFrom(other.SDP_block[l]);
  }
  if (total_judge==SDPA_FAILURE) {
    rError("DenseLinearSpace:: copy miss");
  }
  // for SOCP
#if 0
  if (other.SOCP_nBlock<0) {
    rError("DenseLinearSpace:: SOCP_nBlock is negative");
  }
  if (SOCP_nBlock!=other.SOCP_nBlock) {
    DeleteArray(SOCP_block);
  }
  SOCP_nBlock = other.SOCP_nBlock;
  if ( SOCP_nBlock > 0 && SOCP_block == NULL) {
    NewArray(SOCP_block,DenseMatrix,SOCP_nBlock);
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    total_judge = SOCP_block[l].copyFrom(other.SOCP_block[l]);
  }
  if (total_judge==SDPA_FAILURE) {
    rError("DenseLinearSpace:: copy miss");
  }
#endif

  // for LP
  if (other.LP_nBlock<0) {
    rError("DenseLinearSpace:: LP_nBlock is negative");
  }
  if (LP_nBlock!=other.LP_nBlock) {
    delete[] LP_block;
    LP_block = NULL;
  }
  LP_nBlock = other.LP_nBlock;
  if ((LP_nBlock > 0) && (LP_block == NULL)) {
    LP_block = new double[LP_nBlock];
    if (LP_block==NULL) {
      rError("DenseLinearSpace:: memory exhausted");
    }
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_block[l] = other.LP_block[l];
  }
  return total_judge;
}

void DenseLinearSpace::setElement_SDP(int block, int i, int j, double ele)
{

  // check range
  if (block >= SDP_nBlock){
    rError("out of range in input data");
  }
  if ((i >= SDP_block[block].nRow) || (j >= SDP_block[block].nCol)){
    rError("out of range in input data");
  }

  int nCol = SDP_block[block].nCol;
  SDP_block[block].de_ele[i + j * nCol] = ele;
  SDP_block[block].de_ele[j + i * nCol] = ele;
}

void DenseLinearSpace::setElement_SOCP(int block, int i, int j,
				       double ele)
{
  rError("DenseLinearSpace:: current version does not support SOCP");
}

void DenseLinearSpace::setElement_LP(int block, double ele)
{
  // check range
  if (block >= LP_nBlock){
    rError("out of range in input data");
  }
  LP_block[block] = ele;
}

void DenseLinearSpace::setZero()
{
  // for SDP
  if (SDP_nBlock>0 && SDP_block) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].setZero();
    }
  }

  // for SOCP
#if 0  
  if (SOCP_nBlock>0 && SOCP_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].setZero();
    }
  }
#endif

  // for LP
  if (LP_nBlock>0 && LP_block) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_block[l] = 0.0;
    }
  }

}

void DenseLinearSpace::setIdentity(double scalar)
{
  // for SDP
  if (SDP_nBlock>0 && SDP_block) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_block[l].setIdentity(scalar);
    }
  }
  // for SOCP
#if 0
  if (SOCP_nBlock>0 && SOCP_block) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_block[l].setIdentity(scalar);
    }
  }
#endif
  // for LP
  if (LP_nBlock>0 && LP_block) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_block[l] = scalar;
    }
  }
}

} // end of namespace 'sdpa'

