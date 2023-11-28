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

#ifndef __sdpa_block_h__
#define __sdpa_block_h__

#include "sdpa_include.h"

namespace sdpa {

class BlockStruct
{
public:
  enum BlockType {btSDP,btSOCP,btLP};
  int  nBlock;
  int* blockStruct;
  int* blockNumber;
  BlockType* blockType;
  int  SDP_nBlock;
  int* SDP_blockStruct;
  int  SOCP_nBlock;
  int* SOCP_blockStruct;
  int  LP_nBlock;

  BlockStruct();
  ~BlockStruct();
  void initialize(int nBlock);
  void terminate();
  void makeInternalStructure();
  void display(FILE* fpOut = stdout);

};

} //end of namespace 'sdpa'

#endif // __sdpa_block_h__
