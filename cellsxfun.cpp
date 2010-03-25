//
// cellsxfun.cpp
// Copyright Stuart Yarrow 2010/03/24 (s.yarrow@ed.ac.uk)
// All rights reserved.
//
// Based on mAryCellFcn.ccp by Michael Brost (michaelbrost@yahoo.com).
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include "mex.h"

// STD includes
#include <string>
#include <sstream>
using namespace std;

void usage(void)
{
    mexPrintf("========================================================================================================================\n");
    mexPrintf("| This function is a cell array-based generalization of bsxfun to arbitrary dimensions using cell arrays.              |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| usage: cell_array = cellsxfun(function_handle, cell_array_1, cell_array_2, ..., cell_array_N)                        |\n");
    mexPrintf("|       where function_handle corresponds to a function which accepts N simultaneous cell array contents as inputs and |\n");
    mexPrintf("|       which returns at most one matlab object. The dimensions of each array must either match or be singleton e.g:   |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("|       x = {rand(20), rand(20), rand(20)};             1 x 3 cell array                                               |\n");
    mexPrintf("|       y = {rand(20)};                                 1 x 1 cell array                                               |\n");
    mexPrintf("|       z = {rand(20) rand(20)}';                       2 x 1 cell array                                               |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("|       out = cellsxfun(@plus, x, y, z);                2 x 3 cell array                                               |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("|       Each argument is treated as if singleton dimensions are replicated to the size of the output.                  |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| NOTE: there were limited opportunities to trap errors. I suggest that you defensively program your functions.        |\n");
    mexPrintf("|       Anonymous function handles are supported.                                                                      |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| Copyright Stuart Yarrow 2010/03/24 (s.yarrow@ed.ac.uk)                                                               |\n");
    mexPrintf("| Based on mAryCellFcn.ccp by Michael Brost (michaelbrost@yahoo.com)                                                   |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| This program is free software: you can redistribute it and/or modify                                                 |\n");
    mexPrintf("| it under the terms of the GNU General Public License as published by                                                 |\n");
    mexPrintf("| the Free Software Foundation, either version 3 of the License, or                                                    |\n");
    mexPrintf("| (at your option) any later version.                                                                                  |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| This program is distributed in the hope that it will be useful,                                                      |\n");
    mexPrintf("| but WITHOUT ANY WARRANTY; without even the implied warranty of                                                       |\n");
    mexPrintf("| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                        |\n");
    mexPrintf("| GNU General Public License for more details.                                                                         |\n");
    mexPrintf("|                                                                                                                      |\n");
    mexPrintf("| You should have received a copy of the GNU General Public License                                                    |\n");
    mexPrintf("| along with this program.  If not, see <http://www.gnu.org/licenses/>.                                                |\n");
    mexPrintf("========================================================================================================================\n");
}

// gateway function
void mexFunction( int nOutArgs, mxArray *outPtr[], int nInArgs, const mxArray *inPtr[] )
{
    // string for error messages
    string errMsg;
    
    // usage statement
    if(nInArgs == 0)
    {
        usage();
        return;
    }   
    
    // check input argument counts
    if(nInArgs < 2)
    {
        errMsg = string("function requires at least 2 input args.");
        mexErrMsgTxt(errMsg.c_str());
    }
    
    // check input argument counts
    if(nOutArgs > 1)
    {
        errMsg = string("function supports at most 1 output arg.");
        mexErrMsgTxt(errMsg.c_str());
    }
    
    // check arg 1
    if(mxGetClassID(inPtr[0]) != mxFUNCTION_CLASS)
    {
        errMsg = string("function argument #1 must be a function handle.");
        mexErrMsgTxt(errMsg.c_str());
    }
    
    // check input args - all MUST be cell arrays
    for(int aIndex=1; aIndex<nInArgs; aIndex++)
    {
        if(mxGetClassID(inPtr[aIndex]) != mxCELL_CLASS)
        {
            stringstream ss;
            ss << (aIndex + 1);
            errMsg = "function argument #" + ss.str() + " must be a cell array.";
            mexErrMsgTxt(errMsg.c_str());
        }
    }
    
    // number of cell array inputs
    int nInputs = nInArgs-1;
    
    // calculate output dimensionality
    int nDim = 0;
    for(mwSize argIndex=1; argIndex < nInArgs; argIndex++)
    {
        int nArgDims = mxGetNumberOfDimensions(inPtr[argIndex]);
        if(nArgDims > nDim)
        {
            nDim = nArgDims;
        }
    }
    
    // use matlab memory allocation for vector counter
    mwSize *endVec   = static_cast<mwSize*>(mxCalloc(nDim, sizeof(mwSize)));
    mwSize *startVec = static_cast<mwSize*>(mxCalloc(nDim, sizeof(mwSize)));
    mwSize *thisVec  = static_cast<mwSize*>(mxCalloc(nDim, sizeof(mwSize)));
    mwSize *workVec  = static_cast<mwSize*>(mxCalloc(nDim, sizeof(mwSize)));
    mwSize *argDims  = static_cast<mwSize*>(mxCalloc(nInputs*nDim, sizeof(mwSize)));
    
    // calculate output sizes
    for(mwSize argIndex=0; argIndex < nInputs; argIndex++)
    {
        int nArgDims = mxGetNumberOfDimensions(inPtr[argIndex+1]);
        const mwSize *argD = mxGetDimensions(inPtr[argIndex+1]);
        
        for(mwSize dimIndex=0; dimIndex < nDim; dimIndex++)
        {
            if(dimIndex <= nArgDims)
            {
                argDims[argIndex*nDim + dimIndex] = argD[dimIndex];
            
                if(argD[dimIndex] > endVec[dimIndex])
                {
                    endVec[dimIndex] = argD[dimIndex];
                }
            }
            else
            {
                argDims[argIndex*nDim + dimIndex] = 1;
            }
        }
    }
    
    // check input dims
    for(mwSize argIndex=0; argIndex < nInputs; argIndex++)
    {
        for(mwSize dimIndex=0; dimIndex < nDim; dimIndex++)
        {
            // dims need to either match output dimensionality or be singleton
            if(argDims[argIndex*nDim + dimIndex] != endVec[dimIndex] && argDims[argIndex*nDim + dimIndex] != 1)
            {
                errMsg = "cellsxfun: dimension mismatch\n";
                mexErrMsgTxt(errMsg.c_str());
            }
        }
    }
    
    // create the output cell array
    // normally we'd subtract the startVec but it is all 0s so we just use endVec
    outPtr[0] = mxCreateCellArray(nDim, endVec);
	
    // create an array of pointers to the data
    mxArray **dataPtrArray = (mxArray **)mxCalloc(nInputs+1, sizeof(mxArray *));
    
    // copy and set the function handle in the argument array
    //dataPtrArray[0] = mxDuplicateArray(inPtr[0]);
	dataPtrArray[0] = (mxArray *)(inPtr[0]);
    
    // pointer to the results
    mxArray *resultArray;
    
    // using thisVec, count between startVec and endVec
    for(;;)
    {
        // vector of pointers to cell's contents
        for(mwSize vIndex=0; vIndex < nInputs; vIndex++)
        {
            for(mwSize dIndex=0; dIndex < nDim; dIndex++)
            {
                // calculate which cell we are using from input array
                if(argDims[vIndex*nDim + dIndex] == 1)
                    workVec[dIndex] = 0;
                else
                    workVec[dIndex] = thisVec[dIndex];
            }
            
            // get the contents of the input cell arrays as determined by
            // the vector counter's contents
            mwIndex inputIndex = mxCalcSingleSubscript(inPtr[vIndex+1], nDim, workVec);
            dataPtrArray[vIndex+1] = mxGetCell(inPtr[vIndex+1], inputIndex);
            
            // if the cell array was not initialized correctly, the return pointer
            // will be null.
            if(!(dataPtrArray[vIndex+1]))
            {   
                errMsg = "cellsxfun: could not access data from input cell array\n";
                mexErrMsgTxt(errMsg.c_str());
            }
        }
        
        // calculate the offset to the output cell to update
        mwIndex index = mxCalcSingleSubscript(outPtr[0], nDim, thisVec);
        
        // invoke the user's function here - accept only one output and discrard the others
        int status = false;
        status = mexCallMATLAB(1, &resultArray, (int)nInputs+1, dataPtrArray, "feval");
        
        // was there an error which did not abort directly from mexCallMATLAB first?
        if(status)
        {
            errMsg = "cellsxfun: an error occured while invoking the user's function.\n";
            mexErrMsgTxt(errMsg.c_str());
        }
                    
        // copy the results from calcResult[0] into the cell contents located at offset index
        mxSetCell(outPtr[0], index, mxDuplicateArray(resultArray));
        
        // clear the temporarily allocated memory block
        if(resultArray) mxDestroyArray(resultArray);
        
        // update the counters - we will skip to here if the data was bad (uninitialized ?)
        thisVec[0]++;
        for(mwSize vIndex=0; vIndex < (nDim-1); vIndex++)
        {
            if(thisVec[vIndex] >= endVec[vIndex])
            {
                thisVec[vIndex] = startVec[vIndex];
                thisVec[vIndex+1]++;
            }
        }
        
        // terminal condition - here we break out of the loop
        if(thisVec[nDim-1] >= endVec[nDim-1]) break;
    }
    
    // deallocation of allocated stuff
	//if(dataPtrArray[0]) mxFree(dataPtrArray[0]);
    if(dataPtrArray) mxFree(dataPtrArray);
    if(endVec)       mxFree(endVec);
    if(startVec)     mxFree(startVec);
    if(thisVec)      mxFree(thisVec);
    if(workVec)      mxFree(workVec);
    if(argDims)      mxFree(argDims);
}
