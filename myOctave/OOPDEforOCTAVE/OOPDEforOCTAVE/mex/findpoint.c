// Function to find a point (a,b) in a matrix k x 2.
// Returns the index of the found point or NaN.
// Call with: 'findpoint(point,matrix)'.
// Make sure you've compiled the findpoint.c beforehand on your specific
// system with: 'mex findpoint.c' so MATLAB can use it.
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *newPoints, *point;
  mwSize NrRows, k;
  newPoints  = mxGetPr(prhs[0]);  // pointer to matrix newPoints
  NrRows     = mxGetM(prhs[0]);   // number of rows of newPoints
  point      = mxGetPr(prhs[1]);  // pointer to the point we look for  
  
  for (k = 0; k < NrRows; k++)    // for every row do:
    {
       // The left statement is checked first, if and only if it's correct
       // the second one gets checked aswell.
       // Since the newPoints matrix is read as newPoints(:) in C, we
       // have to add NrRows to k to jump to the next 'column'.
       if (newPoints[k] == point[0] && newPoints[k+NrRows] == point[1])
       {
           // If a matching row is found, return the index k+1, since C
           // starts indexing arrays at 0 instead of 1 like MATLAB.
           plhs[0] = mxCreateDoubleScalar((double) (k + 1));
           return;
       }
    }
  // If for some strange reason the point couldn't be found, return NaN.
  plhs[0] = mxCreateDoubleScalar(mxGetNaN());
}