#include "mex.h"


/* This is a booster function replacement for isPointInAnyTrinagle, a
 * local function in grid2D.isPointInDomain. It tests if a (vector of) 
 * point(s) pt(k) is in any triangle. If yes, b(k) = 1 else b(k).
 *
 * b = isPointInAnyTriangle(pt,p1,x21,x31,y21,y31,J);
 */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* declare variables: nE and nP number of elements/points
     */
    mwSize k, l, nE, nP;
    
    /* local variables 
     */
    double ptx, pty, ptux, ptuy;
    
    /* delare the parameters. all pointers to double matrices/vectors
     */   
    double *b,*pt, *p1,*x21, *x31,*y21,*y31,*J;     
    
    /* some input check
     */
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:nlhs",
                          "One output required.");}
    
    if(nrhs != 7) {
        mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:nrhs",
                          "Seven inputs required.");}
    
    if(mxGetM(prhs[0])!=2) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notMatrix","Input must be a 2 x N matrix.");
     }
      
    if(mxGetM(prhs[1])!=2) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notMatrix","Input must be a 2 x nElements matrix.");
      }
    if(mxGetM(prhs[2])!=1) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notRowVector","Input must be vector.");
      }
    if(mxGetM(prhs[3])!=1) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notRowVector","Input must be vector.");
      }
    if(mxGetM(prhs[4])!=1) {
         mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[5])!=1) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[6])!=1) {
          mexErrMsgIdAndTxt("GRID2D:isPointInAnyTriagle:notRowVector","Input must be a  vector.");
      }
     
    /* get the number of point and elements
     */
    nE = mxGetN(prhs[1]); 
    nP = mxGetN(prhs[0]);  
    
    
    /* create the output vector 
     */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)nP,mxREAL);
    
    /* assign the pointers to the output
     */
    b =  mxGetPr(plhs[0]);
    
    /* assign the pointers to the inputs
     */
    pt  = mxGetPr(prhs[0]);    
    p1  = mxGetPr(prhs[1]);
    x21 = mxGetPr(prhs[2]);
    y21 = mxGetPr(prhs[3]);   
    x31 = mxGetPr(prhs[4]);
    y31 = mxGetPr(prhs[5]);
    J   = mxGetPr(prhs[6]);
    
    for(k=0;k<nP;k++)
    {
        b[k] = 0.0;        
        for (l=0;l<nE;l++)
        {            
            ptx = pt[2*k]-p1[2*l];
            pty = pt[2*k+1]-p1[2*l+1];
            
            ptux =  (y31[l]*ptx-x31[l]*pty)/J[l];
            ptuy =  (x21[l]*pty-y21[l]*ptx)/J[l];
            
            if (((ptux+ptuy)<=1.0) && (ptux>=0.0) && (ptuy>=0.0))                 
            {
                b[k] = 1.0;
                break;
            }   
        }        
    }        
}
