#include "mex.h"


/* This is a booster function replacement for isPointInAnyTrinagle, a
 * local function in GRID3DPR.isPointInDomain. It tests if a (vector of) 
 * point(s) pt(k) is in any triangle. If yes, b(k) = 1 else b(k).
 *
 * b = isPointInAnyTriangle3DPR(pt,p1,xi_x,eta_x,xi_y,...
 *                              eta_y,xi_z,eta_z,zeta_z) 
 */
 
 


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *b;     
    double *pt;    
    double *p1;
    
    double *xi_x;
    double *eta_x;
    
    double *xi_y;
    double *eta_y;
    
    double *xi_z;
    double *eta_z;
    double *zeta_z;
    
    mwSize nE;
    mwSize nP;
    
    mwSize k;
    mwSize l;
    
    double ptx;
    double pty;
    double ptz;
    
    double ptux;
    double ptuy;
    double ptuz;
    
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:nlhs",
                          "One output required.");}
    
    if(nrhs != 9) {
        mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:nrhs",
                          "Eight inputs required.");}
     
    if(mxGetM(prhs[0])!=3) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notMatrix","Input must be a 3 x N matrix.");
     }
      
    if(mxGetM(prhs[1])!=3) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notMatrix","Input must be a 3 x nElements matrix.");
      }
    if(mxGetM(prhs[2])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[3])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector");
      }
    if(mxGetM(prhs[4])!=1) {
         mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[5])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[6])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[7])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      }
    if(mxGetM(prhs[8])!=1) {
          mexErrMsgIdAndTxt("GRID3DPR:isPointInAnyTriagle:notRowVector","Input must be a vector.");
      } 
     
    
    nP = mxGetN(prhs[0]);  
    nE = mxGetN(prhs[1]); 
     
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)nP,mxREAL);
   
    
    /* assign the pointers to the inputs
     * to build the interface
     * b = isPointInAnyTriangle(pt,p1,xi_x,eta_x,xi_y,eta_y,xi_z,eta_z,zeta_z)
     */
    
    /* Left hand side*/
    b =  mxGetPr(plhs[0]); 
    
    /* right hand side */
    pt      =   mxGetPr(prhs[0]);    
    p1      =   mxGetPr(prhs[1]);
    xi_x    =   mxGetPr(prhs[2]);
    eta_x   =   mxGetPr(prhs[3]);   
    xi_y    =   mxGetPr(prhs[4]);
    eta_y   =   mxGetPr(prhs[5]);
    xi_z    =   mxGetPr(prhs[6]);
    eta_z   =   mxGetPr(prhs[7]);
    zeta_z  =   mxGetPr(prhs[8]);
    
    
    
    
    for(k=0;k<nP;k++)
    {            
        b[k] = 0.0;         
        for (l=0;l<nE;l++)
        {   
            ptx = pt[3*k]-p1[3*l];
            pty = pt[3*k+1]-p1[3*l+1];
            ptz = pt[3*k+2]-p1[3*l+2];
            
            ptux = (xi_x[l]*ptx)+(xi_y[l]*pty)+(xi_z[l]*ptz);
            ptuy = (eta_x[l]*ptx)+(eta_y[l]*pty)+(eta_z[l]*ptz);
            ptuz = zeta_z[l]*ptz;
                    
            
            if (((ptux+ptuy)<=1.0)&&(ptux>=0.0)&&(ptuy>=0.0)&&(abs(ptuz)<=1.0)) 
            {
                b[k] = 1.0; 
                break;
            }            
        }                    
    } 
}
