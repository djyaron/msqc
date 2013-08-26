/*
 * G = twoElecFock(P, H2j, H2k)
 *
 * Description:
 *   This function builds the 2-electron components of the Fock matrix (G) from 
 *   the current density matrix (P), the coulomb integrals (H2j), and the
 *   exchange integrals (H2k).
 *
 * Author:
 *   Alex Cappiello (acappiel@andrew.cmu.edu)
 *
 * Date:
 *   7-19-12
 * Updated:
 *   3-21-13
 *
 * Inputs:
 *   prhs[0]: P (nbasis * nbasis) density matrix.
 *   prhs[1]: H2j {nbasis,nbasis} cell array of coulomb integrals.
 *   prhs[2]: H2k {nbasis,nbasis} cell array of exchange integrals.
 *
 * Outputs:
 *   plhs[0]: G (nbasis * nbasis) 2-electron components of Fock matrix.
 *
 * Notes: 
 *   P, G, and each H2j is symmetric. Each H2k is not symmetric.
 *   The code validates the dimensions of the inputs using the CHECK macro.
 *     Overhead appears to be minimal.
 *
 * Compile Line: 
 *   mex twoElecFock.c
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/matlab/programming-interfaces-for-c-c-fortran-com.html
 *
 */

#include "mex.h"
#include "matrix.h"

/* 2D Column-Major-Order array lookup pattern. */
#define CMO2(i, j, I) ((j) * (I) + (i))

/* Input verification macro. */
#define CHECK(A) ((A) ? (void)0 : \
mexErrMsgIdAndTxt("twoElecFock:invalidDim", \
"Invalid dimension found in twoElecFock input."))

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double t1, t2;
    double *G, *P, *H2jj, *H2kk;
    int i, j, k, nbasis;
    mxArray *H2j, *H2k, *mH2jj, *mH2kk;
    
    if (nlhs > 1)
        mexErrMsgIdAndTxt("twoElecFock:nlhs",
                "Invalid number of output variables for twoElecFock.");
    if (nrhs != 3)
        mexErrMsgIdAndTxt("twoElecFock:nrhs",
                "Invalid number of input variables to twoElecFock.");
    
    /* Pull data out of the input arrays. */
    P = mxGetPr(prhs[0]);
    H2j = (mxArray*)prhs[1];
    H2k = (mxArray*)prhs[2];
    nbasis = (int)mxGetM(prhs[0]);
    
    /* Validate input. */
    CHECK(nbasis == mxGetM(prhs[0]));
    CHECK(nbasis == mxGetN(prhs[0]));
    CHECK(nbasis == mxGetM(H2j));
    CHECK(nbasis == mxGetN(H2j));
    CHECK(nbasis == mxGetM(H2k));
    CHECK(nbasis == mxGetN(H2k));

    /* Allocate the return array. */
    G = mxCalloc(nbasis * nbasis, sizeof(double));
    
    /* 
     * Theoretically, this order of loop variables is the most cache efficient,
     * but testing doesn't show any appreciable difference in practice.
     */
    for (j = 0; j < nbasis; j++) {
        for (i = j; i < nbasis; i++) {
            t1 = 0.0;
            t2 = 0.0;
            
            /* Get current arrrays out of the cell arrays and validate input. */
            mH2jj = mxGetCell(H2j, CMO2(i, j, nbasis));
            mH2kk = mxGetCell(H2k, CMO2(i, j, nbasis));
            CHECK(nbasis == mxGetM(mH2jj));
            CHECK(nbasis == mxGetN(mH2jj));
            CHECK(nbasis == mxGetM(mH2kk));
            CHECK(nbasis == mxGetN(mH2kk));
            H2jj = mxGetPr(mH2jj);
            H2kk = mxGetPr(mH2kk);
            
            /* Build the elements of G. */
            for (k = 0; k < nbasis * nbasis; k++) {
                t1 += H2jj[k] * P[k];
                t2 += H2kk[k] * P[k];
            }
            
            G[CMO2(i, j, nbasis)] = t1 - t2 / 2;
            
            /* Taking advantage of symmetry. */
            if (i != j)
                G[CMO2(j, i, nbasis)] = t1 - t2 / 2;
        }
    }
    
    /* Prepare the data to hand back to MATLAB. */
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], G);
    mxSetM(plhs[0], nbasis);
    mxSetN(plhs[0], nbasis);
    
    return;
}
