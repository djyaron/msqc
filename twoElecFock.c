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
 *   In the interest of not adding code overhead, no checks are done
 *     to see that the input data matches this above criteria, so be careful.
 *     Segfaults in MEX functions are not pretty.
 *   P, G, and each H2j is symmetric. Each H2k is not symmetric.
 *
 * Compile Line: 
 *   mex twoElecFock.c
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/techdoc/apiref/bqoqnz0.html
 *
 */

#include "mex.h"
#include "matrix.h"

/* 2D Column-Major-Order array lookup pattern. */
#define CMO2(i, j, I) ((j) * (I) + (i))

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double t1, t2;
    double *G, *P, *H2jj, *H2kk;
    int i, j, k, l, nbasis;
    mxArray *H2j, *H2k;
    
    /* Pull data out of the input arrays. */
    P = mxGetPr(prhs[0]);
    H2j = (mxArray*)prhs[1];
    H2k = (mxArray*)prhs[2];
    nbasis = (int)mxGetM(prhs[0]);
    
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
            
            /* Get current arrrays out of the cell arrays. */
            H2jj = mxGetPr(mxGetCell(H2j, CMO2(i, j, nbasis)));
            H2kk = mxGetPr(mxGetCell(H2k, CMO2(i, j, nbasis)));
            
            /* Build the elements of G. */
            for (l = 0; l < nbasis; l++) {
                for (k = 0; k < nbasis; k++) {
                    t1 += H2jj[CMO2(k, l, nbasis)] * P[CMO2(k, l, nbasis)];
                    t2 += H2kk[CMO2(k, l, nbasis)] * P[CMO2(k, l, nbasis)];
                }
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
