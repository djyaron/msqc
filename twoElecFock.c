/*
 * @author: Alex Cappiello (acappiel)
 * @date: 7-19-12
 *
 * Inputs:
 *   prhs[0]: P nbasis * nbasis density matrix.
 *   prhs[1]: H2j {nbasis,nbasis} cell array of coulomb integrals.
 *   prhs[2]: H2k {nbasis,nbasis} cell array of exchange integrals.
 *
 * Compile Line: mex twoElecFock.c
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/techdoc/apiref/bqoqnz0.html
 *
 */

#include "mex.h"
#include "matrix.h"

#define CMO2(i, j, I) ((j) * (I) + (i))

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double t1, t2;
    double *G, *P, *H2jj, *H2kk;
    int i, j, k, l, nbasis;
    mxArray *H2j, *H2k;
    
    P = mxGetPr(prhs[0]);
    H2j = prhs[1];
    H2k = prhs[2];
    nbasis = mxGetM(prhs[0]);
    G = mxCalloc(nbasis * nbasis, sizeof(double));
    
    for (i = 0; i < nbasis; i++) {
        for (j = i; j < nbasis; j++) {
            t1 = 0.0;
            t2 = 0.0;
            H2jj = mxGetPr(mxGetCell(H2j, CMO2(i, j, nbasis)));
            H2kk = mxGetPr(mxGetCell(H2k, CMO2(i, j, nbasis)));
            for (k = 0; k < nbasis; k++) {
                for (l = 0; l < nbasis; l++) {
                    t1 += H2jj[CMO2(k, l, nbasis)] * P[CMO2(k, l, nbasis)];
                    t2 += H2kk[CMO2(k, l, nbasis)] * P[CMO2(k, l, nbasis)];
                }
            }
            G[CMO2(i, j, nbasis)] = t1 - t2 / 2;
            if (i != j)
                G[CMO2(j, i, nbasis)] = t1 - t2 / 2;
        }
    }
    
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], G);
    mxSetM(plhs[0], nbasis);
    mxSetN(plhs[0], nbasis);
    
    return;
}
