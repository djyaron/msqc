/*
 * C = elementWiseCombine(A, B)
 *
 * Description:
 *   Equivalent to sum(sum(A .* B)) for 2D matrices.
 *   Sometimes, MATLAB isn't very good at handling this, so the same
 *   functionality is offered here.
 *
 * Author:
 *   Alex Cappiello (acappiel@andrew.cmu.edu)
 *
 * Date:
 *   6-4-13
 * Updated:
 *   6-7-13
 *
 * Inputs:
 *   prhs[0]: A (m * n) Matrix.
 *   prhs[1]: B (m * n) Matrix.
 *
 * Outputs:
 *   plhs[0]: C (1 * 1) sum(sum(A .* B)).
 *
 * Notes: 
 *   Sometimes, MATLAB does a better job with with this calculation on its own.
 *   Other times, this function gives a more significant improvement.
 *
 * Compile Line: 
 *   mex elementWiseCombine.c
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/matlab/programming-interfaces-for-c-c-fortran-com.html
 *
 */

#include "mex.h"
#include "matrix.h"

/* Input verification macro. */
#define CHECK(A) ((A) ? (void)0 : \
mexErrMsgIdAndTxt("elementWiseCombine:invalidDim", \
"Invalid dimension found in elementWiseCombine input."))

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *A, *B, *C;
    size_t i, m, n, nelem;
    const mxArray *mA, *mB;
    
    if (nlhs > 1)
        mexErrMsgIdAndTxt("elementWiseCombine:nlhs",
                "Invalid number of output variables for elementWiseCombine.");
    if (nrhs != 2)
        mexErrMsgIdAndTxt("elementWiseCombine:nrhs",
                "Invalid number of input variables to elementWiseCombine.");
    
    /* Pull data out of the input arrays. */
    mA = prhs[0];
    A = mxGetPr(mA);
    mB = prhs[1];
    B = mxGetPr(mB);
    m = mxGetM(mA);
    n = mxGetN(mA);
    nelem = m * n;
    
    /* Validate input. */
    CHECK(mxGetNumberOfDimensions(mA) == 2);
    CHECK(mxGetNumberOfDimensions(mB) == 2);
    CHECK(m == mxGetM(mB));
    CHECK(n == mxGetN(mB));

    /* Allocate the return array. */
    C = mxMalloc(sizeof(double));
    *C = 0.f;
    
    /* Carry out the calculation. */
    for (i = 0; i < nelem; i++)
        *C += A[i] * B[i];
    
    /* Prepare the data to hand back to MATLAB. */
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], C);
    mxSetM(plhs[0], 1);
    mxSetN(plhs[0], 1);
    
    return;
}
