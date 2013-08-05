/*
 * C = mm(S, A, [B])
 *
 * Description:
 *   This function is a general purpose matrix multiply, where the result is
 *   multiplied by a scalar. Combining these two steps gives a marginal
 *   improvement in speed.
 *   If B is not specifed, then A' is used.
 *
 * Author:
 *   Alex Cappiello (acappiel@andrew.cmu.edu)
 *
 * Date:
 *   6-4-13
 * Updated:
 *   6-4-13
 *
 * Inputs:
 *   prhs[0]: S (1 * 1) Scalar.
 *   prhs[1]: A (m * n) Matrix.
 *   prhs[2]: B (n * k) Matrix. Optional.
 *
 * Outputs:
 *   plhs[0]: C (m * k) S * A * B.
 *
 * Notes: 
 *   If S = 1, then you should just multiply normally and not use this.
 *
 * Compile Line: 
 *   mex mm.c "C:\Program Files\MATLAB\R2012b\extern\lib\win64\microsoft\libmwblas.lib"
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/matlab/programming-interfaces-for-c-c-fortran-com.html
 *
 */

#include "mex.h"
#include "matrix.h"
#include "blas.h"

#if !defined(_WIN32)
# define dgemm dgemm_
#endif

/* Input verification macro. */
#define CHECK(A) ((A) ? (void)0 : \
mexErrMsgIdAndTxt("mm:invalidDim", \
"Invalid dimension found in mm input."))

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *S, *A, *B, *C;
    double zero = 0.f;
    ptrdiff_t m, n, k, ldb;
    const mxArray *mS, *mA, *mB;
    char *op;
    
    if (nlhs > 1)
        mexErrMsgIdAndTxt("mm:nlhs",
                "Invalid number of output variables for mm.");
    if (nrhs < 2)
        mexErrMsgIdAndTxt("mm:nrhs",
                "Invalid number of input variables for mm.");
    
    /* Pull data out of the input arrays. */
    mS = prhs[0];
    S = mxGetPr(mS);
    mA = prhs[1];
    A = mxGetPr(mA);
    m = (ptrdiff_t)mxGetM(mA);
    k = (ptrdiff_t)mxGetN(mA);
    
    /* Are we running C = S * A * A' or C = S * A * B ? */
    switch(nrhs) {
        case 2: {
            mB = prhs[1];
            B = A;
            op = "T";
            n = m;
            ldb = n;
            break;
        }
        case 3: {
            mB = prhs[2];
            B = mxGetPr(mB);
            op = "N";
            n = (ptrdiff_t)mxGetN(mB);
            ldb = k;
            CHECK(mxGetN(mA) == mxGetM(mB));
            break;
        }
        default: {
            mexErrMsgIdAndTxt("mm:nrhs",
                    "Invalid number of input variables to mm.");
        }
    }
    
    /* Validate input. */
    CHECK(mxGetNumberOfDimensions(mA) == 2);
    CHECK(mxGetNumberOfDimensions(mB) == 2);
    CHECK(mxGetNumberOfDimensions(mS) == 2);
    CHECK(mxGetM(mS) == 1);
    CHECK(mxGetN(mS) == 1);

    /* Allocate the return array. */
    C = mxMalloc(m * n * sizeof(double));
    
    /* Carry out the calculation. */
    dgemm("N", op, &m, &n, &k, S, A, &m, B, &ldb, &zero, C, &m);
    
    /* Prepare the data to hand back to MATLAB. */
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    mxSetPr(plhs[0], C);
    mxSetM(plhs[0], (int)m);
    mxSetN(plhs[0], (int)m);
    
    return;
}
