/*
 * @author: Alex Cappiello (acappiel)
 * @date: 6-23-12
 *
 * Inputs:
 *   prhs[0]: 1 x nmodel cell array of model3 objects.
 *   prhs[1]: 1 x nmodel cell array of envs to include.
 *   I'll probably assume that nenv is the same for each model.
 *   prhs[2]: nenv x nmodel cell array of obj.H1 data.
 *            each element is nbasis * nbasis.
 *   prhs[3]: nenv x nmodel cell array of obj.H2 data.
 *            each element is 1 * nbasis.
 *
 * Compile Line: mex evalModelsParallel.c -IC:\dev\pthreads-w32\include ...
 *               -LC:\dev\pthreads-w32\lib\x64 -lpthreadVC2 -lmwblas -lmwlapack
 *
 * MEX Library Reference:
 *   http://www.mathworks.com/help/techdoc/apiref/bqoqnz0.html
 *
 */

#include "evalModelsParallel.h"

#ifdef DEBUG
# define dbg_printf(...) fprintf(logfile, __VA_ARGS__);fflush(logfile)
#else
# define dbg_printf(...)
#endif

FILE *logfile;
sem_t semModels, semHF;

bool checkInput(int nrhs, const mxArray *prhs[]);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxArray *model, *envs;
    mwSize imod, nmodel;
    pthread_t *tid;
    void **args;
    
#ifdef DEBUG
logfile = fopen("out.log", "w");
if (!checkInput(nrhs, prhs)) return;
#endif

dbg_printf("Entered mexFunction\n");
/*dbg_printf("Time: %ld\n", time(NULL));*/

sem_init(&semModels, 0, MAXTHREADMODEL);
sem_init(&semHF, 0, MAXTHREADHF);

nmodel = mxGetDimensions(prhs[0])[1];
tid = mxCalloc(nmodel, sizeof(pthread_t));

dbg_printf("nmodel: %d\n", nmodel);
getchar();

for (imod = 0; imod < nmodel; imod++) {
    dbg_printf("imod: %d\n", imod);
    model = mxGetCell(prhs[0], imod);
    envs = mxGetCell(prhs[1], imod);
    args = mxCalloc(5, sizeof(void*));
    args[0] = (void*)model;
    args[1] = (void*)envs;
    args[2] = (void*)prhs[2];
    args[3] = (void*)prhs[3];
    args[4] = VAL(imod);
    sem_wait(&semModels);
    pthread_create(tid+imod, NULL, callSolveHF, (void*)(args));
}
for (imod = 0; imod < nmodel; imod++)
    pthread_join(tid[imod], NULL);
mxFree(tid);

dbg_printf("Leaving mexFunction\n");
return;
}

bool checkInput(int nrhs, const mxArray *prhs[]) {
    dbg_printf("Entered checkInput\n");
    TEST(nrhs == 4, "Wrong # inputs.");
    TEST(mxIsCell(prhs[0]), "Input 1 not cell.");
    TEST(mxIsCell(prhs[1]), "Input 2 not cell.");
    TEST(mxIsCell(prhs[2]), "Input 3 not cell.");
    TEST(mxIsCell(prhs[3]), "Input 4 not cell.");
    
    TEST(mxGetNumberOfDimensions(prhs[0]) == 2 &&
            mxGetDimensions(prhs[0])[0] == 1, "Input 1 has too many dimensions");
    TEST(mxGetNumberOfDimensions(prhs[1]) == 2 &&
            mxGetDimensions(prhs[1])[0] == 1, "Input 2 has too many dimensions");
    TEST(mxGetNumberOfDimensions(prhs[2]) == 2,
            "Input 3 has too many dimensions");
    TEST(mxGetNumberOfDimensions(prhs[3]) == 2,
            "Input 4 has too many dimensions");
    
    TEST(mxGetDimensions(prhs[0])[1] == mxGetDimensions(prhs[1])[1],
            "Input 1 and Input 2 must be same length.");
    /*TEST(mxGetDimensions(prhs[0])[1] == mxGetDimensions(prhs[2])[1] &&
    mxGetDimensions(prhs[1])[1] == mxGetDimensions(prhs[2])[0],
    "Input 3 dimenstions don't match Input 1 and Input 2.");
    TEST(mxGetDimensions(prhs[0])[1] == mxGetDimensions(prhs[3])[1] &&
    mxGetDimensions(prhs[1])[1] == mxGetDimensions(prhs[3])[0],
    "Input 4 dimenstions don't match Input 1 and Input 2.");*/
    dbg_printf("Leaving checkInput\n");
    return true;
}

void *callSolveHF (void *args) {
    mxArray *modelMat, *envs, *H1, *H2;
    model3 modelC;
    int imod;
    
    dbg_printf("Entered callSolveHF\n");
    
    modelMat = (mxArray*)(((void**)args)[0]);
    envs = (mxArray*)(((void**)args)[1]);
    H1 = (mxArray*)(((void**)args)[2]);
    H2 = (mxArray*)(((void**)args)[3]);
    imod = INT((((void**)args)[4]));
    
    mxFree(args);
    
    modelC = makeModelC(modelMat, envs, H1, H2, imod);
    solveHF(modelMat, modelC);
    
    dbg_printf("Leaving callSolveHF\n");
    sem_post(&semModels);    
    return NULL;
}

model3 makeModelC(mxArray *modelMat, mxArray *envs, mxArray *H1, mxArray *H2,
        int imod) {
    model3 modelC;
    mxArray *frag;
    double *denvs;
    int ienv, nenv;

    dbg_printf("Entered makeModelC\n");
    
    modelC = mxCalloc(1, sizeof(struct model3));
    frag = mxGetProperty(modelMat, 0, "frag");
    modelC->frag = frag;
    
    modelC->densitySave = mxGetProperty(modelMat, 0, "densitySave");
    modelC->Hnuc = mxGetPr(mxGetProperty(frag, 0, "Hnuc"))[0];
    modelC->HnucEnv = mxGetPr(mxGetProperty(frag, 0, "HnucEnv"));
    modelC->nbasis = (int)mxGetPr(mxGetProperty(modelMat, 0, "nbasis"))[0];
    modelC->nelec = (int)(mxGetPr(mxGetProperty(frag, 0, "nelec"))[0]);
    modelC->nenvFrag = (int)(mxGetPr(mxGetProperty(modelMat, 0, "nenv"))[0]);
    modelC->X = mxGetPr(mxGetProperty(modelMat, 0, "X"));
    
    nenv = (int)(mxGetDimensions(envs)[1]);
    modelC->nenv = nenv;
    
    denvs = mxGetPr(envs);
    modelC->envs = mxCalloc(nenv, sizeof(int));
    modelC->H1 = mxCalloc(nenv, sizeof(double*));
    modelC->H2 = mxCalloc(nenv, sizeof(double*));
    for (ienv = 0; ienv < nenv; ienv++) {
        modelC->H1[ienv] = mxGetPr(mxGetCell(H1, imod * nenv + ienv));
        modelC->H2[ienv] = mxGetPr(mxGetCell(H2, ((imod) * (nenv) + (ienv))));
        modelC->envs[ienv] = ((int)denvs[ienv]) - 1;
    }        

    dbg_printf("Leaving makeModelC\n");
    return modelC;
}

void solveHF (mxArray *modelMat, model3 modelC) {
    double *orb, *Eorb, *Ehf, *densitySave;
    double *EhfEnv, *EorbEnv, *orbEnv;
    int i, j, k, ienv, nenv = modelC->nenv, nbasis = modelC->nbasis;
    int *envs = modelC->envs;
    mxArray *EhfNew, *EorbNew, *orbNew, *densitySaveNew;
    mwSize dim[3];
    pthread_t *tid;
    void **args, *res;
    
    dbg_printf("Entered solveHF\n");
    
    tid = mxCalloc(nenv, sizeof(pthread_t));
    
    for (ienv = 0; ienv < nenv; ienv++) {
        args = mxCalloc(4, sizeof(void*));
        args[0] = (void*)modelMat;
        args[1] = (void*)modelC;
        args[2] = VAL(ienv);
        sem_wait(&semHF);
        pthread_create(tid+ienv, NULL, callHartreeFock, (void*)(args));
    }
    
    dbg_printf("Done creating HF threads\n");
    
    /* need to memcpy??? */
    EhfEnv  = mxGetPr(mxGetProperty(modelMat, 0, "EhfEnv"));
    EorbEnv = mxGetPr(mxGetProperty(modelMat, 0, "EorbEnv"));
    orbEnv  = mxGetPr(mxGetProperty(modelMat, 0, "orbEnv"));
    
    for (ienv = 0; ienv < nenv; ienv++) {
        pthread_join(tid[ienv], &res);
        orb         = ((void**)res)[0];
        Eorb        = ((void**)res)[1];
        Ehf         = ((void**)res)[2];
        densitySave = ((void**)res)[3];
        
        if (envs[ienv] == 0) {
            EhfNew  = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            EorbNew = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            orbNew  = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
            
            mxSetPr(EhfNew, Ehf);
            mxSetPr(EorbNew, Eorb);
            mxSetPr(orbNew, orb);
            
            dim[0] = 1;
            dim[1] = 1;
            mxSetDimensions(EhfNew, dim, 2);
            dim[0] = modelC->nbasis;
            mxSetDimensions(EorbNew, dim, 2);
            dim[1] = modelC->nbasis;
            mxSetDimensions(orbNew, dim, 2);
            
            mxSetProperty(modelMat, 0, "Ehf", EhfNew);
            mxSetProperty(modelMat, 0, "Eorb", EorbNew);
            mxSetProperty(modelMat, 0, "orb", orbNew);
        }
        else {
            k = (int)envs[ienv];
            EhfEnv[k] = *Ehf;
            for (j = 0; j < nbasis; j++) {
                EorbEnv[CMO2(j, k, nbasis)] = Eorb[j];
                for (i = 0; i < nbasis; i++)
                    orbEnv[CMO3(i, j, k, nbasis, nbasis)] = orb[CMO2(i, j, nbasis)];
            }
        }
        
        densitySaveNew = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
        mxSetPr(densitySaveNew, densitySave);
        dim[0] = nbasis;
        dim[1] = nbasis;
        mxSetDimensions(densitySaveNew, dim, 2);
        mxSetCell(modelC->densitySave, (mwIndex)envs[ienv], densitySaveNew);
        
    }
    
    EhfNew  = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    EorbNew = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    orbNew  = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    
    mxSetPr(EhfNew, EhfEnv);
    mxSetPr(EorbNew, EorbEnv);
    mxSetPr(orbNew, orbEnv);
    
    dim[0] = 1;
    dim[1] = modelC->nenv;
    mxSetDimensions(EhfNew, dim, 2);
    dim[0] = modelC->nbasis;
    mxSetDimensions(EorbNew, dim, 2);
    dim[1] = modelC->nbasis;
    dim[2] = modelC->nenv;
    mxSetDimensions(orbNew, dim, 3);
    
    mxSetProperty(modelMat, 0, "EhfEnv", EhfNew);
    mxSetProperty(modelMat, 0, "EorbEnv", EorbNew);
    mxSetProperty(modelMat, 0, "orbEnv", orbNew);
    
    mxFree(tid);
    
    dbg_printf("Leaving solveHF\n");    
    return;
}

void *callHartreeFock (void *args) {
    mxArray *modelMat = (mxArray*)((void**)args)[0];
    model3 modelC = (model3)((void**)args)[1];
    int ienv = INT(((void**)args)[2]);
    double **res;
    
    mxFree(args);
    
    res = hartreeFock(modelMat, modelC, ienv);
    
    sem_post(&semHF);
    pthread_exit(res);
    
    return res;
}

double **hartreeFock (mxArray *modelMat, model3 modelC, int ienv) {
    bool finished;
    char *norm = "N", *trans = "T", *jobz = "V", *range = "A", *uplo = "U";
    double Ee, *Ehf, t1, t2, eps = 1.0e-8, scalar = 1.0, zero = 0.0;
    double *C, *F, *Ft, *H1, *H2, *H2j, *H2k, *P, *Pn, *Plast, *tmp, **res;
    double *work, *Ct, *Ct1, *e1, *e, *Eorb, *orb, *workspace;
    int *isuppz, *iwork, lwork, liwork, info, *sort, filled;
    int minIter = 5, maxIter = 1000, iter, nbasis = modelC->nbasis, one = 1;
    int inenv = modelC->envs[ienv];
    int i, j, k, l, neig;
    mxArray *densityEnv, *densityNoEnv;
    const mwSize *dimEnv, *dimNoEnv;
    size_t H2Size = nbasis * nbasis * nbasis * nbasis;
    
    dbg_printf("Entered hartreeFock\n");
    
    H1  = modelC->H1[ienv];
    H2  = modelC->H2[ienv];
    H2j = NEWARRAY(H2Size, sizeof(double));
    H2k = NEWARRAY(H2Size, sizeof(double));
    
    for (i = 0; i < nbasis; i++)
        for (j = 0; j < nbasis; j++)
            for (k = 0; k < nbasis; k++)
                for (l = 0; l < nbasis; l++) {
        H2j[RMO4(i, j, k, l, nbasis, nbasis, nbasis)] =
                H2[CMO4(i, j, k, l, nbasis, nbasis, nbasis)];
        H2k[RMO4(i, j, k, l, nbasis, nbasis, nbasis)] =
                H2[CMO4(i, k, l, j, nbasis, nbasis, nbasis)];
                }
    
    densityEnv = mxGetCell(modelC->densitySave, ienv);
    densityNoEnv = mxGetCell(modelC->densitySave, 0);
    dimEnv   = (densityEnv != NULL) ? mxGetDimensions(densityEnv) : NULL;
    dimNoEnv = (densityNoEnv != NULL) ? mxGetDimensions(densityNoEnv): NULL;
    
    P     = NEWARRAY(nbasis * nbasis, sizeof(double));
    Pn    = NEWARRAY(nbasis * nbasis, sizeof(double));
    Plast = NEWARRAY(nbasis * nbasis, sizeof(double));
    
    dbg_printf("Allocated stuff\n");
    
    /* Step 4: Guess at density matrix. */
    if ((dimEnv == NULL || SIZE(dimEnv) == 0) &&
            (dimNoEnv == NULL || SIZE(dimNoEnv) == 0)) {
        dbg_printf("  Choice 1\n");
        fragDensity(modelMat, modelC, inenv, Pn);
    }
    else if (dimEnv == NULL || SIZE(dimEnv) == 0) {
        dbg_printf("  Choice 2\n");
        memcpy(Pn, mxGetPr(mxGetCell(modelC->densitySave, 0)),
                nbasis * nbasis * sizeof(double));
    }
    else {
        dbg_printf("  Choice 3\n");
        memcpy(Pn, mxGetPr(mxGetCell(modelC->densitySave, ienv)),
                nbasis * nbasis * sizeof(double));
    }
    
    dbg_printf("memcpy?\n");
    memcpy(Plast, Pn, nbasis * nbasis * sizeof(double));
    iter = 0;
    finished = false;
    
    dbg_printf("Guessed at density matrix\n");
    
    /* Allocate intermediate matrices. */
    F         = NEWARRAY(nbasis * nbasis, sizeof(double));
    workspace = NEWARRAY(nbasis * nbasis, sizeof(double));
    Ct        = NEWARRAY(nbasis * nbasis, sizeof(double));
    e         = NEWARRAY(nbasis, sizeof(double));
    sort      = NEWARRAY(nbasis, sizeof(int));
    C         = NEWARRAY(nbasis * nbasis, sizeof(double));
    
    /* dsyevr data. */
    Ct1     = NEWARRAY(nbasis * nbasis, sizeof(double));
    e1      = NEWARRAY(nbasis, sizeof(double));
    isuppz  = NEWARRAY(2 * nbasis, sizeof(int));
    work    = NEWARRAY(26 * nbasis, sizeof(double));
    iwork   = NEWARRAY(10 * nbasis, sizeof(int));
    lwork   = 26 * nbasis;
    liwork  = 10 * nbasis;
    
    dbg_printf("Allocated more stuff\n");
    
    while (!finished) {
        /*dbg_printf("  Start loop iteration %d\n", iter);*/
        scalar = 1.0;
        tmp = P; /* Done with this data. Recycle the memory for the new Pn. */
        P   = Pn;
        Pn  = tmp;
        
        if (iter < maxIter / 2)
            for (j = 0; j < nbasis; j++)
                for (i = 0; i < nbasis; i++) {
            k = j * nbasis + i;
            P[k] = 0.5 * Pn[k] + 0.5 * Plast[k];
                }
        /* else P = Pn. */
        
        /* Step 5+6: Build 2-electron components of Fock matrix and the fock matrix. */
        for (i = 0; i < nbasis; i++) {
            for (j = 0; j < nbasis; j++) {
                t1 = 0.0;
                t2 = 0.0;
                for (k = 0; k < nbasis; k++) {
                    for (l = 0; l < nbasis; l++) {
                        /* Accessing P in RMO acts as transposed CMO. */
                        t1 += H2j[RMO4(i, j, k, l, nbasis, nbasis, nbasis)] *
                                P[k * nbasis + l];
                        t2 += H2k[RMO4(i, j, k, l, nbasis, nbasis, nbasis)] *
                                P[k * nbasis + l];
                    }
                }
                F[CMO2(i, j, nbasis)] = H1[CMO2(i, j, nbasis)] + t1 - t2 / 2;
            }
        }
        
        /* Step 7: Calculate the transformed F matrix. */
        /* Ft = (X' * F) * X */
        dgemm(trans, norm, BPTR(&nbasis), BPTR(&nbasis), BPTR(&nbasis), &scalar,
                modelC->X, BPTR(&nbasis), F, BPTR(&nbasis), &zero, workspace,
                BPTR(&nbasis));
        Ft = F; /* Recycle the memory. */
        dgemm(norm, norm, BPTR(&nbasis), BPTR(&nbasis), BPTR(&nbasis), &scalar,
                workspace, BPTR(&nbasis), modelC->X, BPTR(&nbasis), &zero, Ft,
                BPTR(&nbasis));
        
        /* Step 8: Find e and the transformed expansion coeff matrix. */
        dsyevr(jobz, range, uplo, BPTR(&nbasis), Ft, BPTR(&nbasis), NULL, NULL,
                BPTR(&one), BPTR(&nbasis), &eps, BPTR(&neig), e1, Ct1, 
                BPTR(&nbasis), BPTR(isuppz), work, BPTR(&lwork), BPTR(iwork),
                BPTR(&liwork), BPTR(&info));
        
        sorted(e1, sort, nbasis); /* qsort??? */
        for (j = 0; j < nbasis; j++) {
            k = sort[j];
            e[j] = e1[k];
            for (i = 0; i < nbasis; i++) {
                Ct[CMO2(i, j, nbasis)] = Ct1[CMO2(i, k, nbasis)];
            }
        }
        
        /* Step 9: Transform Ct back to C. */
        /* C = X * Ct */
        dgemm(norm, norm, BPTR(&nbasis), BPTR(&nbasis), BPTR(&nbasis), &scalar,
                modelC->X, BPTR(&nbasis), Ct, BPTR(&nbasis), &zero, C, BPTR(&nbasis));
        
        /* Step 10: Calculate the new density matrix. */
        tmp   = Plast;
        Plast = Pn;
        Pn    = tmp;
        /*memset(Pn, 0, nbasis * nbasis);*/
        filled = modelC->nelec / 2;
        scalar = 2.0;
        dgemm(norm, trans, BPTR(&nbasis), BPTR(&nbasis), BPTR(&filled), &scalar,
                C, BPTR(&nbasis), C, BPTR(&nbasis), &zero, Pn, BPTR(&nbasis));
        
        iter++;
        
        if (iter > maxIter)
            finished = true;
        else if (iter > minIter)
            if (maxDiff(P, Pn, nbasis, nbasis) < eps)
                finished = true;
    }
    
    res = NEWARRAY(4, sizeof(double*));
    tmp = P;
    P = Pn;
    
    orb = C;
    Eorb = e;
    
    Ehf = NEWARRAY(1, sizeof(double));
    Ee = 0;
    for (j = 0; j < nbasis; j++)
        for (i = 0; i < nbasis; i++) {
        k = CMO2(i, j, nbasis);
        Ee += P[k] * (H1[k] + F[k]);
        }
    *Ehf = Ee / 2 + ((ienv == 0) ? modelC->Hnuc : modelC->HnucEnv[ienv]);
    
    FREE(H2j);
    FREE(H2k);
    FREE(tmp);
    FREE(Plast);
    FREE(F);
    FREE(workspace);
    FREE(Ct);
    FREE(sort);
    FREE(C);
    FREE(Ct1);
    FREE(e1);
    FREE(isuppz);
    FREE(work);
    FREE(iwork);
    
    res[0] = orb;
    res[1] = Eorb;
    res[2] = Ehf;
    res[3] = P;
    
    if (iter + 1 > maxIter)
        printf("You are living on the edge... hartree fock didn't converge.\n");
    
    dbg_printf("Leaving hartreeFock\n");
    
    return res;
}

void fragDensity (mxArray *modelMat, model3 modelC, int ienv, double *res) {
    int i, j, nocc, nbasis = modelC->nbasis;
    double *orb, *filledOrbs;
    char norm = 'N', trans = 'T';
    double scalar = 2.0, zero = 0.0;
    
    dbg_printf("Entered fragDensity\n");
    
    nocc = modelC->nelec / 2;
    filledOrbs = NEWARRAY(nbasis * nocc, sizeof(double));
    
    if (ienv == 0) {
        orb = mxGetPr(mxGetProperty(modelC->frag, 0, "orb"));
        
        for (j = 0; j < nocc; j++) {
            for (i = 0; i < nbasis; i++) {
                filledOrbs[((j) * (nbasis) + (i))] = orb[((j) * (nbasis) + (i))];
            }
        }
    }
    else {
        orb = mxGetPr(mxGetProperty(modelC->frag, 0, "orbEnv"));
        
        /* replace with memcpy? */
        for (j = 0; j < nocc; j++) {
            for (i = 0; i < nbasis; i++) {
                filledOrbs[((j) * (nbasis) + (i))] =
                        orb[((i) + (nbasis) * ((j) + (nbasis) * (ienv)))];
            }
        }
    }
    
    dgemm(&norm, &trans, BPTR(&nbasis), BPTR(&nbasis), BPTR(&nocc), &scalar,
            filledOrbs, BPTR(&nbasis), filledOrbs, BPTR(&nbasis), &zero, res,
            BPTR(&nbasis));
    FREE(filledOrbs);
    dbg_printf("Leaving fragDensity\n");
}

void sorted (double *in, int *out, int size) {
    int i, j, mini;
    for (i = 0; i < size; i++) {
        mini = 0;
        for (j = 1; j < size; j++) {
            if (in[j] < in[mini])
                mini = j;
        }
        out[i] = mini;
        in[mini] = (double)HUGE_VAL;
    }
}

double maxDiff (double *A, double *B, int row, int col) {
    int i, j;
    double diff, max = 0;
    for (j = 0; j < col; j++)
        for (i = 0; i < row; i++) {
        diff = fabs(A[CMO2(i, j, row)] - B[CMO2(i, j, row)]);
        if (diff > max)
            max = diff;
        }
    return diff;
}
