#include "mex.h"
#include "math.h"
#include "stdlib.h"
//Chi-Hak UY 2017
//Matlab I/O
double *params_laser, *params_temp, *PAST, *RESULTS, *X, *Xd, *dX, *K;


//Parameters
double ka, eta, omega, alpha, beta, T, P, tau_R, theta, h, horizon; // eta, alpha, tau_R, P
int steps, delay;


//RK coefficients and parameters
double a[4] =
{       0.0,            1.0/3.0,          2.0/3.0,           1.0};
double b[4][4] = {
{            0.0,           0.0,              0.0,           0.0},
{        1.0/3.0,           0.0,              0.0,           0.0},
{       -1.0/3.0,           1.0,              0.0,           0.0},
{            1.0,          -1.0,              1.0,           0.0}};
double c[4] =
{       1.0/8.0,        3.0/8.0,          3.0/8.0,       1.0/8.0};


int DIM = 5, S = 4;


//Equations of the system
void eval(double *dX, double *X, double *Xd) {
    dX[0] = X[0]*X[4] + eta*X[2]*cos(X[3] - X[1]);
    dX[1] = alpha*X[4] + eta*X[2]/X[0]*sin(X[3] - X[1]);
    dX[2] = (Xd[0]*cos(Xd[1] + X[3]) - X[2])/tau_R;
    dX[3] = -Xd[0]*sin(Xd[1] + X[3])/(tau_R*X[2]);
    dX[4] = (P - X[4] - (1.00 + 2.0*X[4])*X[0]*X[0])/T;
} // eta, alpha, tau_R, P, T




//Runge-Kutta
void simu() {
    int i, j, s, n;
    X = (double *)malloc(DIM *sizeof(double));     //input of eval
    Xd = (double *)malloc(DIM *sizeof(double));    //delayed input of eval
    dX = (double *)malloc(DIM *sizeof(double));    //output of eval
    K = (double *)malloc(S*DIM *sizeof(double));   //for intermediate values in RK
    
    //RESULTS initialized with PAST
    for (i=0; i<delay+1; i++){
        for (n=0; n<DIM; n++){
            RESULTS[n+DIM*i] = PAST[n+DIM*i];
        }
    }
    
    //K initialized to zero
    for (s=0; s<S; s++){
        for (n=0; n<DIM; n++){
            K[n+DIM*s] = 0;
        }
    }
    
    for (i=delay; i<steps+delay+1; i++){
        for (s=0; s<S; s++){
            for (n=0; n<DIM; n++){
                X[n] = RESULTS[n+DIM*i];
                Xd[n] = RESULTS[n+DIM*(i-delay)]*(1-a[s]) + RESULTS[n+DIM*(i-delay+1)]*a[s];
            }
            
            for (j=0; j<s; j++){
                for (n=0; n<DIM; n++)
                    X[n] += h*b[s][j]*K[n+DIM*j];
            }
            
            eval(dX, X, Xd);
            
            for (n=0; n<DIM; n++)
                K[n+DIM*s] = dX[n];
        }
        
        for (n=0; n<DIM; n++){
            RESULTS[n+DIM*(i+1)] = RESULTS[n+DIM*i];
            for (s=0; s<S; s++)
                RESULTS[n+DIM*(i+1)] += h*c[s]*K[n+DIM*s];
        }
    }
    
    free(dX);
    free(X);
    free(Xd);
    free(K);
}


//Adapt Matlab and C I/O
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
    //Check number of inputs and outputs
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of outputs!");
    if (nrhs !=3)
        mexErrMsgTxt("Wrong number of inputs!");
    
    //Get inputs
    params_laser = mxGetPr(prhs[0]);
    params_temp = mxGetPr(prhs[1]);
    PAST = mxGetPr(prhs[2]);
    
    //Parameters [eta, omega, alpha, beta, ka, T, P, theta, tau_R];
    //           [eta, alpha, tau_R, P, T]
    eta = params_laser[0];
    omega = params_laser[1];
    alpha = params_laser[2];
    beta = params_laser[3];
    ka = params_laser[4];
    T = params_laser[5];
    P = params_laser[6];
    theta = params_laser[7];
    tau_R = params_laser[8];
    
    h = params_temp[0];
    horizon = params_temp[1];
    
    steps = floor(horizon/h);
    delay = floor(theta/h);
    
    //Prepare OUTPUT matrix
    plhs[0] = mxCreateDoubleMatrix(DIM*(steps+delay+1), 1, mxREAL);
    RESULTS = mxGetPr(plhs[0]);
    
    //MAIN
    simu();
}
