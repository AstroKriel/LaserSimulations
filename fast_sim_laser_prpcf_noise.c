// Authors: Chi-Hak UY (2017)
//          Neco Kriel (2019)
#include "mex.h"
#include "math.h"
#include "stdlib.h"

// Initialise Matlab I/O
double *params_laser, *params_temp, *PAST, *RESULTS, *X, *Xd, *dX, *K, *NOISE;

// Initialise Parameters Used
double P, T, theta, eta, beta, ka, alpha, tau_R, omega, R, h, horizon;
int steps, delay;

// RK coefficients and parameters
double a[4] =
{       0.0,            1.0/3.0,          2.0/3.0,           1.0};
double b[4][4] = {
{            0.0,           0.0,              0.0,           0.0},
{        1.0/3.0,           0.0,              0.0,           0.0},
{       -1.0/3.0,           1.0,              0.0,           0.0},
{            1.0,          -1.0,              1.0,           0.0}};
double c[4] =
{       1.0/8.0,        3.0/8.0,          3.0/8.0,       1.0/8.0};

int DIM = 7, S = 4;

// Equations of the system
void eval(double *dX, double *X, double *Xd) {
    dX[0] = X[0]*X[6]; // E_1
    dX[1] = alpha*X[6]; // phi_1
    dX[2] = ka*X[2]*(X[6]-beta) + eta*sqrt(ka)*X[4]*cos(X[5]-X[3]); // E_2
    dX[3] = alpha*ka*(X[6]-beta) + eta*sqrt(ka)*X[4]*sin(X[5]-X[3])/X[2] - omega; // phi_2
    dX[4] = (Xd[0]*cos(Xd[1]+X[5]) - X[4])/tau_R; // E_F
    dX[5] = (-Xd[0]*sin(Xd[1]+X[5]))/(X[4]*tau_R); // phi_F
    dX[6] = P - X[6] - (1+2*X[6])*(X[0]*X[0] + X[2]*X[2]); // Z
}

// Runge-Kutta
void simu() {
    int i, j, s, n;
    X   = (double *)malloc(DIM *sizeof(double));     // input of eval
    Xd  = (double *)malloc(DIM *sizeof(double));     // delayed input of eval
    dX  = (double *)malloc(DIM *sizeof(double));     // output of eval
    K   = (double *)malloc(S*DIM *sizeof(double));   // for intermediate values in RK
    // RESULTS initialized with PAST
    for (i=0; i<delay+1; i++){
        for (n=0; n<DIM; n++){
            RESULTS[n+DIM*i] = PAST[n+DIM*i];
        }
    }
    // K initialized to zero
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
        RESULTS[0+DIM*(i+1)] += sqrt(h*R)*NOISE[0+2*(i-delay)];
        RESULTS[2+DIM*(i+1)] += sqrt(h*R)*NOISE[1+2*(i-delay)];
    }
    free(dX);
    free(X);
    free(Xd);
    free(K);
}

// Adapt Matlab and C I/O
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
    // Check number of inputs and outputs
    if (nlhs != 1)
        mexErrMsgTxt("Wrong number of outputs!");
    if (nrhs !=4)
        mexErrMsgTxt("Wrong number of inputs!");
    
    // Get inputs
    params_laser = mxGetPr(prhs[0]);
    params_temp  = mxGetPr(prhs[1]);
    PAST         = mxGetPr(prhs[2]);
    NOISE        = mxGetPr(prhs[3]);
    
    // Parameters Received: [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R]
    // Parameters Used:     [P, T, theta, eta, beta, ka, alpha, tau_R, omega,  R]
    P       = params_laser[0];
    T       = params_laser[1];
    theta   = params_laser[2];
    eta     = params_laser[3];
    beta    = params_laser[4];
    ka      = params_laser[5];
    alpha   = params_laser[6];
    tau_R   = params_laser[7];
    omega   = params_laser[8];
    R       = params_laser[9];
    
    h       = params_temp[0];
    horizon = params_temp[1];
    
    steps   = floor(horizon/h);
    delay   = floor(theta/h);
    
    // Prepare OUTPUT matrix
    plhs[0] = mxCreateDoubleMatrix(DIM*(steps+delay+1), 1, mxREAL);
    RESULTS = mxGetPr(plhs[0]);
    
    // MAIN
    simu();
}
