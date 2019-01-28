#include "mex.h"
#include "math.h"
#include "stdlib.h"

double *DATA, *params, *EXTREMA, *return_count;
double epsilon;
int num_extrema, length, count;

void detect() {

	int inc, n, s, save;
	double current, previous, next;
	inc = 1;
	count = 0;
	
	for (s = 0; s < num_extrema; s++) {
		EXTREMA[s] = 0;
	}

	while ((inc < (length-1)) && (count < num_extrema)) {
        save = 1;
		current = DATA[inc];
		previous = current-DATA[inc-1];
		next = DATA[inc+1]-current;
	
		if (previous*next < 0) {
            n = 0;
			while ((n < count) && (save == 1)) {
				if (fabs(EXTREMA[n]-current) < epsilon) {
                    save = 0;
				}
                n++;
			}
            if (save == 1) {
                EXTREMA[count] = current;
                count++;
            }
		}
		inc++;
	}

	if (count == 0) {
		EXTREMA[count] = DATA[inc];
		count++;
	}
    
    return_count[0] = count;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray * prhs[]) {
// Check number of inputs and outputs
	if (nlhs != 2)
		mexErrMsgTxt("Wrong number of outputs!");
	if (nrhs != 2)
		mexErrMsgTxt("Wrong number of inputs!");
	
// Get inputs
	DATA = mxGetPr(prhs[0]);
	params = mxGetPr(prhs[1]);
	num_extrema = params[0];
	epsilon = params[1];
	length = params[2];
	
// Prepare outputs
	plhs[0] = mxCreateDoubleMatrix(num_extrema, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
	EXTREMA = mxGetPr(plhs[0]);
    return_count = mxGetPr(plhs[1]);
	
// MAIN
	detect();
}
