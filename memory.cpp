#include "functions.h"


void mem_allocate()
{
	NEW1 (dx,		double,		NX);
	NEW1 (x1,		double,		NX);
	NEW1 (xc,		double,		NX);
	NEW1 (x2,		double,		NX);

    NEW2 (xgp, double, NX, NGP);
    NEW1 (jgp, double, NX);
    NEW1 (wgp, double, NGP);

	NEW2 (ro,		double,		NX,		NF);
	NEW2 (ru,		double,		NX,		NF);
	NEW2 (re,		double,		NX,		NF);

	NEW2 (ro_,		double,		NX,		NF);
	NEW2 (ru_,		double,		NX,		NF);
	NEW2 (re_,		double,		NX,		NF);

	NEW2 (ro_1,		double,		NX,		NF);
	NEW2 (ru_1,		double,		NX,		NF);
	NEW2 (re_1,		double,		NX,		NF);

	NEW2 (ro_2,		double,		NX,		NF);
	NEW2 (ru_2,		double,		NX,		NF);
	NEW2 (re_2,		double,		NX,		NF);

	NEW2 (ro_half,		double,		NX,		NF);
	NEW2 (ru_half,		double,		NX,		NF);
	NEW2 (re_half,		double,		NX,		NF);

	NEW2 (ro_int,	double,		NX,		NF);
	NEW2 (ru_int,	double,		NX,		NF);
	NEW2 (re_int,	double,		NX,		NF);

	A = static_cast<double ***>(malloc(sizeof(double **) * NX));
	for (int i = 0; i < NX; i++) {
		A[i] = static_cast<double **>(malloc(sizeof(double *) * NF));
		for (int j = 0; j < NF; j++) {
			A[i][j] = static_cast<double *>(malloc(sizeof(double) * NF));
		}
	}
}

void mem_deallocate()
{
	DELETE1 (dx);
	DELETE1 (x1);
	DELETE1 (xc);
	DELETE1 (x2);

    DELETE2 (xgp, NX);
    DELETE1 (jgp);
    DELETE1 (wgp);

	DELETE2 (ro,		NX);
	DELETE2 (ru,		NX);
	DELETE2 (re,		NX);

	DELETE2 (ro_,		NX);
	DELETE2 (ru_,		NX);
	DELETE2 (re_,		NX);

	DELETE2 (ro_1,		NX);
	DELETE2 (ru_1,		NX);
	DELETE2 (re_1,		NX);

	DELETE2 (ro_2,		NX);
	DELETE2 (ru_2,		NX);
	DELETE2 (re_2,		NX);

	DELETE2 (ro_half,		NX);
	DELETE2 (ru_half,		NX);
	DELETE2 (re_half,		NX);

	DELETE2 (ro_int,	NX);
	DELETE2 (ru_int,	NX);
	DELETE2 (re_int,	NX);

	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NF; j++) {
			free( A[i][j] );
		}
		free( A[i] );
	}
	free( A );
}
