#include  "functions.h"
#include <cassert>

double get_f(int i_cell, int i_func, double x) {
	switch (i_func) {
		case 0:
			return 1.0;
		case 1:
			return (x-xc[i_cell])/dx[i_cell];
		case 2:
			return (x-xc[i_cell])*(x-xc[i_cell])/dx[i_cell]/dx[i_cell];
		default:
			assert(false);
	}
}

double get_df(const int i_cell, const int i_func, double x) {
	switch (i_func)	{
		case 0:
			return 0.0;
		case 1:
			return 1.0/dx[i_cell];
		case 2:
			return 2.0*(x-xc[i_cell])/dx[i_cell]/dx[i_cell];
		default:
			assert(false);
	}
}

double get_field(const double x, int i, fid_t field_id) {
	double * field;
	switch (field_id) {
		case ID_RO:		field = ro[i];		break;
		case ID_RU:		field = ru[i];		break;
		case ID_RE:		field = re[i];		break;
		default:		assert(false);
	}
	double result = 0.0;
	for (int i_func = 0; i_func < NF; i_func++)
		result += field[i_func]*get_f(i, i_func, x);
	return result;
}

double get_field(double **fld, int i_cell, double x)
{
    return fld[i_cell][0] + fld[i_cell][1] * get_f(i_cell, 1, x) +  fld[i_cell][2] * get_f(i_cell, 2, x);
}


double get_field_avg(double **fld, int i_cell)
{
    return fld[i_cell][0] + fld[i_cell][2];// 12.;
}


