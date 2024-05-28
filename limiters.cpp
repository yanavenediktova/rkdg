#include "functions.h"


inline double min_(double x, double y, double z)
{
	if (x <= y && x <= z) return x;
	if (y <= x && y <= z) return y;
	return z;
}

inline double sign_(double x) {if (x < 0.0) return -1.0; if (x == 0.0) return 0.0; return 1.0;}

double minmod(double x, double y, double z)
{
    if ( sign_(x) == sign_(y) && sign_(y) == sign_(z) )
    {
        return sign_(x) * min_(fabs(x), fabs(y), fabs(z));
    }
    return 0.0;
}


double minmod_(double x, double y, double z)
{
    if ( fabs(x) <= M_LIM*HX*HX )
    {
        return x;
    }
    return minmod(x, y, z);
}


#define AVG(fld, cell) ((fld)[(cell)][0]+(fld)[(cell)][2]/12.)


void calc_limiters_()
{
	for (int i = 1; i < NX-1; i++)
	{
		double ro_l0 = AVG(ro,i);
		double ro_l1 = ro[i][1];
		double ru_l0 = AVG(ru,i);
		double ru_l1 = ru[i][1];
		double re_l0 = AVG(re,i);
		double re_l1 = re[i][1];

		ro_l1 = minmod_(ro_l1, LIM_ALPHA*(ro_l0-AVG(ro,i-1)), LIM_ALPHA*(AVG(ro,i+1)-ro_l0));
		ru_l1 = minmod_(ru_l1, LIM_ALPHA*(ru_l0-AVG(ru,i-1)), LIM_ALPHA*(AVG(ru,i+1)-ru_l0));
		re_l1 = minmod_(re_l1, LIM_ALPHA*(re_l0-AVG(re,i-1)), LIM_ALPHA*(AVG(re,i+1)-re_l0));

		if (ro_l1 != ro[i][1])
		{
			ro[i][0] = ro_l0;
			ro[i][1] = ro_l1;
			ro[i][2] = 0.0;
		}
		if (ru_l1 != ru[i][1])
		{
			ru[i][0] = ru_l0;
			ru[i][1] = ru_l1;
			ru[i][2] = 0.0;
		}
		if (re_l1 != re[i][1])
		{
			re[i][0] = re_l0;
			re[i][1] = re_l1;
			re[i][2] = 0.0;
		}
	}

	{
		int i = 0;
		double ro_l0 = ro[i][0]+ro[i][2]/12.0;
		double ro_l1 = ro[i][1];
		double ru_l0 = ru[i][0]+ru[i][2]/12.0;
		double ru_l1 = ru[i][1];
		double re_l0 = re[i][0]+re[i][2]/12.0;
		double re_l1 = re[i][1];

		ro_l1 = minmod_(ro_l1, LIM_ALPHA*(ro[i][0]-ro[i][0])/HX, LIM_ALPHA*(ro[i+1][0]-ro[i][0])/HX);
		ru_l1 = minmod_(ru_l1, LIM_ALPHA*(ru[i][0]-ru[i][0])/HX, LIM_ALPHA*(ru[i+1][0]-ru[i][0])/HX);
		re_l1 = minmod_(re_l1, LIM_ALPHA*(re[i][0]-re[i][0])/HX, LIM_ALPHA*(re[i+1][0]-re[i][0])/HX);

		if (ro_l1 != ro[i][1])
		{
			ro[i][0] = ro_l0;
			ro[i][1] = ro_l1;
			ro[i][2] = 0.0;
		}
		if (ru_l1 != ru[i][1])
		{
			ru[i][0] = ru_l0;
			ru[i][1] = ru_l1;
			ru[i][2] = 0.0;
		}
		if (re_l1 != re[i][1])
		{
			re[i][0] = re_l0;
			re[i][1] = re_l1;
			re[i][2] = 0.0;
		}
	}

	{
		int i = NX-1;
		double ro_l0 = ro[i][0]+ro[i][2]/12.0;
		double ro_l1 = ro[i][1];
		double ru_l0 = ru[i][0]+ru[i][2]/12.0;
		double ru_l1 = ru[i][1];
		double re_l0 = re[i][0]+re[i][2]/12.0;
		double re_l1 = re[i][1];

		ro_l1 = minmod_(ro_l1, LIM_ALPHA*(ro[i][0]-ro[i-1][0])/HX, LIM_ALPHA*(ro[i][0]-ro[i][0])/HX);
		ru_l1 = minmod_(ru_l1, LIM_ALPHA*(ru[i][0]-ru[i-1][0])/HX, LIM_ALPHA*(ru[i][0]-ru[i][0])/HX);
		re_l1 = minmod_(re_l1, LIM_ALPHA*(re[i][0]-re[i-1][0])/HX, LIM_ALPHA*(re[i][0]-re[i][0])/HX);

		if (ro_l1 != ro[i][1])
		{
			ro[i][0] = ro_l0;
			ro[i][1] = ro_l1;
			ro[i][2] = 0.0;
		}
		if (ru_l1 != ru[i][1])
		{
			ru[i][0] = ru_l0;
			ru[i][1] = ru_l1;
			ru[i][2] = 0.0;
		}
		if (re_l1 != re[i][1])
		{
			re[i][0] = re_l0;
			re[i][1] = re_l1;
			re[i][2] = 0.0;
		}
	}

}

static double L[3][3];
static double R[3][3];

void calc_limiters_weno()
{
    double H = 1.0;
    double c2_, u_;
    double AGAM = GAM-1.;


    double ***fld_arr = new double**[3];
    fld_arr[0] = ro_;
    fld_arr[1] = ru_;
    fld_arr[2] = re_;

    double ***fld_arr_1 = new double**[3];
    fld_arr_1[0] = ro_1;
    fld_arr_1[1] = ru_1;
    fld_arr_1[2] = re_1;

    double psi = 1.e-3;
    double gam[3] = {psi, 1.-psi-psi, psi};
    double beta[3], omega[3];

    for (int i = 0; i < NX; i++) {
        int ind[3] = {i-1, i, i+1};
        if (i == 0) ind[0] = 0;
        if (i == NX-1) ind[2] = NX-1;

        int& i0 = ind[0];
        int& i1 = ind[1];
        int& i2 = ind[2];

        u_ = AVG(ru,i1)/AVG(ro,i1);
        c2_ = (AVG(re,i1) / AVG(ro,i1) - u_*u_*0.5)*GAM*AGAM;
        calc_matr_L(c2_, u_, GAM, L);
        calc_matr_R(c2_, u_, GAM, R);

        for (int i_ = 0; i_ < 3; i_++) {
            int ii = ind[i_];
            for (int j = 0; j < NF; j++) {
                ro_[ii][j] = L[0][0] * ro[ii][j] + L[0][1] * ru[ii][j] + L[0][2] * re[ii][j];
                ru_[ii][j] = L[1][0] * ro[ii][j] + L[1][1] * ru[ii][j] + L[1][2] * re[ii][j];
                re_[ii][j] = L[2][0] * ro[ii][j] + L[2][1] * ru[ii][j] + L[2][2] * re[ii][j];
            }
        }

        memcpy(ro_1[i1], ro_[i1], sizeof(double)*NF);
        memcpy(ru_1[i1], ru_[i1], sizeof(double)*NF);
        memcpy(re_1[i1], re_[i1], sizeof(double)*NF);

        for (int f = 0; f < 3; f++) {
            double **fld = fld_arr[f];
            double **fld_1 = fld_arr_1[f];

            double xp12 = xc[i1]+0.5*HX;
            double xm12 = xc[i1]-0.5*HX;
            double uj = AVG(fld,i1);
            double uj1 = get_field(fld, i1, xp12)-uj;
            double uj2 = uj-get_field(fld, i1, xm12);
            double dp_uj = get_field_avg(fld, i2)-uj;
            double dm_uj = uj-get_field_avg(fld, i0);
            double uj1_mod = minmod_(uj1, dp_uj, dm_uj);
            double uj2_mod = minmod_(uj2, dp_uj, dm_uj);
            if ((fabs(uj1_mod-uj1) > EPS_LIM) || (fabs(uj2_mod-uj2) > EPS_LIM)) {
                double s = 0.0;
                for (int j = 0; j < 3; j++) {
                    int jj = ind[j];
                    beta[j] = fld[jj][1] * fld[jj][1] + (5.0 * fld[jj][2] * fld[jj][2]);
                    omega[j] = gam[j]/_SQR_(1.e-6 + beta[j]);
                    s += omega[j];
                }
                for (int j = 0; j < 3; j++)
                    omega[j] /= s;
                double p0 = AVG(fld,i0);
                double p1 = AVG(fld,i1);
                double p2 = AVG(fld,i2);
                fld_1[i1][0] =
                        omega[0] * (fld[i0][0]-p0+p1) +
                        omega[1] * (fld[i1][0]) +
                        omega[2] * (fld[i2][0]-p2+p1);
                fld_1[i1][1] =
                        omega[0] * (fld[i0][1]) +
                        omega[1] * (fld[i1][1]) +
                        omega[2] * (fld[i2][1]);
                fld_1[i1][2] =
                        omega[0] * (fld[i0][2]) +
                        omega[1] * (fld[i1][2]) +
                        omega[2] * (fld[i2][2]);
            }
        }
        for (int j = 0; j < NF; j++) {
            ro_2[i1][j] = R[0][0] * ro_1[i1][j] + R[0][1] * ru_1[i1][j] + R[0][2] * re_1[i1][j];
            ru_2[i1][j] = R[1][0] * ro_1[i1][j] + R[1][1] * ru_1[i1][j] + R[1][2] * re_1[i1][j];
            re_2[i1][j] = R[2][0] * ro_1[i1][j] + R[2][1] * ru_1[i1][j] + R[2][2] * re_1[i1][j];
        }
    }

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NF; j++) {
            ro[i][j] = ro_2[i][j];
            ru[i][j] = ru_2[i][j];
            re[i][j] = re_2[i][j];
        }
    }

    delete[] fld_arr;
    delete[] fld_arr_1;

} //calcLimiterWENO()
