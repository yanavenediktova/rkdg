#include "functions.h"

double t_glob;

double rho_0(double x) {
    if (fabs(x) < 0.5) {
        return 1;
    }
    return 0.125;
}

double p_0(double x) {
    if (fabs(x) < 0.5) {
        return 1;
    }
    return 0.1;
}

void init()
{
	double r, u, e, p;
	mem_allocate();
	init_arrays();
	for (int i = 0; i < NX; i++) {
		r = rho_0(xc[i]);
		p = p_0(xc[i]);
		e = p/(r*(GAM-1));
		u = 0;
		ro[i][0] = r;
		ru[i][0] = r*u;
		re[i][0] = r*e + r*u*u*0.5;
		for (int j = 1; j < NF; j++) {
			ro[i][j] = 0.0;
			ru[i][j] = 0.0;
			re[i][j] = 0.0;
		}
	}
}

double calc_energy(double p, double rho_0, double gamma) {
	return p / (rho_0 * (gamma - 1.0));
}

void inverse_matr(double** a_src, double **am) {
	double **a = a_src;
	double detA = a[0][0]*a[1][1]*a[2][2] + a[1][0]*a[2][1]*a[0][2] + a[0][1]*a[1][2]*a[2][0]
	            - a[2][0]*a[1][1]*a[0][2] - a[1][0]*a[0][1]*a[2][2] - a[0][0]*a[2][1]*a[1][2];
	double m[3][3];
	m[0][0] = a[1][1]*a[2][2]-a[2][1]*a[1][2];
	m[0][1] = a[2][0]*a[1][2]-a[1][0]*a[2][2];
	m[0][2] = a[1][0]*a[2][1]-a[2][0]*a[1][1];
	m[1][0] = a[2][1]*a[0][2]-a[0][1]*a[2][2];
	m[1][1] = a[0][0]*a[2][2]-a[2][0]*a[0][2];
	m[1][2] = a[2][0]*a[0][1]-a[0][0]*a[2][1];
	m[2][0] = a[0][1]*a[1][2]-a[1][1]*a[0][2];
	m[2][1] = a[1][0]*a[0][2]-a[0][0]*a[1][2];
	m[2][2] = a[0][0]*a[1][1]-a[1][0]*a[0][1];

	am[0][0] = m[0][0]/detA;
	am[0][1] = m[1][0]/detA;
	am[0][2] = m[2][0]/detA;
	am[1][0] = m[0][1]/detA;
	am[1][1] = m[1][1]/detA;
	am[1][2] = m[2][1]/detA;
	am[2][0] = m[0][2]/detA;
	am[2][1] = m[1][2]/detA;
	am[2][2] = m[2][2]/detA;
}

static void calc_matrix(double** A, int i_cell) {
	double ** tmp_a;
    double sqrt3 = 1./sqrt(3.);
    NEW2(tmp_a, double, NF, NF);

    double xx1 = 0.5*(x1[i_cell]+x2[i_cell]-(x2[i_cell]-x1[i_cell])*sqrt3);
    double xx2 = 0.5*(x1[i_cell]+x2[i_cell]+(x2[i_cell]-x1[i_cell])*sqrt3);
    double J   = 0.5*(x2[i_cell]-x1[i_cell]);
	for (int i = 0; i < NF; i++) {
		for (int j = 0; j < NF; j++) {
            tmp_a[i][j] = 0.;
            for (int igp = 0; igp < NGP; igp++) {
                tmp_a[i][j] += wgp[igp]*get_f(i_cell, i, xgp[i_cell][igp])*get_f(i_cell, j, xgp[i_cell][igp]);
            }
            tmp_a[i][j] *= jgp[i_cell];
		}
	}

	inverse_matr(tmp_a, A);
    DELETE2(tmp_a, NF);
}

void init_arrays() {
    double gp5_x[] = {0.0000000000000000, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640};
    double gp5_w[] = {0.5688888888888889,  0.4786286704993665, 0.4786286704993665,  0.2369268850561891, 0.2369268850561891};

    double hx_half = HX*0.5;
    wgp[0] = wgp[2] = 5./9.;
    wgp[1] = 8./9.;
    for (int i = 0; i < NX; i++) {
		dx[i] = HX;
		x1[i] = XMIN+i*HX;
		xc[i] = x1[i]+hx_half;
		x2[i] = x1[i]+HX;

        for (int j = 0; j < NGP; j++) {
            xgp[i][j] = 0.5*(x1[i]+x2[i]+gp5_x[j]*(x2[i]-x1[i]));
            wgp[j] = gp5_w[j];
        }

        jgp[i] = 0.5*(x2[i]-x1[i]);

        calc_matrix(A[i], i);
	}
}

double calc_rho_toch(double x, double t) {
    const double eps = 1.e-8;
    double a = -4.*sqrt(10.)/3.;
    double x0, x1, xl, xr, u0, c0, r0;
    if (x < -0.2+t*a || x > 0.2+t*a) {
        return 1.;
    }
    xl = -0.2;
    xr =  0.2;
    do {
        x0 = 0.5*(xl+xr);
        r0 = rho_0(x0);
        c0 = sqrt(pow(r0, GAM-1.)*(GAM-1.)*GAM);
        u0 = (-2.*c0)/(GAM-1.);
        x1 = x0+t*(u0-c0);
        if (x1 < x) {
            xl = x0;
        }
        else {
            xr = x0;
        }
    } while (fabs(x1-x) > eps);
    return r0;
}

void save_csv(int step) {
	char fname[100];
	sprintf(fname, "res_%010d.csv", step);
	FILE *fp = fopen(fname, "w");
	fprintf(fp, "x,Rho,U,P,E\n");
	for (int i = 0; i < NX; i++) {
		double ror = ro[i][0] + ro[i][2]/ 12.;
		double rur = ru[i][0] + ru[i][2]/ 12.;
		double rer = re[i][0] + re[i][2]/ 12.;
		double ur = rur / ror;
		double pr = ror * (GAM - 1.0) * (rer / ror - ur * ur * 0.5);
		double energy = calc_energy(pr, ror, GAM); // Вычисление энергии
		//double rho_toch = calc_rho_toch(xc[i], t_glob);
		//double rho_ep = pow(rho_toch,GAM -1);
		//double rho_u = -2*sqrt(rho_ep*(GAM-1)*GAM)/(GAM -1);
		//double rho_e = rho_toch*rho_ep+rho_toch*(pow(rho_u,2)/2);
		fprintf(fp, "%25.16e,%25.16e,%25.16e,%25.16e,%25.16e\n", xc[i], ror, ur, pr, energy);
	}
	fclose(fp);
	printf("File '%s' is saved...\n", fname);
}

void calc_matr_L(double c2, double u, double GAM, double L[3][3])
{
	double cz, kk, deta;

	cz = sqrt(c2);
	kk = GAM - 1.0;
	deta = 0.5*kk*u*u;

	L[0][0] = -c2 + deta;
	L[0][1] = -kk*u;
	L[0][2] = kk;

	L[1][0] = -cz*u + deta;
	L[1][1] = cz - kk*u;
	L[1][2] = kk;

	L[2][0] = cz*u + deta;
	L[2][1] = -cz - kk*u;
	L[2][2] = kk;

}

void calc_matr_R(double c2, double u, double GAM, double R[3][3])
{
	double cz, kk, deta;

	cz = sqrt(c2);
	kk = GAM - 1.0;
	R[0][0] = -1.0 / c2;
	R[0][1] = 0.5 / c2;
	R[0][2] = 0.5 / c2;

	R[1][0] = -u / c2;
	R[1][1] = 0.5*u / c2 + 0.5 / cz;
	R[1][2] = 0.5*u / c2 - 0.5 / cz;

	R[2][0] = -0.5*u*u / c2;
	R[2][1] = 0.25*u*u / c2 + 0.5 / kk + 0.5*u / cz;
	R[2][2] = 0.25*u*u / c2 + 0.5 / kk - 0.5*u / cz;

}

void multMatrVec3(double A[3][3], double x[3], double result[3]) {
    int i, j;
    for (i = 0; i < 3; i++) {
        result[i] = 0.;
        for (j = 0; j < 3; j++) {
            result[i] += A[i][j]*x[j];
        }
    }
}

double calc_err_L2() {
    double e = 0.;
    for (int i = 0; i < NX; i++) {
        double s = 0.;
        for (int j = 0; j < NGP; j++) {
            double r = calc_rho_toch(xgp[i][j], t_glob) - get_field(xgp[i][j], i, ID_RO);
            s += r * r * wgp[j];
        }
        e += s*jgp[i];
    }
    return sqrt(e);
}