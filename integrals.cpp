#include "functions.h"
#include "iostream"

using namespace std;

void rim_orig(double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
              double  RB,             double  PB, double  UB, double  VB, double  WB,
              double  RE,             double  PE, double  UE, double  VE, double  WE, double GAM);


void calc_integral_rhs() {
    double ro_[NGP], ru_[NGP], re_[NGP], u_[NGP], p_[NGP], fr_[NGP], fu_[NGP], fe_[NGP];
    for (int i = 0; i < NX; i++) {
		memset(ro_int[i], 0, sizeof(double)*NF);
		memset(ru_int[i], 0, sizeof(double)*NF);
		memset(re_int[i], 0, sizeof(double)*NF);
	}

	// интегралы
	for (int i = 0; i < NX; i++) {

        for (int j = 0; j < NGP; j++) {
            ro_[j] = get_field(xgp[i][j], i, ID_RO);
            ru_[j] = get_field(xgp[i][j], i, ID_RU);
            re_[j] = get_field(xgp[i][j], i, ID_RE);
            u_[j] = ru_[j]/ro_[j];
            p_[j] = ro_[j]*(GAM-1.0)*(re_[j]/ro_[j]-u_[j]*u_[j]*0.5);

            fr_[j] = ru_[j];
            fu_[j] = fr_[j]*u_[j]+p_[j];
            fe_[j] = (re_[j]+p_[j])*u_[j];
        }
        for (int j = 0; j < NF; j++) {
			double int1 = 0.0, int2 = 0.0, int3 = 0.0;

            for (int jg = 0; jg < NGP; jg++) {
                int1 += fr_[jg] * get_df(i, j, xgp[i][jg]) * wgp[jg];
                int2 += fu_[jg] * get_df(i, j, xgp[i][jg]) * wgp[jg];
                int3 += fe_[jg] * get_df(i, j, xgp[i][jg]) * wgp[jg];
            }

			int1 *= jgp[i];
			int2 *= jgp[i];
			int3 *= jgp[i];

			ro_int[i][j] -= int1;
			ru_int[i][j] -= int2;
			re_int[i][j] -= int3;
		}
	}

	{	int i = 0;
		double ror = get_field(x1[i],i, ID_RO);
		double rur = get_field(x1[i],i, ID_RU);
		double rer = get_field(x1[i],i, ID_RE);
		double ur = rur/ror;
		double pr = ror*(GAM-1.0)*(rer/ror-ur*ur*0.5);
		double cr = sqrt(pr*GAM/ror);

		double rol =  ror;
		double rul = rur;
		double rel =  rer;
		double ul = rul/rol;
		double pl = rol*(GAM-1.0)*(rel/rol-ul*ul*0.5);
		double cl = sqrt(pl*GAM/rol);

        double ri, ei, pi, ui, vi, wi;
        rim_orig(&ri, &ei, &pi, &ui, &vi, &wi, rol, pl, ul, 0., 0., ror, pr, ur, 0., 0., GAM);
        double fr = ri*ui;
        double fu = fr*ui+pi;
        double fe = (ri*(ei+ui*ui*0.5)+pi)*ui;

//		double alpha = max(fabs(ur)+cr, fabs(ul)+cl);
//		double fr = 0.5*(rur            + rul          -  alpha*(ror-rol));
//		double fu = 0.5*(rur*ur+pr      + rul*ul+pl    -  alpha*(rur-rul));
//		double fe = 0.5*((rer+pr)*ur    + (rel+pl)*ul  -  alpha*(rer-rel));

		for (int j = 0; j < NF; j++) {
			ro_int[i][j] -= fr*get_f(i, j, x1[i]);
			ru_int[i][j] -= fu*get_f(i, j, x1[i]);
			re_int[i][j] -= fe*get_f(i, j, x1[i]);
		}
	}

	for (int i = 1; i < NX; i++) {
		double rol = get_field(x2[i-1],i-1, ID_RO);
		double rul = get_field(x2[i-1],i-1, ID_RU);
		double rel = get_field(x2[i-1],i-1, ID_RE);
		double ul = rul/rol;
		double pl = rol*(GAM-1.0)*(rel/rol-ul*ul*0.5);
		double cl = sqrt(pl*GAM/rol);

		double ror = get_field(x1[i],i, ID_RO);
		double rur = get_field(x1[i],i, ID_RU);
		double rer = get_field(x1[i],i, ID_RE);
		double ur = rur/ror;
		double pr = ror*(GAM-1.0)*(rer/ror-ur*ur*0.5);
		double cr = sqrt(pr*GAM/ror);

        double ri, ei, pi, ui, vi, wi;
        rim_orig(&ri, &ei, &pi, &ui, &vi, &wi, rol, pl, ul, 0., 0., ror, pr, ur, 0., 0., GAM);
        double fr = ri*ui;
        double fu = fr*ui+pi;
        double fe = (ri*(ei+ui*ui*0.5)+pi)*ui;

//		double alpha = max(fabs(ur)+cr, fabs(ul)+cl);
//		double fr = 0.5*(rur            + rul          -  alpha*(ror-rol));
//		double fu = 0.5*(rur*ur+pr      + rul*ul+pl    -  alpha*(rur-rul));
//		double fe = 0.5*((rer+pr)*ur    + (rel+pl)*ul  -  alpha*(rer-rel));

		for (int j = 0; j < NF; j++) {
			ro_int[i-1][j] += fr*get_f(i-1, j, x2[i-1]);
			ru_int[i-1][j] += fu*get_f(i-1, j, x2[i-1]);
			re_int[i-1][j] += fe*get_f(i-1, j, x2[i-1]);

			ro_int[i][j]   -= fr*get_f(  i, j, x1[i]  );
			ru_int[i][j]   -= fu*get_f(  i, j, x1[i]  );
			re_int[i][j]   -= fe*get_f(  i, j, x1[i]  );
		}
	}

	{	int i = NX;
		double rol = get_field(x2[i-1],i-1, ID_RO);
		double rul = get_field(x2[i-1],i-1, ID_RU);
		double rel = get_field(x2[i-1],i-1, ID_RE);
		double ul = rul/rol;
		double pl = rol*(GAM-1.0)*(rel/rol-ul*ul*0.5);
		double cl = sqrt(pl*GAM/rol);

		double ror =  rol;
		double rur = rul;
		double rer =  rel;
		double ur = rur/ror;
		double pr = ror*(GAM-1.0)*(rer/ror-ur*ur*0.5);
		double cr = sqrt(pr*GAM/ror);

        double ri, ei, pi, ui, vi, wi;
        rim_orig(&ri, &ei, &pi, &ui, &vi, &wi, rol, pl, ul, 0., 0., ror, pr, ur, 0., 0., GAM);
        double fr = ri*ui;
        double fu = fr*ui+pi;
        double fe = (ri*(ei+ui*ui*0.5)+pi)*ui;

//		double alpha = max(fabs(ur)+cr, fabs(ul)+cl);
//		double fr = 0.5*(rur            + rul          -  alpha*(ror-rol));
//		double fu = 0.5*(rur*ur+pr      + rul*ul+pl    -  alpha*(rur-rul));
//		double fe = 0.5*((rer+pr)*ur    + (rel+pl)*ul  -  alpha*(rer-rel));

		for (int j = 0; j < NF; j++) {
			ro_int[i-1][j] += fr*get_f(i-1, j, x2[i-1]);
			ru_int[i-1][j] += fu*get_f(i-1, j, x2[i-1]);
			re_int[i-1][j] += fe*get_f(i-1, j, x2[i-1]);
		}
	}
}

void calc_new() {
	for (int i = 0; i < NX; i++) {
		for (int j = 0; j < NF; j++) {
			double fr = 0.0;
			double fu = 0.0;
			double fe = 0.0;
			for (int k = 0; k < NF; k++) {
				fr += A[i][j][k]*ro_int[i][k];
				fu += A[i][j][k]*ru_int[i][k];
				fe += A[i][j][k]*re_int[i][k];
			}
			ro[i][j] -= TAU*fr;
			ru[i][j] -= TAU*fu;
			re[i][j] -= TAU*fe;
		}
	}
}


void copy_to_half() {
    for (int i = 0; i < NX; i++) {
        memcpy(ro_half[i], ro[i], sizeof(double)*NF);
        memcpy(ru_half[i], ru[i], sizeof(double)*NF);
        memcpy(re_half[i], re[i], sizeof(double)*NF);
    }
}


void calc_new_with_half() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NF; j++) {
            ro[i][j] += ro_half[i][j];
            ru[i][j] += ru_half[i][j];
            re[i][j] += re_half[i][j];

            ro[i][j] *= 0.5;
            ru[i][j] *= 0.5;
            re[i][j] *= 0.5;
        }
    }
}


void calc_new_1() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NF; j++) {
            ro[i][j] += 3.*ro_half[i][j];
            ru[i][j] += 3.*ru_half[i][j];
            re[i][j] += 3.*re_half[i][j];

            ro[i][j] *= 0.25;
            ru[i][j] *= 0.25;
            re[i][j] *= 0.25;
        }
    }
}

void calc_new_2() {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NF; j++) {
            ro[i][j] *= 2.;
            ru[i][j] *= 2.;
            re[i][j] *= 2.;

            ro[i][j] += ro_half[i][j];
            ru[i][j] += ru_half[i][j];
            re[i][j] += re_half[i][j];

            ro[i][j] /= 3.;
            ru[i][j] /= 3.;
            re[i][j] /= 3.;
        }
    }
}


void rim_orig(  double* RI, double* EI, double* PI, double* UI, double* VI, double* WI,
                double RB, double PB, double UB, double VB, double WB,
                double RE, double PE, double UE, double VE, double WE, double gam) {
    int    step;
    double AGAM = (gam - 1.0);
    double BGAM = (2.0 * sqrt(gam / AGAM));
    double CGAM = (1.0 / gam);
    double DGAM = (2.0 / AGAM);
    double EGAM = (AGAM / (gam + 1.0));
    double GGAM = (sqrt(gam * AGAM));
    double HGAM = (AGAM / 2.0);
    double FGAM = (3.0 * gam - 1.0);
    double OGAM = (AGAM / (2.0 * gam));
    double QGAM = (gam + 1.0);
    double PGAM = (QGAM / (2.0 * gam));
    double RGAM = (4.0 * gam);
    double SGAM = (gam * AGAM);
    double TGAM = (QGAM / 2.0);
    double UGAM = (sqrt(AGAM / gam));

    double RF, RS, EF, ES, SBL, SFL, SSL, SEL, D, FS1, F1, ZNB, PKB, ZFB, PKE, ZNE, F2, FS2, ZFE, DP, UBD, RUBD, UF, UED, RUED, US, PPE, PPB, P;
    double eps = 1.e-7;
    double CB = sqrt(gam * PB / RB);
    double CE = sqrt(gam * PE / RE);
    double EB = CB * CB / SGAM;
    double EE = CE * CE / SGAM;
    double RCB = RB * CB;
    double RCE = RE * CE;
    double DU = UB - UE;
    if (DU < -2.0 * (CB + CE) / AGAM) {
        printf ("%s\n", " RIEMANN PROBLEM SOLVER: ATTENTION!!!  VACUUM!!!");
        RF = 0.0;
        RS = 0.0;
        EF = 0.0;
        ES = 0.0;
        SBL = UB - CB;
        SFL = UB + 2.0 * CB / AGAM;
        SSL = UE - 2.0 * CE / AGAM;
        SEL = UE + CE;
    } else {
        P = (PB * RCE + PE * RCB + DU * RCB * RCE) / (RCB + RCE);
        step = 0;
        do {

            if (P < eps) P = eps;

            PPB = P / PB;
            if (PB <= P) {
                PKB = PGAM * PPB + OGAM;
                ZNB = RCB * sqrt(PKB);
                F1 = (P - PB) / ZNB;
                FS1 = (QGAM * PPB + FGAM) / (RGAM * ZNB * PKB);
            }
            else {
                ZFB = CB * exp(log(PPB) * OGAM);
                F1 = DGAM * (ZFB - CB);
                FS1 = ZFB / (gam * P);
            }
            PPE = P / PE;
            if (PE <= P) {
                PKE = PGAM * PPE + OGAM;
                ZNE = RCE * sqrt(PKE);
                F2 = (P - PE) / ZNE;
                FS2 = (QGAM * PPE + FGAM) / (RGAM * ZNE * PKE);
            }
            else {
                ZFE = CE * exp(log(PPE) * OGAM);
                F2 = DGAM * (ZFE - CE);
                FS2 = ZFE / (gam * P);
            }
            DP = (DU - F1 - F2) / (FS1 + FS2);
            P = P + DP;
        } while ((fabs(DU - F1 - F2) > eps) && (++step < 5000));


        PPB = P / PB;
        PPE = P / PE;

//       ZFB=CB*PPB**OGAM;
//       ZFE=CE*PPE**OGAM;
        ZFB = CB * exp(log(PPB) * OGAM);
        ZFE = CE * exp(log(PPE) * OGAM);
        if (PB <= P) {
            D = UB - sqrt((TGAM * P + HGAM * PB) / RB);
            UBD = UB - D;
            RUBD = RB * UBD;
            RF = RUBD * RUBD / (PB - P + RUBD * UBD);
            UF = D + RUBD / RF;
            EF = P / (AGAM * RF);
            SBL = D;
            SFL = D;
        }
        else {
            EF = ZFB * ZFB / SGAM;
            UF = UB + DGAM * (CB - ZFB);
            RF = P / (AGAM * EF);
            SBL = UB - CB;
            SFL = UF - ZFB;
        }
        if (PE <= P) {
            D = UE + sqrt((TGAM * P + HGAM * PE) / RE);
            UED = UE - D;
            RUED = RE * UED;
            RS = RUED * RUED / (PE - P + RUED * UED);
            US = D + RUED / RS;
            ES = P / (AGAM * RS);
            SEL = D;
            SSL = D;
        }
        else {
            ES = ZFE * ZFE / SGAM;
            US = UE - DGAM * (CE - ZFE);
            RS = P / (AGAM * ES);
            SSL = US + ZFE;
            SEL = UE + CE;
        }
    }
//
// C     compute the interpolation value
    if (SEL<=0.0) {
        *RI= RE;
        *EI= EE;
        *UI= UE;
        *VI= VE;
        *WI= WE;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SBL>=0.0) {
        *RI= RB;
        *EI= EB;
        *UI= UB;
        *VI= VB;
        *WI= WB;
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if ((SSL>=0.0)&&(SFL<=0.0)) {
        if (US>=0.0) {
            *RI= RF;
            *EI= EF;
            *UI= UF;
            *VI= VB;
            *WI= WB;
        } else {
            *RI= RS;
            *EI= ES;
            *UI= US;
            *VI= VE;
            *WI= WE;
        }
        *PI= AGAM*(*EI)*(*RI);
        return;
    }

    if (SFL>0.0) {
        *UI= (UB+DGAM*GGAM*sqrt(EB))/(1+DGAM);
        *VI= VB;
        *WI= WB;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RB*exp(log((*EI)/EB)*(1/AGAM));
    } else {
        *UI= (UE-DGAM*GGAM*sqrt(EE))/(1+DGAM);
        *VI= VE;
        *WI= WE;
        *EI= ((*UI)*(*UI))/SGAM;
        *RI= RE*exp(log((*EI)/EE)*(1/AGAM)) ;
    }

    *PI= AGAM*(*EI)*(*RI);

    return;
}
