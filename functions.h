#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//5.677714e-01 50
//5.690845e-01 100
//5.698114e-01 200
enum fid_t {ID_RO, ID_RU, ID_RE};

/*Константы*/
const auto NX = 100; // Количество ячеек
const double XMIN = 0.0;
const double XMAX = 1.0 ;
const double HX = (XMAX - XMIN) / NX; // Шаг по пространственной координате

const int NF = 3; // Число функций базиса
const int NGP = 5; // Число точек Гаусса

const double CFL = 1.e-2; // Число Куранта-Фридрихса-Леви
const double TMAX = 0.2; // Максимальное время
const double TAU = CFL * HX; // Шаг по времени

const double GAM = 7.0 / 5.0; // Показатель адиабаты

const double LIM_ALPHA = 2.0;
const double EPS_LIM = 1.e-6;
const double M_LIM = 1.e-3;

const int STEP_OUTPUT = 1000; // Интервал вывода данных

void WENO(double* UU, double& Um, double& Up);
void calc_flux_hllc(
    double rl, double ul, double pl, double gaml,
    double rr, double ur, double pr, double gamr,
    double& qr, double& qu, double& qe);

extern double t_glob;

/*Массивы*/
extern double ** ro;
extern double ** ru;
extern double ** re;

extern double ** ro_;
extern double ** ru_;
extern double ** re_;

extern double ** ro_1;
extern double ** ru_1;
extern double ** re_1;

extern double ** ro_2;
extern double ** ru_2;
extern double ** re_2;

extern double ** ro_half;
extern double ** ru_half;
extern double ** re_half;

extern double ** ro_int;
extern double ** ru_int;
extern double ** re_int;

extern double * dx;
extern double * xc;
extern double * x1;
extern double * x2;
extern double ** xgp;
extern double * wgp;
extern double * jgp;

extern double ***A;

extern int step;

void mem_allocate();
void mem_deallocate();
void init_arrays();

double get_f	(int i_cell, int i_func, double x);
double get_df	(int i_cell, int i_func, double x);
double get_field(double x, int i, fid_t field_id);
double get_field(double **fld, int i_cell, double x);
double get_field_avg(double **fld, int i_cell);

void calc_integral_rhs();
void calc_limiters_weno();
void calc_new();
void save_csv(int step);
void init();
void copy_to_half();
void calc_new_1();
void calc_new_2();

inline double max(double x, double y) {if (x > y) return x; return y;}
void calc_matr_L(double c2, double u, double GAM, double L[3][3]);
void calc_matr_R(double c2, double u, double GAM, double R[3][3]);

#define NEW1(_X_, _TYPE_, _N_) {_X_ = (_TYPE_*)malloc(sizeof(_TYPE_)*_N_);}
#define DELETE1(_X_) {free(_X_);}

#define NEW2(_X_, _TYPE_, _N1_, _N2_) {NEW1(_X_, _TYPE_*, _N1_); for(int i = 0; i < _N1_;i++) NEW1(_X_[i], _TYPE_, _N2_);}
#define DELETE2(_X_, _N_) {for (int i = 0; i < _N_; i++) {DELETE1(_X_[i])}; DELETE1(_X_);}

#define _SQR_(_X_) ((_X_)*(_X_))

#define calc_limiters() calc_limiters_weno()
//#define calc_limiters()

double calc_err_L2();

#endif
