//
// Created by Alexey Prikhod'ko on 17.02.2021.
//

#ifndef HEAT_3D_SOLVER_3D_H
#define HEAT_3D_SOLVER_3D_H


#pragma once

#include <stdio.h>

#include "Special_functions.h"
#include "Solution.h"

class solver_3D
{
public:
    solver_3D(Solution solut, double t_F, int Nt, double A, int Nx, double B, int Ny,
              double h, int Nz, double noise);

    solver_3D(double t_F, double Nt, double A_b, double Nx,
              double B_b, double Ny, double h_b,double Nz,
              double alpha, double lambda, double q_V, double Tinf,
              double rho, double C_p,
              vec T1);
    ~solver_3D();

    Solution sol;

    vec get_solution();

    //init data
    double T_0(double x, double y, double z);

    //boundary data
    // x
    double g_0(double t, double y, double z);

    double g_1(double t, double y, double z);

    // y
    double h_0(double t, double x, double z);

    double h_1(double t, double x, double z);
    // z
    //double T_1(double t, double x, double y);

    double T_a(double t, double x, double y);

    double q_(double t, double x, double y);

    // DATA
    double f_(double t, double x, double y);

//    vec get_q();

    //coef
    double a_(double x, double y, double z);

    double G_source(double t, double x, double y, double z);

    //solver
    vec solve_3D(vec& q_n);

    vec solve_3D_conj(vec& U_z, vec& Ut);

    vec solve_3D_inverse_explicit_t();

    vec solve_3D_inverse_explicit_z();

    vec get_U_x_y(vec&  U, int j, int k);

    double integrate_omega_z(vec& arg);

    double get_norm(vec& arr);

    double get_J(vec& U);

    vec get_grad(vec& a_psi);

    vec get_Uz_0(vec& U);

    vec get_Uz_h(vec&  U);

    vec get_Uh(vec& U);

    vec get_Ut_T(vec& U);

    vec get_f(double noise);


    vec grad_method(vec& q0, int num_of_iters, double alph);

    int get_Nt();
    int get_Nx();
    int get_Ny();
    int get_Nz();

    double get_tau();
    double get_hx();
    double get_hy();
    double get_hz();
    long long get_q_size();

//    void set_q(vec& q);
private:
    double T_b; int Nt; double tau;
    double A_b; int Nx; double hx;
    double B_b; int Ny; double hy;
    double h_b; int Nz; double hz;
    double sigma;

    double alpha;
    double lambda;
    double q_V, T_inf;
    double rho, C_p;
    bool flag_simulated;
    vec T_1_arr;
    vec f_arr;
    //vec q_arr;
    long long q_size;

    const double c0 = -1.5, c1 = 2, c2 = -0.5;
    const double b0 = -c0, b1 = -c1, b2 = -c2;
    //vec solve_3D_explicit();
    //vec solve_3D_pardiso();
};




#endif //HEAT_3D_SOLVER_3D_H
