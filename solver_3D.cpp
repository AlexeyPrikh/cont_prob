#include "solver_3D.h"


solver_3D::solver_3D(Solution solut, double t_F, int nt,
                     double a, int nx, double b, int ny,
                     double c, int nz, double noise) {
    flag_simulated = true;

    sol = solut;
    T_b = t_F; Nt =nt ; tau = t_F / (Nt - 1);
    A_b = a; Nx = nx; hx = A_b / (Nx - 1);
    B_b = b; Ny = ny; hy = B_b / (Ny - 1);
    h_b = c; Nz = nz; hz = h_b / (Nz - 1);
    q_size = Nt*Nx*Ny;

    sigma = 1;
    f_arr = get_f(noise);
}

solver_3D::solver_3D(double t_F, double nt, double a, double nx,
          double b,double ny,double c,double nz,
          double alph, double lambd, double qV, double Tinf,
          double rh, double Cp,
          vec T1) {
    flag_simulated = false;

    T_b = t_F; Nt =nt ; tau = t_F / (Nt - 1);
    A_b = a; Nx = nx; hx = A_b / (Nx - 1);
    B_b = b; Ny = ny; hy = B_b / (Ny - 1);
    h_b = c; Nz = nz; hz = h_b / (Nz - 1);
    q_size = Nt*Nx*Ny;
    alpha = alph; lambda = lambd;
    q_V = qV; T_inf = Tinf;
    rho = rh; C_p = Cp;
    sigma = 1;
    T_1_arr = T1;
}


solver_3D::~solver_3D() {}


vec solver_3D::get_solution() {
    vec U;
    long long full_size = ((long long)Nt*Nx*Ny*Nz);
    U = getZeros(full_size);
// #pragma omp parallel for
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                for (int s = 0; s < Nz; s++) {
                    U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = sol.u_(n*tau,j*hx, k*hy, s*hz);
                }
            }
        }
    }
    return U;
}

//init data
double solver_3D::T_0(double x, double y, double z) {
    if (flag_simulated) {
        return sol.u_(0, x, y, z);
    } else {
        int j = (int)x/hx;
        int k = (int)y/hy;
        int s = (int)z/hz;
        return T_1_arr[Nx*Ny*0+Ny*j+k];
    }
}

//boundary data

// x
double solver_3D::g_0(double t, double y, double z) {
    if (flag_simulated) {
        return sol.u_x(t, 0, y, z);
    }else {
        return 0;
    }
}

double solver_3D::g_1(double t, double y, double z) {
    if (flag_simulated) {
        return sol.u_x(t, A_b, y, z);
    }else {
        return 0;
    }
}

// y
double solver_3D::h_0(double t, double x, double z) {
    if (flag_simulated) {
        return sol.u_y(t, x, 0, z);
    }else {
        return 0;
    }
}

double solver_3D::h_1(double t, double x, double z) {
    if (flag_simulated) {
        return sol.u_y(t, x, B_b, z);
    }else {
        return 0;
    }
}

// z
//double solver_3D::T_1(double t, double x, double y) {
//    if (flag_simulated) {
//        return sol.u_(t,x,y,0);
//    }else {
//        int n = (int)t/tau;
//        int j = (int)x/hx;
//        int k = (int)y/hy;
//        return T_1_arr[Nx*Ny*n+Ny*j+k];
//    }
//}

double solver_3D::T_a(double t, double x, double y) {
    if (flag_simulated) {
        return sol.u_(t,x,y,0);
    }else {
        int j = (int)x/hx;
        int k = (int)y/hy;
        if (t<=1) {
            int n = (int)t/0.1;
            return T_1_arr[Nx*Ny*n+Ny*j+k];
        } else {
            return T_1_arr[Nx*Ny*(Nt-1)+Ny*j+k];
        }
    }
}

//data    u_z |z=0
double solver_3D::f_(double t, double x, double y) {
    if (flag_simulated) {
        //return sol.u_z(t, x, y, 0);
        int n = (int)t/tau;
        int j = (int)x/hx;
        int k = (int)y/hy;
        return f_arr[Nx*Ny*n+Ny*j+k];
    } else {
        return alpha/lambda*(T_a(t,x,y)-T_inf);
    }
}

//data    u_z |z=h
//double solver_3D::q_(double t, double x, double y) {
//        //return sol.u_z(t, x, y, h_b);
//        int n = (int)t/tau;
//        int j = (int)x/hx;
//        int k = (int)y/hy;
//        return q_arr[Nx*Ny*n+Ny*j+k];
//}
//coef
double solver_3D::a_(double x, double y, double z) {
    if (flag_simulated){
        return 1.;
    } else {
        return 1. / (rho * Cp);
    }

}

double solver_3D::G_source(double t, double x, double y, double z) {
    if (flag_simulated) {
        return sol.u_t(t, x, y, z) - a_(x, y, z) * sol.laplace_u(t, x, y, z);
    } else {
        return q_V;
    }
}

//solver
vec solver_3D::solve_3D(vec& q_n) {
    clock_t t0 = clock();

    vec U, U_tmp, U_tmp2;
    long long full_size = ((long long)Nt*Nx*Ny*Nz);
    U = getZeros(full_size);
    U_tmp = getZeros(Nx*Ny*Nz);
    U_tmp2 = getZeros(Nx*Ny*Nz);

// #pragma omp parallel for
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int s = 0; s < Nz; s++) {
                U[0 * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = T_0(j*hx, k*hy, s*hz);
            }
        }
    }
    double beta, rx, ry, rz;
    beta = 1 - sigma;
    rx = tau / (hx*hx);
    ry = tau / (hy*hy);
    rz = tau / (hz*hz);
    clock_t t1;


//	int n_regul = Nt/8;
    for (int n = 0; n < Nt - 1; n++) {
        if (n == Nt - 2)t1 = clock();

        ////////////////////// 1st sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double rxn, f1;
            vec Ax, Bx, Cx, Fx;

            Ax = getZeros(Nx);
            Bx = getZeros(Nx);
            Cx = getZeros(Nx);
            Fx = getZeros(Nx);

            for (int k = 0; k < Ny; k++) {
                for (int j = 1; j < Nx - 1; j++) {
                    rxn = a_(j*hx, k*hy, s*hz)*rx;

                    Ax[j] = -sigma * rxn;
                    Cx[j] = 1 + 2 * sigma*rxn;
                    Bx[j] = -sigma * rxn;
                    f1 = U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + beta * rxn
                                                                  *(U[n*Nx*Ny*Nz + (j - 1)*Ny*Nz + k * Nz + s] - 2 * U[n*Nx*Ny*Nz + (j)*Ny*Nz + k * Nz + s] + U[n*Nx*Ny*Nz + (j + 1)*Ny*Nz + k * Nz + s]);
                    Fx[j] = f1 + tau * G_source(n*tau, j*hx, k*hy, s*hz);
                }
                Cx[0] = c0 - c2 / Bx[1] * Ax[1];
                Bx[0] = c1 - c2 / Bx[1] * Cx[1];
                Fx[0] = g_0((n + 1 / 3)*tau, k*hy, s*hz)*hx - c2 / Bx[1] * Fx[1];

                Ax[Nx - 1] = b1 - b2 / Ax[Nx - 2] * Cx[Nx - 2];
                Cx[Nx - 1] = b0 - b2 / Ax[Nx - 2] * Bx[Nx - 2];
                Fx[Nx - 1] = g_1((n + 1 / 3)*tau, k*hy, s*hz)*hx - b2 / Ax[Nx - 2] * Fx[Nx - 2];

                Fx = progonka(Nx, Ax, Bx, Cx, Fx);

                for (int j = 0; j < Nx; j++)U_tmp[j*Ny*Nz + k * Nz + s] = Fx[j];

            }
        }

        /////////////	2nd sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double ryn;

            vec Ay, By, Cy, Fy;

            Ay = getZeros(Ny);
            By = getZeros(Ny);
            Cy = getZeros(Ny);
            Fy = getZeros(Ny);

            for (int j = 0; j < Nx; j++) {
                for (int k = 1; k < Ny - 1; k++) {
                    ryn = a_(j*hx, k*hy, s*hz)*ry;

                    Ay[k] = -sigma * ryn;
                    Cy[k] = 1 + 2 * sigma*ryn;
                    By[k] = -sigma * ryn;
                    Fy[k] = U_tmp[j*Ny*Nz + k * Nz + s] + beta * ryn*(U_tmp[j*Ny*Nz + (k - 1)*Nz + s] - 2 * U_tmp[j*Ny*Nz + (k)*Nz + s] + U_tmp[j*Ny*Nz + (k + 1)*Nz + s]);
                }
                Cy[0] = c0 - c2 / By[1] * Ay[1];
                By[0] = c1 - c2 / By[1] * Cy[1];
                Fy[0] = h_0((n + 2 / 3)*tau, j*hx, s*hz)*hy - c2 / By[1] * Fy[1];

                Ay[Ny - 1] = b1 - b2 / Ay[Ny - 2] * Cy[Ny - 2];
                Cy[Ny - 1] = b0 - b2 / Ay[Ny - 2] * By[Ny - 2];
                Fy[Ny - 1] = h_1((n + 2 / 3)*tau, j*hx, s*hz)*hy - b2 / Ay[Ny - 2] * Fy[Ny - 2];

                Fy = progonka(Ny, Ay, By, Cy, Fy);

                for (int k = 0; k < Ny; k++)U_tmp2[j*Ny*Nz + k * Nz + s] = Fy[k];
            }
        }

        //////////////// 3rd sublayer
// #pragma omp parallel for
        for (int k = 0; k < Ny; k++) {
            double rzn;
            
            vec Az, Bz, Cz, Fz;

            Az = getZeros(Nz);
            Bz = getZeros(Nz);
            Cz = getZeros(Nz);
            Fz = getZeros(Nz);

            for (int j = 0; j < Nx; j++) {
                for (int s = 1; s < Nz - 1; s++) {
                    rzn = a_(j*hx, k*hy, s*hz)*rz;

                    Az[s] = -sigma * rzn;
                    Cz[s] = 1 + 2 * sigma*rzn;
                    Bz[s] = -sigma * rzn;
                    Fz[s] = U_tmp2[j*Ny*Nz + k * Nz + s] + beta * rzn*(U_tmp2[j*Ny*Nz + k * Nz + s - 1] - 2 * U_tmp2[j*Ny*Nz + k * Nz + s] + U_tmp2[j*Ny*Nz + k * Nz + s + 1]);
                }
                Cz[0] = 1.;
                Bz[0] = 0.;
                Fz[0] = T_a((n + 1)*tau, j*hx, k*hy);

                Az[Nz - 1] = b1 - b2 / Az[Nz - 2] * Cz[Nz - 2];
                Cz[Nz - 1] = b0 - b2 / Az[Nz - 2] * Bz[Nz - 2];
                Fz[Nz - 1] = q_n[(n + 1)*Nx*Ny + j*Ny + k]*hz - b2 / Az[Nz - 2] * Fz[Nz - 2]; // CHECK!!!

                Fz = progonka(Nz, Az, Bz, Cz, Fz);
                for (int s = 0; s < Nz; s++)U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = Fz[s];
            }
        }
    }

    //cout << "Calculation time: " << ((float)(clock() - t0)) / CLOCKS_PER_SEC << endl;
    //cout << "Calculation of one time layer: " << ((float)(clock() - t1)) / CLOCKS_PER_SEC << endl;

    return U;
}

//solver
vec solver_3D::solve_3D_conj(vec& U_z, vec& Ut) {
    clock_t t0 = clock();

    vec psi, psi_tmp, psi_tmp2;
    long long full_size = ((long long)Nt*Nx*Ny*Nz);
    psi = getZeros(full_size);
    psi_tmp = getZeros(Nx*Ny*Nz);
    psi_tmp2 = getZeros(Nx*Ny*Nz);


    //initial data psi=0| t=t_F
// #pragma omp parallel for
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int s = 0; s < Nz; s++) {
                psi[(Nt - 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = 0;
                //regularization
                //psi[(Nt - 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = Ut[j * Ny*Nz + k * Nz + s];
            }
        }
    }
    double beta, rx, ry, rz;
    beta = 1 - sigma;
    rx = tau / (hx*hx);
    ry = tau / (hy*hy);
    rz = tau / (hz*hz);
    clock_t t1;

#define APSI_DERIV
    for (int n = Nt - 1; n > 0; n--) {
        if (n == 1)t1 = clock();
        ////////////////////// 1st sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double rxn;
            vec Ax, Bx, Cx, Fx;

            Ax = getZeros(Nx);
            Bx = getZeros(Nx);
            Cx = getZeros(Nx);
            Fx = getZeros(Nx);

            for (int k = 0; k < Ny; k++) {
                for (int j = 1; j < Nx - 1; j++) {

                    rxn = a_(j*hx, k*hy, s*hz)*rx;

                    Ax[j] = -sigma * rxn;
                    Cx[j] = 1 + 2 * sigma*rxn;
                    Bx[j] = -sigma * rxn;
                    Fx[j] = psi[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + beta * rxn*(psi[n*Nx*Ny*Nz + (j - 1)*Ny*Nz + k * Nz + s] - 2 * psi[n*Nx*Ny*Nz + (j)*Ny*Nz + k * Nz + s] + psi[n*Nx*Ny*Nz + (j + 1)*Ny*Nz + k * Nz + s]);
                }

#ifdef APSI_DERIV
                Cx[0] = c0 - c2 / Bx[1] * Ax[1];
                Bx[0] = c1 - c2 / Bx[1] * Cx[1];
                Fx[0] = 0 - c2 / Bx[1] * Fx[1];

                Ax[Nx - 1] = b1 - b2 / Ax[Nx - 2] * Cx[Nx - 2];
                Cx[Nx - 1] = b0 - b2 / Ax[Nx - 2] * Bx[Nx - 2];
                Fx[Nx - 1] = 0 - b2 / Ax[Nx - 2] * Fx[Nx - 2];
#else
                Cx[0] = 1;
				Bx[0] = 0;
				Fx[0] = 0;

				Ax[Nx - 1] = 0;
				Cx[Nx - 1] = 1;
				Fx[Nx - 1] = 0;
#endif
                Fx = progonka(Nx, Ax, Bx, Cx, Fx);

                for (int j = 0; j < Nx; j++)psi_tmp[j*Ny*Nz + k * Nz + s] = Fx[j];

            }
        }

        ///////////// 2nd sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double ryn;
            vec Ay, By, Cy, Fy;

            Ay = getZeros(Ny);
            By = getZeros(Ny);
            Cy = getZeros(Ny);
            Fy = getZeros(Ny);

            for (int j = 0; j < Nx; j++) {
                for (int k = 1; k < Ny - 1; k++) {
                    ryn = a_(j*hx, k*hy, s*hz)*ry;

                    Ay[k] = -sigma * ryn;
                    Cy[k] = 1 + 2 * sigma*ryn;
                    By[k] = -sigma * ryn;
                    Fy[k] = psi_tmp[j*Ny*Nz + k * Nz + s] + beta * ryn*(psi_tmp[j*Ny*Nz + (k - 1)*Nz + s] - 2 * psi_tmp[j*Ny*Nz + (k)*Nz + s] + psi_tmp[j*Ny*Nz + (k + 1)*Nz + s]);
                }

#ifdef APSI_DERIV
                Cy[0] = c0 - c2 / By[1] * Ay[1];
                By[0] = c1 - c2 / By[1] * Cy[1];
                Fy[0] = 0 - c2 / By[1] * Fy[1];

                Ay[Ny - 1] = b1 - b2 / Ay[Ny - 2] * Cy[Ny - 2];
                Cy[Ny - 1] = b0 - b2 / Ay[Ny - 2] * By[Ny - 2];
                Fy[Ny - 1] = 0 - b2 / Ay[Ny - 2] * Fy[Ny - 2];
#else
                Cy[0] = 1;
				By[0] = 0;
				Fy[0] = 0;

				Ay[Ny - 1] = 0;
				Cy[Ny - 1] = 1;
				Fy[Ny - 1] = 0;
#endif
                Fy = progonka(Ny, Ay, By, Cy, Fy);

                for (int k = 0; k < Ny; k++)psi_tmp2[j*Ny*Nz + k * Nz + s] = Fy[k];
            }
        }

        //////////////// 3rd sublayer
// #pragma omp parallel for
        for (int k = 0; k < Ny; k++) {
            double rzn;
            vec Az, Bz, Cz, Fz;

            Az = getZeros(Nz);
            Bz = getZeros(Nz);
            Cz = getZeros(Nz);
            Fz = getZeros(Nz);

            for (int j = 0; j < Nx; j++) {
                for (int s = 1; s < Nz - 1; s++) {
                    rzn = a_(j*hx, k*hy, s*hz)*rz;


                    Az[s] = -sigma * rzn;
                    Cz[s] = 1 + 2 * sigma*rzn;
                    Bz[s] = -sigma * rzn;
                    Fz[s] = psi_tmp2[j*Ny*Nz + k * Nz + s] + beta * rzn
                                                             *(psi_tmp2[j*Ny*Nz + k * Nz + s - 1] - 2 * psi_tmp2[j*Ny*Nz + k * Nz + s] + psi_tmp2[j*Ny*Nz + k * Nz + s + 1]);
                }

                Cz[0] = 1.;
                Bz[0] = 0.;


                Fz[0] = 2 * (U_z[(n-1) * Nx*Ny + j * Ny + k] - f_((n - 1)*tau,j*hx, k*hy)); // CHECK
                //Fz[0] = sin((n - 1)*tau)*cos(j*hx)*cos(k*hy)*cos(0);


                Az[Nz - 1] = b1 - b2 / Az[Nz - 2] * Cz[Nz - 2];
                Cz[Nz - 1] = b0 - b2 / Az[Nz - 2] * Bz[Nz - 2];
                Fz[Nz - 1] = 0 - b2 / Az[Nz - 2] * Fz[Nz - 2];

                Fz = progonka(Nz, Az, Bz, Cz, Fz);
                for (int s = 0; s < Nz; s++)psi[(n - 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = Fz[s];
            }
        }
    }

    //cout << "Calculation time: " << ((float)(clock() - t0)) / CLOCKS_PER_SEC << endl;
    //cout << "Calculation of one time layer: " << ((float)(clock() - t1)) / CLOCKS_PER_SEC << endl;

    return psi;
}

vec solver_3D::solve_3D_inverse_explicit_t()
{
    clock_t t0 = clock();

    vec U, U_tmp, U_tmp2;
    long long full_size = ((long long)Nt*Nx*Ny*Nz);
    U = getZeros(full_size);
    U_tmp = getZeros(Nx*Ny*Nz);
    U_tmp2 = getZeros(Nx*Ny*Nz);

// #pragma omp parallel for
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int s = 0; s < Nz; s++) {
                U[0 * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = T_0(j*hx, k*hy, s*hz);
            }
        }
    }
    double beta, rx, ry, rz;
    beta = 1 - sigma;
    rx = tau / (hx*hx);
    ry = tau / (hy*hy);
    rz = tau / (hz*hz);
    clock_t t1;


    //	int n_regul = Nt/8;
    for (int n = 0; n < Nt - 1; n++) {
        if (n == Nt - 2)t1 = clock();

        //REGULIZATION
        /*
        if (n >= Nt-n_regul) {
// #pragma omp parallel for
            for (int j = 0; j < Nx; j++) {
                for (int k = 0; k < Ny; k++) {
                    for (int s = 0; s < Nz; s++) {
                        U[(n + 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = U[(Nt - n_regul)* Nx*Ny*Nz + j * Ny*Nz + k * Nz + s];
                    }
                }
            }
            continue;
        }
        */

        ////////////////////// 1st sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double rxn, f1;
            vec Ax, Bx, Cx, Fx;

            Ax = getZeros(Nx);
            Bx = getZeros(Nx);
            Cx = getZeros(Nx);
            Fx = getZeros(Nx);

            for (int k = 0; k < Ny; k++) {
                for (int j = 1; j < Nx - 1; j++) {
                    rxn = a_(j*hx, k*hy, s*hz)*rx;

                    Ax[j] = -sigma * rxn;
                    Cx[j] = 1 + 2 * sigma*rxn;
                    Bx[j] = -sigma * rxn;
                    f1 = U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + beta * rxn
                                                                  *(U[n*Nx*Ny*Nz + (j - 1)*Ny*Nz + k * Nz + s] - 2 * U[n*Nx*Ny*Nz + (j)*Ny*Nz + k * Nz + s] + U[n*Nx*Ny*Nz + (j + 1)*Ny*Nz + k * Nz + s]);
                    Fx[j] = f1 + tau * G_source(n*tau, j*hx, k*hy, s*hz);
                }
                Cx[0] = c0 - c2 / Bx[1] * Ax[1];
                Bx[0] = c1 - c2 / Bx[1] * Cx[1];
                Fx[0] = g_0((n + 1 / 3)*tau, k*hy, s*hz)*hx - c2 / Bx[1] * Fx[1];

                Ax[Nx - 1] = b1 - b2 / Ax[Nx - 2] * Cx[Nx - 2];
                Cx[Nx - 1] = b0 - b2 / Ax[Nx - 2] * Bx[Nx - 2];
                Fx[Nx - 1] = g_1((n + 1 / 3)*tau, k*hy, s*hz)*hx - b2 / Ax[Nx - 2] * Fx[Nx - 2];

                Fx = progonka(Nx, Ax, Bx, Cx, Fx);

                for (int j = 0; j < Nx; j++)U_tmp[j*Ny*Nz + k * Nz + s] = Fx[j];

            }
        }

        /////////////	2nd sublayer
// #pragma omp parallel for
        for (int s = 0; s < Nz; s++) {
            double ryn;
            vec Ay, By, Cy, Fy;

            Ay = getZeros(Ny);
            By = getZeros(Ny);
            Cy = getZeros(Ny);
            Fy = getZeros(Ny);

            for (int j = 0; j < Nx; j++) {
                for (int k = 1; k < Ny - 1; k++) {
                    ryn = a_(j*hx, k*hy, s*hz)*ry;

                    Ay[k] = -sigma * ryn;
                    Cy[k] = 1 + 2 * sigma*ryn;
                    By[k] = -sigma * ryn;
                    Fy[k] = U_tmp[j*Ny*Nz + k * Nz + s] + beta * ryn*(U_tmp[j*Ny*Nz + (k - 1)*Nz + s] - 2 * U_tmp[j*Ny*Nz + (k)*Nz + s] + U_tmp[j*Ny*Nz + (k + 1)*Nz + s]);
                }
                Cy[0] = c0 - c2 / By[1] * Ay[1];
                By[0] = c1 - c2 / By[1] * Cy[1];
                Fy[0] = h_0((n + 2 / 3)*tau, j*hx, s*hz)*hy - c2 / By[1] * Fy[1];

                Ay[Ny - 1] = b1 - b2 / Ay[Ny - 2] * Cy[Ny - 2];
                Cy[Ny - 1] = b0 - b2 / Ay[Ny - 2] * By[Ny - 2];
                Fy[Ny - 1] = h_1((n + 2 / 3)*tau, j*hx, s*hz)*hy - b2 / Ay[Ny - 2] * Fy[Ny - 2];

                Fy = progonka(Ny, Ay, By, Cy, Fy);

                for (int k = 0; k < Ny; k++)U_tmp2[j*Ny*Nz + k * Nz + s] = Fy[k];
            }
        }

        //////////////// 3rd sublayer
// #pragma omp parallel for
        for (int k = 0; k < Ny; k++) {
            double rzn;
            for (int j = 0; j < Nx; j++) {
                U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 0] = T_a((n+1)*tau,j*hx,k*hy);
                rzn = a_(j*hx, k*hy, 1*hz)*rz;
                double G= U_tmp2[j*Ny*Nz + k * Nz + 2] + beta * rzn*(U_tmp2[j*Ny*Nz + k * Nz + 2 - 1] - 2 * U_tmp2[j*Ny*Nz + k * Nz + 2] + U_tmp2[j*Ny*Nz + k * Nz + 2 + 1]);
                U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 1] = ( hz*f_((n+1)*tau, j*hx, k*hy) - c2*hz*hz*G - U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 0]*(c0-c2))/ ( c1+2*c2 );
                for (int s = 2; s < Nz; s++) {
                    rzn = a_(j*hx, k*hy, s*hz)*rz;
                    G = U_tmp2[j*Ny*Nz + k * Nz + s] + beta * rzn*(U_tmp2[j*Ny*Nz + k * Nz + s - 1] - 2 * U_tmp2[j*Ny*Nz + k * Nz + s] + U_tmp2[j*Ny*Nz + k * Nz + s + 1]);
                    U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = hz * hz*G - U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s - 2] + 2 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s - 1];
                }
            }
        }
    }


    //cout << "Calculation time: " << ((float)(clock() - t0)) / CLOCKS_PER_SEC << endl;
    //cout << "Calculation of one time layer: " << ((float)(clock() - t1)) / CLOCKS_PER_SEC << endl;

    return U;
}


vec solver_3D::solve_3D_inverse_explicit_z()
{
    cout << "Difference circuit inversion method" << endl;
    cout << "Kurant condition: " << hz << " <= " << a_(0,0,0)*pow(1. / tau / tau + 1. / hx / hx + 1. / hy / hy, -1) << endl;
    clock_t t0 = clock();

    vec U;
    long long full_size = ((long long)Nt*Nx*Ny*Nz);
    U = getZeros(full_size);

    // U|z=0
    // Uz|z=0
// #pragma omp parallel for
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int n = 0; n < Nt; n++) {
                U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + 0] = T_a(n*tau,j*hx, k*hy);
                U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + 1] = T_a(n*tau, j*hx, k*hy) +
                        hz * f_(n*tau, j*hx, k*hy);
            }
        }
    }

    //clock_t t1;

    /*
    int s = 1;
    double g0;
    a = a_(j*hx, k*hy, s*hz);
    ut = (U[(n + 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] - U[(n - 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s]) / (2 * tau);
    uxx = (U[n * Nx*Ny*Nz + (j + 1) * Ny*Nz + k * Nz + s] - 2 * U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + U[n * Nx*Ny*Nz + (j - 1) * Ny*Nz + k * Nz + s]) / (hx*hx);
    uyy = (U[n * Nx*Ny*Nz + j * Ny*Nz + (k + 1) * Nz + s] - 2 * U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + U[n * Nx*Ny*Nz + j * Ny*Nz + (k - 1) * Nz + s]) / (hy*hy);
    G = G_source(n*tau, j*hx, k*hy, s*hz);
    g0 = ut / a - uxx - uyy - G / a;
    U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + 1] = (U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + 0]*(1-c0/c2) - (hz/c2) * f_data[n * Nx*Ny + j * Ny + k] + (hz*hz)*g0)/(2+c1/c2);
    */


    for (int s = 1; s < Nz-1; s++) {
// #pragma omp parallel for
        for (int n = 1; n < Nt-1; n++) {/////////////////////
            double ut, uxx, uyy, G, a;
            for (int j = 1; j < Nx-1; j++) {
                for (int k = 1; k < Ny-1; k++) {
                    a = a_(j*hx, k*hy, s*hz);
                    ut = (U[(n+1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s]- U[(n-1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s])/(2*tau);
                    //ut = (U[(n) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] - U[(n - 1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s]) / (tau);

                    uxx = (U[n * Nx*Ny*Nz + (j+1) * Ny*Nz + k * Nz + s]- 2*U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s]+ U[n * Nx*Ny*Nz + (j-1) * Ny*Nz + k * Nz + s])/(hx*hx);
                    uyy = (U[n * Nx*Ny*Nz + j * Ny*Nz + (k+1) * Nz + s]- 2*U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s]+ U[n * Nx*Ny*Nz + j * Ny*Nz + (k-1) * Nz + s])/(hy*hy);
                    G = G_source(n*tau,j*hx,k*hy,s*hz);
                    U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s+1] = 2 * U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] - U[n * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s - 1] + hz * hz *(ut/a - uxx - uyy - G/a);
                }
            }
        }
        //use data on x=0 x = A
        int s1 = s + 1;
// #pragma omp parallel for
        for (int n = 1; n < Nt-1; n++) {
            for (int k = 1; k < Ny-1; k++) {
                U[n * Nx*Ny*Nz + 0 * Ny*Nz + k * Nz + s1] = (g_0(n*tau, k*hy, s1*hz)*hx - c1 * U[n * Nx*Ny*Nz + 1 * Ny*Nz + k * Nz + s1] - c2 * U[n * Nx*Ny*Nz + 2 * Ny*Nz + k * Nz + s1]) / c0;
                U[n * Nx*Ny*Nz + (Nx-1) * Ny*Nz + k * Nz + s1] = (g_1(n*tau, k*hy, s1*hz)*hx - b1 * U[n * Nx*Ny*Nz + (Nx-2) * Ny*Nz + k * Nz + s1] - b2 * U[n * Nx*Ny*Nz + (Nx-3) * Ny*Nz + k * Nz + s1]) / b0;
            }
        }
        //use data in y=0, y =B
// #pragma omp parallel for
        for (int n = 1; n < Nt-1; n++) {
            for (int j = 0; j < Nx; j++) {
                U[n * Nx*Ny*Nz + j * Ny*Nz + 0 * Nz + s1] = (h_0(n*tau, j*hx, s1*hz)*hy - c1 * U[n * Nx*Ny*Nz + j * Ny*Nz + 1 * Nz + s1] - c2 * U[n * Nx*Ny*Nz + j * Ny*Nz + 2 * Nz + s1]) / c0;
                U[n * Nx*Ny*Nz + j * Ny*Nz + (Ny-1) * Nz + s1] = (h_1(n*tau, j*hx, s1*hz)*hy - b1 * U[n * Nx*Ny*Nz + j * Ny*Nz + (Ny-2) * Nz + s1] - b2 * U[n * Nx*Ny*Nz + j * Ny*Nz + (Ny-3) * Nz + s1]) / b0;
            }
        }

        //use data in t=0; t=T
// #pragma omp parallel for
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                U[0 * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s1] = T_0(j*hx, k*Ny, s1*hz);
                U[(Nt-1) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s1] = (0*tau - b1 * U[(Nt-2) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s1] - b2 * U[(Nt-3) * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s1]) / b0;
            }
        }

    }

    //cout << "Calculation time: " << ((float)(clock() - t0)) / CLOCKS_PER_SEC << endl;
    //cout << "Calculation of one time layer: " << ((float)(clock() - t1)) / CLOCKS_PER_SEC << endl;

    return U;
}

double solver_3D::integrate_omega_z(vec& arg) {
    long long n = Nt * Nx*Ny;
    double sum = 0;
    //#pragma omp parallel for
    for (long long i = 0; i < n; i++)sum += pow(arg[i], 2);
    return sum*tau*hx*hy;
}

double solver_3D::get_norm(vec& arr) {
    long long n = Nt * Nx*Ny*Nz;
    double sum = 0;
    //#pragma omp parallel for
    for (long long i = 0; i < n; i++)sum += pow(arr[i], 2);
    return sqrt(sum*tau*hx*hy*hz);
}

double solver_3D::get_J(vec& U_z) {
    double J = 0;
//#pragma omp parallel for
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                J += pow((U_z[n * Nx*Ny + j * Ny + k] - f_(n*tau, j*hx, k*hy)), 2)* tau*hx*hy;
            }
        }
    }
    return J;

}

//solver
vec solver_3D::get_grad(vec& a_psi) {
    vec grad;
    grad = getZeros(q_size);

//#pragma omp parallel for
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                grad[n * Nx*Ny + j * Ny + k] = a_psi[n * Nx*Ny*Nz + j * Ny*Nz + k*Nz + (Nz-1)];
            }
        }
    }
    return grad;
}
//
//vec solver_3D::get_q() {
//    vec q_arr;
//    q_arr = getZeros(get_q_size());
//// #pragma omp parallel for
//    for (int n = 0; n < Nt; n++) {
//        for (int j = 0; j < Nx; j++) {
//            for (int k = 0; k < Ny; k++) {
//                q_arr[n * Nx*Ny + j * Ny + k] = q_(n*tau, j*hx, k*hy);
//            }
//        }
//    }
//    return q_arr;
//}

vec solver_3D::get_U_x_y(vec& U, int j, int k) {
    vec Uxy; Uxy = getZeros(Nt*Nz);
    for (int n = 0; n < Nt; n++) {
        for (int s = 0; s < Nz; s++) {
            Uxy[n*Nz+s] = U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s];
        }
    }
    return Uxy;
}



vec solver_3D::get_Uz_0(vec& U) {
    vec U_z;U_z= getZeros(q_size);
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                U_z[(n)*Nx*Ny + j * Ny + k] = (c0*U[(n )*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 0] + c1 * U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 1] + c2 * U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + 2]) / hz;
            }
        }
    }
    return U_z;
}

vec solver_3D::get_Uz_h(vec& U)
{
    vec U_z; U_z = getZeros(q_size);
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                U_z[(n)*Nx*Ny + j * Ny + k] = (b0*U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + Nz-1] + b1 * U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + Nz-2] + b2 * U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + Nz-3]) / hz;
            }
        }
    }
    return U_z;
}

vec solver_3D::get_Uh(vec& U) {
    vec Uh; Uh = getZeros(q_size);
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                Uh[(n)*Nx*Ny + j * Ny + k] = U[(n)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (Nz-1)];
            }
        }
    }
    return Uh;
}

vec solver_3D::get_Ut_T(vec& U) {
    vec Ut; Ut = getZeros(Nx*Ny*Nz);
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int s = 0; s < Nz; s++) {
                Ut[j*Ny*Nz + k * Nz + s] = (b0*U[(Nt-1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + b1*U[(Nt -2)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + b2*U[(Nt-3)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s])/tau;
            }
        }
    }
    return Ut;
}

vec solver_3D::get_f(double noise) {

    const double mean = 0.0;
    const double stddev = 1.;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    vec f;
    f = getZeros(q_size);
// #pragma omp parallel for
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                f[(n)*Nx*Ny + j * Ny + k] = f_(n*tau, j*hx, k*hy)+ noise*dist(generator);
            }
        }
    }
    return f;
}

vec solver_3D::grad_method(vec& q0, int num_of_iters, double alph){

    //double *q_sol; q_sol = get_q();
    vec q;

    vec U,a_psi, A_grad_J;
    vec grad_J, U_z;
    double norm_grad, norm_A_grad;
    double alpha_n=0;
    double J_val = 0;
    double J_val_prev = 0;
    cout << "Gradient method started"<<endl;
    //cout << "for n=0: || q_sol - q_n || : "<<get_MSE(q_sol,q,q_size)<< endl;

    vec Ut;
    alpha_n = alph;
    int k_point=2;

    string sss = "";
    //double epsilon = pow(10, -4);
    //double coef = 0.9;
    for (int it = 0; it < num_of_iters; it++) {
        J_val_prev = J_val;
        U = solve_3D(q);
        U_z = get_Uz_0(U);

        Ut = get_Ut_T(U);

        a_psi = solve_3D_conj(U_z,  Ut);

        grad_J = get_grad(a_psi);
        norm_grad = integrate_omega_z(grad_J);

        A_grad_J = solve_3D(grad_J);
        norm_A_grad = integrate_omega_z(A_grad_J);

        //alpha_n = 0.5 * norm_grad/norm_A_grad;

        for (int i=0; i<  q_size; i++)
            q[i] -=alpha_n*grad_J[i];

        J_val = get_J(U_z);
        //if (J_val > (J_val_prev))alpha_n *= coef;
        //else if(J_val < (J_val_prev - epsilon)) alpha_n /= coef;
        cout << "it : " << it << "; J(qn) : " << J_val << "; alph_n : " << alpha_n << endl;
        //cout << "||J'(qn)|| : " << norm_grad << "; ||q_sol-q_n|| : " << get_MSE(q_sol, q, q_size) << endl;
        //cout << "alpha_n : "<< alpha_n << "; max(q_n) : " << get_max(q, q_size) << endl;
        if ((it % 5 == 0)&& (it/5 <5)){
            savetofile("q_n/q_" + to_string(it/5), q, Nt, Nx, Ny, k_point);
            //savetofile("q_n/q_" + sss, q, Nt, Nx, Ny, k_point);
            //sss += 'a';
        }
    }

    cout <<"Gradient method ended"<<endl;

    return q;
}

int solver_3D::get_Nt() {
    return Nt;
}
int solver_3D::get_Nx() {
    return Nx;
}
int solver_3D::get_Ny() {
    return Ny;
}
int solver_3D::get_Nz() {
    return Nz;
}

double solver_3D::get_tau() {
    return tau;
}
double solver_3D::get_hx() {
    return hx;
}
double solver_3D::get_hy() {
    return hy;
}
double solver_3D::get_hz() {
    return hz;
}

long long solver_3D::get_q_size()
{
    return q_size;
}

//void solver_3D::set_q(vec &q) {
//    q_arr = q;
//}


// EXPlicit	and pardiso
/*



double* solver_3D::solve_3D_pardiso() {
	clock_t t0 = clock();

	cout << "Solving forward 3D problem ... " << endl;

	cout << "tau: " << tau << ", hx: " << hx << ", hy: " << hy << ", hz: " << hz << endl;
	double *U, *U_tmp, *U_tmp2;
	//U = (double *)mkl_malloc((Nt*Nx*Ny*Nz) * sizeof(double), 64);
	U = (double *)malloc((Nt*Nx*Ny*Nz) * sizeof(double));
	U_tmp = (double *)mkl_malloc((Nx*Ny*Nz) * sizeof(double), 64);
	U_tmp2 = (double *)mkl_malloc((Nx*Ny*Nz) * sizeof(double), 64);

	// memory
   //omp_set_dynamic(0);
   //omp_set_num_threads(4);
	int nthr = omp_get_max_threads();
	cout << "Threads : " << nthr << endl;
	cout << "Memory : " << (Nt*Nx*Ny*Nz + 2 * Nx*Ny*Nz + nthr * 4 * (Nx + Ny + Nz)) * 8.0 / pow(10, 9) << endl;

	//#define MEMORY_INSIDE
		// outside Nt=30 -> 86.277 / layer -> 2.946 ->0.00065041
		// inside Nt=30 -> 84.993 / layer -> 2.892
#ifndef MEMORY_INSIDE
	double *Ax_gl, *Bx_gl, *Cx_gl, *Fx_gl;
	double *Ay_gl, *By_gl, *Cy_gl, *Fy_gl;
	double *Az_gl, *Bz_gl, *Cz_gl, *Fz_gl;

	Ax_gl = (double *)mkl_malloc((Nx*nthr) * sizeof(double), 64);
	Bx_gl = (double *)mkl_malloc((Nx*nthr) * sizeof(double), 64);
	Cx_gl = (double *)mkl_malloc(Nx*nthr * sizeof(double), 64);
	Fx_gl = (double *)mkl_malloc(Nx*nthr * sizeof(double), 64);

	Ay_gl = (double *)mkl_malloc((Ny*nthr) * sizeof(double), 64);
	By_gl = (double *)mkl_malloc((Ny*nthr) * sizeof(double), 64);
	Cy_gl = (double *)mkl_malloc(Ny*nthr * sizeof(double), 64);
	Fy_gl = (double *)mkl_malloc(Ny*nthr * sizeof(double), 64);

	Az_gl = (double *)mkl_malloc((Nz*nthr) * sizeof(double), 64);
	Bz_gl = (double *)mkl_malloc((Nz*nthr) * sizeof(double), 64);
	Cz_gl = (double *)mkl_malloc(Nz*nthr * sizeof(double), 64);
	Fz_gl = (double *)mkl_malloc(Nz*nthr * sizeof(double), 64);
#endif
	double c0 = -1.5, c1 = 2, c2 = -0.5;
	double b0 = -c0, b1 = -c1, b2 = -c2;
// #pragma omp parallel for
	for (int j = 0; j < Nx; j++) {
		for (int k = 0; k < Ny; k++) {
			for (int s = 0; s < Nz; s++) {
				U[0 * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = T_0(j*hx, k*hy, s*hz);
			}
		}

	}
	double beta, rx, ry, rz;
	beta = 1 - sigma;
	rx = tau / (hx*hx);
	ry = tau / (hy*hy);
	rz = tau / (hz*hz);
	clock_t t1;

	for (int n = 0; n < Nt - 1; n++) {
		if (n == Nt - 2)t1 = clock();
		//////////////////////
// #pragma omp parallel for //num_threads(nthr) //private(Ax_gl,Bx_gl,Cx_gl,Fx_gl)
		for (int s = 0; s < Nz; s++) {
			double rxn1, rxn, f1;
			double *Ax, *Bx, *Cx, *Fx;

#ifdef MEMORY_INSIDE
			Ax = (double *)mkl_malloc((Nx) * sizeof(double), 64);
			Bx = (double *)mkl_malloc((Nx) * sizeof(double), 64);
			Cx = (double *)mkl_malloc(Nx * sizeof(double), 64);
			Fx = (double *)mkl_malloc(Nx * sizeof(double), 64);
#else
			int tid = omp_get_thread_num();
			Ax = Ax_gl + tid * Nx;
			Bx = Bx_gl + tid * Nx;
			Cx = Cx_gl + tid * Nx;
			Fx = Fx_gl + tid * Nx;

#endif
			for (int k = 0; k < Ny; k++) {
				for (int j = 1; j < Nx - 1; j++) {
					rxn1 = a_((n + 1 / 3)*tau, j*hx, k*hy, s*hz)*rx;
					rxn = a_(n*tau, j*hx, k*hy, s*hz)*rx;

					Ax[j] = -sigma * rxn1;
					Cx[j] = 1 + 2 * sigma*rxn1;
					Bx[j] = -sigma * rxn1;
					f1 = U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + beta * rxn*(U[n*Nx*Ny*Nz + (j - 1)*Ny*Nz + k * Nz + s] - 2 * U[n*Nx*Ny*Nz + (j)*Ny*Nz + k * Nz + s] + U[n*Nx*Ny*Nz + (j + 1)*Ny*Nz + k * Nz + s]);
					Fx[j] = f1 + tau * G_source(n*tau, j*hx, k*hy, s*hz);
				}
				Cx[0] = c0 - c2 / Bx[1] * Ax[1];
				Bx[0] = c1 - c2 / Bx[1] * Cx[1];
				Fx[0] = g_0((n + 1 / 3)*tau, k*hy, s*hz)*hx - c2 / Bx[1] * Fx[1];

				Ax[Nx - 1] = b1 - b2 / Ax[Nx - 2] * Cx[Nx - 2];
				Cx[Nx - 1] = b0 - b2 / Ax[Nx - 2] * Bx[Nx - 2];
				Fx[Nx - 1] = g_1((n + 1 / 3)*tau, k*hy, s*hz)*hx - b2 / Ax[Nx - 2] * Fx[Nx - 2];

				progonka(Nx, Ax, Bx, Cx, Fx);

				//pardiso(pt,maxfct,mnum,mtype,phase,n,a,ia,ja,perm,nrhs,iparm,msglvl,b,x,error);

				for (int j = 0; j < Nx; j++)U_tmp[j*Ny*Nz + k * Nz + s] = Fx[j];

			}
#ifdef MEMORY_INSIDE
			mkl_free(Ax); mkl_free(Bx); mkl_free(Cx); mkl_free(Fx);
#endif
		}
		/////////////

// #pragma omp parallel for //num_threads(nthr) //private(Ay_gl,By_gl,Cy_gl,Fy_gl)
		for (int s = 0; s < Nz; s++) {
			double ryn1, ryn2;
			double *Ay, *By, *Cy, *Fy;
#ifdef MEMORY_INSIDE
			Ay = (double *)mkl_malloc((Ny) * sizeof(double), 64);
			By = (double *)mkl_malloc((Ny) * sizeof(double), 64);
			Cy = (double *)mkl_malloc(Ny * sizeof(double), 64);
			Fy = (double *)mkl_malloc(Ny * sizeof(double), 64);
#else
			int tid = omp_get_thread_num();
			Ay = Ay_gl + tid * Ny;
			By = By_gl + tid * Ny;
			Cy = Cy_gl + tid * Ny;
			Fy = Fy_gl + tid * Ny;
#endif
			for (int j = 0; j < Nx; j++) {
				for (int k = 1; k < Ny - 1; k++) {
					ryn1 = a_((n + 1. / 3)*tau, j*hx, k*hy, s*hz)*ry;
					ryn2 = a_((n + 2. / 3)*tau, j*hx, k*hy, s*hz)*ry;

					Ay[k] = -sigma * ryn2;
					Cy[k] = 1 + 2 * sigma*ryn2;
					By[k] = -sigma * ryn2;
					Fy[k] = U_tmp[j*Ny*Nz + k * Nz + s] + beta * ryn1*(U_tmp[j*Ny*Nz + (k - 1)*Nz + s] - 2 * U_tmp[j*Ny*Nz + (k)*Nz + s] + U_tmp[j*Ny*Nz + (k + 1)*Nz + s]);
				}
				Cy[0] = c0 - c2 / By[1] * Ay[1];
				By[0] = c1 - c2 / By[1] * Cy[1];
				Fy[0] = h_0((n + 2 / 3)*tau, j*hx, s*hz)*hy - c2 / By[1] * Fy[1];

				Ay[Ny - 1] = b1 - b2 / Ay[Ny - 2] * Cy[Ny - 2];
				Cy[Ny - 1] = b0 - b2 / Ay[Ny - 2] * By[Ny - 2];
				Fy[Ny - 1] = h_1((n + 2 / 3)*tau, j*hx, s*hz)*hy - b2 / Ay[Ny - 2] * Fy[Ny - 2];

				progonka(Ny, Ay, By, Cy, Fy);

				for (int k = 0; k < Ny; k++)U_tmp2[j*Ny*Nz + k * Nz + s] = Fy[k];
			}
#ifdef MEMORY_INSIDE
			mkl_free(Ay); mkl_free(By); mkl_free(Cy); mkl_free(Fy);
#endif
		}

		////////////////
// #pragma omp parallel for //num_threads(nthr) //private(Az_gl,Bz_gl,Cz_gl,Fz_gl)
		for (int k = 0; k < Ny; k++) {
			double rzn2, rzn3;
			double *Az, *Bz, *Cz, *Fz;
#ifdef MEMORY_INSIDE
			Az = (double *)mkl_malloc((Nz) * sizeof(double), 64);
			Bz = (double *)mkl_malloc((Nz) * sizeof(double), 64);
			Cz = (double *)mkl_malloc(Nz * sizeof(double), 64);
			Fz = (double *)mkl_malloc(Nz * sizeof(double), 64);
#else
			int tid = omp_get_thread_num();
			Az = Az_gl + tid * Nz;
			Bz = Bz_gl + tid * Nz;
			Cz = Cz_gl + tid * Nz;
			Fz = Fz_gl + tid * Nz;
#endif
			for (int j = 0; j < Nx; j++) {
				for (int s = 1; s < Nz - 1; s++) {
					rzn2 = a_((n + 2 / 3)*tau, j*hx, k*hy, s*hz)*rz;
					rzn3 = a_((n + 1)*tau, j*hx, k*hy, s*hz)*rz;

					Az[s] = -sigma * rzn3;
					Cz[s] = 1 + 2 * sigma*rzn3;
					Bz[s] = -sigma * rzn3;
					Fz[s] = U_tmp2[j*Ny*Nz + k * Nz + s] + beta * rzn2*(U_tmp2[j*Ny*Nz + k * Nz + s - 1] - 2 * U_tmp2[j*Ny*Nz + k * Nz + s] + U_tmp2[j*Ny*Nz + k * Nz + s + 1]);
				}
				Cz[0] = 1.;
				Bz[0] = 0.;
				Fz[0] = T_1((n + 1)*tau, j*hx, k*hy);

				Az[Nz - 1] = b1 - b2 / Az[Nz - 2] * Cz[Nz - 2];
				Cz[Nz - 1] = b0 - b2 / Az[Nz - 2] * Bz[Nz - 2];
				Fz[Nz - 1] = q_((n + 1)*tau, j*hx, k*hy)*hz - b2 / Az[Nz - 2] * Fz[Nz - 2];

				progonka(Nz, Az, Bz, Cz, Fz);
				for (int s = 0; s < Nz; s++)U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = Fz[s];
			}
#ifdef MEMORY_INSIDE
			mkl_free(Az); mkl_free(Bz); mkl_free(Cz); mkl_free(Fz);
#endif
		}
	}
#ifndef MEMORY_INSIDE
	mkl_free(Ax_gl); mkl_free(Bx_gl); mkl_free(Cx_gl); mkl_free(Fx_gl);
	mkl_free(Ay_gl); mkl_free(By_gl); mkl_free(Cy_gl); mkl_free(Fy_gl);
	mkl_free(Az_gl); mkl_free(Bz_gl); mkl_free(Cz_gl); mkl_free(Fz_gl);
#endif
	mkl_free(U_tmp); mkl_free(U_tmp2);

	cout << "Calculation time: " << ((float)(clock() - t0)) / CLOCKS_PER_SEC << endl;
	cout << "Calculation of one time layer: " << ((float)(clock() - t1)) / CLOCKS_PER_SEC << endl;

	return U;
}



double* solver_3D::solve_3D_explicit() {
	clock_t t;
	//t = clock();

	cout << "tau: " << tau << ", hx: " << hx << endl;
	double *U;
	U = getZeros(Nt*Nx*Ny*Nz);

	double c0 = -1.5, c1 = 2, c2 = -0.5;
	double b0 = -c0, b1 = -c1, b2 = -c2;
// #pragma omp parallel for
	for (int j = 0; j < Nx; j++) {
		for (int k = 0; k < Ny; k++) {
			for (int s = 0; s < Nz; s++) {
				U[0 * Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = T_0(j*hx, k*hy, s*hz);
			}
		}

	}
	double rx, ry, rz;
	rx = 1. / (hx*hx);
	ry = 1. / (hy*hy);
	rz = 1. / (hz*hz);

	t = clock();
	for (int n = 0; n < Nt - 1; n++) {
		//////////////////////
		for (int s = 1; s < Nz - 1; s++) {
			for (int k = 1; k < Ny - 1; k++) {
				for (int j = 1; j < Nx - 1; j++) {
					U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] = U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + s] + tau * a_(n*tau, j*hx, k*hy, s*hz)* (rx*(U[n*Nx*Ny*Nz + (j - 1)*Ny*Nz + k * Nz + s] - 2 * U[n*Nx*Ny*Nz + (j)*Ny*Nz + k * Nz + s] + U[n*Nx*Ny*Nz + (j + 1)*Ny*Nz + k * Nz + s]) + ry * (U[n*Nx*Ny*Nz + j * Ny*Nz + (k - 1)*Nz + s] - 2 * U[n*Nx*Ny*Nz + j * Ny*Nz + (k)*Nz + s] + U[n*Nx*Ny*Nz + j * Ny*Nz + (k + 1)*Nz + s]) + rz * (U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (s - 1)] - 2 * U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (s)] + U[n*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (s + 1)]) + G_source(n*tau, j*hx, k*hy, s*hz));
				}
			}
		}
		// use boundary conditions
		for (int j = 0; j < Nx; j++) {
			for (int k = 0; k < Ny; k++) {
				// z=0
				U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (0)] = T_1((n + 1)*tau, j*hx, k*hy);
				// z=h
				U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (Nz - 1)] = (hz*q_((n + 1)*tau, j*hx, k*hy) - b1 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (Nz - 2)] - b2 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + k * Nz + (Nz - 3)]) / b0;
			}
		}


		for (int k = 0; k < Ny; k++) {
			for (int s = 1; s < Nz - 1; s++) {
				// x=0
				U[(n + 1)*Nx*Ny*Nz + (0)*Ny*Nz + k * Nz + s] = (hx*g_0((n + 1)*tau, k*hy, s*hz) - c1 * U[(n + 1)*Nx*Ny*Nz + (1)*Ny*Nz + k * Nz + s] - c2 * U[(n + 1)*Nx*Ny*Nz + (2)*Ny*Nz + k * Nz + s]) / c0;
				// x=A
				U[(n + 1)*Nx*Ny*Nz + (Nx - 1)*Ny*Nz + k * Nz + s] = (hx*g_1((n + 1)*tau, k*hy, s*hz) - b1 * U[(n + 1)*Nx*Ny*Nz + (Nx - 2)*Ny*Nz + k * Nz + s] - b2 * U[(n + 1)*Nx*Ny*Nz + (Nx - 3)*Ny*Nz + k * Nz + s]) / b0;
			}
		}
		for (int j = 1; j < Nx - 1; j++) {
			for (int s = 1; s < Nz - 1; s++) {
				// y=0
				U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (0)*Nz + s] = (hy*h_0((n + 1)*tau, j*hx, s*hz) - c1 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (1)*Nz + s] - c2 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (2)*Nz + s]) / c0;
				// y=B
				U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (Ny - 1)*Nz + s] = (hy*h_1((n + 1)*tau, j*hx, s*hz) - b1 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (Ny - 2)*Nz + s] - b2 * U[(n + 1)*Nx*Ny*Nz + j * Ny*Nz + (Ny - 3)*Nz + s]) / b0;
			}
		}

	}

	t = clock() - t;
	cout << "Calculation time: " << ((float)t) / CLOCKS_PER_SEC << endl;

	return U;
}

*/

