#include "solver_3D.h"

// [n,j,k,s] -> [n*Nx*Ny*Nz + j*Ny*Nz + k*Nz + s]

/*
solver_3D prepare_data() {

    double t_F = 0.1; // sec
    double A_b = 36 * pow(10, -6);//90 * pow(10, -3); //m
    double B_b = 36 * pow(10, -6);//36 * pow(10, -3);
    double h_b = 25 * pow(10, -9);//25 * pow(10, -6);

    int Nt = 11;
    int Nx = 170;//392 - 222;
    int Ny = 170;//337 - 167;
    int Nz = 40;//200;

    // Experement

    double alpha = 23.765;
    double lambda = 23;
    double q_V = 1140740;//11407.407; // W / m3

    double T_a = 30.84;

    double rho = 8900; // kg / m3
    double C_p = 410; // J / kg K

    vec T0;
    vec T1;

    T1 = getZeros(Nt*Nx*Ny);
    T0 = getZeros(Nx*Ny*Nz);

    // get T1
    double mean_T = 0;
    string name = "";
    for (int n = 0; n < Nt; n++) {
        ifstream file;
        name = "data_T/T_" + to_string(n) + ".txt";
        file.open(name);
        cout << name << endl;

        for (int j = 0; j < Nx; j++) {
            for (int k = 0; k < Ny; k++) {
                file >> T1[n*Nx*Ny + j * Ny + k];

                if (n == 0) {
                    //cout << T1[n*Nx*Ny + j * Ny + k] << " ";
                    mean_T += T1[n*Nx*Ny + j * Ny + k];
                }
            }
            //if (n == 0) cout << endl;
        }
        file.close();
    }
    mean_T /= Nx * Ny;
    cout << "T_0 = " << mean_T;


    //get T1
    for (int j = 0; j < Nx; j++) {
        for (int k = 0; k < Ny; k++) {
            for (int s = 0; s < Nz; s++) {
                T0[j*Ny*Nz + k * Nz + s] = mean_T;
            }
        }
    }

    return solver_3D(t_F, Nt, A_b, Nx,
                     B_b, Ny, h_b, Nz,
                     alpha, lambda, q_V, T_a,
                     rho, C_p,
            //T0,
                     T1);
}
*/

solver_3D prepare_data_simulated(double noise) {

    double t_F = M_PI/2; // sec
    double A_b = M_PI*2; //m
    double B_b = M_PI*2;//36 * pow(10, -3);
    double h_b = M_PI/10;//25 * pow(10, -6);

    int Nt = 10;
    int Nx = 20;//392 - 222;
    int Ny = 20;//337 - 167;
    int Nz = 40;//200;

    Solution solut;
    return solver_3D(solut, t_F, Nt, A_b, Nx, B_b, Ny,
             h_b, Nz, noise);
}

int main(int argc, const char * argv[]) {

    double noise = 0.0;
    solver_3D solver = prepare_data_simulated(noise);

    cout << "Solver created" << endl;

    // Grid
    int Nt = solver.get_Nt();
    int Nx = solver.get_Nx(); int Ny = solver.get_Ny(); int Nz = solver.get_Nz();

    cout << "T_0 = " << solver.T_0(Nx / 2, Ny / 2, Nz / 2) << endl;
    cout << "T_a = " << solver.T_a(Nt / 2, Nx / 2, Ny / 2) << endl;
    //test forward and conj
    //Check_forward_problem(solver);
    //Check_conjugated_problem(solver);



    int k_point = Ny / 2;
    int j_point = Ny / 2;

    vec q_sol = getZeros(solver.get_q_size());
    //double *U = solver.get_solution();
    //double *Uxy_sol = solver.get_U_x_y(U, k_point, j_point);

    //savetofile("Uh_sol", solver.get_Uh(U), Nt, Nx, Ny, k_point);
    //mkl_free(U);
    //double *U_fdsi = solver.solve_3D_inverse_explicit_z();

    //gradient!
    int num_of_iter = 50;
    double alph = pow(10, -2);
    //double *q_sol = solver.get_q();

    auto q_size = solver.get_q_size();
    vec q_gm, q0;

    // INIT APPROXIMATION
    q0 = getZeros(q_size);

    //double noise = 0.1;

    q_gm = solver.grad_method(q0, num_of_iter, alph);

    //double *q_fdsi = solver.get_Uz_h(U_fdsi);
    //savetofile("U_fdsi", U_fdsi, Nt*Nx*Ny*Nz, 1);


    //cout << "Save to file q_sol and q_fdsi." << endl;
    //savetofile("q_sol", q_sol, Nt, Nx, Ny, k_point);
    //savetofile("q_fdsi", q_fdsi, Nt, Nx, Ny, k_point);
    //savetofile("q_gm", q_gm, Nt, Nx, Ny, k_point);

    /*
    double *Uxy_fdsi = solver.get_U_x_y(U_fdsi, j_point, k_point);

    cout << "Save to file Uxy_sol and Uxy_fdsi." << endl;
    savetofile("Uxy_fdsi", Uxy_fdsi, Nt, Nz);
    //savetofile("Uxy_sol", Uxy_sol, Nt, Nz);

    cout << "Save to file Uh_sol and Uh_fdsi." << endl;

    savetofile("Uh_fdsi", solver.get_Uh(U_fdsi), Nt, Nx, Ny, k_point);
    */

    /*
    //gradient method
    int num_of_iter = 100;
    double alph = pow(10, 0);
    double *q_sol = solver.get_q();

    auto q_size = solver.get_q_size();
    double* q_gm,*q0;

    // INIT APPROXIMATION
    q0 = (double *)mkl_malloc((q_size) * sizeof(double), 64);
    for (int i = 0; i < q_size; i++) q0[i] = 0;

    double noise = 0.1;
    double *f_data; f_data = solver.get_f(noise);
    q_gm =solver.grad_method(q0, f_data, num_of_iter, alph);
    cout << "Mean_Sq_Err(q_sol - q_gm) " << get_MSE(q_sol, q_gm, q_size) << endl;
    cout << "Max(q_sol - q_gm) " << get_max(q_sol, q_gm, q_size) << endl;

    double *U,*U_gm, *Uz, *Uz_gm;
    U = solver.solve_3D(q_sol);
    Uz = solver.get_Uz_0(U);

    U_gm = solver.solve_3D(q_gm);
    Uz_gm = solver.get_Uz_0(U_gm);

    cout << "Save to file Uz_sol and Uz_gm." << endl;
    int k_point = 0;

    savetofile("Uz_sol", f_data, Nt, Nx, Ny, k_point);
    savetofile("Uz_gm", Uz_gm, Nt, Nx, Ny, k_point);

    cout<< "Save to file q_sol and q_gm."<<endl;
    savetofile("q_sol", q_sol, Nt, Nx, Ny, k_point);
    savetofile("q_gm", q_gm, Nt, Nx, Ny, k_point);

    double *Uh, *Uh_gm;
    Uh = solver.get_Uh(U);
    Uh_gm = solver.get_Uh(U_gm);

    cout << "Save to file Uh_sol and Uh_gm." << endl;
    savetofile("Uh_sol", Uh, Nt, Nx, Ny, k_point);
    savetofile("Uh_gm", Uh_gm, Nt, Nx, Ny, k_point);


    mkl_free(q_sol); mkl_free(q_gm); mkl_free(q0);
    */


    system("pause");
    return 0;
}



