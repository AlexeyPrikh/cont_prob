//
// Created by Alexey Prikhod'ko on 17.02.2021.
//

#include "Special_functions.h"

double get_MSE(vec& a, vec& b, long long n) {
    double sum = 0;
    for (long long i = 0; i < n; i++)sum += pow(a[i]-b[i], 2);
    return sqrt(sum) / n;
}

double get_max(vec& arr, long long n) {
    double max = 0;
    for (long long i = 0; i < n; i++) {
        if (abs(arr[i]) > max) max = abs(arr[i]);
    }
    return max;
}

double get_max(vec& arr1, vec& arr2, long long n) {
    double max = 0;
    for (long long i = 0; i < n; i++) {
        if (abs(arr1[i] - arr2[i]) > max) max = abs(arr1[i]-arr2[i]);
    }
    return max;
}

vec getZeros(int N) {
    vec res;
    return res;
}

vec progonka(int N, vec& A, vec& B, vec& C, vec& F) {
    vec alph, beta, x;
    alph = getZeros(N);
    beta = getZeros(N);
    x = getZeros(N);

    alph[1] = -B[0] / C[0];
    beta[1] = F[0] / C[0];

    for (int i = 1; i < N - 1; i++) {
        alph[i + 1] = -B[i] / (A[i] * alph[i] + C[i]);
        beta[i + 1] = (F[i] - A[i] * beta[i]) / (A[i] * alph[i] + C[i]);
    }
    x[N - 1] = (F[N - 1] - A[N - 1] * beta[N - 1]) / (C[N - 1] + A[N - 1] * alph[N - 1]);
    int s = 0;
    for (int i = 2; i < N + 1; i++) {
        s = N - i;
        x[s] = alph[s + 1] * x[s + 1] + beta[s + 1];
    }
    return x;
}


void check_progonka() {
    vec A(4,0);
    vec B(4,0);
    vec C(4,0);
    vec F(4,0);
    A[0] = 0; A[1] = 1; A[2] = 1; A[3] = 1;
    C[0] = 2; C[1] = 10; C[2] = -5; C[3] = 4;
    B[0] = 1; B[1] = -5; B[2] = 2; B[3] = 0;
    F[0] = -5; F[1] = -18; F[2] = -40; F[3] = -27;

    cout << "Solution is: -3, 1, 5 ,-8";

    auto x = progonka(4, A, B, C, F);
    for (int i = 0; i < 4; i++)cout << x[i] << " ";
    cout << endl;
}


void savetofile(string name, vec& array, int Nt, int Nx) {
    ofstream myfile;
    myfile.open("Result/" + name + ".txt");
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            myfile << array[n*Nx + j] << " ";
        }
        myfile << endl;
    }
    myfile.close();
}

void savetofile(string& name, vec& array, int Nt, int Nx, int Ny,int k) {
    ofstream myfile;
    myfile.open("Result/"+name+".txt");
    for (int n = 0; n < Nt; n++) {
        for (int j = 0; j < Nx; j++) {
            myfile << array[n*Nx*Ny + j*Ny + k] << " ";
        }
        myfile << endl;
    }
    myfile.close();
}


/*
void check_daxpy() {
    vec  q = new double[3];
    q[0] = 10; q[1] = 20; q[2] = 30;

    vec  q1 = new double[3];
    q1[0] = 1; q1[1] = 1; q1[2] = 1;

    double alpha = 0.5;

    cblas_daxpy(3, -alpha * 10, q1, 1, q, 1);

    cout << q[0] << ", " << q[1] << ", " << q[2] << endl;
}
*/
/*
vec progonka_mkl(int N, vec A, vec B, vec C, vec F) {
    lapack_int nrhs = 1;
    int ierr;
    ierr = LAPACKE_dgtsv(LAPACK_ROW_MAJOR, N, nrhs,
                         A+1, C, B, F, nrhs);
    return F;
}

void progonka_mkl_2(int N, vec A, vec B, vec C, vec F) {
    int nrhs = 1;
    int ldb = N;
    int info;
    ddtsvb(&N, &nrhs, A + 1, C, B, F, &ldb, &info);
}
*/

/*
void progonka_pardiso(int N, vec A, vec B, vec C, vec F) {
	int nrhs = 1;
	int ldb = N;
	int info;
	ddtsvb(&N, &nrhs, A + 1, C, B, F, &ldb, &info);
}
*/

/*

 PARDISO();

 lapack_int LAPACKE_dgbsv( int matrix_layout, lapack_int n, lapack_int kl,
 lapack_int ku, lapack_int nrhs, vec ab,
 lapack_int ldab, lapack_int* ipiv, vec b,
 lapack_int ldb )

 lapack_int LAPACKE_dgtsv( int matrix_layout, lapack_int n, lapack_int nrhs,
 <datatype>* dl, <datatype>* d, <datatype>* du, <datatype>* b, lapack_int ldb );
 */
