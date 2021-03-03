//
// Created by Alexey Prikhod'ko on 17.02.2021.
//

#ifndef HEAT_3D_SPECIAL_FUNCTIONS_H
#define HEAT_3D_SPECIAL_FUNCTIONS_H

#endif //HEAT_3D_SPECIAL_FUNCTIONS_H


#pragma once

#ifndef Special_functions_hpp
#define Special_functions_hpp

#define _USE_MATH_DEFINES
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <random>
#include <cmath>
#include <string>

//#include <omp.h>
//#include <mkl.h>

#endif /* Special_functions_hpp */

using namespace std;

typedef vector<double> vec;
const double rho = 1.;
const double Cp = 1.;

double get_MSE(vec& a, vec& b, long long n);

double get_max(vec& arr, long long n);

double get_max(vec& arr1, vec& arr2, long long n);

vec getZeros(int N);

vec progonka(int N, vec& A, vec& B, vec& C, vec& F);

vec progonka_mkl(int N, vec& A, vec& B, vec& C, vec& F);

void progonka_mkl_2(int N, vec& A, vec& B, vec& C, vec& F);

void savetofile(string name, vec& array, int Nt, int Nx);

void savetofile(string name, vec& array, int Nt, int Nx, int Ny, int k);

void check_daxpy();

void check_progonka();