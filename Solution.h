//
// Created by Alexey Prikhod'ko on 17.02.2021.
//

#ifndef HEAT_3D_SOLUTION_H
#define HEAT_3D_SOLUTION_H

#include "Special_functions.h"

class Solution {
public:
    Solution();
    ~Solution();

    double u_(double t, double x, double y, double z);

    double u_t(double t, double x, double y, double z);

    double u_x(double t, double x, double y, double z);

    double u_y(double t, double x, double y, double z);

    double u_z(double t, double x, double y, double z);

    double laplace_u(double t, double x, double y, double z);
};


#endif //HEAT_3D_SOLUTION_H
