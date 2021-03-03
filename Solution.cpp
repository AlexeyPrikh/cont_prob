//
// Created by Alexey Prikhod'ko on 17.02.2021.
//

#include "Solution.h"


Solution::Solution()
{
}


Solution::~Solution()
{
}

double Solution::u_(double t, double x, double y, double z) {
    //return pow(0.5 + pow(t - 0.5, 2), -1)*z + cos(x)*cos(y)*z;
    //return pow(0.5 + pow(t - 0.5, 2), -1) * sin(x + M_PI/6)*sin(y + M_PI / 6)*cos(z + M_PI / 6);
    return sin(t)*cos(x)*cos(y)*cos(z);
    //return sin(t)*1*cos(z);
}

double Solution::u_t(double t, double x, double y, double z) {
    //return -2*t*pow(0.5 + pow(t - 0.5, 2), -2)*z;
    //return -2*t*pow(0.5 + pow(t - 0.5, 2), -2) * sin(x + M_PI / 6)*sin(y + M_PI / 6)*cos(z + M_PI / 6);
    return cos(t)*cos(x)*cos(y)*cos(z);
    //return cos(t)*1*cos(z);
}

double Solution::u_x(double t, double x, double y, double z) {
    //return -sin(x)*cos(y)*z;
    //return pow(0.5 + pow(t - 0.5, 2), -1) * cos(x + M_PI / 6)*sin(y + M_PI / 6)*cos(z + M_PI / 6);
    return -sin(t)*sin(x)*cos(y)*cos(z);
    //return 0;
}

double Solution::u_y(double t, double x, double y, double z) {
    //return -cos(x)*sin(y)*z;
    //return pow(0.5 + pow(t - 0.5, 2), -1) * sin(x + M_PI / 6)*cos(y + M_PI / 6)*cos(z + M_PI / 6);
    return -sin(t)*cos(x)*sin(y)*cos(z);
    //return 0;
}

double Solution::u_z(double t, double x, double y, double z) {
    //return pow(0.5 + pow(t - 0.5, 2), -1) + cos(x)*cos(y);
    //return -pow(0.5 + pow(t - 0.5, 2), -1) * sin(x + M_PI / 6)*sin(y + M_PI / 6)*sin(z + M_PI / 6);
    return -sin(t)*cos(x)*cos(y)*sin(z);
    //return -sin(t)*1*sin(z);
}

double Solution::laplace_u(double t, double x, double y, double z) {
    //return -2 * cos(x)*cos(y)*z;
    //return pow(0.5 + pow(t - 0.5, 2), -1) * (-3* sin(x + M_PI / 6)*sin(y + M_PI / 6)*cos(z + M_PI / 6) );
    return -3*sin(t)*cos(x)*cos(y)*cos(z);
    //return -sin(t)*1*cos(z);
}