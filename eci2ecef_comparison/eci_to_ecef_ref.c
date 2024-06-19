#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Constants
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

double calculate_gmst(double JulianDate) {
    double T = (JulianDate - 2451545.0) / 36525.0;
    double gmst = fmod(280.46061837 + 360.98564736629 * (JulianDate - 2451545.0) + 0.000387933 * T * T - T * T * T / 38710000.0, 360.0);
    
    return gmst * M_PI / 180.0; // Convert to radians
}

void compute_quaternion(double JulianDate, double* q) {
    double psi = calculate_gmst(JulianDate);
    
    q[0] = 0;
    q[1] = 0;
    q[2] = sin(psi / 2.0);
    q[3] = cos(psi / 2.0);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <JulianDate>\n", argv[0]);
        return 1;
    }
    double JulianDate = atof(argv[1]);
    double q[4];
    compute_quaternion(JulianDate, q);
    printf("%lf %lf %lf %lf\n", q[0], q[1], q[2], q[3]);
    return 0;
}
