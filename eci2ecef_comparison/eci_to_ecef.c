#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Constants
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

// Function 1
void eci2ecef(double JulianDate, double* q) {
    double N_r;
    double T_n;
    float psi;
    // Seconds on a day
    double daysec = 86400.0;
    // Earth precession [rad/sec]
    double omega = 7.292115855377075e-05;
    // Sidereal epoch (raf's)
    double s = (6.23018356e+04 / daysec + 2450448.5);
    // Time for a revolution
    double T_r = (2.0 * M_PI / omega);

    // Number of revolutions since epoch
    N_r = roundf((float)(((JulianDate - s) * daysec) / T_r));

    // Time into a revolution
    T_n = (JulianDate - s) * daysec - N_r * T_r;

    // Calculate angle with 2pi-boundary
    psi = (float) fmod(omega * T_n, 2 * M_PI);

    q[0] = 0;
    q[1] = 0;
    q[2] = sin(psi / 2.0);
    q[3] = cos(psi / 2.0);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input>\n", argv[0]);
        return 1;
    }
    double JulianDate = atof(argv[1]);
    double q[4];
    eci2ecef(JulianDate, q);
    printf("%lf %lf %lf %lf\n", q[0], q[1], q[2], q[3]);
    return 0;
}