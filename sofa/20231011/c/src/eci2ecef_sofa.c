#include <stdio.h>
#include "sofa.h"

#define arcsec2rad 4.848136811095360e-06;

double JD_UTC_f1, JD_UTC_f2, JD_TT_f1, JD_TT_f2, JD_UTC1_f1, JD_UTC1_f2, xp, yp, dut1, rc2t[3][3];

// Function to print a 3x3 matrix
void printMatrix(double matrix[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {

    xp = 0.073751*arcsec2rad;
    yp = 0.475647*arcsec2rad;

    dut1 = -0.0077275; // seconds

    JD_UTC_f1 = 2460485.5; // Fraction of JD UTC of the date part - 2024-06-24 00:00
    JD_UTC_f2 = 0; // Fraction of JD UTC of the time part

    JD_TT_f1 = JD_UTC_f1;
    JD_TT_f2 = JD_UTC_f2 + 59.184/86400.0; // offset from UTC to TT -> (TT = TAI + 32.184) & (TAI = UTC + 37 but the offset is not fixed) // Also (TAI = GPS + 19) this one should be used

    JD_UTC1_f1 = JD_UTC_f1;
    JD_UTC1_f2 = JD_UTC_f2 + dut1/86400.0; // offset from UTC to TT -> (TT = TAI + 32.184) & (TAI = UTC + 33)


    iauC2t06a(JD_TT_f1, JD_TT_f2, JD_UTC1_f1, JD_UTC1_f2, xp, yp, rc2t);

    // Print the matrix
    printf("3x3 Matrix:\n");
    printMatrix(rc2t);

    return 0;
}

  