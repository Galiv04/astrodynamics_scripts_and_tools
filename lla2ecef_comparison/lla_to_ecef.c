#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>

#define GS_ACDS_deg2rad 0.01745329252

void lla_to_ecef(float lat, float lon, float alt, float *x, float *y, float *z);

void lla_to_ecef(float lat, float lon, float alt, float *x, float *y, float *z){
    
    float R_oblat;
    lat = lat * GS_ACDS_deg2rad;
    lon = lon * GS_ACDS_deg2rad;
    
    // WGS84 ellipsoid constants:
    float a = 6378137.0;
    float e = 8.1819190842622e-2;

    // intermediate calculation
    // (prime vertical radius of curvature)
    R_oblat = a / sqrt(1.0 - powf(e, 2.0) * powf(sin(lat), 2.0));
    
    *x = (R_oblat + alt) * cos(lat)*cos(lon);
    *y = (R_oblat + alt) * cos(lat)*sin(lon);
    *z = ((1.0 - powf(e, 2.0)) * R_oblat + alt) * sin(lat);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <latitude> <longitude> <altitude>\n", argv[0]);
        return 1;
    }

    float lat = atof(argv[1]);
    float lon = atof(argv[2]);
    float alt = atof(argv[3]);
    float x, y, z;
    
    lla_to_ecef(lat, lon, alt, &x, &y, &z);
    
    printf("{\"x\": %.6f, \"y\": %.6f, \"z\": %.6f}\n", x, y, z);
    return 0;
}