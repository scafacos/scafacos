#ifdef akosoaksoaks



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "triclinic.h"

#define unit_volume(x,y,z) sqrt((1-cos(x)*cos(x)-cos(y)*cos(y)-cos(z)*cos(z)+2*cos(x)*cos(z)*cos(y)))//volume of triclinic cell in cartesian units for sidelengths 1 (in cartesian units)


//pseudo init




/*transforms vec[3] from triclinic coordinates to cartesian coordiantes
 *TODO: get angles and lengths automatically from box vectors
 *
 *
 *
 *
 */
  
void cartFROMtric(float *vec[]) {

    vec[0] = BASIS_A * vec[0] + BASIS_B * cos(GAMMA) * vec[1]
            + BASIS_C * cos(BETA) * vec[2]; //a*m+b*cos(gamma)*n+c*cos(beta)*l
    vec[1] = BASIS_B * sin(GAMMA) * vec[1]
            + BASIS_C * (cos(ALPHA) - cos(BETA) * cos(GAMMA)) / sin(GAMMA)
            * vec[2]; //b*sin(gamma)*n+c*(cos(alpha)-cos(beta)*cos(gamma))/(sin(gamma))*l
    vec[2] = BASIS_C * unit_volume(ALPHA, BETA, GAMMA) / sin(GAMMA) * vec[2]; //c*v/sin(gamma)*l
}

/*transforms vec[3] from cartesians to triclinic coordiantes
 *
 *
 *
 */
void tricFROMcart(float *vec[]) {

    vec[0] = vec[0] / BASIS_A - cos(GAMMA) / (BASIS_A * sin(GAMMA)) * vec[1]
            + (cos(ALPHA) * cos(GAMMA) - cos(BETA))
            / (BASIS_A * unit_volume(ALPHA, BETA, GAMMA) * sin(GAMMA))
            * vec[2];
    vec[1] = 1 / (BASIS_B * sin(GAMMA)) * vec[1]
            + (cos(BETA) * cos(GAMMA) - cos(ALPHA)) / sin(GAMMA) / BASIS_B
            / unit_volume(ALPHA, BETA, GAMMA) * vec[2];
    vec[2] = 1 / (BASIS_C * unit_volume(ALPHA, BETA, GAMMA)) * sin(GAMMA)
            * vec[2];
}

/*calculates and returns the angle between to vectors vec_a and vec_b.
 *the vectors have to be given in cartesian coordinates.
 *
 *@param *vec_a pointer to a float array that contains the vector in cartesian coordinates.
 *@param *vec_b pointer to a float array that contains the vector in cartesian coordinates.
 *@return the angle between vec_a and vec_b in radians as a float
 */
float angleFROMvecs(float *vec_a, float *vec_b) {
    return acos(
            (vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] + vec_b[2])
            / (absVec(vec_a) * absVec(vec_b)));
}




/*Calculate the Euclidean norm of a cartesian vector 
 *@param *a pointer to a float array that contains a cartesian vector
 *
 *@return the Euclidean norm of the vector a
 */
float absVec(float *a) {
    return (float) sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}
int box_a[3] = {1, 0, 0};
int box_b[3] = {0, 1, 0};
int box_c[3] = {0, 0, 1};


/*

int main(void) {


    //calc fundamental variables



    //real program


    float *vec = NULL;
    vec = malloc(3 * sizeof (float));

    vec[0] = 3.5;
    vec[1] = 1.5;
    vec[2] = 0;
    int i;
    tricFROMcart(vec);

    for (i = 0; i < 3; i++) {

        printf("%d th component: %f \n", i, vec[i]);

    }

    cartFROMtric(vec);
    for (i = 0; i < 3; i++) {

        printf("%d th component: %f \n", i, vec[i]);

    }
    free(vec);

    printf("sinus Gamma is %f\n", sin(GAMMA));
    printf("unitvolume is :%f \n",
            unit_volume(ALPHA, BETA, GAMMA));
    return 0;
}

*/
#endif