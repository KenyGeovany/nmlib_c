#ifndef NMLIB_H
#define NMLIB_H

/* COMMON LIBRARIES. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* CONSTANTS. */

#ifndef PI_constant
#define PI_constant 3.1415926535897932
#endif // PI_constant

#ifndef E_constant
#define E_constant 2,71828182845904523
#endif // E_constant

/*FUNCTIONS*/

/*1. Matrix tools. */
double **crearMatriz(int rows, int cols);
void imprimirMatriz(double ** matriz,int rows, int cols);
double ** copiarMatriz(double ** matriz, int rows, int cols, int err_matriz);
double ** crearMatrizTXT(char * archivo, int * rows, int * cols);
void writeMatrizTXT(double ** matriz, int rows, int cols, char * nombre_archivo);
/*---------------------*/
double ** prodMatrices(double ** M, double ** N, int rows_M, int cols_M,int rows_N, int cols_N, int err_M, int err_N);
double ** matrizTransp(double ** matriz, int rows, int cols,int err_matriz);
/*---------------------*/
void DiagFillSup(double ** matriz, int rows, int cols, int k, int numero);
void DiagFillInf(double ** matriz, int rows, int cols, int k, int numero);
void encontrarMayor(double ** matriz, int dim, double * mayor, int *r_mayor, int *c_mayor);

/*2. Vector tools. */
double * crearVector(int rows);
void imprimirVector(double * vect,int rows);
double * copiarVector(double * vect, int rows,int err_vect);
double * crearVectorTXT(char * archivo, int * rows);
void writeVectorTXT(double * vect, int rows, char * nombre_archivo);
/*---------------------*/
double * sumaVectores(double * vector1,double * vector2,int dim, int err_vector1, int err_vector2);
double productoPunto(double * vector1, double * vector2,int dim);
double moduloVector(double * vect,int dim);
double * normalizarVector(double * vect, int dim ,int err_vect);
double vectorDistance(double * V1, double * V2, int rows);
double * sumaRealVectorCW(double * vect, double numero,int rows, int err_vect);
double * multiplicaRealVectorCW(double * vect, double numero,int rows, int err_vect);

/*3. Lineal algebra tools. */
double * prodMatVec(double ** M, double * V, int rows_M, int cols_M,int rows_V);
double ** GramSchmidt(double ** matriz, int dim, int err_matriz);

/*4. Lineal system equations solvers.*/
double * resolvDiag(double ** matriz,double * V, int dim);
double * resolvTrianInf(double ** matriz,double * vect,int dim);
double * resolvTrianSup(double ** matriz,double * vect,int dim);
void factorLUCrout(double ** matriz, int dim, double ** L, double ** U);
double * resolvLUCrout(double ** matriz, double * vect, int dim);
void factorLUDoolittle(double ** matriz, int dim, double ** L, double ** U);
double * resolvLUDoolittle(double ** matriz, double * vect, int dim);
double * resolvLUCroutOpt(double ** matriz, double * vect, int dim);
double ** factorCholesky(double ** matriz, int dim);
double * resolvCholesky(double ** matriz, double * vect, int dim);
void factorCholeskyMejorado(double ** matriz, int dim, double ** L,double ** D);
double * resolvCholeskyMejorado(double ** matriz, double * V,int dim);
double * resolvJacobi(double ** matriz, double * vect, double * X_init, double tol, int dim);
double * resolvGaussSeidel(double ** matriz, double * vect, double * X_init, double tol, int dim);
double * resolvElimGauss(double ** matriz, double * vect, int dim);

/*5. Eigenvalue problem solvers.*/
double * metodoPotencia(double ** matriz, double * V_init, int dim, double lambda_init, double tol,double * valor_propio);
double * metodoPotenciaInversa(double ** matriz, double * V_init, double dim, double lambda_init, double tol,double * valor_propio);
double * metodoPotenciaConDeflacion(double ** matriz, double * V_init, int dim, int numVProp , int numIter, double tol,double ** PHI);
double * metodoPotenciaConDeflacionInversa(double ** matriz, double * V_init, int dim, int numVProp , int numIter, double tol,double ** PHI);
double * metodoJacobiVPropios(double ** matriz, int dim, int numIter, double tol,double ** PHI);
double metodoRaylegh(double ** matriz,int dim, double lambda_init,double * phi_init,int numIter, double tol, double * phi);
double * metodoSubespacios(double ** matriz,int dim,int s,double ** PHI_init,int numIter, double tol,double ** PHI);
double ** factorQR(double ** matriz,int dim,double ** Q);
double * metodoQR(double ** matriz,int dim, int numIter, double tol, double **PHI);

/*6. Root finding methods. */
double newton_method_modified(double(*f)(double),double(*f_der1)(double),double(*f_der2)(double),double x0, double tol_x,double tol_f);


#endif /* NMLIB_H */
