#include "nmlib.h"

/*Realiza el producto de una matriz por un vector. La función devuelve una matriz creada dinámicamente.*/
double * prodMatVec(double ** M, double * V, int rows_M, int cols_M,int rows_V){
    int i=0,k=0;
    /*Verificamos que el producto se pueda hacer.*/
    if(cols_M != rows_V){
        printf("El número de columnas de la matriz no coincide con el número de renglones del vector.");
        return NULL;
    }
    double sum=0;
    /*Creamos el vector dinámico*/
    double * producto;
    producto = crearVector(rows_M);
    /*Relizamos el producto guardando el resultado en cada iteración.*/
    for(i=0;i<rows_M;i++){
        sum=0;
        for(k=0;k<cols_M;k++){
            sum = sum + M[i][k]*V[k];
        }
        producto[i] = sum;
    }
    /*Se retorna producto*/
    return producto;
}

/*Algoritmo de Ortonormalización de Gram-Schmidt*/
double ** GramSchmidt(double ** matriz, int dim, int err_matriz){
    /*Variables importantes*/
    double ** E;
    double **E_T, **M_T;
    /*Variables auxiliares*/
    int i,j;
    double * sum, * V;

    /*Definimos las transpuestas de matriz y E*/
    E_T = (double**)calloc(dim,sizeof(double*));
    M_T = matrizTransp(matriz,dim,dim,0);

    /*Procedemos a ortonormalizar*/
    for(i=0;i<dim;i++){
        sum = copiarVector(M_T[i],dim,0);
        for(j=0;j<i;j++){
            /*Restamos las contribuciones de los vectores ya ortnormalizados*/
            V = multiplicaRealVectorCW(E_T[j],-productoPunto(E_T[j],M_T[i],dim),dim,0);
            sum = sumaVectores(sum,V,dim,1,1);
        }
        /*Normalizamos y guardamos el resultado en el siguiente renglón de E_T*/
        E_T[i] = normalizarVector(sum,dim,1);

    }

    /*Transponemos E_T para convertir los vectores ortnormales en renglones a columnas*/
    E = matrizTransp(E_T,dim,dim,1);

    /*Opción de borrar matriz*/
    if(err_matriz == 1){
        for(i=0;i<dim;i++) free(matriz[i]); free(matriz);
    }

    //for(i=0;i<dim;i++) free(E_T[i]); free(E_T);
    for(i=0;i<dim;i++) free(M_T[i]); free(M_T);

    /*Retornamos E*/
    return E;
}

