#include "nmlib.h"

/*Crea matriz de ceros.*/
double **crearMatriz(int rows, int cols){
    double ** matriz;
    matriz = (double **)calloc(rows, sizeof(double *));
    for(int i=0; i<rows;i++){
        matriz[i] = (double *)calloc(cols, sizeof(double));
    }
    return matriz;
}

void imprimirMatriz(double ** matriz,int rows, int cols){
    for(int i=0; i<rows;i++){
        for(int j=0; j<cols;j++){
            printf("%lf ",matriz[i][j]);
        }
        printf("\n");
    }
}

/*Esta función realiza una copia de una matriz sobre un puntero doble*/
double ** copiarMatriz(double ** matriz, int rows, int cols, int err_matriz){
    /*Se crea una matriz dinámicamente*/
    double ** copia = crearMatriz(rows,cols);
    /*Se copian los valores de la matriz dada en la nueva matriz*/
    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            copia[i][j] = matriz[i][j];
        }
    }
    /*Borramos la matriz ingresada si err_* == 1*/
    if(err_matriz == 1){
        for(int i=0;i<rows;i++){
            free(matriz[i]);
        }
        free(matriz);
    }
    /*Retornamos copia*/
    return copia;
}

/*Crea una matriz a partir de un archivo txt*/
double ** crearMatrizTXT(char * archivo, int * rows, int * cols){
    int i=0,j=0;
    int dim1=0, dim2=0;
    /*Se abre el archivo*/
    FILE * file = fopen(archivo,"r");
    /*Se verifica que se haya abierto correctamente*/
    if(file == NULL){
        perror("\nError al abrir archivo.\n");
        return NULL;
    }
    /*Escaneamos las dimensiones de la matriz*/
    fscanf(file,"%d %d",&dim1,&dim2);
    /*Actualizamos rows y cols*/
    * rows = dim1;
    * cols = dim2;
    /*Creamos una matriz dinámica para guardar la matriz*/
    double ** matriz = (double**)calloc(dim1,sizeof(double *));
    /*Verificamos que se haya asignado memoria correctamente*/
    if(matriz == NULL){
        printf("\nError al crear matriz.\n");
        return NULL;
    }
    /*Asignamos memoria a los renglones*/
    for(int i=0; i<dim1;i++){
        matriz[i] = (double *)calloc(dim2, sizeof(double));
    }
    /*Guardamos la matriz del txt en la matriz dinámica*/
    for(i=0;i<dim1;i++){
        for(j=0;j<dim2;j++){
            fscanf(file,"%lf",&matriz[i][j]);
        }
    }
    /*Cerramos el archivo*/
    fclose(file);
    /*Retornamos la matriz*/
    return matriz;
}

/*Esta función escribe en un archivo txt los valores de una matriz dada*/
void writeMatrizTXT(double ** matriz, int rows, int cols, char * nombre_archivo){
    int i=0,j=0;
    /*Abrimos el archivo*/
    FILE * salida = fopen(nombre_archivo,"w");
    /*Guardamos las dimensiones de la matriz*/
    fprintf(salida,"%d%c%d%c",rows,' ',cols,'\n');
    /*Guardamos la matriz*/
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            fprintf(salida, "%.12lf ", matriz[i][j]);
        }
        fprintf(salida, "\n");
    }
    /*Cerramos el archivo*/
    fclose(salida);
}

/*---------------------*/

/*Realiza el producto de dos matrices*/
double ** prodMatrices(double ** M, double ** N, int rows_M, int cols_M,int rows_N, int cols_N, int err_M, int err_N){
    int i=0,j=0,k=0;
    double sum=0;
    double ** producto;
    /*Se verficia que el producto se pueda hacer*/
    if(cols_M != rows_N){
        printf("El número de columnas de la primera matriz no coincide con el número de renglones de la segunda matriz.");
        return NULL;
    }
    /*Se crea una matriz dinámica producto.*/
    producto = crearMatriz(rows_M,cols_N);
    /*Se guarda el producto de las matrices en producto.*/
    for(i=0;i<rows_M;i++){
        for(j=0;j<cols_N;j++){
            sum=0;
            for(k=0;k<cols_M;k++){
                sum = sum + M[i][k]*N[k][j];
            }
            producto[i][j] = sum;
        }
    }
    /*Borramos la matriz ingresada si err_* == 1*/
    if( err_M == 1){ for(i=0;i<rows_M;i++){ free(M[i]); } free(M);}
    if( err_N == 1){ for(i=0;i<rows_N;i++){ free(N[i]); } free(N);}
    /*Retornamos producto*/
    return producto;
}

/*Calcula la matriz transpuesta de una matriz dada*/
double ** matrizTransp(double ** matriz, int rows, int cols,int err_matriz){
    double ** matTrans;
    int i=0,j=0;
    /*Creamos dinámicamente la matriz*/
    matTrans = crearMatriz(cols,rows);
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++){
            /*Guardamos la transpuesta*/
            matTrans[j][i] = matriz[i][j];
        }
    }
    /*Borramos el vector ingresado si err_* == 1*/
    if( err_matriz == 1){ for(i=0;i<rows;i++){ free(matriz[i]); } free(matriz);}
    /*Retornamos matTrans*/
    return matTrans;
}

/*---------------------*/

/*Rellena la k-ésima diagonal superior con un número dado*/
void DiagFillSup(double ** matriz, int rows, int cols, int k, int numero){
    //Verificar que k este entre 0 y cols-1
    int i=0, j=0;
    for(i=0;i<rows;i++){
        if(k+i<cols){
            matriz[i][k+i] = numero;
        }
    }
}

/*Rellena la k-ésima diagonal inferior con un número dado*/
void DiagFillInf(double ** matriz, int rows, int cols, int k, int numero){
    int i=0, j=0;
    for(i=0;i<cols;i++){
        if(k+i<rows){
            matriz[k+i][i] = numero;
        }
    }
}

/*Esta función halla la entrada mayor en valor absoluto de la parte triangular superior de una matriz.*/
void encontrarMayor(double ** matriz, int dim, double * mayor, int *r_mayor, int *c_mayor){
    int i,j;
    *mayor = matriz[0][1];
    for(i=0;i<dim;i++){
        for(j=i+1;j<dim;j++){
            if(fabs(matriz[i][j]) > fabs(*mayor)){
                *mayor = matriz[i][j]; //elemento mayor
                *r_mayor = i; //renglon del elemento mayor
                *c_mayor = j; //columna del elemento mayor
            }
        }
    }
}

