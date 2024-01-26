#include "nmlib.h"

/*Crea vector de ceros.*/
double * crearVector(int rows){
    double * vector;
    vector = (double *)calloc(rows, sizeof(double));
    return vector;
}

/*Imprime un vector*/
void imprimirVector(double * vect,int rows){
    for(int i=0; i<rows;i++){
        printf("%lf \n",vect[i]);
    }
}

/*Copia un vector, en una nueva matriz creada dinámicamente*/
double * copiarVector(double * vect, int rows,int err_vect){
    double * copia = crearVector(rows);
    for(int i=0;i<rows;i++){
        copia[i] = vect[i];
    }
    if(err_vect == 1) free(vect); /*Borramos el vector ingresado si err_* == 1*/
    return copia;
}

/*Crea un vector a partir de un archivo txt*/
double * crearVectorTXT(char * archivo, int * rows){
    int i=0,j=0;
    int dim=0;
    /*Abrimos el archivo txt*/
    FILE * file = fopen(archivo,"r");
    /*Verificamos que se halla abierto correctamente*/
    if(file == NULL){
        perror("\nError al abrir archivo.\n");
        return NULL;
    }
    /*Guardamos las dimensiones de la matriz*/
    fscanf(file,"%d",&dim);
    * rows = dim;
    /*Creamos un vector dinámico*/
    double * vector = (double*)calloc(dim,sizeof(double));
    /*Verificamos que se haya asignado memoria correctamente*/
    if(vector == NULL){
        printf("\nError al crear matriz.\n");
        return NULL;
    }
    /*Guardamos los valores en el vector*/
    for(i=0;i<dim;i++){
        fscanf(file,"%lf",&vector[i]);
    }
    /*Cerramos el archivo*/
    fclose(file);
    /*Retornamos vector*/
    return vector;
}

/*Esta función escribe un vector en un archivo txt*/
void writeVectorTXT(double * vect, int rows, char * nombre_archivo){
    int i=0,j=0;
    /*Abrimos el archivo*/
    FILE * salida = fopen(nombre_archivo,"w");
    /*Guardamos la dimension del vector*/
    fprintf(salida,"%d%c",rows,'\n');
    /*Guardamos los valores del vector*/
    for(i=0;i<rows;i++){
        fprintf(salida, "%.12lf%c", vect[i],'\n');
    }
    /*Cerramos el archivo*/
    fclose(salida);
}

/*---------------------*/

/*Realiza la suma de dos vectores*/
double * sumaVectores(double * vector1,double * vector2,int dim, int err_vector1, int err_vector2){
    double * resultado = (double *)calloc(dim,sizeof(double));
    for(int i=0;i<dim;i++){
        resultado[i] = vector1[i]+vector2[i];
    }
    /*Borramos el vector ingresado si err_* == 1*/
    if(err_vector1 == 1){ free(vector1); }
    if(err_vector2 == 1){ free(vector2); }
    /*Retornamos resultado*/
    return resultado;
}

/*Realiza el producto punto de dos vectores*/
double productoPunto(double * vector1, double * vector2,int dim){
    double sum = 0;
    for(int i=0;i<dim;i++){
        sum = sum + vector1[i]*vector2[i];
    }
    return sum;
}

/*Calcula el módulo de un vector dado*/
double moduloVector(double * vect,int dim){
    double sum = 0;
    for(int i=0;i<dim;i++){
        sum = sum + vect[i]*vect[i];
    }
    sum = sqrt(sum);
    return sum;
}

/*Normaliza un vector dado*/
double * normalizarVector(double * vect, int dim ,int err_vect){
    double * vectNorm = (double *)calloc(dim,sizeof(double));
    double modulo = moduloVector(vect,dim);
    for(int i=0;i<dim;i++){
        vectNorm[i] = vect[i]/modulo;
    }
    /*Borramos el vector ingresado si err_* == 1*/
    if(err_vect == 1){ free(vect); }
    /*Retornamos vectNorm*/
    return vectNorm;
}

/*Calcula la distancia euclidiana entre dos vectores dados*/
double vectorDistance(double * V1, double * V2, int rows){
    double distance = 0;
    double sum=0;
    for(int i=0;i<rows;i++){
        sum = sum + (V1[i]-V2[i])*(V1[i]-V2[i]);
    }
    /*Aquí asemos uso de la librería math.h*/
    distance = sqrt(sum);
    /*Retornamos distance*/
    return distance;
}

/*Suma un real a cada entrada de un vector*/
double * sumaRealVectorCW(double * vect, double numero,int rows, int err_vect){
    //err_vect == 1 borrar vect antes de salir de la función.
    //Esta opción evita posibles perdidas de memoria como por ejemplo en autoreferencias: v = f(v,...)
    double * resultado = crearVector(rows);
    for(int i=0;i<rows;i++){
        resultado[i] = vect[i]+numero;
    }
    /*Borramos el vector ingresado si err_* == 1*/
    if(err_vect == 1){ free(vect); }
    /*Retornamos resultado*/
    return resultado;
}

/*Multiplicar un escalara a cada componente de un vector */
double * multiplicaRealVectorCW(double * vect, double numero,int rows, int err_vect){
    double * resultado = crearVector(rows);
    for(int i=0;i<rows;i++){
        resultado[i] = vect[i]*numero;
    }
    /*Borramos el vector ingresado si err_* == 1*/
    if(err_vect == 1){ free(vect); }
    /*Retornamos resultado*/
    return resultado;
}

