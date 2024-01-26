#include "nmlib.h"

/*Métodos para resolver sistemas de ecuaciones.*/
double * resolvDiag(double ** matriz,double * V, int dim){
    double * X;
    X = crearVector(dim);
    for(int i=0;i<dim;i++){
        /*Simplemente se realiza una división por cada renglón.*/
        X[i] = V[i]/matriz[i][i];
    }
    return X;
}

double * resolvTrianInf(double ** matriz,double * vect,int dim){
    double * X;
    int i=0, j=0;
    double sum=0;
    X = crearVector(dim);
    for(i=0;i<dim;i++){
        sum=0;
        for(j=0;j<i;j++){
            sum = sum + matriz[i][j]*X[j];
        }
        X[i] = (vect[i]-sum)/(matriz[i][i]) ;
    }
    return X;
}

double * resolvTrianSup(double ** matriz,double * vect,int dim){
    double * X;
    int i=0, j=0;
    double sum=0;
    X = crearVector(dim);
    for(i=dim-1;i>=0;i--){
        sum=0;
        for(j=i;j<dim;j++){
            sum = sum + matriz[i][j]*X[j];
        }
        X[i] = (vect[i]-sum)/(matriz[i][i]) ;
    }
    return X;
}

void factorLUCrout(double ** matriz, int dim, double ** L, double ** U){
    //No olvide inicializar los parmetros L y U con ceros antes de ingresarlos a esta función.
    int i=0,j=0,k=0;
    double sum=0;
    for(i=0;i<dim;i++){
        U[i][i]=1;
    }
    for(i=0;i<dim;i++){
        for(j=0;j<i;j++){//j = columna; i = renglon
            if(j==0){
                L[i][j] = matriz[i][j];
            }
            else{
                sum=0;
                for(k=0;k<j;k++){
                    sum = sum + L[i][k]*U[k][j];
                }
                L[i][j] = matriz[i][j] - sum;
            }
        }
        for(j=0;j<i;j++){//j = renglon; i = columna
            if(j==0){
                U[j][i] = matriz[j][i]/L[j][j];
            }
            else{
                sum=0;
                for(k=0;k<i;k++){
                    sum = sum + L[j][k]*U[k][i];
                }
                U[j][i] = (matriz[j][i] - sum)/L[j][j];
            }
        }
        sum=0;
        for(k=0;k<i;k++){
            sum = sum + L[i][k]*U[k][i];
        }
        L[i][i] = matriz[i][i] - sum;
    }
}

double * resolvLUCrout(double ** matriz, double * vect, int dim){
    double * X;
    double * Y;
    double ** L;
    double ** U;
    L = crearMatriz(dim,dim);
    U = crearMatriz(dim,dim);
    factorLUCrout(matriz,dim,L,U);
    Y = resolvTrianInf(L,vect,dim);
    X = resolvTrianSup(U,Y,dim);
    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(Y);
    return X;
}

void factorLUDoolittle(double ** matriz, int dim, double ** L, double ** U){
    //No olvide inicializar los parmetros L y U con ceros antes de ingresarlos a esta función.
    int i=0,j=0,k=0;
    double sum=0;
    for(i=0;i<dim;i++){
        L[i][i]=1;
    }
    for(i=0;i<dim;i++){
        for(j=0;j<i;j++){//j = columna; i = renglon
            if(j==0){
                L[i][j] = matriz[i][j]/U[j][j];
            }
            else{
                sum=0;
                for(k=0;k<j;k++){
                    sum = sum + L[i][k]*U[k][j];
                }
                L[i][j] = (matriz[i][j] - sum)/U[j][j];
            }
        }
        for(j=0;j<i;j++){//j = renglon; i = columna
            if(j==0){
                U[j][i] = matriz[j][i];
            }
            else{
                sum=0;
                for(k=0;k<i;k++){
                    sum = sum + L[j][k]*U[k][i];
                }
                U[j][i] = matriz[j][i] - sum;
            }
        }
        sum=0;
        for(k=0;k<i;k++){
            sum = sum + L[i][k]*U[k][i];
        }
        U[i][i] = matriz[i][i] - sum;
    }
}

double * resolvLUDoolittle(double ** matriz, double * vect, int dim){
    double * X;
    double * Y;
    double ** L;
    double ** U;
    L = crearMatriz(dim,dim);
    U = crearMatriz(dim,dim);
    factorLUDoolittle(matriz,dim,L,U);
    Y = resolvTrianInf(L,vect,dim);
    X = resolvTrianSup(U,Y,dim);
    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(Y);
    return X;
}

double * resolvLUCroutOpt(double ** matriz, double * vect, int dim){
    double * X;
    double * Y;
    double ** L;

    int i=0,j=0,k=0;
    double sum=0, sum1=0, sum2=0;

    L = crearMatriz(dim,dim);
    for(i=0;i<dim;i++){
        for(j=0;j<i;j++){//j = columna; i = renglon
            sum1=0;
            sum2=0;
            for(k=0;k<j;k++){
                sum1 = sum1 + L[i][k]*L[k][j];
                sum2 = sum2 + L[j][k]*L[k][i];
            }
            L[i][j] = matriz[i][j] - sum1;
            L[j][i] = (matriz[j][i] - sum2)/L[j][j];
        }
        sum=0;
        for(k=0;k<i;k++){
            sum = sum + L[i][k]*L[k][i];
        }
        L[i][i] = matriz[i][i] - sum;
    }

    Y = crearVector(dim);
    for(i=0;i<dim;i++){
        sum=0;
        for(j=0;j<i;j++){
            sum = sum + L[i][j]*Y[j];
        }
        Y[i] = (vect[i]-sum)/(L[i][i]);
    }

    X = crearVector(dim);
    for(i=dim-1;i>=0;i--){
        sum=0;
        for(j=i;j<dim;j++){
            sum = sum + L[i][j]*X[j];
        }
        X[i] = (Y[i]-sum) ;
    }

    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(L[i]);
    }
    free(L);
    free(Y);
    return X;
}

double ** factorCholesky(double ** matriz, int dim){
    double ** L;
    L = crearMatriz(dim,dim);
    int i=0,j=0, k=0;
    double sum=0;
    for(i=0;i<dim;i++){
        for(j=0;j<i;j++){
            sum = 0;
            for(k=0;k<j;k++){
                sum = sum + L[i][k]*L[j][k];
            }
            L[i][j] = (matriz[i][j] - sum)/L[j][j];
        }
        sum=0;
        for(k=0;k<i;k++){
            sum = sum + L[i][k]*L[i][k];
        }
        L[i][i] = sqrt(matriz[i][i]-sum);
    }
    return L;
}

double * resolvCholesky(double ** matriz, double * vect, int dim){
    double ** L;
    double ** LT;
    double * Y;
    double * X;
    L = factorCholesky(matriz,dim);
    LT = matrizTransp(L,dim,dim,0);
    Y = resolvTrianInf(L,vect,dim);
    X = resolvTrianSup(LT,Y,dim);
    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(L[i]);
        free(LT[i]);
    }
    free(L);
    free(LT);
    free(Y);
    return X;
}

void factorCholeskyMejorado(double ** matriz, int dim, double ** L,double ** D){
    /*No olvidar inicializar L y D antes de ingresarla a la función.*/
    int i=0,j=0, k=0;
    double sum=0;

    for(i=0;i<dim;i++){
        L[i][i] = 1;
    }

    D[0][0] = matriz[0][0];
    for(i=1;i<dim;i++){
        for(j=0;j<i;j++){
            if(j==0){
                L[i][j] = matriz[i][j]/D[j][j];
            }
            else{
                sum=0;
                for(k=0;k<j;k++){
                    sum = sum + D[k][k]*L[i][k]*L[j][k];
                }
                L[i][j] = (matriz[i][j]-sum)/D[j][j];
            }
        }
        sum=0;
        for(k=0;k<i;k++){
            sum = sum + D[k][k]*L[i][k]*L[i][k];
        }
        D[i][i] = matriz[i][i] - sum;
    }
}

double * resolvCholeskyMejorado(double ** matriz, double * V,int dim){
    double ** L;
    double ** D;
    double ** LT;
    double * X;
    double * Y;
    double * Z;

    L = crearMatriz(dim,dim);
    D = crearMatriz(dim,dim);
    factorCholeskyMejorado(matriz,dim,L,D);
    LT = matrizTransp(L,dim,dim,0);
    Z = resolvTrianInf(L,V,dim);
    Y = resolvDiag(D,Z,dim);
    X = resolvTrianSup(LT,Y,dim);
    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(L[i]);
        free(D[i]);
        free(LT[i]);
    }
    free(L);
    free(D);
    free(LT);
    free(Y);
    free(Z);
    return X;
}

double * resolvJacobi(double ** matriz, double * vect, double * X_init, double tol, int dim){
    double * X_old;
    double * X_new;
    double vector_dist = tol+1;
    int i=0,j=0,k=0;
    double sum=0;
    double * zero_vector;
    X_old = copiarVector(X_init,dim,0);
    X_new = crearVector(dim);
    zero_vector = crearVector(dim);

    while(vector_dist > tol){
        for(i=0;i<dim;i++){
            sum=0;
            for(j=0;j<dim;j++){
                if(i!=j){
                    sum = sum + matriz[i][j]*X_old[j];
                }
            }
            X_new[i] = (vect[i]-sum)/matriz[i][i];
        }
        vector_dist = vectorDistance(X_old,X_new,dim)/vectorDistance(X_new,zero_vector,dim);
        free(X_old);//Liberamos memoria antes reiniciar el vector
        X_old =  copiarVector(X_new,dim,0);
    }
    /*Liberar memoria*/
    free(X_old);
    free(zero_vector);

    return X_new;
}

double * resolvGaussSeidel(double ** matriz, double * vect, double * X_init, double tol, int dim){
    double * X_old;
    double * X_new;
    double vector_dist = tol+1;
    int i=0,j=0,k=0;
    double sum=0;
    double * zero_vector;
    X_old = copiarVector(X_init,dim,0);
    X_new = crearVector(dim);
    zero_vector = crearVector(dim);

    while(vector_dist > tol){
        for(i=0;i<dim;i++){
            sum=0;
            for(j=0;j<dim;j++){
                if(i!=j){
                    sum = sum + matriz[i][j]*X_new[j];
                }
            }
            X_new[i] = (vect[i]-sum)/matriz[i][i];
        }
        vector_dist = vectorDistance(X_old,X_new,dim)/vectorDistance(X_new,zero_vector,dim);
        for(i=0;i<dim;i++){
            X_old[i] =  X_new[i];
        }
    }

    /*Liberar memoria*/
    free(X_old);
    free(zero_vector);

    return X_new;
}

double * resolvElimGauss(double ** matriz, double * vect, int dim){
    int sigma[dim]; //Guarda la permutación de los renglones de X, cuando se mueven las columnas.
    double pivote =0;//Variable que guarda el pivote.
    int r_pivote = 0;//Variable que guarda el renglón donde esta ubicado el pivote.
    int c_pivote = 0;//Variable que guarada la columna donde esta ubicado el pivote.
    double aux=0;
    int i=0,j=0,k=0;
    double sum=0;
    double * X;
    X = crearVector(dim);
    double ** M;
    M = copiarMatriz(matriz,dim,dim,0);
    double * V;
    V = copiarVector(vect,dim,0);
    //Inicializamos X_sol
    for(k=0;k<dim;k++){
        sigma[k] = k;
    }
    for(k=0;k<dim-1;k++){
        //Pivoteo
        pivote = 0;
        r_pivote = k;
        c_pivote = k;
        //Buscamos el pivote. Lo guardamos junto con sus coordenadas.
        for(i=k;i<dim;i++){
            for(j=k;j<dim;j++){
                if(fabs(pivote)<fabs(M[i][j])){
                    pivote = M[i][j];
                    r_pivote = i;
                    c_pivote = j;
                }
            }
        }
        //Movemos el pivote a la primera posición de la matriz.
        //Primero movemos los renglones.
        aux=0;
        for(j=0;j<dim;j++){
            aux = M[k][j];
            M[k][j] = M[r_pivote][j] ;
            M[r_pivote][j] = aux;
        }
        aux = V[k];
        V[k] = V[r_pivote];
        V[r_pivote] = aux;
        //Luego movemos las columnas.
        aux=0;
        for(j=0;j<dim;j++){
            aux = M[j][k];
            M[j][k] = M[j][c_pivote] ;
            M[j][c_pivote] = aux;
        }
        aux = sigma[k];
        sigma[k] = sigma[c_pivote];
        sigma[c_pivote] = aux;
        //Gauss
        pivote = M[k][k];
        for(i=k;i<dim;i++){
            M[k][i] = M[k][i]/pivote;
        }
        V[k] = V[k]/pivote;

        for(i=k+1;i<dim;i++){
            pivote = M[i][k];
            for(j=k;j<dim;j++){
                M[i][j] = M[i][j]-pivote*M[k][j];
            }
            V[i] = V[i]-pivote*V[k];
        }
    }
    //Solución de la triangular superior.
        for(i=dim-1;i>=0;i--){
            sum=0;
            for(j=i+1;j<dim;j++){
                sum = sum + M[i][j]*X[sigma[j]];
            }
            X[sigma[i]] = (V[i]-sum)/(M[i][i]) ;
        }

    /*Liberar memoria*/
    for(int i=0;i<dim;i++){
        free(M[i]);
    }
    free(M);
    free(V);

    return X;
}


