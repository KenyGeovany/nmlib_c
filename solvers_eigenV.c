#include "nmlib.h"

double * metodoPotencia(double ** matriz, double * V_init, int dim, double lambda_init, double tol,double * valor_propio){
    /*Variables*/
    double * V0;
    double * V1;
    double numerador, denominador, lambda_new, lambda_old;
    int i=0;
    double sum=0;

    /*Inicializamos V0 como el vector V_init*/
    V0 = copiarVector(V_init,dim,0);
    /*Inicializamos lambda_old y lambda_new*/
    lambda_old = lambda_init;
    lambda_new = tol+1;

    for(i=0;(fabs(lambda_new-lambda_old))>tol;i++){
        /*En cada iteración hacemos el producto MV0 -> V1*/
        V1 = prodMatVec(matriz,V0,dim,dim,dim);
        numerador = productoPunto(V1,V0,dim);
        denominador = productoPunto(V0,V0,dim);
        /*Guardamos el anterior lambda en lambda_old*/
        lambda_old = lambda_new;
        /*Actualizamos lambda_new*/
        lambda_new = numerador/denominador;
        free(V0);
        /*Normalizamos V0*/
        V0 = normalizarVector(V1,dim,1);
    }

    /*Guardamos el valor final de lambda_new en valor_propio*/
    *valor_propio = lambda_new;

    /*Retornamos V0 que es el vector propio correspondiente al valor propio hallado.*/
    return V0;
}

double * metodoPotenciaInversa(double ** matriz, double * V_init, double dim, double lambda_init, double tol,double * valor_propio){
    /*Variables*/
    double * V0;
    double * V1;
    double numerador, denominador, lambda_new, lambda_old;
    int i=0;
    double sum=0;

    /*Inicializamos V0 como el vector V_init*/
    V0 = copiarVector(V_init,dim,0);
    /*Inicializamos lambda_old y lambda_new*/
    lambda_old = 1000000000000;
    lambda_new = 0;

    for(i=0;(fabs(lambda_new-lambda_old))>tol;i++){
        /*En cada iteración resolvemos el sistema MV1 = V0*/
        V1 = resolvElimGauss(matriz,V0,dim);
        numerador = productoPunto(V0,V0,dim);
        denominador = productoPunto(V1,V0,dim);
        /*Guardamos el anterior lambda en lambda_old*/
        lambda_old = lambda_new;
        /*Actualizamos lambda_new*/
        lambda_new = numerador/denominador;
        free(V0);
        /*Normalizamos V0*/
        V0 = normalizarVector(V1,dim,1);
    }

    /*Guardamos el valor final de lambda_new en valor_propio*/
    *valor_propio = lambda_new;

    /*Retornamos V0 que es el vector propio correspondiente al valor propio hallado.*/
    return V0;
}

double * metodoPotenciaConDeflacion(double ** matriz, double * V_init, int dim, int numVProp , int numIter, double tol,double ** PHI){
   /*Inicializar previamente PHI con ceros.*/
    //Variables
    double * v;
    double lambda;
    double a;
    double * V1;
    double numerador, denominador, lambda_new, lambda_old;
    double * V_copy;
    //Variables Auxiliares
    int i=0,j=0,k=0;

    /*Creamos un arreglo que contendra los valores propios hallados.*/
    double * LAMBDA;
    LAMBDA = (double*)calloc(dim,sizeof(double));
    //------------------------------------------

    for(i=0;i<numVProp;i++){
        printf("\n%d\n",i);
        /*En cada iteracion reiniciamos v como el vector inicial V_init*/
        v = copiarVector(V_init  ,dim,0);
        /*Inicializamos vaiables*/
        lambda_old = 0;
        lambda_new = tol+1;
        /*Ejecutamos una versión modificada del método de la potencia.*/
        k=0;
        while((k<numIter)&&(fabs(lambda_new-lambda_old)>tol)){
            /* En cada iteración del método de al potencia iniciamos el while
             * con con vector inicial el último vector propio hallado. A este vector
             * se le restan las contribuciones de los anteriores vectores propios.
             */
            V_copy = copiarVector(v,dim,0);
            for(j=0;j<i;j++){
                a = productoPunto(PHI[j],V_copy  ,dim);
                v = sumaVectores(v,multiplicaRealVectorCW(PHI[j],-a,dim,0),dim,1,1) ;
            }
            /* Se normaliza el vector. El 1 indica que antes de guardar el vector normalizado
             * en v, se debe liberar la memoria dde v en la función antes de salir.
             */
            v = normalizarVector(v,dim,1);

            /*Hacemos el producto MV -> V1*/
            V1 = prodMatVec(matriz,v,dim,dim,dim);
            numerador = productoPunto(V1,v,dim);
            denominador = productoPunto(v,v,dim);
            /*Guardamos el anterior lambda en lambda_old*/
            lambda_old = lambda_new;
            /*Actualizamos lambda_new*/
            lambda_new = numerador/denominador;
            free(v);
            /*Normalizamos V0*/
            v = normalizarVector(V1,dim,1);
            /*Liberamos memoria y sumamos k en uno*/
            free(V_copy);
            k++;
        }
        /*Al salir del while guardamos el valor de lambda_new en LAMBDA[i]
         *y el valor de v en PHI[i]*/
        LAMBDA[i] = lambda_new;
        PHI[i] = copiarVector(v,dim,0);
        /*Liberamos memoria*/
        free(v);
    }

    /*Retornamos PHI, una matriz de tamaño numIter*dim, es decir con numIter renglones.*/
    /*Retornamos LAMBDA que es el vector de valores propios de tamaño numIter*/
    return LAMBDA;

}

/*Método de la potencia Inversa con Deflación.*/
double * metodoPotenciaConDeflacionInversa(double ** matriz, double * V_init, int dim, int numVProp , int numIter, double tol,double ** PHI){
    //Variables
    double a;
    double * v;
    double lambda;
    double * V1;
    double numerador, denominador, lambda_new, lambda_old;
    double * V_copy;
    //Auxiliares
    int i=0,j=0,k=0;

    /*Creamos un arreglo que contendra los valores propios hallados.*/
    double * LAMBDA;
    LAMBDA = (double*)calloc(dim,sizeof(double));
    //------------------------------------------

    for(i=0;i<numVProp;i++){
        /*En cada iteracion reiniciamos v como el vector inicial V_init*/
        v = copiarVector(V_init  ,dim,0);
        /*Inicializamos vaiables*/
        lambda_old = 0;
        lambda_new = tol+1;
        k=0;
        /*Ejecutamos una versión modificada del método de la potencia.*/
        while((k<numIter)&&(fabs(lambda_new-lambda_old)>tol)){
            /* En cada iteración del método de al potencia iniciamos el while
             * con con vector inicial el último vector propio hallado. A este vector
             * se le restan las contribuciones de los anteriores vectores propios.
             */
            V_copy = copiarVector(v,dim,0);
            for(j=0;j<i;j++){
                a = productoPunto(PHI[j],V_copy  ,dim);
                v = sumaVectores(v,multiplicaRealVectorCW(PHI[j],-a,dim,0),dim,1,1) ;
            }
            /* Se normaliza el vector. El 1 indica que antes de guardar el vector normalizado
             * en v, se debe liberar la memoria dde v en la función antes de salir.
             */
            v = normalizarVector(v,dim,1);

            /*Hacemos el producto MV -> V1*/
            V1 = resolvElimGauss(matriz,v,dim);
            numerador = productoPunto(v,v,dim);
            denominador = productoPunto(V1,v,dim);
            /*Guardamos el anterior lambda en lambda_old*/
            lambda_old = lambda_new;
            /*Actualizamos lambda_new*/
            lambda_new = numerador/denominador;
            free(v);
            /*Normalizamos V0*/
            v = normalizarVector(V1,dim,1);
            /*Liberamos memoria y sumamos k en uno*/
            free(V_copy);
            k++;
        }
        /*Al salir del while guardamos el valor de lambda_new en LAMBDA[i]
         *y el valor de v en PHI[i]*/
        LAMBDA[i] = lambda_new;
        PHI[i] = copiarVector(v,dim,0);
        /*Liberamos memoria*/
        free(v);
    }

    /*Retornamos LAMBDA que es el vector de valores propios*/
    return LAMBDA;

}

/*Metodo de Jacobi*/
double * metodoJacobiVPropios(double ** matriz, int dim, int numIter, double tol,double ** PHI){
    //Nota importante: PHI debe semi-inicializarse antes, PHI = (double**)calloc(dim,sizeof(double*))
    //Variables
    int i,j,k;
    double ** M;
    double ** P;
    double ** Q;
    double * LAMBDA;
    double mayor;
    int r_mayor,c_mayor;
    double theta;
    double C,S;

    /*Inicializamos las matrices y vectores*/
    /*Copiamos la matriz en M, por que este ultimo va a ser modificado*/
    M = copiarMatriz(matriz,dim,dim,0);
    Q = (double **)calloc(dim,sizeof(double *));
    P = (double **)calloc(dim,sizeof(double *));
    for(i=0;i<dim;i++){
        Q[i] = (double *)calloc(dim,sizeof(double));
        P[i] = (double *)calloc(dim,sizeof(double));
    }
    LAMBDA = (double *)calloc(dim,sizeof(double));


    /*Hacemos Q igual a la matriz identidad*/
    for(i=0;i<dim;i++){
        Q[i][i] = 1;
    }

    /*Inicializamos algunas variables*/
    k=0;
    mayor = tol+1;
    encontrarMayor(M,dim,&mayor,&r_mayor,&c_mayor);

    /*Ejecutamos el algoritmo con dos criterios de paro.*/
    while((k<numIter)&&(fabs(mayor) > tol)){
        if(M[r_mayor][r_mayor] == M[c_mayor][c_mayor]){
              theta = PI_constant/4;
        }
        else{
            theta = (0.5)*atan((2*M[r_mayor][c_mayor])/M[r_mayor][r_mayor]-M[c_mayor][c_mayor]);
        }
        C = cos(theta);
        S = sin(theta);

        /*Construimos la matriz de rotación. Empezamos iniciando P como la identidad*/
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
                if(i==j){
                    P[i][j] = 1;
                }
                else{
                    P[i][j] = 0;
                }
            }
        }
        /*Luego modificamos las cuatro posiciones de la matriz*/
        P[r_mayor][r_mayor] = C;
        P[c_mayor][c_mayor] = C;
        P[r_mayor][c_mayor] = -S;
        P[c_mayor][r_mayor] = S;

        /*En cada iteración multiplicamos Q por la nueva matriz P y lo guardamos en Q.*/
        Q = prodMatrices(Q,P,dim,dim,dim,dim,1,0);

        /*Actualizamos M como M <- Pt*M*P*/
        M = prodMatrices(M,P,dim,dim,dim,dim,1,0);
        P = matrizTransp(P,dim,dim,1);
        M = prodMatrices(P,M,dim,dim,dim,dim,0,1);

        /*Encontramos nuevamente el elemento mayor*/
        encontrarMayor(M,dim,&mayor,&r_mayor,&c_mayor);
        k++;
    }

    /*Guardamos el vector LAMBDA como la diagonal del último valor de M*/
    for(i=0;i<dim;i++){
        LAMBDA[i] = M[i][i];
    }

    /*Guardamos PHI como él último valor de Q*/
    for(i=0;i<dim;i++){
        PHI[i] = copiarVector(Q[i],dim,0);
    }

    /*Liberamos memoria*/
    for(i=0;i<dim;i++){
        free(M[i]);
        free(Q[i]);
        free(P[i]);
    }
    free(M);
    free(Q);
    free(P);

    /*Retornamos LAMBDA que es la matriz de Eigenvalores.*/
    return LAMBDA;

}

double metodoRaylegh(double ** matriz,int dim, double lambda_init,double * phi_init,int numIter, double tol, double * phi){
    /*Variables importantes*/
    double L0 = lambda_init;
    double L1;
    double * V0 = copiarVector(phi_init,dim,0);
    double * V1;
    /* Variables auxiliares */
    int i=0, j=0, k=0;
    double suma=0;
    double ** W;
    double numerador;
    double denominador;
    double dif=tol+1;
    double * MV0;
    /* Output */
    double lambda;

    i=0;
    while((i<numIter)&&(fabs(dif)>tol)){
        /*Construimos la matriz (M-L0*I)*/
        W = copiarMatriz(matriz,dim,dim,0);
        for(j=0;j<dim;j++){
            W[j][j] = W[j][j] - L0;
        }
        /*Obtenemos V1*/
        V1 = resolvElimGauss(W,V0,dim);
        V1 = normalizarVector(V1,dim,1);
        /*Obtenemos L1*/
        MV0 = prodMatVec(matriz,V0,dim,dim,dim);
        numerador = productoPunto(V0,MV0,dim);
        denominador = productoPunto(V0,V0,dim);
        L1 = numerador/denominador;
        /*Hacemos la diferencia entre el lambda anterior y el actual*/
        dif = L0-L1;
        /*Actualizamos V0*/
        free(V0);
        V0 = copiarVector(V1,dim,1);
        /*Actualizamos L0*/
        L0 = L1;
        /*Actualizamos el contador y liberamos memoria*/
        i++;
        for(j=0;j<dim;j++) free(W[j]); free(W);
        free(MV0);
    }

    /*Guardamos el valor final de V0 en phi*/
    for(i=0;i<dim;i++) phi[i] = V0[i];
    /*Guardamos el valor final de L0 en lambda*/
    lambda = L0;

    /*Liberamos memoria*/
    free(V0);

    /*Retornamos lambda*/
    return lambda;
}

double * metodoSubespacios(double ** matriz,int dim,int s,double ** PHI_init,int numIter, double tol,double ** PHI){
    //Nota imporante, inicializar PHI antes: PHI = crearMatriz(dim,s);
    /*Variables importantes*/
    double **R, *L;
    double ** PHI_T ;
    double ** PHI_;
    /*Parmetros de paro*/
    int numIterJacobi = 10000000;
    double tolJacobi = 0.001;
    /* Variables auxiliares */
    int i=0, j=0, k=0, l=0;
    double * V0, * V, a;
    double mayor=tol+1; int r,c;
    /*Output*/
    double ** Q;
    double * LAMBDA = crearVector(s);

    /*Inicializamos PHI_ como el PHI_init*/
    PHI_ = copiarMatriz(PHI_init,dim,s,0);

    i=0;
    while((i < numIter)&&(fabs(mayor) > tol)){
        /* Para el Gram Schmidt usamos la transpuesta de PHI por que es más facil
         * manejar renglones como apuntadores que columnas.
         */
        PHI_T = matrizTransp(PHI_,dim,s,0);
        for(j=0;j<s;j++){
            /*Inicializamos valores*/
            V0 = copiarVector(PHI_T[j],dim,0);
            V = copiarVector(PHI_T[j],dim,0);
            for(k=0;k<j;k++){
                /*Restamos las contribuciones*/
                a = productoPunto(PHI_T[k],V0,dim);
                V  = sumaVectores(V,multiplicaRealVectorCW(PHI_T[k],-a,dim,0),dim,1,1) ;
            }
            /*Resolvemos el sistema para hallar el nuevo valor de PHI_T[j]*/
            free(PHI_T[j]);
            PHI_T[j] = resolvElimGauss(matriz,V,dim);
            /*Normalizamos*/
            PHI_T[j] = normalizarVector(PHI_T[j],dim,1);
            /*Liberamos memoria*/
            free(V0); free(V);
        }

        /*Actualizamos PHI_ como la transpuesta del ya actualizado PHI_T*/
        for(l=0;l<dim;l++) free(PHI_[l]); free(PHI_);
        PHI_ = matrizTransp(PHI_T,s,dim,0);

        /*Realizamos el producto PHI_T*M*PHI = Q*/
        if(i>0){ for(l=0;l<s;l++) free(Q[l]); free(Q); }
        Q = prodMatrices(matriz,PHI_,dim,dim,dim,s,0,0);
        Q = prodMatrices(PHI_T,Q,s,dim,dim,s,1,1);

        /*Hallamos la matriz de rotación con el método de Jacobi*/
        R = (double **)calloc(s,sizeof(double*));
        L = metodoJacobiVPropios(Q,s,numIterJacobi,tolJacobi,R);

        /*Rotamos PHI_*/
        PHI_ = prodMatrices(PHI_,R,dim,s,s,s,1,0);
        for(l=0;l<s;l++) free(R[l]); free(R);

        /*Verificamos si Q es diagonal*/
        encontrarMayor(Q,s,&mayor,&r,&c);

        /*Actualizamos el contador y liberamos memoria*/
        i++ ; printf("\n%d\n",i);
        free(L);
    }

    /*Guardamos el valor de PHI_ en PHI*/
    for(i=0;i<dim;i++){
        for(j=0;j<s;j++){
            PHI[i][j] = PHI_[i][j];
        }
    }

    /*Guardamos el valor de la diagonal de Q en LAMBDA*/
    for(i=0;i<s;i++){
        LAMBDA[i] = Q[i][i];
    }

    /*Liberamos memoria*/
    for(l=0;l<dim;l++) free(PHI_[l]); free(PHI_);
    for(l=0;l<s;l++) free(Q[l]); free(Q);

    /*Retornamos*/
    return LAMBDA;
}

/*Algoritmo de factorizacion QR*/
double ** factorQR(double ** matriz,int dim,double ** Q){
    /*Variables*/
    double ** E;
    double ** R;
    double ** E_T;
    int i,j;

    /*Ortonormalizamos matriz y lo guardamos en E. Luego transponemos E.*/
    E = GramSchmidt(matriz,dim,0);
    E_T = matrizTransp(E,dim,dim,0);

    /*Hallamos R realizando el producto E_T*matriz.*/
    R = prodMatrices(E_T,matriz,dim,dim,dim,dim,0,0);

    /*Guardamos la matriz ortonormal E en Q*/
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            Q[i][j] = E[i][j];
        }
    }

    for(i=0;i<dim;i++) free(E[i]); free(E);
    for(i=0;i<dim;i++) free(E_T[i]); free(E_T);

    /*Retornamos R*/
    return R;
}

/*Algoritmo que halla los eigenvalores de una matriz dada.*/
double * metodoQR(double ** matriz,int dim, int numIter, double tol, double **PHI){
    /*Inicializar antes PHI con ceros.*/
    /*Variables importantes*/
    double **A;
    double **Q;
    double **R;
    double **_PHI;
    /* Variables auxiliares */
    int i=0, j=0,k=0;
    double mayor;
    /* Output */
    double * LAMBDA;

    /*Inicializamos las matrices y vectores*/
    LAMBDA = crearVector(dim);
    Q = crearMatriz(dim,dim);
    A = copiarMatriz(matriz,dim,dim,0);

    /*Inicializamos _PHI como la matriz identidad*/
    _PHI = crearMatriz(dim,dim);
    for(i=0;i<dim;i++) _PHI[i][i] = 1;

    /*Ejecutamos la recursion*/
    mayor = tol+1;
    i=0;
    while((i<numIter)&&(mayor > tol)){
        /*Factorizamos A en QR*/
        R = factorQR(A,dim,Q);
        /*Actualizamos _PHI*/
        _PHI = prodMatrices(_PHI,Q,dim,dim,dim,dim,1,0);
        /*Actualizamos A como RQ*/
        for(j=0;j<dim;j++) free(A[j]); free(A);
        A = prodMatrices(R,Q,dim,dim,dim,dim,1,0);

        /*Determinamos la mayor entrada fuera de la diagonal de A en valor absoluto.*/
        mayor=0;
        for(j=0;j<dim;j++){
            for(k=0;k<dim;k++){
                if((j!=k)&&(mayor<fabs(A[j][k]))){
                    mayor = fabs(A[j][k]);
                }
            }
        }

        i++;
    }
    printf("El número de iteraciones requerido es: %i\n",i);

    /*Guardamos LAMBDA como la diagonal de A*/
    for(i=0;i<dim;i++){
        LAMBDA[i] = A[i][i];
    }

    /*Guardamos los valores de _PHI en PHI*/
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            PHI[i][j] = _PHI[i][j];
        }
    }

    /*Liberamos la memoria*/
    for(i=0;i<dim;i++) free(_PHI[i]); free(_PHI);
    for(i=0;i<dim;i++) free(A[i]); free(A);
    for(i=0;i<dim;i++) free(Q[i]); free(Q);

    /*Retornamos los eigenvalores */
    return LAMBDA;
}



