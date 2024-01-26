#include "nmlib.h"

/*Método de Newton modificado.*/
double newton_method_modified(double(*f)(double),double(*f_der1)(double),double(*f_der2)(double),double x0, double tol_x,double tol_f){
    double x_old=x0; /*x_old = x_k  y x_new = x_k+1.*/
    double x_new = x_old+tol_x+1;
    double y=0;
    int k=1;
    printf("Iteración k |    x_k   |   |f(x_k)| \n");
    while( (fabs(x_new - x_old)> tol_x)&&(fabs(f(x_new))> tol_f) ){
        /*Se guarda x_old como el valor de la iteración anterior*/
        x_old = x_new;
        y = (f_der1(x_old)*f_der1(x_old))-(f_der2(x_old)*f(x_old)); /*Denominador*/
        if( y==0 ) return x_old;
        else{
            /*Se actualiza x_new*/
            x_new = x_old-(f_der1(x_old)*f(x_old))/y;
            printf("      %d     | %lf | %lf \n",k,x_new,fabs(f(x_new)));
        }
        k++;
    }
    printf("Una raíz usando el Metodo de Newton Modificado es: %lf\n",x_new);
}
