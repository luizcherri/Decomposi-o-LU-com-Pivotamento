/// Decomposição LU com pivotamento
/**
    Entrada: A entrada é dada pelo arquivo matriz.txt cujo o formato é:

    ordem da matriz
    Matriz
    vetor b

    exemplo: seja o sistema
    |1 2| |x| =  |5|
    |3 4| |y| =  |6|
    então a entrana para o problema neste algoritmo eh
    2      <-- ordem da matriz
    1 2    <-- matriz
    3 3
    5      <-- vetor b
    6
**/


#include <stdio.h>
#include <math.h>

/// Este procedimento atribui valores nulos a uma matriz de ordem n
void limpa (int n,double **A){
     for (int i = 0 ; i<n ; i++){
            for (int j = 0 ; j<n ; j++){
                A[i][j] = 0;
            }
     }
}

/// Este procedimento imprime uma matriz de ordem n
void Impm (int n, double **A) {
     for (int i = 0 ; i<n ; i++){
            printf ("\n");
            for (int j = 0 ; j<n ; j++){
                printf ("%.3f\t", A[i][j]);
            }
        }
}

/// Este procedimento encontra a solução de um sistema de equaçoes cuja a matriz é triangular inferior
void SistTriInf(int n,double **A,double *b,double *y) {
     float soma;
     for (int i = 0 ; i < n ; i++) {
         soma = 0;
         for (int j = 0 ; j < i ; j++)
             soma = soma + y[j]*A[i][j];
         y[i]=(b[i]-soma)/A[i][i];
     }
}


/// Este procedimento encontra a solução de um sistema de equaçoes cuja a matriz é triangular superior
void SistTriSup(int n,double **A,double *y,double *x) {
     float soma;
     for (int i = n-1 ; i >= 0;i--) {
         soma=0;
         for (int j=i+1 ; j < n ; j++)
             soma = soma + x[j]*A[i][j];
         x[i]=(y[i] - soma)/A[i][i];
     }
}


/// Faz a decomposição LU com pivotamento de uma matriz A
void DecLU(int n,double **U,double *b,double *x){

     double **L = new double*[n];
     for (int i = 0 ; i < n ; i++){
         L[i] = new double[n];
     }
     limpa (n,L);

     double *y = new double[n];
     double aux = 0;
     int aux_pos;
     for (int i=0 ; i < n ; i++){
         // Pivotamento
         aux = U[i][i];
         aux_pos = i;
         for (int j=i+1 ; j<n; j++){
            if (U[j][i] > aux){
                aux = U[j][i];
                aux_pos = j;
            }
         }
         for (int k = i ; k<n ; k++){
             aux = U[i][k];
             U[i][k] = U[aux_pos][k];
             U[aux_pos][k] = aux;
         }
         aux = b[i];
         b[i] = b[aux_pos];
         b[aux_pos] = aux;

         // Eliminação
         for (int j = 0 ; j<i ; j++) {
             U[i][j]= 0;
         }
         for (int j=i ; j < n ; j++){
             L[i][i]=1;
             for (int k=0 ; k < i; k++) {
                 U[i][j]=U[i][j] - (L[i][k]*U[k][j]);
             }
         }

         for (int m=i+1 ; m<n ; m++){
             L[m][i] = U[m][i];
             for (int k=0 ; k < i ; k++) {
                 L[m][i]=U[m][i] - (L[m][k]*U[k][i]);
             }
             L[m][i]=L[m][i]/U[i][i];
         }
     }

     // Resolução do problema.

     printf ("\n\nA matriz L eh: \n");
     Impm (n, L);

     printf ("\n\n");
     printf ("A matriz U eh: \n");
     Impm (n, U);

     printf ("\n\n");

     SistTriInf(n,L,b,y);
     SistTriSup(n,U,y,x);

}



int main (){
     /// LEITURA DO PROBLEMA
     int n;         // Dimensão da matriz
     double aux;
     FILE *arq;
     if (!(arq = fopen ("matriz.txt" , "r"))){
        return (0);
     }
     fscanf (arq , "%d" , &n);

     double **A = new double*[n];
     for (int i = 0 ; i < n ; i++){
         A[i] = new double[n];
     }
     double *b = new double[n];

     for (int i =0 ; i < n ; i++){
           for (int j = 0 ; j < n ; j++){
               fscanf (arq , "%lf" , &aux);
               A[i][j] = aux;
           }
     }

     for (int i = 0 ; i<n ; i++){
           fscanf (arq , "%lf" , &aux);
           b[i] = aux;
     }
     fclose (arq);
     /// FIM DA LEITURA DO PROBLEMA

     /// Unidade de teste, imprime a matriz
     //     Impm ( n, A);
     //     printf ("\n\n");
     //     for (int i = 0 ; i<n ; i++){
     //           printf ("%f  " , b[i]);
     //     }

     double *x = new double[n]; // Vetor solução
     DecLU(n,A,b,x);            // Resolve o sistema Ax=b por decomposição LU com pivotamento


     /// Impresão da solução encontrada
     printf("\nSolucao encontrada:\n");
     for (int i = 0 ; i < n ; i++){
        printf ("%lf\n", x[i]);
     }
     /// Fim da impresão



}
