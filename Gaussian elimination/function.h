//
//  function.h
//  Gaussian elimination
//
//  Created by Vika Granadzer on 14.03.2023.
//
using namespace std;
#include <cmath>

#ifndef function_h
#define function_h

int N = 7;

void set_Identity_matrix(double **P){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                P[i][j] = 1;
            }
            else{
                P[i][j] = 0;
            }
        }
    }
}

void set_Null_matrix(double **mult){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N+1; j++){
            mult[i][j] = 0;
        }
    }
}

void change_Identity_matrix(double **P, int col_number, int row_number){
    double *pass = P[col_number];
    P[col_number] = P[row_number];
    P[row_number] = pass;
}

void change_M_matrix(double **M, double** A, int col_num){
    for (int i = col_num; i < N; i++){
        double num = A[i][col_num];
        double mid_num = A[col_num][col_num];
        if (i != col_num){
            if (num==0){
                M[i][col_num] = 0;
            }
            else{
                if (mid_num!=0){
                    M[i][col_num] = -num/mid_num;
                }
            }
        }
        else{
            if (mid_num!=0){
                M[i][col_num] = 1/mid_num;
            }
        }
    }
}

void findMax(double **A, int col_index, int& row_index, double& max){
    for (int i = col_index; i < N; i++){
        if (abs(A[i][col_index]) > max){
            max = A[i][col_index];
            row_index = i;
        }
    }
}

double** matrix_NxN1_creation(){
    double **creation;
    creation = new double *[N];
    for(int i = 0; i <N; i++)
        creation[i] = new double[N+1];
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N+1; j++){
            creation[i][j] = 0;
        }
    }
    
    return creation;
}

void print_E_matrix(double **result){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            cout << setw(20) << setprecision(15) << result[i][j];
        }
        cout << endl;
    }
}

void print_with_zeros(double **result){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N+1; j++){
            cout << setw(20) << setprecision(13) << result[i][j];
        }
        cout << endl;
    }
}

void print_Augmented_matrix(double **result){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N+1; j++){
            cout << setw(20) << setprecision(13) << result[i][j];
        }
        cout << endl;
    }
}

double** mult(double **Augmented_matrix, double** P){
    double **res;
    res = new double *[N];
    for(int i = 0; i <N; i++)
        res[i] = new double[N+1];
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N+1; j++){
            res[i][j] = 0;
        }
    }
    
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N+1; ++j)
            for(int k = 0; k < N; ++k)
            {
                if (P[i][k]!=0 && Augmented_matrix[k][j]!=0){

                    res[i][j] += P[i][k] * Augmented_matrix[k][j];
                }
            }
    return res;

}

double* residual(double **A, double *X, double **B){
    static double *res;
    res = new double [N];
    
    for(int i = 0; i < N; ++i)
        for(int k = 0; k < N; ++k)
              {
                  if (X[k]!=0 && A[i][k]!=0){
                      res[i] += A[i][k] * X[k];
                  }
              }
    
    for (int i = 0; i < N; i++){
        res[i] = res[i] - B[i][0];
    }
    
    return res;
}



double** matrix_NxN_creation(){
    double **creation;
    creation = new double *[N];
    for(int i = 0; i <N; i++)
        creation[i] = new double[N];
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            creation[i][j] = 0;
        }
    }
    
    return creation;
}

double* result_fill(double** Augmented_matrix){
    static double *X;
    X = new double [N];
    for (int i = N-1; i>=0; i--){
        double main = Augmented_matrix[i][N];
        X[i] = main;
        for (int m = i+1; m < N; m++){
            if (Augmented_matrix[i][m]!=0){
                X[i] =  main - Augmented_matrix[i][m]*X[m];
                main -= Augmented_matrix[i][m]*X[m];
            }
        }
    }
    return X;
}

void result_print(double* X){
    cout << endl << "(";
    for (int i = 0; i < N-1; i++){
        cout << setprecision(13) << X[i] << ", ";
    }
    cout << setprecision(13) << X[N-1] << ")" << endl;
}

double** inverse_matrix(double **A, double **E){
    double** X = matrix_NxN_creation();
    for (int j = 0; j < N; j++){
        X[N-1][j] = E[N-1][j];
    }
    
    for (int i = N-2; i >= 0; i--){
        for (int j = 0; j < N; j++){
            double main = E[i][j];
            for (int n = i; n < N; n++){
                if (A[i][n]!= 0 && X[n][j]!=0){
                    X[i][j] = main - A[i][n]*X[n][j];
                    main -= A[i][n]*X[n][j];
                }
            }
        }
    }
    
    return X;
}

double determinant(double** P){
    double det = (-1)*P[0][0];
    for (int  i =1; i < N; i++){
        det = P[i][i]* det;
    }
    return det;
}

double condition_num(double** A, double **A_rev){
    double first_norm_mas[N];
    for (int i = 0; i < N; i++){
        first_norm_mas[i] = 0;
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            first_norm_mas[i] += abs(A[i][j]);
        }
    }
    double first_norm = first_norm_mas[0];
    for (int i=0; i <N; i++){
        if (first_norm_mas[i] > first_norm){
            first_norm = first_norm_mas[i];
        }
    }
    
    double second_norm_mas[N];
    for (int i = 0; i < N; i++){
        second_norm_mas[i] = 0;
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            second_norm_mas[i] += abs(A_rev[i][j]);
        }
    }
    double second_norm = second_norm_mas[0];
    for (int i=0; i <N; i++){
        if (second_norm_mas[i] > second_norm){
            second_norm = second_norm_mas[i];
        }
    }
    
    return first_norm*second_norm;
}

double** fill_Augmented_matrix(double** A, double** B){
    double **creation = matrix_NxN1_creation();
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            creation[i][j] = A[i][j];
        }
        creation[i][N] = B[i][0];
    }
    return creation;
}

double** create_A(){
    double** creation = matrix_NxN_creation();
    
    for (int i = 1; i <= N; i++){
        for (int j = 1; j <= N; j++){
            creation[i-1][j-1] = 10000/(2*pow(i*j, 2) + pow(i+j, 3));
        }
    }
    return creation;
}

double** create_B(){
    double **creation;
    creation = new double *[N];
    for(int i = 0; i <N; i++)
        creation[i] = new double[1];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            creation[i][0] = 0;
        }
    }
    
    for (int i = 0; i < N; i++){
        creation[i][0] = 7;
    }
    
    return creation;
}


#endif /* function_h */
