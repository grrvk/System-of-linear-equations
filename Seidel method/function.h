//
//  function.h
//  Seidel method
//
//  Created by Vika Granadzer on 16.03.2023.
//
using namespace std;
#include <cmath>

int N = 7;

void change_rows(double **P, int col_number, int row_number){
    double *pass = P[col_number-1];
    P[col_number-1] = P[row_number-1];
    P[row_number-1] = pass;
}

void change_cols(double **P,int n, int m)
{
    for(int i=0;i<N;i++){
        swap(P[i][m-1],P[i][n-1]); 
    }
}

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

void calc_q(double** A){
    for (int j = 0; j < N; j++){
        double summ = 0;
        for (int i = j+1; i < N; i++){
            summ += abs(A[i][j]);
        }
        cout << "q = " << A[j][j]/summ << endl;
    }
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
void change_M_matrix(double **M, double** A, int col_num){
    for (int i = col_num; i < N; i++){
        if (i != col_num){
            M[i][col_num] = -A[i][col_num]/A[col_num][col_num];
        }
        else{
            M[i][col_num] = 1/A[col_num][col_num];
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

void change_Identity_matrix(double **P, int col_number, int row_number){
    double *pass = P[col_number];
    P[col_number] = P[row_number];
    P[row_number] = pass;
}

double** mult(double **Augmented_matrix, double** P){
    double **res;
    res = new double *[N];
    for(int i = 0; i <N; i++)
        res[i] = new double[N+1];
    
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

double** create_B(){
    double **creation;
    creation = new double *[N];
    for(int i = 0; i <N; i++)
        creation[i] = new double[1];
    
    for (int i = 0; i < N; i++){
        creation[i][0] = 7;
    }
    
    return creation;
}


void summ_rows(double** A){
    for (int i = 2; i < N-1; i++){
        for (int k = 0; k < N+1; k++){
            A[i+1][k]+=A[i][k];
        }
    }
}

void dec_rows(double** A){
    for (int i = 3; i < N; i++){
        for (int k = 0; k < N+1; k++){
            A[i-1][k]-=A[i][k];
        }
    }
}

void summ_two_rows(double** A, int first, int second){
    for (int k = 0; k < N+1; k++){
        A[first][k]+=A[second][k];
    }
}

void summ_with_par_two_rows(double** A, int first, int second, double par){
    for (int k = 0; k < N+1; k++){
        A[first][k]+=A[second][k]*par;
    }
}

void dec_two_rows(double** A, int first, int second){
    for (int k = 0; k < N+1; k++){
        A[first][k]-=A[second][k];
    }
}

double** create_custom(int M){
    double **creation;
    creation = new double *[M];
    for(int i = 0; i <M; i++)
        creation[i] = new double[1];
    
    return creation;
}

double first_Summ(double** Augmented_matrix, double* x, int i){
    double sum = 0;
    for (int j = 0; j < i; j++){
        if (Augmented_matrix[i][j]!=0 && x[j]!=0){
            sum += (Augmented_matrix[i][j]*x[j])/Augmented_matrix[i][i];
        }
    }
    return (-1)*sum;
}

double second_Summ(double** Augmented_matrix, double* x, int i){
    double sum = 0;
    for (int j = i+1; j < N; j++){
        if (Augmented_matrix[i][j]!=0 && x[j]!=0){
            sum += (Augmented_matrix[i][j]*x[j])/Augmented_matrix[i][i];
        }
    }
    return (-1)*sum;
}

void calculate(double** Augmented_matrix, double* x){
    static double *new_x;
    new_x = new double [N];
    
    for (int i = 0; i < N; i++){
        double third_sum = 0;
        if (Augmented_matrix[i][N]==0){
            third_sum = 0;
        }
        else{
            third_sum = Augmented_matrix[i][N]/Augmented_matrix[i][i];
        }
        x[i] = first_Summ(Augmented_matrix, x, i) + second_Summ(Augmented_matrix, x, i) + third_sum;
    }
    
}

double get_norm(double* new_x, double *x){
    static double *norm_mas;
    norm_mas = new double [N];
    
    for (int i = 0; i < N; i++){
        norm_mas[i] = new_x[i] - x[i];
    }
    
    double norm = abs(norm_mas[0]);
    for (int i = 0; i < N; i++){
        if (abs(norm_mas[i]) > norm){
            norm = abs(norm_mas[i]);
        }
    }
    
    return norm;
}

void set_equal(double* x, double* previous){
    for (int i = 0; i < N; i++){
        previous[i] = x[i];
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

double **submatrix(double **matrix, int n, int x, int y) {
    double **submatrix = new double *[n - 1];
    int subi = 0;
    for (int i = 0; i < n; i++) {
        submatrix[subi] = new double[n - 1];
        int subj = 0;
        if (i == y) {
            continue;
        }
        for (int j = 0; j < n; j++) {
            if (j == x) {
                continue;
            }
            submatrix[subi][subj] = matrix[i][j];
            subj++;
        }
        subi++;
    }
    return submatrix;
}

double determinant(double **matrix, int n) {
    double det = 0;
    if (n==1){
        return matrix[0][0];
    }
    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
    }
    for (int x = 0; x < n; ++x) {
        if (matrix[0][x]!=0){
            det += ((x % 2 == 0 ? 1 : -1) * matrix[0][x] * determinant(submatrix(matrix, n, x, 0), n - 1));
        }
    }

    return det;
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
