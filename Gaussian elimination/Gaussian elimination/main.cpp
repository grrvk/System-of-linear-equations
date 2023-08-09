//
//  main.cpp
//  Gaussian elimination
//
//  Created by Vika Granadzer on 10.03.2023.
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include "function.h"

using namespace std;

int main(int argc, const char * argv[]) {
    
    double **A = create_A();
    double **B = create_B();
    
    
    static double *x;
    x = new double [N];
    
    static double *r;
    r = new double [N];
    
    double **Augmented_matrix = matrix_NxN1_creation();
    
    double **E = matrix_NxN_creation();
    set_Identity_matrix(E);
    
    double **A_inverse = matrix_NxN_creation();
    
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Augmented_matrix[i][j] = A[i][j];
        }
        Augmented_matrix[i][N] = B[i][0];
    }
    
    double **P = matrix_NxN_creation();
    double **P_global = matrix_NxN_creation();
    
    double **M = matrix_NxN_creation();
    for (int i = 0; i<N; i++){
        cout << setw(17) << "A[" << i << "]";
    }
    cout << setw(17) << "B" << endl;
    print_Augmented_matrix(Augmented_matrix);
    cout << endl << "Процес" << endl;
    
    for (int i = 0; i<N; i++){
        cout << setw(17) << "A[" << i << "]";
    }
    cout << setw(17) << "B" << endl;
    
    for (int i = 0; i < N; i++){
        set_Identity_matrix(P);
        set_Identity_matrix(M);
        int row_index = 0;
        double max = 0;
        
        findMax(Augmented_matrix, i, row_index, max);
        change_Identity_matrix(P, i, row_index);
        
        Augmented_matrix = mult(Augmented_matrix, P);
        E = mult(E, P);
        P_global[i][i] = Augmented_matrix[i][i];
        
        change_M_matrix(M, Augmented_matrix, i);
        
        Augmented_matrix = mult(Augmented_matrix, M);
        E = mult(E, M);
        
        print_Augmented_matrix(Augmented_matrix);
        cout << endl;
    }
    
    cout << endl << "Розвʼязок системи";
    x = result_fill(Augmented_matrix);
    result_print(x);
    
    cout << endl << "Невʼязка";
    r = residual(A, x, B);
    result_print(r);
    
    cout << endl << "Визначник";
    double det = determinant(P_global);
    cout << endl << det << endl;
    
    cout <<endl << "Обернена матриця" << endl;
    A_inverse = inverse_matrix(Augmented_matrix, E);
    print_E_matrix(A_inverse);
    
    cout << endl;
    
    cout << "Число обумовленості" << endl;
    double condition_n = condition_num(A, A_inverse);
    cout << setprecision(13) << condition_n << endl;
    
    cout << endl << "Добуток А*А^(-1)" << endl;
    double **m = matrix_NxN_creation();
    m = mult(A_inverse, A);
    print_with_zeros(m);
    
    return 0;
}
