//
//  main.cpp
//  Seidel method
//
//  Created by Vika Granadzer on 16.03.2023.
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
    for (int i = 0; i < N; i++){
        x[i] = 0;
    }
    
    static double *r;
    r = new double [N];
    
    static double *previous_x;
    previous_x = new double [N];
    for (int i = 0; i < N; i++){
        previous_x[i] = 0;
    }
    
    double **Augmented_matrix = matrix_NxN1_creation();
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Augmented_matrix[i][j] = A[i][j];
        }
        Augmented_matrix[i][N] = B[i][0];
    }
    cout << "Початкова розширена матриця" << endl;
    for (int i = 0; i<N; i++){
        cout << setw(17) << "A[" << i << "]";
    }
    cout << setw(17) << "B" << endl;
    print_Augmented_matrix(Augmented_matrix);
    cout << endl;
    
    cout << "Головні мінори початкової розширеної матриці" << endl;
    for (int i = 0; i < N; i++){
        cout << "det " << i+1 << "x" << i+1 << " " << setw(20) << setprecision(13) <<
        determinant(Augmented_matrix, i+1) << endl;
    }
    cout << endl;
    int k = 1;
    
    double **P = matrix_NxN_creation();
    
    double **M = matrix_NxN_creation();
    
    double **E = matrix_NxN_creation();
    set_Identity_matrix(E);
    
    for (int i = 0; i < N; i++){
        set_Identity_matrix(P);
        set_Identity_matrix(M);
        int row_index = 0;
        double max = 0;
        
        findMax(Augmented_matrix, i, row_index, max);
        change_Identity_matrix(P, i, row_index);
        
        Augmented_matrix = mult(Augmented_matrix, P);
        E = mult(E, P);
        
        change_M_matrix(M, Augmented_matrix, i);
        
        Augmented_matrix = mult(Augmented_matrix, M);
        E = mult(E, M);
    }
    cout << "Матриця після модифікації" << endl;
    for (int i = 0; i<N; i++){
        cout << setw(17) << "A[" << i << "]";
    }
    cout << setw(17) << "B" << endl;
    print_Augmented_matrix(Augmented_matrix);
    cout << endl;
    dec_rows(Augmented_matrix);
    summ_two_rows(Augmented_matrix, 2, 4);
    dec_two_rows(Augmented_matrix, 3, 4);
    dec_two_rows(Augmented_matrix, 5, 6);
    dec_two_rows(Augmented_matrix, 4, 5);
    summ_with_par_two_rows(Augmented_matrix, 3, 5, 0.61);
    summ_with_par_two_rows(Augmented_matrix, 2, 6, 1.5);
    summ_with_par_two_rows(Augmented_matrix, 4, 6, 0.34);
    print_Augmented_matrix(Augmented_matrix);
    cout << endl;
    
    cout << "Головні мінори модифікованої матриці" << endl;
    for (int i = 0; i < N; i++){
        cout << "det " << i+1 << "x" << i+1 << " " << setw(20) << setprecision(13) <<
        determinant(Augmented_matrix, i+1) << endl;
    }
    cout << endl;
    
    cout << "Ітераційний процес" << endl << endl;
    do{
        set_equal(x, previous_x);
        calculate(Augmented_matrix, x);
        
        cout << "Ітерація " << k << " ";
        for (int i = 0; i < N; i++){
            cout << setw(20) << setprecision(13) << x[i] << " ";
        }
        cout << endl;
        k++;
    }
    while (get_norm(x, previous_x) > pow(10, -6));
    
    cout << endl << "Розвʼязок системи";
    result_print(x);
    
    cout << endl << "Невʼязка";
    r = residual(A, x, B);
    result_print(r);
    
    //get_norm(x, previous_x) < 0.5
    return 0;
}

