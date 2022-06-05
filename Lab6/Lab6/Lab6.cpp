#include <chrono>
#include <iostream>
#include <random>

#include "Matrix.h"

using namespace std;
using namespace std::chrono;

#pragma region Wolne

double **createMatrix(int rows, int columns)
{
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++)
        matrix[i] = new double[columns];
    return matrix;
}


void printMatrix(double ** matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
            cout << matrix[i][j] << "\t";
        cout << endl;
    }
    cout << endl;
}

void matrixCopy(double** nmatrix, double** matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++) {
        memcpy_s(nmatrix[i], columns * sizeof(matrix[0][0]), matrix[i], columns * sizeof(matrix[0][0]));
    }
}

void fillMatrix(double** matrix, int rows, int columns, double value)
{
    for (int i = 0; i < rows; i++)
        fill_n(matrix[i], columns, value);
}

void fillIdentity(double **matrix, int rows, int columns)
{
    if (rows != columns) return;
	for (int i = 0; i < rows; i++)
	{
        fill_n(matrix[i], columns, 0);
        matrix[i][i] = 1;
	}
}

void subMatrixCopy(double **nmatrix, double** matrix, int startRow, int startCol, int endRow, int endCol)
{
	int rows = endRow - startRow;
    int cols = endCol - startCol;
    for (int i = startRow, j = 0; i < endRow; i++, j++) {
        //for (int k = 0, l = startCol; k < cols; k++, l++)
        memcpy_s(nmatrix[j], cols * sizeof(matrix[0][0]), matrix[i] + startCol, cols * sizeof(matrix[0][0]));
            //nmatrix[i][k] = matrix[j][l];
    }
}

void dealloc(double** matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        delete[] matrix[i];
    delete[] matrix;
}

#pragma endregion 


int main()
{
    std::random_device rd; // obtain a random number from hardware
    std::mt19937_64 gen(rd()); // seed the generator

    //srand(time(nullptr));
    /*
    int rows = 5;
    int columns = 4;
    double** A = createMatrix(rows, columns);
    double** L = createMatrix(rows - 1, columns);

    for (int i = 0; i < rows; i++)
    {
	    for (int j = 0; j < columns; j++)
	    {
            A[i][j] = rand() % 100;
	    }
    }

    //fillIdentity(L, rows - 1, columns);
    fillMatrix(L, rows - 1, columns, 0);

    printMatrix(A, rows, columns);
    printMatrix(L, rows - 1, columns);

    double** C = createMatrix(rows, columns);
    matrixCopy(C, A, rows, columns);

    int startRow = 1, endRow = 4;
    int startCol = 2, endCol = 4;

    int nrows = endRow - startRow;
    int ncols = endCol - startCol;

    double** D = createMatrix(nrows, ncols);
	subMatrixCopy(D, A, startRow, startCol, endRow, endCol);

    dealloc(A, rows);
    dealloc(L, rows - 1);

    printMatrix(C, rows, columns);
    printMatrix(D, nrows, ncols);

    dealloc(C, rows);
    dealloc(D, nrows);
    */

    try
    {
        //int rows = 5;
        //int cols = 4;

        //Matrix<double> m(rows, cols);
        //m.fillValue(5);
        ////m.fillIdentity();
        //m.printMatrix();

        //Matrix<double> m1(m);
        //m1.printMatrix();

        //Matrix<double> m2(5, 3);

        //m2 = m1;

        //m2.fillValue(12);
        //m1.printMatrix();
        //m2.printMatrix();

        //Matrix<double> m3;

        //m.fillRandom(0, 100, gen);
        //m.printMatrix();

        //m3 = m.subMatrix(1, 2, 2, 3);

        //m3.printMatrix();

        //Matrix<double> m4 = m3 * 2;

        //m4.printMatrix();

        //Matrix<double> m5 = m4 + m3;

        //m5.printMatrix();

        //Matrix<double> m6 = m.transpose();

        //m6.printMatrix();

        //Matrix<double> m(1000, 1000);
        //m.fillRandom(0, 10, gen);

        //gen.seed(rd());

        //Matrix<double> m2(1000, 1000);
        //m2.fillRandom(0, 10, gen);

        //clock_t t1 = clock();
        //Matrix<double> res = m * m2;
        //cout << (clock() - t1) / CLOCKS_PER_SEC << endl;

        //m.printMatrix();
        //m2.printMatrix();
        //res.printMatrix();

        //Matrix<double> m(4);

        
        //m[0][0] = 1;
        //m[0][1] = 1;
        //m[0][2] = 1;
        //m[1][0] = 4;
        //m[1][1] = 3;
        //m[1][2] = -1;
        //m[2][0] = 3;
        //m[2][1] = 5;
        //m[2][2] = 3;

        //m[0][0] = 0;
        //m[0][1] = 3;
        //m[0][2] = 7;
        //m[0][3] = 4;

        //m[1][0] = 9;
        //m[1][1] = 2;
        //m[1][2] = 2;
        //m[1][3] = 1;

        //m[2][0] = 3;
        //m[2][1] = 6;
        //m[2][2] = 2;
        //m[2][3] = 8;

        //m[3][0] = 9;
        //m[3][1] = 4;
        //m[3][2] = -2;
        //m[3][3] = -1;

        //m[0][0] = 5;
        //m[0][1] = 3;
        //m[0][2] = 7;
        //m[0][3] = 4;

        //m[1][0] = 2;
        //m[1][1] = 9;
        //m[1][2] = 2;
        //m[1][3] = 2;

        //m[2][0] = 1;
        //m[2][1] = 1;
        //m[2][2] = 3;
        //m[2][3] = 6;

        //m[3][0] = 2;
        //m[3][1] = 8;
        //m[3][2] = 9;
        //m[3][3] = 9;

        //Matrix<int> b(3, 1);

        //b[0][0] = 1;
        //b[1][0] = 6;
        //b[2][0] = 4;


        //m.fillRandom(0, 100, gen);

        //Matrix<double> upper;
        //Matrix<double> lower;

        //clock_t t1 = clock();
        ////cout << (clock() - t1) / CLOCKS_PER_SEC << endl;

        //m.decompose(upper, lower);

        //upper.printMatrix();
        //lower.printMatrix();

        //cout << m.det() << endl << endl;

        //m.decompose();

        //m.pivotDecompose();

        //m.printMatrix();

        //Matrix<double> res = m.solve(b);

        //res.printMatrix();


		//Matrix<double> m(3);

  //      m[0][0] = 5;
  //      m[0][1] = -2;
  //      m[0][2] = 3;
  //      m[1][0] = -2;
  //      m[1][1] = 3;
  //      m[1][2] = 1;
  //      m[2][0] = -1;
  //      m[2][1] = 2;
  //      m[2][2] = 3;


  //      Matrix<double> b(3, 1);

  //      b[0][0] = 21;
  //      b[1][0] = -4;
  //      b[2][0] = 5;

  //      m.printMatrix();

        //Matrix<double> result = m.solve(b);
        //result.printMatrix();

        //double d = m.det();

        //m[0][0] = 21;
        //m[1][0] = -4;
        //m[2][0] = 5;

        //double x1 = m.det();

        //m[0][0] = 5;
        //m[1][0] = -2;
        //m[2][0] = -1;

        //m[0][1] = 21;
        //m[1][1] = -4;
        //m[2][1] = 5;

        //double x2 = m.det();

        //m[0][1] = -2;
        //m[1][1] = 3;
        //m[2][1] = 2;

        //m[0][2] = 21;
        //m[1][2] = -4;
        //m[2][2] = 5;

        //double x3 = m.det();

        //cout << x1 / d << " " << x2 / d << " " << x3 / d << endl;

        //auto start = high_resolution_clock::now();
        //Matrix<double> res =  m.Cramer(b);
        //auto stop = high_resolution_clock::now();
        //auto duration = duration_cast<microseconds>(stop - start);
        //res.printMatrix();
        //cout << duration.count() << " microseconds" << endl << endl;


        //auto start1 = high_resolution_clock::now();
        //Matrix<double> res1 = m.Cramer2(b);
        //auto stop1 = high_resolution_clock::now();
        //auto duration1 = duration_cast<microseconds>(stop1 - start1);
        //res1.printMatrix();
        //cout << duration1.count() << " microseconds" << endl << endl;

    	Matrix<double> D(5);


		D[0][0] = 5;
		D[0][1] = 3;
		D[0][2] = 7;
		D[0][3] = 4;
		D[0][4] = 2;

		D[1][0] = 9;
		D[1][1] = 2;
		D[1][2] = 2;
		D[1][3] = 1;
		D[1][4] = 1;

		D[2][0] = 3;
		D[2][1] = 6;
		D[2][2] = 2;
		D[2][3] = 9;
		D[2][4] = 8;

		D[3][0] = 9;
		D[3][1] = 4;
		D[3][2] = -2;
		D[3][3] = -1;
		D[3][4] = -3;

		D[4][0] = 0;
		D[4][1] = 5;
		D[4][2] = 3;
		D[4][3] = -6;
		D[4][4] = -11;

        Matrix<double> L, U;

        D.decompose(U, L);

        cout << "Zadanie 4" << endl;
        L.printMatrix();
        U.printMatrix();
    }
    catch (exception& e)
    {
        cout << e.what() << endl;
    }

}
