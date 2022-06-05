#include <iostream>
#include <string>

using namespace std;

void laplaceMatrix(double** matrix, double** temp, int exRow, int exCol, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) 
    {
        for (int col = 0; col < n; col++) 
        {
            if (row != exRow && col != exCol) 
            {
                temp[i][j++] = matrix[row][col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(double** matrix, int n) {
    double det = 0;
    if (n == 1)
        return matrix[0][0];

    if (n == 2)
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);

    double** temp = new double* [n - 1];

    for (int i = 0; i < n - 1; i++)
    {
        temp[i] = new double[n - 1];
    }

	int sign = 1;

    for (int i = 0; i < n; i++) {
        laplaceMatrix(matrix, temp, 0, i, n);
        det += sign * matrix[0][i] * determinant(temp, n - 1);
        sign = -sign;
    }

    delete[] temp;
    return det;
}

double* GaussianElimination(double** matrix, int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double f = matrix[j][i] / matrix[i][i];
            for (int k = 0; k < n + 1; k++)
            {
                matrix[j][k] = matrix[j][k] - f * matrix[i][k];
            }
        }
    }

    double* res = new double[n];

    for (int i = n - 1; i >= 0; i--)
    {
        res[i] = matrix[i][n];

        for (int j = i + 1; j < n; j++)
        {
            if (i != j)
            {
                res[i] = res[i] - matrix[i][j] * res[j];
            }
        }
        res[i] = res[i] / matrix[i][i];
    }
    return res;
}

int main()
{
    /*/
    int n = 3;
    double** tab = new double*[n];
    double** tempMatrix = new double* [n];
    for (int i = 0; i < n; i++)
    {
        tab[i] = new double[n];
        tempMatrix[i] = new double[n];
    }
    

    tab[0][0] = 5;
    tab[0][1] = -2;
    tab[0][2] = 3;
    
    tab[1][0] = -2;
    tab[1][1] = 3;
    tab[1][2] = 1;
    
    tab[2][0] = -1;
    tab[2][1] = 2;
    tab[2][2] = 3;

    //tab[0][0] = 1;
    //tab[0][1] = 7;
    //tab[0][2] = 3;
    //tab[0][3] = 5;

    //tab[1][0] = 8;
    //tab[1][1] = 4;
    //tab[1][2] = 6;
    //tab[1][3] = 2;

    //tab[2][0] = 2;
    //tab[2][1] = 6;
    //tab[2][2] = 4;
    //tab[2][3] = 8;

    //tab[3][0] = 5;
    //tab[3][1] = 3;
    //tab[3][2] = 7;
    //tab[3][3] = 1;

    double* outcome = new double[n];

    outcome[0] = 21;
    outcome[1] = -4;
    outcome[2] = 5;

    //outcome[0] = 16;
    //outcome[1] = -16;
    //outcome[2] = 16;
    //outcome[3] = -16;


    double* dets = new double[n];

    for (int k = 0; k < n; k++) 
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (k == i)
                {
                    tempMatrix[j][i] = outcome[j];
                }
                else 
				{
                    tempMatrix[j][i] = tab[j][i];
                }
            }
        }
        dets[k] = determinant(tempMatrix, n);
    }

    double res = determinant(tab, n);
    //cout << res << endl;

    for (int i = 0; i < n; i++)
    {
        outcome[i] = dets[i] / res;
        cout << outcome[i] << endl;
    }

    for (int i = 0; i < n; i++)
    {
        delete[] tab[i];
        delete[] tempMatrix[i];
    }
    delete[] tab;
    delete[] tempMatrix;
    delete[] outcome;
    delete[] dets;
    */

    int n = 3;

    double** tab = new double* [n];
    for (int i = 0; i < n; i++)
    {
        tab[i] = new double[n + 1];
    }

    //tab[0][0] = 5;
    //tab[0][1] = -2;
    //tab[0][2] = 3;


    //tab[1][0] = -2;
    //tab[1][1] = 3;
    //tab[1][2] = 1;


    //tab[2][0] = -1;
    //tab[2][1] = 2;
    //tab[2][2] = 3;

    //tab[0][3] = 21;
    //tab[1][3] = -4;
    //tab[2][3] = 5;

    //tab[0][0] = 2;
    //tab[0][1] = 2;
    //tab[0][2] = 8;

    //tab[1][0] = 2;
    //tab[1][1] = 1;
    //tab[1][2] = -2;


    double* res = GaussianElimination(tab, n);
    double r = determinant(tab, n);
    cout << r << endl;

    for (int i = 0; i < n; i++)
		cout << res[i] << endl;


    for (int i = 0; i < n; i++)
    {
        delete[] tab[i];
    }
    delete[] tab;
    delete[] res;
}

