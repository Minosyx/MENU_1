#include <iomanip>
#include <iostream>

using namespace std;

struct Data
{
	int x, y;
};

double interpolate(Data f[], int xi, int n)
{
	double res = 0;

	for (int i = 0; i < n; i++)
	{
		double term = f[i].y;

		for (int j = 0; j < n; j++)
		{
			if (j != i)
			{
				term = term * (xi - f[j].x) / static_cast<double>(f[i].x - f[j].x);
			}
		}
		res += term;
	}
	return res;
}

double findPoly(double x[], double y[], int n, double x0)
{
	double res = 0;
	for (int i = 0; i < n ;i++)
	{
		double val = y[i];
		for (int j = 0; j < n; j++)
		{
			if (j != i)
			{
				val = val * (x0 - x[j]) / (x[i] - x[j]);
			}
		}
		res += val;
	}
	return res;
}

double NewtonInterpolation(double *x, double **y, int n, double x0)
{
	for(int i = 1; i < n; i++)
	{
		for (int j = 0; j < n - i; j++)
		{
			y[j][i] = (y[j + 1][i - 1] - y[j][i - 1]) / (x[j + i] - x[j]);
		}
	}

	for (int i = 0; i < n; i++)
	{
		cout << x[i];
		for (int j = 0; j < n - i; j++)
		{
			cout << "\t" << y[i][j];
		}
		cout << endl;
	}

	double sum = y[0][0];
	for (int i = 1; i < n; i++)
	{
		double part = y[0][i];
		for (int j = 0; j <= i - 1; j++)
		{
			part *= x0 - x[j];
		}
		sum += part;
	}

	return sum;
}

int main()
{
	//Data f[] = { {-2,0}, {-1,1}, {0,1}, {2,2} };
	//cout << "Value of f(2) is : " << interpolate(f, 2, 3) << endl;

	//double x[] = { -2,-1,0,2 };
	//double y[] = { 0,1,1,2 };
	//cout << "Value of f(1) is : " << findPoly(x, y, 3, 2) << endl;

	int n = 4;
	double x[] = {100,121,144,169};
	double** y = new double* [n];
	for (int i = 0; i < n ; i++)
	{
		y[i] = new double[n];
	}

	y[0][0] = 10;
	y[1][0] = 11;
	y[2][0] = 12;
	y[3][0] = 13;

	double x0 = 122;
	cout << "Value of f(" << x0 << ") : " << NewtonInterpolation(x, y, 4, x0) << endl;

	for (int i = 0; i < n; i++)
	{
		delete[] y[i];
	}
	delete[] y;
}

