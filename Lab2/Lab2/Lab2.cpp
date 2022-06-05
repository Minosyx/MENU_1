#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

using Funkcja = double(*)(double);

double bisection(Funkcja f, double a, double b, double eps)
{
	//double eps = 1e-5;
	if (f(a) * f(b) >= 0)
	{
		throw "There are possibly no roots!";
	}
	int it = 0;

	double c = a;

	double temp = c;

	while (fabs(b - a) >= eps)
	{
		it++;
		c = (a + b) / 2;
		double funcVal = f(c);
		if (fabs(funcVal) < eps)
			break;
		if (funcVal * f(a) < 0)
			b = c;
		else
			a = c;

		cout << "eps[" << it << "]: " << abs(c - temp) << endl;
		temp = c;
	}

	cout << "Bisekcja - Liczba iteracji: " << it << endl;
	return c;
}

double funct(double x) {
	return x - cbrt(x) - 2;
}

double functImproved(double x)
{
	return cbrt(x) + 2;
}

double f1(double x1, double x2, double x3)
{
	return 3 * x1 + 2 * x2 - x3;
}

double f2(double x1, double x2, double x3)
{
	return 2 * x1 - 2 * x2 + 4 * x3;
}

double f3(double x1, double x2, double x3)
{
	return -x1 + 0.5 * x2 - x3;
}

double test(double x)
{
	return x * x * x + 4 * x * x - 3 * x + 9;
}

double derivative(Funkcja f, double x)
{
	double h = 1e-5;
	return (f(x + h) - f(x - h)) / 2 / h;
}

double Fderivative(Funkcja f, double x)
{
	double h = 1e-5;
	return (f(x + h) - f(x)) / h;
}

double Bderivative(Funkcja f, double x)
{
	double h = 1e-5;
	return (f(x) - f(x - h)) / h;
}

double Newton(Funkcja f, double x)
{
	double eps = 1e-5;
	double check, h;
	do
	{
		h = f(x) / derivative(f, x);
		x -= h;
		check = f(x);
	} while (fabs(check) >= eps && fabs(h) >= eps);
	return x;
}

double NewtonRaphson(Funkcja f, double x, double eps)
{
	//double eps = 1e-5;
	double h = f(x) / derivative(f, x);
	int it = 0;
	while (fabs(h) >= eps)
	{
		it++;
		h = f(x) / derivative(f, x);
		x -= h;
		cout << "eps[" << it << "]: " << fabs(h) << endl;
	}
	cout << "Newton - Liczba iteracji: " << it << endl;
	return x;
}

double sieczne(Funkcja f, double x1, double x2, double eps)
{
	//double eps = 1e-8;
	double f1, f2, f3, x3;
	f1 = f(x1);
	f2 = f(x2);
	int it = 0;

	double temp = 0;

	while (fabs(x1 - x2) > eps)
	{
		it++;
		if (fabs(f1 - f2) < eps)
		{
			throw "Wybrano zle punkty poczatkowe";
		}
		x3 = x2 - f2 * (x2 - x1) / (f2 - f1);
		f3 = f(x3);

		cout << "eps[" << it << "]: " << fabs(x3 - temp) << endl;
		temp = x3;

		if (fabs(f3) <= eps) break;
		x1 = x2;
		f1 = f2;
		x2 = x3;
		f2 = f3;
	}

	cout << "Sieczne - Liczba iteracji: " << it << endl;
	return x3;
}

double simpleIt(Funkcja f, double x0)
{
	double eps = 1e-5, maxIt = 10000, i = 0;
	double x1;
	while(i <= maxIt)
	{
		x1 = f(x0);
		if (fabs(x0 - x1) < eps)
		{
			return x1;
		}

		i++;
		x0 = x1;

		if (i > maxIt)
		{
			throw "Nie bylo mozliwe znalezienie rozwiazania w zadanej liczbie krokow";
		}
	}
	return x1;
}

double regulaFalsi(Funkcja f, double a, double b)
{
	double eps = 1e-5, maxIt = 5000;
	if (f(a) * f(b) >= 0)
	{
		throw "You have not assumed right a and b";
	}

	double c = a;

	for (int i = 0; i < maxIt; i++)
	{
		c = (a * f(b) - b * f(a)) / (f(b) - f(a));

		if (fabs(f(c)) < eps)
			break;

		if (f(c) * f(a) < 0)
			b = c;
		else
			a = c;
	}
	return c;
}

double ff(double x)
{
	return x*x*x * (x + sin(x*x - 1) - 1) - 3;
}

double fff(double x)
{
	return -2 * (x - 1) * (x - 1) * (x - 1) * (x - 1) + x * x * x * (x + sin(0.25 * x * x + 1));
}

double newtonTest(double x)
{
	return x * exp(-x);
}

using FSystem = double(*)(double*);

double partialDerivative(FSystem f, double *x, const int index, int dim)
{
	double h = 1e-8;
	double *fh = new double[dim];
	double *fn = new double[dim];
	for (int i = 0; i < dim; i++)
	{
		if (i == index)
			fh[i] = x[i] + h;
		else fh[i] = x[i];
	}
	for (int i = 0; i < dim; i++)
	{
		if (i == index)
			fn[i] = x[i] - h;
		else fn[i] = x[i];
	}

	double res = (f(fh) - f(fn)) / 2 / h;
	delete[] fh;
	delete[] fn;
	return res;
}

double *GaussianElimination(double **matrix, int n)
{
	for(int i = 0; i < n - 1; i++)
    {
        for(int j = i + 1; j < n; j++)
        {
            double f = matrix[j][i] / matrix[i][i];
            for(int k = 0; k < n + 1; k++)
            {
            	matrix[j][k] = matrix[j][k] - f * matrix[i][k];
            }
        }
    }

	double *res = new double[n];

    for(int i = n - 1; i >= 0; i--)          
    {                     
        res[i] = matrix[i][n];
                    
        for(int j = i + 1; j < n; j++)
        {
        	if(i!=j)
        	{
        	    res[i] = res[i] - matrix[i][j] * res[j];
        	}          
        }
		res[i] = res[i] / matrix[i][i];  
	}
	return res;
}

double *multiDimensionalEq(FSystem *f, double *x, int size)
{
	double eps = 1e-8;
	double h = 1e-8;
	double *funcs = new double[size];
	double **jacob = new double*[size];
	double len = 0;
	double *res;
	for (int i = 0; i < size; i++)
	{
		jacob[i] = new double[size + 1];
	}

	do {
		for (int i = 0; i < size; i++)
		{
			funcs[i] = f[i](x);
			for (int j = 0; j < size + 1; j++)
			{
				if (j < size)
					jacob[i][j] = partialDerivative(f[i], x, j, size);
				else jacob[i][j] = -funcs[i];
			}
		}

		len = 0;
		res = GaussianElimination(jacob, size);

		for (int i = 0; i < size; i++)
		{
			x[i] += res[i];
			len += res[i] * res[i];
		}

		len = sqrt(len);
	} while (len > eps);

	delete[] funcs;
	for (int i = 0; i < size; i++)
	{
		delete[] jacob[i];
	}
	delete[] jacob;
	delete[] res;

	return x;
}

/// <summary>
/// Set number 1
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s1(double *x)
{
	return x[0] * x[0] + x[1] * x[1] - 10; 
}

double s2(double *x)
{
	return 2 * x[0] + x[1] - 1;
}


/// <summary>
/// Set number 2
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s3(double *x)
{
	return 2*x[0]*x[0] + x[1]*x[1] - 24;
}

double s4(double *x)
{
	return x[0] * x[0] - x[1] * x[1] + 12;
}

/// <summary>
/// Set number 3
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s5(double *x)
{
	return 2 * x[0] + x[1] + 2 * x[2] * x[2] - 5;
}

double s6(double *x)
{
	return x[1] * x[1] * x[1] + 4 * x[2] - 4;
}

double s7(double *x)
{
	return x[0] * x[1] + x[2] - exp(x[2]);
}

/// <summary>
/// Set number 4
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s8(double *x)
{
	return x[1] - exp(x[0]);
}

double s9(double *x)
{
	return x[1] - 4 * x[0] * x[0] + 1;
}

/// <summary>
/// Set number 5
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s10(double *x)
{
	return x[0] * x[0] - x[0] * x[1] + x[1] * x[1] - 21;
}

double s11(double *x)
{
	return x[0] * x[0] + 2 * x[0] * x[1] - 8 * x[1] * x[1];
}

/// <summary>
/// Set number 6
/// </summary>
/// <param name="x"></param>
/// <returns></returns>
double s12(double *x)
{
	return x[0] * x[0] * x[1] * x[1] * x[2] * x[3] + 2 * x[0] * x[0] * x[1] * x[2] * x[2] * x[3] - 1;
}

double s13(double* x)
{
	return x[0] * x[0] * x[1] * x[1] * x[2] * x[3] + x[0] * x[1] * x[1] * x[2] * x[2] - 1;
}

double s14(double* x)
{
	return 2 * x[0] * x[0] * x[1] * x[2] * x[2] * x[3] + x[0] * x[1] * x[1] * x[2] * x[2] - 1;
}

double s15(double* x)
{
	return x[0] * x[1] + x[0] * x[2] - 4;
}

int main()
{
	try {
		//cout << bisection(funct, 0, 100) << endl;
		//cout << NewtonRaphson(funct, 50) << endl;
		//cout << Newton(funct, 50) << endl;
		//cout << Newton(ff, 0.75) << endl;
		//cout << sieczne(funct, 5, 10) << endl;
		//cout << sieczne(ff, -5, -2) << endl;
		//cout << sieczne(ff, 5, 4) << endl << endl;
		//cout << sieczne(fff, 5, 6) << endl;
		//cout << sieczne(fff, 6.5, 7) << endl;
		//cout << sieczne(fff, 7.5, 8) << endl;
		//cout << sieczne(fff, -2, 0) << endl;
		//cout << simpleIt(functImproved, 1) << endl;
		//cout << regulaFalsi(funct, 0, 100) << endl;


		//cout << bisection(fff, -2, 1, 1e-8) << endl << endl;
		//cout << NewtonRaphson(fff, 1, 1e-8) << endl << endl;
		//cout << sieczne(fff, 1, 2, 1e-8) << endl << endl;

		//cout << bisection(ff, -2, 0, 1e-8) << endl << endl;
		//cout << NewtonRaphson(ff, -1, 1e-8) << endl << endl;
		//cout << sieczne(ff, -2, -1, 1e-8) << endl << endl;


		FSystem fs[] = { s12, s13, s14, s15 };
		double x[] = { 1,1,1,1 };
		//for (auto f : fs)
		//{
		//	cout << f(x) << endl;
		//}
		int n = size(fs);

		double *res = multiDimensionalEq(fs, x, n);

		for (int i = 0; i < n; i++)
		{
			cout << res[i] << endl;
		}
		////cout << partialDerivative(s3, x, 1, size(x)) << endl;
	}
	catch (const char* msg) {
		cout << msg << endl;
	}
}