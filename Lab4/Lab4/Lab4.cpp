#include <iostream>

using namespace std;

using Funkcja = double(*)(double);

double square(Funkcja f, double start, double end, double h)
{
	double sum = 0;
	for (double x = start; x < end; x+= h)
		sum += f((x + x + h) / 2) * h;
	return sum;
}

double trapezium(Funkcja f, double start, double end, double h)
{
	double sum = 0;
	for (double x = start; x < end; x += h)
		sum += (f(x) + f(x + h)) * h / 2;
	return sum;
}

double Simpson(Funkcja f, double start, double end, double h)
{
	double sum = 0;
	for (double x = start; x < end; x += h)
		sum += (f(x) + 4 * f((x + x + h) / 2) + f(x + h)) * h / 6;
	return sum;
}

double f(double x)
{
	return sin(x);
}

int main()
{
	cout << square(f, 0, 3.14, 1e-3) << endl;
	cout << trapezium(f, 0, 3.14, 1e-3) << endl;
	cout << Simpson(f, 0, 3.14, 1e-3) << endl;
}

