#pragma once
#include <random>
#include <iostream>
#include <stdexcept>
#include <stdexcept>

template <typename T>
class Matrix
{
	int size;
	T* data;
	int rows, columns;
	int swapped = 1;

	void swap(Matrix& rhs)
	{
		using std::swap;
		swap(size, rhs.size);
		swap(data, rhs.data);
		swap(rows, rhs.rows);
		swap(columns, rhs.columns);
	}
public:
	Matrix(): size(0), data(nullptr), rows(0), columns(0)  {}

	Matrix(const int n) : size(n* n), data(new T[size]), rows(n), columns(n) {}

	Matrix(const int rows, const int columns): size(rows * columns), data(new T[size]), rows(rows), columns(columns) {}

	Matrix(const Matrix& matrix): size(matrix.size), data(new T[size]), rows(matrix.rows), columns(matrix.columns)
	{
		memcpy_s(data, size * sizeof(T), matrix.data, size * sizeof(T));
	}

	~Matrix()
	{
		delete[] data;
	}

	T* operator[](const int r)
	{
		return data + (r * columns);
	}

	T* operator[](const int r) const
	{
		return data + (r * columns);
	}
	
	Matrix& operator=(const Matrix& other)
	{
		Matrix tmp(other);
		swap(tmp);
		return *this;
	}

	Matrix operator+(const Matrix& other)
	{
		if (rows != other.rows || columns != other.columns)
			throw std::invalid_argument("Nierowne macierze");
		Matrix tmp(rows, columns);
		for (int i = 0; i < size; i++)
			tmp.data[i] = data[i] + other.data[i];
		return tmp;
	}

	Matrix operator*(const double scalar)
	{
		Matrix tmp(rows, columns);
		for (int i = 0; i < size; i++)
			tmp.data[i] = data[i] * scalar;
		return tmp;
	}

	Matrix operator*(const Matrix& other)
	{
		if (columns != other.rows)
			throw std::invalid_argument("Nieprawidlowe macierze wejsciowe");

		Matrix tmp(rows, other.columns);
		tmp.fillValue();

		Matrix tmp2 = other.transpose();

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < other.columns; j++)
			{
				for (int k = 0; k < columns; k++)
					tmp[i][j] += (*this)[i][k] * tmp2[j][k];
			}
		return tmp;
	}

	Matrix transpose() const
	{
		Matrix tmp(columns, rows);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < columns; j++)
				//tmp[j][i] = data[i * rows + j];
				tmp[j][i] =  (*this)[i][j];
		return tmp;
	}

	void fillValue(int value = 0)
	{
		std::fill_n(data, size, value);
	}

	void fillIdentity()
	{
		if (rows != columns) return;
		fillValue();
		for (int i = 0; i < size; i += columns + 1)
			*(data + i) = 1;
	}

	void fillRandom(int from, int to, std::mt19937_64 generator)
	{
		std::uniform_int_distribution<> distr(from, to);
		for (int i = 0; i < size; i++)
			data[i] = distr(generator);
	}

	Matrix subMatrix(int startRow, int endRow, int startCol, int endCol)
	{
		int rows = endRow - startRow + 1;
		int cols = endCol - startCol + 1;

		Matrix temp(rows, cols);
		for (int i = startRow, j = 0; i <= endRow; i++, j++)
			memcpy_s(temp[j], cols * sizeof(T), (*this)[i] + startCol, cols * sizeof(T));
		return temp;
	}

	void decompose(Matrix<double>& upper, Matrix <double>& lower)
	{
		if (rows != columns) return;

		upper = Matrix<double>(rows);
		lower = Matrix<double>(rows);

		upper.fillValue();
		lower.fillValue();

		for (int i = 0; i < rows; i++)
		{
			lower[i][i] = 1;
			for (int j = 0; j <= i; j++)
			{
				double sum = 0;
				for (int k = 0; k < j; k++)
					sum += lower[j][k] * upper[k][i];

				upper[j][i] = (*this)[j][i] - sum;
			}

			for (int j = i + 1; j < rows; j++)
			{
				double sum = 0;
				for (int k = 0; k < i; k++)
					sum += lower[j][k] * upper[k][i];

				lower[j][i] = ((*this)[j][i] - sum) / upper[i][i];
			}
		}
	}

	void decompose()
	{
		for (int i = 0; i < rows - 1; i++)
		{
			double divider = (*this)[i][i];

			for (int j = i + 1; j < rows; j++)
				(*this)[j][i] /= divider;

			for (int j = i + 1; j < rows; j++)
				for (int k = i + 1; k < rows; k++)
					(*this)[j][k] -= (*this)[j][i] * (*this)[i][k];
		}
	}

	double pivotDecompose()
	{
		int *W = new int[rows];

		for (int i = 0; i < rows; i++) W[i] = i;

		for (int k = 0; k < rows; k++)
		{
			int maxw = k;
			double maxe = fabs((*this)[W[k]][k]);

			for (int i = k + 1; i < rows; i++)
				if ((*this)[W[i]][k] > maxe)
				{
					maxw = i;
					maxe = (*this)[W[i]][k];
				}

			if (maxw != k)
			{
				swapped = -swapped;
				std::swap(W[k], W[maxw]);
			}

			double divider = (*this)[W[k]][k];

			for (int i = k + 1; i < rows; i++) 
				(*this)[W[i]][k] /= divider;

			for (int i = k + 1; i < rows; i++)
				for (int j = k + 1; j < rows; j++) 
					(*this)[W[i]][j] -= (*this)[W[i]][k] * (*this)[W[k]][j];
		}

		double sum = swapped * (*this)[W[0]][0];
		for (int k = 1; k < rows; k++) sum *= (*this)[W[k]][k];
		delete[] W;
		return sum;
	}

	Matrix<double> solve(Matrix &b)
	{
		if (!(rows == columns == b.rows && b.columns == 1) == false) return NULL;
		Matrix<double> upper, lower;
		decompose(upper, lower);

		Matrix<double> z(b.rows, 1);
		Matrix<double> x(b.rows, 1);

		z[0][0] = b[0][0];
		for (int i = 1; i < rows; i++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
				sum += z[j][0] * lower[i][j];
			z[i][0] = b[i][0] - sum;
		}

		x[rows - 1][0] = z[rows - 1][0] / upper[rows - 1][rows - 1];
		for (int i = rows - 2; i >= 0; i--)
		{
			double sum = 0;
			for (int j = rows - 1; j > i; j--)
				sum += x[j][0] * upper[i][j];
			x[i][0] = (z[i][0] - sum) / upper[i][i];
		}
		return x;
	}

	double det()
	{
		if (rows != columns)
			throw std::invalid_argument("Macierz nie jest kwadratowa");
		Matrix<double> u, l;
		decompose(u, l);

		double mul = 1;
		for (int i = 0; i < rows; i++)
		{
			mul *= u[i][i];
		}

		return mul;
	}

	Matrix<double> Cramer(Matrix<double> b)
	{
		Matrix cp(*this);

		double a = cp.pivotDecompose();
		Matrix<double> res(rows, 1);
		for (int j = 0; j < rows; ++j) {
			cp = *this;
			cp.swapped = 1;
			for (int i = 0; i < rows; ++i)
			{
				cp[i][j] = b[i][0];
			}
			res[j][0] = cp.pivotDecompose() / a;
		}

		return res;
	}

	Matrix<double> Cramer2(Matrix<double> b)
	{
		Matrix X(*this);
		double detA = X.pivotDecompose();
		Matrix<double> res(rows, 1);

		for (int k = 0; k < rows; k++)
		{
			X.swapped = 1;
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < rows; j++)
					if (j == k) 
						X[i][j] = b[i][0];
					else 
						X[i][j] = (*this)[i][j];
			res[k][0] = X.pivotDecompose() / detA;
		}
		return res;
	}

	void printMatrix()
	{
		for (int i = 0; i < rows * columns; i++)
		{
			std::cout << *(data + i) << '\t';
			if ((i + 1) % columns == 0)
				std::cout << std::endl;
		}
		std::cout << std::endl;
	}
};

