#ifndef _MATRIX_H
#define _MATRIX_H

#include <cstdio>
#include <algorithm>
#include <assert.h>
#include "parallel.h"

#define FOREACH for (int i = 0; i < n; i++) for(int j = 0; j < m; j++)
#define ASSERT  assert
#define sqr(x) ((x)*(x))

class Matrix
{
	public:
	int n, m;
	double *Data;

	Matrix()
	{
		n = m = 0;
		Data = NULL;
	}
	
	Matrix(int n, int m)
	{
		this->n = n;
		this->m = m;	
		Data = (double *)malloc(sizeof(*Data) * n * m);
	}

	Matrix(const Matrix &x) 
	{
		this->n = x.n;
		this->m = x.m;
		Data = (double *)malloc(sizeof(*Data) * n * m);	
		memcpy(Data, x.Data, sizeof(*Data) * n * m);
	}
	
	static Matrix&  ReadFromFile(FILE *f)
	{
		int n, m;
		fscanf(f, "%d%d", &n, &m);
		Matrix &x = *new Matrix(n, m);
		FOREACH fscanf(f, "%lf", &x[i][j]);			
		return x;
	}
	
	Matrix operator * (const Matrix &x)
	{
		ASSERT(m == x.n);
		return ParallelMatrixMultiplication::Instance()->Mul(*this, x);
	}
	
	Matrix operator + (const Matrix &x)
	{
		ASSERT(n == x.n && m == x.m);
		Matrix s(n, m);	
		FOREACH	
			s[i][j] = (*this)[i][j] + x[i][j];
		return s;
	}
	
	Matrix operator - (const Matrix &x)
	{
		ASSERT(n == x.n && m == x.m);
		Matrix s(n, m);	
		FOREACH	
			s[i][j] = (*this)[i][j] - x[i][j];

		return s;
	}

	Matrix operator * (double x)
	{
		Matrix s(n, m);	
		FOREACH	
			s[i][j] = (*this)[i][j] * x;
		return s;
	}


	
	Matrix Inverse()
	{	
		ASSERT(n == m);
		Matrix E(n,n);
		Matrix T(n,n);
		T = *this;
		
		FOREACH	
			E[i][j] = (i == j);
		
		
		for (int i = 0; i < n; i++)
		{
			int t = i;
			while (t < n && T[t][i] == 0) t++;
			ASSERT(t < n);
			
			T.SwapRows(i, t);
			E.SwapRows(i, t);

			double k = T[i][i];
			for (int j = 0; j < m; j++) 
					T[i][j] /= k, E[i][j] /= k;
			
			for (int t = i + 1; t < n; t++)
			{
				double k = T[t][i];
				for(int j = 0; j < m; j++)
				{
					T[t][j] -= T[i][j] * k;
					E[t][j] -= E[i][j] * k;
				}
			}
		}	
		
		for (int i = n - 1; i >= 0; i--)
			for (int j = 0; j < m; j++)
				for (int t = i + 1; t < n; t++)
					E[i][j] -= E[t][j] * T[i][t];
		
		return E;			
	}
	
	void PrintToFile(FILE *f)
	{
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
				fprintf(f, "%lf ", (*this)[i][j]);
			fprintf(f, "\n");
		}
	}

	// should be called AppendRight
	Matrix operator & (const Matrix &x)
	{
		Matrix s(n,  m + x.m);
		
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < m; j++)
				s[i][j] = (*this)[i][j];
			for(int j = 0; j < x.m; j++)
				s[i][j + m] = x[i][j];
		}
		return s;
	}
	
	Matrix& operator = (const Matrix &x)
	{
		if (this == &x) return *this;
		free(Data);
		this->n = x.n;
		this->m = x.m;
		Data = (double *)malloc(sizeof(*Data) * n * m);	
		memcpy(Data, x.Data, sizeof(*Data) * n * m);		
		return *this;
	}
	
	double* operator [] (int x) const 
	{
		return &Data[x * this->m];
	}
	
	// call given method for each element
	void foreach( void (*f)(double *value) )
	{
		FOREACH f( (*this)[i] + j );
	}

	Matrix Transpose()
	{
		Matrix res(m, n);
		FOREACH 
		{
			res[j][i] = (*this)[i][j];
		}
		return res;
	}

	// should be replaced with GetNorm()
	double SqrDifference( Matrix x )
	{
		ASSERT(n == x.n && m == x.m);
		double t = 0;
		FOREACH
			t += sqr((*this)[i][j] - x[i][j]);
		
		return t;
	}

	~Matrix()
	{
		free(Data);
	}
private: 

	void SwapRows(int i, int j)
	{
		double *p1 = (*this)[i];
		double *p2 = (*this)[j];
		for (int i = 0; i < this->m; i++)
			std::swap(*(p1++), *(p2++));
	}
};


#endif

