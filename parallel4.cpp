#include <cstdio>
#include <algorithm>
#include <assert.h>
#include "Matrix.h"

#define MAX_ERROR 0.00001
#define STEP 0.001

int main(int argc, char *argv[])
{
	FILE *f = fopen("t.txt", "r");
	ASSERT(f != NULL);
	
	Matrix Ax, Ay, Sx, Sy, R, Kx, Ky, H;
	
	Ax = Matrix::ReadFromFile(f);
	Ay = Matrix::ReadFromFile(f);
	
	Sx = Matrix::ReadFromFile(f);
	Sy = Matrix::ReadFromFile(f);
	
	Kx = Matrix::ReadFromFile(f);
	Ky = Matrix::ReadFromFile(f);
	
	R = Matrix::ReadFromFile(f);
	H = Matrix::ReadFromFile(f);
		
	fclose(f);
	
	Matrix W  = Ax.Inverse() * Ay;	
	Matrix TP = (Sy * Ky - Sx * Kx * W).Inverse() * (Sx & Sy);	
	Matrix Y(TP.n, 1);
	Matrix X(W.n,  1);
	
	Y.foreach( [](double &x) { x = 0;} );
	X.foreach( [](double &x) { x = 0;} );

	Matrix LastZ = (X.Transpose() & Y.Transpose()).Transpose();
	Matrix Z(LastZ);

	do 
	{
		LastZ = Z;

		Y = Y + TP * STEP * (H - R * Z);
		X = W * Y * -1;

		Z = (X.Transpose() & Y.Transpose()).Transpose();
		Z.foreach([](double &x){ x = x * fabs(x); });

	} while ( Z.SqrDifference( LastZ ) > sqr( MAX_ERROR ) );

	f = fopen("out2.txt", "w");
	ASSERT( f );
	
	fprintf(f, "W(%d, %d) = \n", W.n, W.m);
	W.PrintToFile( f );
	fprintf(f, "TP(%d, %d) = \n", TP.n, TP.m);
	TP.PrintToFile( f );

	for (int i = 0; i < X.n; i++)
		fprintf(f, "X[ %d ] = %lf\n", i, X[i][0]);
	
	for (int i = 0; i < Y.n; i++)
		fprintf(f, "Y[ %d ] = %lf\n", i, Y[i][0]);

	fclose(f);	
	return 0;
}
