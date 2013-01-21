#include <cstdio>
#include <algorithm>
#include <assert.h>
#include "Matrix.h"

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
	
	f = fopen("out.txt", "w");
	ASSERT( f );
	
	fprintf(f, "W(%d, %d) = \n", W.n, W.m);
	W.PrintToFile( f );
	fprintf(f, "TP(%d, %d) = \n", TP.n, TP.m);
	TP.PrintToFile( f );
	fclose(f);	
	return 0;
}
