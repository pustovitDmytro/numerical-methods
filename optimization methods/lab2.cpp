// lab2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>

#define eps 0.00001

using namespace std;
void outMat(ostream &stm, double ** A, int n){
	int i, j;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++)
			stm << A[i][j] << "\t";
		stm << endl;
	}
};
double delivery(double* a, int k, double f(double*)){
	double f0, f1 = NAN, fp, fm, dx = 0.0000001;
	do{
		f0 = f1;
		a[k] += dx;
		fp = f(a);
		a[k] -= 2 * dx;
		fm = f(a);
		f1 = (fp - fm) / (2 * dx);
		dx = dx / 2;
	} while (fabs(f0 - f1) > eps / 100);
	return f1;
}
double second_delivery(double* a, int k ,int m, double f(double*)){
	//(f(x + hx, y + hy) - f(x + hx, y - hy) - f(x - hx, y + hy) + f(x - hx, y - hy)) / (4 * hx*hy);
	double f1,f2,f3,f4,d = 0.00001;
	a[k] += d;	a[m] += d;
	f1 = f(a);
	a[m] -= 2*d;
	f2 = f(a);
	a[k] -= 2*d; a[m] += 2*d;
	f3 = f(a);
	a[m] -= 2*d;
	f4 = f(a);
	return (f1-f2-f3+f4) / (4*pow(d, 2));
}
void grad(double*a, double f(double*), double*g, int n){
	for (int i = 0; i < n; i++)	g[i] = delivery(a, i, f);
}
double myfunc(double* x){
	return //pow(x[0], 2) + 8 * pow(x[1], 2) + 0.001*x[0] * x[1] - x[0] - x[1];
		pow(x[0] - 4.4, 20) + 30 * pow(x[1] + 3, 2);// +2 * (x[1] + 3)*(x[0] - 4.4);
}
double funk25(double* x){
	return 100 * pow(x[1] - pow(x[0], 3), 2) + 100 * pow(1 - x[0], 2);
}
void Gauss(double** A, double** B , int n){
	int k, i, j, t = 0;
	double L;
	for (k = 0; k < n; k++){
		for (i = 0; i < n; i++){
			if (i == k) ++i;
			if (i == n) break;
			if (A[k][k] == 0) {
				cout << "Error, matrix is wrong" << endl;
				system("Pause");
				exit(-1);
			}
			L = A[i][k] / A[k][k];
			for (j = 0; j < n; j++)
				A[i][j] = A[i][j] - L*A[k][j];
			for (j = 0; j < n; j++)
				B[i][j] = B[i][j] - L*B[k][j];
		}
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			B[i][j] = B[i][j] / A[i][i];
};
double normvec(double* x, double* y, int n){
	double s = 0;
	//	for (int i = 0; i < n; i++) s += pow(x[i] - y[i], 2);
	//	s=pow(s, 0.5);
	for (int i = 0; i < n; i++)
		if (fabs(x[i] - y[i])>s) s = fabs(x[i] - y[i]);
	cout << "||Xk - Xk+1|| = " << s << endl;
	return s;
}
int Newton(double *x1, double f(double*), int n){
	int k = 1;
	double **A = new double*[n], *x = new double[n], *g = new double[n], *h = new double[n], **B = new double*[n] , alfa, e = 0.1;
	for (int row = 0; row < n; row++){
		A[row] = new double[n];
		B[row] = new double[n];
	}
	for (int i = 0; i < n; i++)	x[i] = x1[i];
	do{
		for (int i = 0; i < n; i++)	 x1[i] = x[i];
		for (int i = 0; i < n; i++){
			x[i] = x1[i];
			for (int j = 0; j <= i; j++)
				A[i][j] = A[j][i] = second_delivery(x, i, j, f);
			for (int j = 0; j < n; j++) B[i][j] = 0;
			B[i][i] = 1;			
		}
		Gauss(A, B, n);
		grad(x, f, g, n);
		cout << k << " iteration: " << endl; //<<"(f'')^-1:\n";
	//	outMat(cout, B, n);
	//	cout << "grad = (" << g[0] << " , " << g[1] << ") \n solution: ";
		for (int i = 0; i < n; i++){
			h[i] = 0;
			for (int j = 0; j < n; j++)
				h[i] += -B[i][j] * g[j];
		}
		double s = 0, alfa = 2;
		for (int i = 0; i < n; i++) s += g[i] * h[i];
		do{
			alfa = alfa / 2;
			for (int i = 0; i < n; i++) x[i] = x1[i] + alfa*h[i];
		} while (fabs(f(x1) - f(x)) < e*alfa*s);		
		
		for (int i = 0; i < n; i++) cout << x[i] << " ";
		k++;
	} while (normvec(x, x1, n)>eps);
	cout << "f* = f(x*) = " << f(x1) << endl;
	delete[]x;
	delete[]g;
	delete[]h;
	for (int i = 0; i < n; i++) delete A[i];
	for (int i = 0; i < n; i++) delete B[i];
	delete A;
	delete B;
	return k - 1;
}

int _tmain(int argc, _TCHAR* argv[])
{	
	int n = 2;
	double *x = new double[n];
	x[0] = x[1] = 5;
	Newton(x, myfunc, n);
//	Newton(x, funk25, n);
	system("Pause");
	return 0;
}

