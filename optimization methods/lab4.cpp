// lab4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <iostream>
#include <fstream>
#define alfamin 0
#define alfamax 1
#define eps 0.00001

using namespace std;
double normvec(double* x, double* y, int n){
	double s = 0;
	//	for (int i = 0; i < n; i++) s += pow(x[i] - y[i], 2);
	//	s=pow(s, 0.5);
	for (int i = 0; i < n; i++)
		if (fabs(x[i] - y[i])>s) s = fabs(x[i] - y[i]);
	cout << " ||Xk - Xk-1|| = " << s << endl;
	return s;
}
double f(double*x, double**A, double*b, int n){
	double *Ax = new double[n] , s=0;
	for (int i = 0; i < n; i++) {
		s += x[i] * b[i];
		Ax[i] = 0;
		for (int j = 0; j < n; j++) Ax[i] += A[i][j] * x[j];
	}
	for (int i = 0; i < n; i++) s += Ax[i] * x[i]/2;
	delete[] Ax;
	return s;
}
double golden_section(double t1, double t2, double *x, double* g, int n, double **A, double*b){
	double * tmp = new double[n], x1, x2, f1, f2, fi = (1 + pow(5, 0.5)) / 2;
	int i, k = 0;
	x1 = t2 - (t2 - t1) / fi;
	x2 = t1 + (t2 - t1) / fi;
	for (int i = 0; i < n; i++) tmp[i] = x[i] + x1*g[i]; f1 = f(tmp,A,b,n);
	for (int i = 0; i < n; i++) tmp[i] = x[i] + x2*g[i]; f2 = f(tmp, A, b, n);
	do{
		if (f1 < f2){
			t2 = x2;
			x2 = x1;
			f2 = f1;
			x1 = t2 - (t2 - t1) / fi;
			for (int i = 0; i < n; i++) tmp[i] = x[i] + x1*g[i]; f1 = f(tmp, A, b, n);
		}
		else{
			t1 = x1;
			x1 = x2;
			f1 = f2;
			x2 = t1 + (t2 - t1) / fi;
			for (int i = 0; i < n; i++) tmp[i] = x[i] + x2*g[i]; f2 = f(tmp, A, b, n);
		}
	} while (fabs(t1 - t2) > eps/1000000000);
	delete[] tmp;
	return (t1 + t2) / 2.0;
}
void kvadr_grad(double*x, double*b, double**A, int n, double* g){
	for (int i = 0; i < n; i++){
		g[i] = b[i];
		for (int j = 0; j < n; j++)
			g[i] += A[i][j] * x[j];
	}
}
double find_beta(double* g, double* h, double**A, int n){
	double* Ah = new double[n], s1 = 0,s2=0;
	for (int i = 0; i < n; i++){
		Ah[i] = 0;
		for (int j = 0; j < n; j++)
			Ah[i] += A[i][j] * h[j];
	}
	for (int i = 0; i < n; i++){
		s1 += g[i] * Ah[i];
		s2 += h[i] * Ah[i];
	}
	delete[] Ah;
	return s1 / s2;
}
int Sprazen_grad(double*x,double*b,double**A, int n){
	double *hk = new double[n], alfa,beta, *y = new double[n], *g = new double[n];
	kvadr_grad(x, b, A, n, hk);
	for (int i = 0; i < n; i++) hk[i] = -hk[i];
	for (int k = 0; k <=n+1; k++){
		alfa = golden_section(alfamin, alfamax, x, hk, n, A, b);
		for (int i = 0; i < n; i++) y[i] = x[i] + alfa*hk[i];
		cout << "iteration " << k << ": x[k] = ( ";	for (int i = 0; i < n; i++) cout << y[i] << "  "; cout << ") ;"<<"f="<<f(y,A,b,n);
		if (normvec(x, y, n) < eps) break;
		for (int i = 0; i < n; i++) x[i] = y[i];
		kvadr_grad(x, b, A, n, g);
		beta = find_beta(g, hk, A, n);
			for (int i = 0; i < n; i++) hk[i] = -g[i] + beta*hk[i];
	}
	delete[]hk;
	delete[]y;
	delete[]g;
	return 0;
}


int _tmain(int argc, _TCHAR* argv[])
{
	int n = 3;
	double *x = new double[n], ** A = new double*[n], *b = new double[n];
	for (int row = 0; row < n; row++) A[row] = new double[n];
	ifstream stm("input.txt");
		int i, j;
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++)
				stm >> A[i][j];
			stm >> b[i] >> x[i];
		}
		stm.close();
		Sprazen_grad(x, b, A, n);
		system("Pause");
	return 0;
}

