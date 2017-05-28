#include "stdafx.h"
#include <math.h>
#include <iostream>
#define eps 0.00001
#define alfamin 0
#define alfamax 10
using namespace std;

double f(double* x){
	return 1*pow(x[0], 2) + 1 * pow(x[1]-1, 2) + 1 * pow(x[2], 2);
}
double skal(double* x, double* p, int n){
	double s = 0;
	for (int i = 0; i < n; i++) s += x[i] * p[i];
	return s;
}
void PiX(double* x, double* p, double d, int n, double* a){
	for (int i = 0; i < n; i++)	a[i] = x[i] + (d - skal(x, p, n))*p[i] / skal(p, p, n);
}
double normvec(double* x, double* y, int n){
	double s = 0;
	for (int i = 0; i < n; i++) s += pow(x[i] - y[i], 2);
	s = pow(s, 0.5);
	cout << "||Xk - Xk+1|| = " << s << endl;
	return s;
}
void Oper(double **A, double*p, double d, int n){
	double *x = new double[n], *a = new double[n];
	for (int k = 0; k < n; k++){
		x[0] = 0;
		for (int i = 0; i < n; i++) x[i] = 0;
		x[k] = 1;
		PiX(x, p, d, n, a);
		for (int i = 0; i < n; i++) A[i][k] = a[i];
	}
	delete[]x;
	delete[]a;
}
double delivery(double* a, int k, double f(double*)){
	double f0, f1 = NAN, fp, fm, dx = 0.00001;
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
void grad(double*a, double f(double*), double*g, int n){
	for (int i = 0; i < n; i++)	g[i] = delivery(a, i, f);
}
double golden_section(double t1, double t2, double *x, double* g, int n, double f(double*)){
	double * tmp = new double[n], x1, x2, f1, f2, fi = (1 + pow(5, 0.5)) / 2;
	int i, k = 0;
	x1 = t2 - (t2 - t1) / fi;
	x2 = t1 + (t2 - t1) / fi;
	for (int i = 0; i < n; i++) tmp[i] = x[i] - x1*g[i]; f1 = f(tmp);
	for (int i = 0; i < n; i++) tmp[i] = x[i] - x2*g[i]; f2 = f(tmp);
	do{
		if (f1 < f2){
			t2 = x2;
			x2 = x1;
			f2 = f1;
			x1 = t2 - (t2 - t1) / fi;
			for (int i = 0; i < n; i++) tmp[i] = x[i] - x1*g[i]; f1 = f(tmp);
		}
		else{
			t1 = x1;
			x1 = x2;
			f1 = f2;
			x2 = t1 + (t2 - t1) / fi;
			for (int i = 0; i < n; i++) tmp[i] = x[i] - x2*g[i]; f2 = f(tmp);
		}
	} while (fabs(t1 - t2) > eps);
	delete[] tmp;
	return (t1 + t2) / 2.0;
}
double met_drob(double *x, double* g, int n, double f(double*)){
	double al = 10, lambda = 0.4, *t = new double[n];
	for (int i = 0; i<n; i++) t[i] = x[i] - al*g[i];
	while (f(t) > f(x)){
		al = lambda*al;
		for (int i = 0; i < n; i++) t[i] = x[i] - al*g[i];
	}
	return al;
}
int PrGrad(double *x1, double f(double*), int n, double* p){
	int k = 1;
	double *x = new double[n], alfa, *g = new double[n], *a = new double[n], s;
	do{
		for (int i = 0; i < n; i++) x[i] = x1[i];
		grad(x, f, g, n);
		alfa = golden_section(alfamin, alfamax, x, g, n, f);
		//	alfa = met_drob(x, g, n, f);
		cout << "grad = " << g[0] << " , " << g[1] << " ;  alfa = " << alfa << "  f(xk) = " << f(x1) << endl;
		cout << k << " iteration: ( ";
		for (int i = 0; i < n; i++) a[i] = x[i] - alfa*g[i];
		PiX(a, p, 1, n, x1);
		PiX(x, p, 1, n, a);
		for (int i = 0; i < n - 1; i++) cout << x1[i] << " , ";
		cout << x1[n - 1] << " );" << endl;
		k++;
		if (k == 500){
			cout << " error, too much iterations";
			return -1;
		}
	} while (normvec(a, x1, n)>eps);
	cout << "f* = f(x*) = " << f(x1) << endl;
	return k - 1;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 3;
	double *x = new double[n], *p = new double[n];
	x[0] = x[1] = x[2] = 5;
	p[0] = 1; p[1] = 1; p[2] = 1;
	PrGrad(x, f, n, p);
	system("Pause");
	return 0;
}

