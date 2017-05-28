// lab1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <iostream>
#define eps 0.00001
#define alfamin 0
#define alfamax 1.5
using namespace std;

double myfunc(double* x){
	return //pow(x[0], 2) + 8*pow(x[1], 2)+0.001*x[0]*x[1]-x[0]-x[1];
	exp(pow(x[0] - 4.4, 2) + 3 * pow(x[1] + 3, 2) + 2 * (x[1] + 3)*(x[0] - 4.4));// 0.001*x[0] * x[1] - x[0] - x[1];
}
double funk25(double* x){
	return 100 * pow(x[1] - pow(x[0], 3), 2) + 100 * pow(1 - x[0], 2);
}
double delivery(double* a, int k, double f(double*)){
	double f0, f1 = NAN, fp, fm, dx = 0.00001;
	do{
		f0 = f1;
		a[k] += dx;
		fp = f(a);
		a[k] -= 2 * dx;
		fm = f(a);
		f1 = (fp - fm)/(2*dx);
		dx = dx / 2;
	} while (fabs(f0 - f1) > eps/100);
	return f1;
}

void grad( double*a, double f(double*), double*g, int n){
	for (int i = 0; i < n; i++)	g[i] = delivery(a, i, f);
}

double golden_section(double t1, double t2, double *x, double* g, int n, double f(double*)){
	double * tmp = new double[n], x1, x2,f1,f2, fi = (1 + pow(5, 0.5)) / 2;
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

double normvec(double* x, double* y, int n){
	double s = 0;
//	for (int i = 0; i < n; i++) s += pow(x[i] - y[i], 2);
//	s=pow(s, 0.5);
	for (int i = 0; i < n; i++)
		if (fabs(x[i] - y[i])>s) s = fabs(x[i] - y[i]); 
	cout << "||Xk - Xk+1|| = " << s << endl;
	return s;
}

int grad_spusk(double *a, double f(double*), int n){
	double *a0 = new double[n], alfa, *g = new double[n];
	int k = 1;
	do{
		for (int i = 0; i < n; i++) a0[i] = a[i];
		grad(a0, f, g, n);
		alfa = golden_section(alfamin, alfamax, a0, g, n, f);
		cout << "grad = " << g[0] << " , " << g[1] << ";  alfa = " << alfa << "  f(xk) = " << f(a) << endl;
		cout << k << " iteration: ( ";
		for (int i = 0; i < n; i++) a[i] = a0[i] - alfa*g[i];
		for (int i = 0; i < n-1; i++) cout << a[i] << " , ";
		cout << a[n-1] << " );" << endl;
		k++; 
		if (k == 50000000){
			cout << " error, too much iterations";
			return -1;
		}
	} while (normvec(a, a0, n)>eps);
	cout << "f* = f(x*) = " << f(a) << endl;
	return k - 1;
	delete[] a0;
	delete[] g;
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 2;
	double* x = new double[n];
	x[0] = 4.5;
	x[1] = -3.5;
	grad_spusk(x, myfunc, n);
	//grad_spusk(x, funk25, n);
	system("Pause");
	return 0;
}