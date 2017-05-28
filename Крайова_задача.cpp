// boundary_value_problem.cpp : Defines the entry point for the console application.
//



#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

#define ca 1.135
#define cb -2.643
#define cc -1.41
#define cd -2.767
#define ce 1.777
#define cf1	2.09303
#define	cf2 -3.3906

#define calp0	0
#define calp1	1
#define cbet0	0.4
#define cbet1	-1
#define cx1 2
#define cx2 2.3


using namespace std;

inline double answ(double x){
	return ca*pow(x, 2) + cb*x + cc + 1.0 / (cd*x + ce);
}
inline double dansw(double x){
	return 2 * ca * x + cb - cd / pow(ce + cd*x, 2);
}
inline double ddansw(double x){
	return 2 * ca + 2 * cd * cd / pow(ce + cd*x, 3);
}
void findparams(){
	cout << "cf1=" << dansw(cx1);
	cout << "\ncf2=" << cbet0*answ(cx2) + cbet1*dansw(cx2)<<endl;
}

double err(vector<double> y, double a, double b, int n){
	ofstream fout("error.txt");
	fout << "помилка методу кінцевих різниць:" << endl;
	double max=0,x = a, e, h = (b - a) / (n - 1);
	for (int i = 0; i < n; i++){
		e = fabs(answ(x) - y[i]);
		if (e > max) max = e;
		fout << " err[ " << x << " ] = " << e << endl;
		x += h;
	}
	return max;
}

void simple_progonka(vector<double> &y, vector<double> m, vector<double> n, vector<double> f, double h, int N){
	// y[i+2] + m[i]*y[i+1]+n[i]*y[i]=f[i]*h^2 i=0,N-2	
	//прямий хід
	vector<double> c(N), d(N);
	c[0] = (calp1 - calp0*h) / (m[0] * (calp1 - calp0*h) + n[0] * calp1);
	d[0] = (n[0] * cf1*h) / (calp1 - calp0*h) + f[0] * h*h;
	for (int i = 1; i < N - 1; i++){
		c[i] = 1.0 / (m[i] - n[i] * c[i - 1]);
		d[i] = f[i] * h*h - n[i] * c[i - 1] * d[i - 1];
	}
	//зворотній хід
	y[N] = (cbet1*c[N-2]*d[N-2]+cf1*h) / (cbet1*(1+c[N-2]+cbet0*h));
	for (int i = N; i > 1; i--) 	y[i - 1] = c[i - 2] * (d[i - 2] - y[i]);
	y[0] = (calp0*y[1] - cf1*h) / (calp1 - calp0*h);
}

void progonka(vector<double> &y, vector<double> m, vector<double> n, vector<double> p,vector<double> f, double h, int N){
	vector<double> c(N+1), d(N+1);
	c[1] = (calp1 - calp0*h) / (m[1] * (calp1 - calp0*h) + n[1] * calp1);
	d[1] = (2 * f[1] * h*h) / (2 + p[1] * h) + n[1] * cf1*h / (calp1 - calp0*h);
	for (int i = 2; i <= N; i++){
		c[i] = 1.0 / (m[i] - n[i] * c[i - 1]);
		d[i] = 2 * f[i] * h*h / (2 + h*p[i]) - n[i] * c[i - 1] * d[i - 1];
	}
	y[N] = (2 * cf2*h - cbet1*(d[N] - c[N - 1] * d[N - 1])) / (2 * cbet0*h + cbet1*(c[N - 1] - 1.0 / c[N]));
	for (int i = N - 1; i > 0; i--) y[i] = c[i] * (d[i] - y[i + 1]);
	y[0] = (cf1*h - calp1*y[1]) / (calp0*h - calp1);
}


void kin_rizn(double a, double b, int n, vector<double> &y){
	double h = (b - a) / (n-1);
	double x=a, q=-2.0;
	vector<double> p(n), f(n), t1(n),t2(n);
	for (int i = 0; i < n; i++){
		p[i] = 2 * x;
		t1[i] = (2 * q*h*h - 4) / (2+h*p[i]);	//m
		t2[i] = (2-h*p[i])/(2+h*p[i]);			 //n
		f[i] = ddansw(x) + p[i] * dansw(x) + q * answ(x);
		x += h;
	}
	progonka(y,t1, t2, p,f, h, n-1);
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n = 10;
	vector<double> y(n);
		//findparams();
	kin_rizn(cx1, cx2, n, y);
	cout << "Max error value = " << err(y, cx1, cx2, n) << endl;
	system("Pause");
	return 0;
}

