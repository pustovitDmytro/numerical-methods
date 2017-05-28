// ZadKoshi.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#define X0 1
#define Y0 0

inline double f(double x, double y){
	return 3 * pow(1 + pow(y / x, 2), 1 / 2.0) + y / x;
}

inline double ansv(double x){
	return (pow(x, 6) - 1) / (2 * x * x);
}

void err(double *y, double h, int n,double x0){
	ofstream stm("error.txt", ios_base::app);
	double t=abs(ansv(x0) - y[0]),s=t;
	stm << "err(" << x0 << ") = " << t << endl;
	for (int i = 1; i < n; i++){
		x0 += h;
		t = abs(ansv(x0) - y[i]);
		if (t > s) s = t;
		stm << "err(" <<x0 << ") = "<<t<<endl;
	}
	cout << "maxerr =" << s<<endl;
}
void print(char* s,double *y, double h, int n, double x0){
	ofstream stm(s);
	stm << "x = " << x0 << " y= " << ansv(x0)<<" `y= " << *y << endl;
	for (int i = 1; i < n; i++){
		x0 += h;
		stm << "x = " << x0 << " y= " << ansv(x0) << " `y= " << y[i] << endl;
	}
}

void runge_kut(double* y, double h,double x,int n){
	double k1, k2, k3, k4;
	for (int i = 0; i < n-1;i++){
		k1 = f(x, y[i]);
		k2 = f(x + h / 2.0, y[i] + h* k1 / 2.0 );
		k3 = f(x + h / 2.0, y[i] + h* k2 / 2.0 );
		k4 = f(x + h, y[i] + h*k3);	
		x += h;
		y[i + 1] = y[i] + h*(k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
	}
}
void adams(double* y, double h, double x, int n){
	for (int i = 0; i < n-3; i++){
		y[i + 3] = y[i+2]+h*(3*f(x+3*h,y[i+3])/8.0+
			19*f(x+2*h,y[i+2])/24.0-5*f(x+h,y[i+1])/24.0+
			f(x,y[i])/24.0);
		x += h;
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	int n;
	double h, x = X0;
	cin >> n >> h;
	//n = 10;
	double *y = new double[n];
	y[0] = Y0;
	
	//for (h = 0.05; h < 0.5; h += 0.05){
		//cout << endl<<"h="<<h << endl;
		runge_kut(y, h, x, n);
		print("runge_kut.txt", y, h, n, x);
		err(y, h, n, x);
		adams(y, h, x, n);
		print("adams.txt", y, h, n, x);
		err(y, h, n, x);
	//}
	system("Pause");
	return 0;
}

