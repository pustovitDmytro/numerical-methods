// lab_9.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#define L 1.0
#define T 5.0
#define p M_PI

using namespace std;

double u0(double x){
	return (1 - x*x)*cos(p*x);
}
double u1(double t){
	return cos(p*t);
}
double u2(double t){
	return 0;
}
double F(double x, double t){
	return 2 * cos(p*t)*cos(p*x)-4*p*x*cos(p*t)*sin(p*x);
}
double ansv(double x, double t){
	return u0(x)*cos(p*t);
}

int _tmain(int argc, _TCHAR* argv[])
{
	int N = 100, M=1000;
	double dt = T / (M-1), dx = L /(N-1);
	double sigma = 0.00001;
	vector<vector<double>> y;
	y.resize(N); // к-ть строк
	for (int i = 0; i < y.size(); i++) y[i].resize(M);
	for (int i = 0; i < N; i++) y[i][0] = u0(dx*i);
	for (int i = 0; i < M; i++){
		y[0][i] = u1(dt*i);
		y[N-1][i] = u2(dt*i);
	}
	for (int i = 1; i < N-1; i++)	{
		y[i][1] = y[i][0] + (dt*dt/2)*( ( (y[i + 1][0] - 2 * y[i][0] + y[i - 1][0]) / (dx*dx) ) + F(i*dx, 0) );
	}	
	vector<double> a(N), b(N), z(N), w(N), q(N), s(N);
	int k = 1;
	while (k < M-1)
	{
		for (int i = 1; i < N-1; i++)
		{
			a[i] = z[i] = (sigma*dt*dt) / (dx*dx);
			b[i] = -(1 + 2 * ((sigma*dt*dt) / (dx*dx)));
			s[i] = (2 * y[i][k] - y[i][k - 1] + (((1 - 2 * sigma)*dt*dt) / (dx*dx))*(y[i + 1][k] - 2 * y[i][k] + y[i - 1][k])
				+ ((sigma*dt*dt) / (dx*dx))*(y[i + 1][k - 1] - 2 * y[i][k - 1] + y[i - 1][k - 1])
				+ dt*dt*(sigma*F(i*dx, (k + 1)*dt) + (1 - 2 * sigma)*F(dx*i, (k)*dt) + sigma*F(dx*i, (k - 1)*dt)));
			//
		}
		a[1] = 0;
		z[N - 2] = 0;
		w[1] = z[1] / b[1];
		q[1] = -s[1] / b[1];
		for (int i = 2; i < N - 2; i++)
		{
			w[i] = z[i] / (b[i] - a[i] * w[i - 1]);
			q[i] = (a[i] * q[i - 1] - s[i]) / (b[i] - a[i] * w[i - 1]);
		}

		y[N - 2][k + 1] = (a[N - 2] * q[N - 3] - s[N - 2]) / (b[N - 2] - a[N - 2] * w[N - 3]);

		for (int i = N - 2; i >= 1; i--)
			y[i][k + 1] = w[i] * y[i + 1][k + 1] + q[i];
		y[0][k + 1] = u1((k + 1)*dt);
		y[N-1][k + 1] = u2((k + 1)*dt);
		k++;
	}

	ofstream fout("out.txt");
	double norm = 0, tmp,t,maxy=0;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++){
			if (abs(y[i][j])>maxy) maxy = abs(y[i][j]);
			tmp = abs(y[i][j] - ansv(dx*i, dt*j));
			if (tmp > norm) norm = tmp;
		}
		fout << "|| y || = " << maxy << " || u - y || = " << norm << endl;
		for (int i = 0; i < N; i++){
			fout << " yk[" << dx*i << "] = " << y[i][M - 1] <<
				"  \t(yk - u)[" << dx*i << "] = " << y[i][M - 1] - ansv(dx*i, T) << "      \tu =  " << ansv(dx*i, T) << endl;
		}
		cout << norm << endl;
	system("Pause");
	return 0;
}

