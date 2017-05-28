// OwnNumbers.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <math.h>

#define n 4
#define eps 0.00001
#define zn 3

using namespace std;

typedef double Mat[n][n];
typedef double vec[n];

void fileMat(ofstream &);
void inMat(char*, Mat*);
void outMat(char*,char*,Mat*);
void outVec(char*, char*, vec*, bool);
void QR(Mat&, Mat&);
void Matprod(Mat, Mat, Mat &);
void Hausholder(Mat &, vec);
void findEigen(Mat&, vec&, Mat&);
void correct(Mat, vec, Mat);
double step(Mat);

int _tmain(int argc, _TCHAR* argv[])
{
	Mat A,B;
	vec b;
	inMat("inMat.txt", &A);
	ofstream fout("ans.txt");
	fout << "Розвязок часткової проблеми (Max eigen): " << step(A) << endl;
	findEigen(A, b, B);
	outVec("ans.txt", "\nвласні числа:", &b, false);
	outMat("ans.txt", "власні вектори:", &B);
	inMat("inMat.txt", &A);
	correct(A,b,B);
	system("Pause");
	return 0;
}

double step(Mat A){
	vec x, x1;
	for (int i = 0; i < n; i++) x[i] = 1;
	bool flag;
	double t;
	do
	{
		for (int i = 0; i < n; i++) x1[i] = x[i];
		for (int i = 0; i < n; i++){ //x=Ax
			x[i] = 0;
			for (int j = 0; j < n; j++)
				x[i] += A[i][j] * x1[j];
		}
		t = x[0] / x1[0];
		flag = true;
		for (int i = 1; i < n; i++)
			if ( abs(x[i]/x1[i]-t)>eps ){
				flag = false;
				break;
			}
	} while (!flag);
	double s1=0, s2=0;
	for (int i = 0; i < n; i++){
		s1 += x[i] * x[i];	// ||x||^2
		s2 += x1[i] * x1[i];	// ||x1||^2
	}
	return pow(s1 / s2, 0.5);
}
void Matprod(Mat A, Mat B, Mat &AB){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++){
			AB[i][j] = 0;
			for (int k = 0; k < n; k++) AB[i][j] += A[i][k] * B[k][j];
		}
}
void Hausholder(Mat &P, vec v){
	double s=0;
	for (int i = 0; i < n; i++) s += v[i] * v[i];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			P[i][j] = -2 * v[i] * v[j] / s;
	for (int i = 0; i < n; i++) P[i][i] += 1;
}
void QR(Mat& A, Mat &Q){
	double s, t;
	vec v;
	Mat H,R,I;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Q[i][j] = 0; 
	for (int i = 0; i < n; i++) Q[i][i] = 1;
	// Q=1
	for (int k = 0; k < n - 1; k++){
		t = 0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++){
				I[i][j] = Q[i][j];
				R[i][j] = A[i][j];
		}
		//I=Q; R=A;
		for (int j = k; j < n; j++) t += A[j][k]*A[j][k];
		s = sqrtf(t);
		if (A[k][k] < 0) s = -s;
		for (int j = 0; j < k; j++) v[j] = 0;
		v[k] = A[k][k] + s;
		for (int j = k + 1; j < n; j++) v[j] = A[j][k];
		Hausholder(H, v);
		Matprod(I, H, Q);	//Q=I*H
		Matprod(H, R, A);	//A=H*R
	}
}
void findEigen(Mat& A, vec& b, Mat& B){
	Mat Q, P ,R;
	int p=0;
	double s;
	char* str = new char[20];
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++)
			if(i!=j) P[i][j] = 0;
		P[i][i] = 1;
	}
	//P=1
	ofstream fout("debug.txt");
	fout << "Матриці А в QR методі:"<<endl;
	while (true){
		p++;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				R[i][j] = A[i][j];
		//R=A
		QR(R, Q);
		Matprod(P,Q, B);
		Matprod(R, Q, A);
		//виведення проміжних данних
		_itoa(p, str, 10);
		strcat(str," iteration :");
		outMat("debug.txt",str, &A);
		//перевірка виходу
		s = 0;
		for (int k = 0; k < n - 1; k++)
			for (int i = k + 1; i < n; i++)	//s=||x||^2 
				s += A[i][k] * A[i][k];
		if (s < eps*eps) break;
		for (int i = 0; i < n; i++)	//P=B
			for (int j = 0; j < n; j++)
				P[i][j] = B[i][j];
	}
	for (int i = 0; i < n; i++) b[i] = A[i][i];
}
void correct(Mat A, vec b, Mat B){
	vec x;
	char * c = new char[10];
	for (int i = 0; i < n; i++){
		_itoa(i + 1, c, 10);
		strcat(c, "-а Невязка:\n");
		for (int j = 0; j < n; j++){
			x[j] = 0;
			for (int k = 0; k < n; k++)
				x[j] += A[j][k] * B[k][i];
			x[j] = x[j] - b[i] * B[j][i];
		}
		outVec("ans.txt", c, &x, true);
	}
}
void fileMat(ofstream &fout){
	fout.width(zn + 10);
	fout.precision(zn);
	fout.setf(ios::right);
	fout.setf(ios::fixed);
}
void inMat(char* s, Mat* A){
	ifstream stm(s);
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			stm >> (*A)[i][j];
	stm.close();
};
void outMat(char* s,char* msg, Mat* A){
	ofstream stm(s,ios_base::app);
	fileMat(stm);
	stm << msg << "\n";
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++)
			stm << (*A)[i][j] << "\t";
		stm << endl;
	}
	stm << "----------------------------\n";
};
void outVec(char* s, char* msg, vec* x, bool b){
	ofstream stm(s,ios_base::app);
	stm << msg << "\n";
	if (b) stm.setf(ios::scientific);
	int i;
	for (i = 0; i < n; i++)
		stm << (*x)[i] << endl;
	stm << endl;
}