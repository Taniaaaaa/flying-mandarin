// conf_inter.cpp: определяет точку входа для консольного приложения.
//

//Task 1

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;
#define PI 3.14159265
double sigma_kv, m;

//const double z_hamma = 1.96;//hamma=0,05

/*double generator()// ~ N(0, 1)
{
int k = 1000;
double S = 0.;
for (int i = 0; i < k; ++i)
S += rand() % 10 + 1;
S = (double)(S - 5.5*k) / (sqrt(8.25*k));
return S ;
}
*/

//double *z0;
const int nSize = 30;
double Arr[nSize];
//double Arr[] = { 0.3781, 0.69, 1.28, 0.38, 0.63, 0.69, 1.28, 1.41, 0.97, 1.28, 0.97, 1.41, 1.57 };
//double Arr[] = { 0.34, -1.42, -1.26, 0.34, -1.30, 1.78, 0.086, 1.28, -0.76, -1.37, -0.13, 1.36, -0.49, -1.38, 0.26, -0.68, 0.1, -0.75, 0.89, 0.63 };
//double Arr[] = { 0.34, -1.42, -1.26, 0.34, -1.30, 1.78, 0.086, 1.28, -0.76, -1.37, -0.13, 1.36, -0.49 };
double *pAr1 = new double[nSize];

double Random()// ~ U([0;1])
{
	static const unsigned int A = 1686629717;
	static const unsigned int C = 907633385;
	static unsigned int n = GetTickCount();

	return ((double)(n = A*n + C)) / 0xFFFFFFFF;
}

double BoxMuller(bool standart)//if standart=1, then number ~ N(0,1), else ~N(a, sigma^2)
{
	//double Arr[nSize];
	double z0, z1, r, f;
	r = Random();
	f = Random();
    z0 = cos(2 * PI*f)*sqrt(2 * log(1 / r));
	* pAr1 = z0;
	//z1 = sin(2 * PI*f)*sqrt(2 * log(1 / r));
	if (standart == 1) return * pAr1;
	

}

double* task_a(double* x, int n, double sigma)//unknown parametr is a, sigma^2 == 1
{//build a trusting interval for a
	double x_ser = 0;
	double a[2];
	const double z_hamma = 1.96;//hamma=0,05

	for (int i = 0; i < n; i++)
	{
		x_ser += x[i];
	}
	x_ser /= n;

	a[0] = x_ser - (z_hamma*sigma) / (sqrt(n));//lower border
	a[1] = x_ser + (z_hamma*sigma) / (sqrt(n));//higher border

	cout << a[0] << " < a <" << a[1] << endl;
	return a;
}




double* conf_interval(double * pAr, int nSize)//unknown parametr is a, sigma^2 == 1
{//build a trusting interval for a
	double sig_n = 0, x_ser = 0;
	double dispersion[2];
	//const double z_hamma = 1.96;//hamma=0,05
	//if (n == 30) 
	const double z_hamma1 = 16.791;
	const double z_hamma2 = 46.979;
			
	for (int i = 0; i < nSize; i++)
	{
		Arr[i] = *pAr1;
		x_ser += Arr[i];
		sig_n += (Arr[i] - x_ser)*(Arr[i] - x_ser);
		
	}
	x_ser /= nSize;
	sig_n /= nSize;


	//for (int i = 0; i < nSize; i++)
	//{
	//	//Arr[i] = *pAr1;
	//	
	//}
	//

	dispersion[0] = double(sig_n) / z_hamma2;//lower border 2
	dispersion[1] = double(sig_n) / z_hamma1;//higher border 1

	//dispersion[0] = x_ser - (z_hamma*sigma) / (sqrt(nSize));//lower border
	//dispersion[1] = x_ser + (z_hamma*sigma) / (sqrt(nSize));//higher border

	cout << dispersion[0] << " < sigma^2 <" << dispersion[1] << endl;
	return dispersion;
}


int main()
{
	for (int i = 0; i < nSize; i++)
		cout << rand() << "  ";
	cout << endl;
	for (int i = 0; i < nSize; ++i)
		cout << BoxMuller(1) << ' ';
		
	cout << endl;
    conf_interval(pAr1, nSize);


	system("pause");
	return 0;
}


