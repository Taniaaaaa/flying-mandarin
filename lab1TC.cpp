// lab1TC.cpp: определ€ет точку входа дл€ консольного приложени€.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <windows.h>

using namespace std;
#define PI 3.14159265
double sigma_kv, m;


double Random()// ~ U([0;1])
{
	static const unsigned int A = 1686629717;
	static const unsigned int C = 907633385;
	static unsigned int n = GetTickCount();

	return ((double)(n = A*n + C)) / 0xFFFFFFFF;
}

double BoxMuller(bool standart, double expected_value, double dispersion)//if standart=1, then number ~ N(0,1), else ~N(a, sigma^2)
{
	double z0, z1, r, f;
	r = Random();
	f = Random();
	z0 = cos(2 * PI*f)*sqrt(2 * log(1 / r));
	//z1 = sin(2 * PI*f)*sqrt(2 * log(1 / r));

	if (standart == 1) return z0;
		
	return (expected_value + z0*sqrt(dispersion));
}

double* task1_a(double* x, int n, double sigma)//unknown parametr is a, sigma^2 == 1
{//build a confidence interval for a
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

	cout << a[0] << " < expected_value <" << a[1] << endl;
	return a;
}

double* task1_b(double* x, int n)//unknown parametr is a
{//build a confidence interval for a

	double sig_n = 0, x_ser = 0;
	double a[2];
	double z_hamma;//hamma=0,05
	if (n == 30) z_hamma = 1.697;
	if (n == 60) z_hamma = 1.671;

	for (int i = 0; i < n; i++)
	{
		x_ser += x[i];
	}
	x_ser /= n;

	for (int i = 0; i < n; i++)
	{
		sig_n += (x[i] - x_ser)*(x[i] - x_ser);
	}
	sig_n /= n;


	a[0] = x_ser - (z_hamma*sig_n) / (sqrt(n - 1));//lower border
	a[1] = x_ser + (z_hamma*sig_n) / (sqrt(n - 1));//higher border

	cout << a[0] << " < expected_value <" << a[1] << endl;
	return a;
}

double* task1_c(double* x, int n)//build a confidence interval for a
{
	double x_ser = 0, sig_n = 0;
	double a[2];
	const double z_hamma = 1.96;//hamma=0,05

	for (int i = 0; i < n; i++)
	{
		x_ser += x[i];
	}
	x_ser /= n;

	for (int i = 0; i < n; i++)
	{
		sig_n += (x[i] - x_ser)*(x[i] - x_ser);
	}
	sig_n /= n;

	a[0] = x_ser - (z_hamma*sig_n) / (sqrt(n));//lower border
	a[1] = x_ser + (z_hamma*sig_n) / (sqrt(n));//higher border

	cout << a[0] << " < expected_value <" << a[1] << endl;
	return a;
}

/*
double* task_2(double* x, int n)//confidence interval for sigma^2 in ~N(0,1)
{
double sig_n = 0, x_ser = 0;
double dispersion[2];
double z_hamma1, z_hamma2;//hamma=0,05
//if (n == 30)
z_hamma1 = 16.791;
z_hamma2 = 46.979;

for (int i = 0; i < n; i++)
{
x_ser += x[i];
}
x_ser /= n;

for (int i = 0; i < n; i++)
{
sig_n += (x[i] - x_ser)*(x[i] - x_ser);
}

dispersion[0] = double(sig_n)/z_hamma2;//lower border
dispersion[1] = double(sig_n) / z_hamma1;//higher border

cout << dispersion[0] << " < dispersion <" << dispersion[1] << endl;
return dispersion;
}
*/


double* conf_interval(int nSize)// sigma^2 == 1
{//build a trusting interval for sigma^2
	double sig_n = 0, x_ser = 0;
	double dispersion[2];
	double Arr[30];

	//const double z_hamma = 1.96;//hamma=0,05
	//if (n == 30) 
	const double z_hamma1 = 16.791;
	const double z_hamma2 = 46.979;
	for (int i = 0; i < 30; ++i)
	{

		Arr[i] = BoxMuller(1, 0, 1);
		x_ser += Arr[i];
		
	}
	x_ser /= nSize;

	for (int i = 0; i < nSize; i++)
	{
		sig_n += (Arr[i] - x_ser)*(Arr[i] - x_ser);

	}
	
		
	dispersion[0] = double(sig_n) / z_hamma2;//lower border 2
	dispersion[1] = double(sig_n) / z_hamma1;//higher border 1

   	cout << dispersion[0] << " < sigma^2 <" << dispersion[1] << endl;
	return dispersion;
}

double Weibull(double alpha, double beta)
{
	double result, w;
	w = Random();
	result = pow(-log(w), 1 / beta) / alpha;
	return result;
}

double F(double u, double alpha)//u >= 0
{
	if (u < 0)
	{
		cout << "Error. u must be positive!";
		return -1;
	}
	else
	{
		double result;
		result = (1 - exp(-(alpha*u)*(alpha*u)));
		return result;
	}
}

double G(double u)//u >= 0
{
	if (u < 0)
	{
		cout << "Error. u must be positive!";
		return -1;
	}
	else
	{
		double result;
		result = (1 - exp(-u*u*u));
		return result;
	}
}

double I_firstWay(int n)
{
	double result = 0;
	double* mas = new double[n];
	for (int i = 0; i < n; i++)
	{
		mas[i] = Weibull(1, 3);
		mas[i] = F(mas[i], 0.1);
	}
	for (int i = 0; i < n; i++)
	{
		result += mas[i];
	}
	result = result / n;
	return result;
}

double I_secondWay(int n)
{
	double result = 0;
	double* mas = new double[n];
	for (int i = 0; i < n; i++)
	{
		mas[i] = Weibull(0.1, 2);
		mas[i] = G(mas[i]);
		mas[i] = 1 - mas[i];
	}
	for (int i = 0; i < n; i++)
	{
		result += mas[i];
	}
	result = result / n;
	return result;
}



int main()
{
	setlocale(LC_ALL, "Rus");
	/*cout << "¬ычисление первым способом:" << endl;
	for (int i = 0; i < 10; i++)
	cout << I_firstWay(100000) << " ";
	cout << endl;
	cout << "¬ычисление вторым способом:" << endl;
	for (int i = 0; i < 10; i++)
	cout << I_secondWay(100000) << " ";
	cout << endl;
	*/

	for (int i = 0; i < 30; ++i)
		cout << BoxMuller(1, 0.0, 1.0) << ' ';

	cout << endl;
	conf_interval(30);

	//delete[] Arr; 

	//int n = 30;
	//int l = 2;
	//double* x_array = new double[n];

	//for (int i = 0; i < 30; i++)
	//{
	//	x_array[i] = BoxMuller(0, 3., 4.4);
	//	//cout << x_array[i] << "  ";
	//}
	////cout << endl;

	//double* result = new double[l];
	//for (int i = 0; i < l; i++)
	//	result[i] = 0;

	//result = task1_c(x_array, n);

	////cout << result[0] << endl << result[1] << endl;
	////delete[] result;
	//delete[] x_array;

	system("pause");
	return 0;

}

