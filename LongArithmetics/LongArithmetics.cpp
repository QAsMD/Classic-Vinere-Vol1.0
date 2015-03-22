// LongArithmetics.cpp: определяет точку входа для консольного приложения.
//

#pragma once
#include "stdafx.h"

#define E_LOWER_VALUE "200000000"
#define E_UPPER_VALUE "1000000000000000000000000000000"
const int n = 50;

using namespace std;

///
/// Function podh_drob
///
/// This function calculates continued fraction's 
/// denominator for specified E and N.
///
/// The answer is stored two LINT array Q[n]
///
/// Arguments:
///
///		E - public key
///		
///		N - ...
///
///		*Q - pointer to denominator array
///		
/// Return value: none
///

void podh_drob(LINT E, LINT N, LINT* Q)
{
	LINT V[n];
	LINT H[n];
	LINT Z[n];

	for (int k = 0; k < n; k++) //занулили числители и знаменатели подходящих дробей
	{
		H[k] = Z[k] = V[k] = Q[k] = 0;
	}

	H[0] = E;
	Z[0] = N;

	for (int i = 0; i<n; i++)
	{
		if (Z[i] == 0) { break; }
		else
		{
			V[i] = H[i] / Z[i];


			Z[i + 1] = H[i] - Z[i] * V[i];
			H[i + 1] = Z[i];
		}
	}
	Q[0] = 1;
	Q[1] = V[1];

	for (int i = 2; i < n - 1; i++)
	{
	
		//cout << endl << "V[" << i << "]= " << V[i].decstr();
		if ((V[i]) == 0) break;
		Q[i] = Q[i - 1] * V[i] + Q[i - 2];
		cout << "Q[i]: " << Q[i].decstr() << endl;
	}
}

void extended_euclid(
	__in	LINT a, 
	__in	LINT b, 
	__out	LINT *x, 
	__out	LINT *y, 
	__out	LINT *d)


/* calculates a * *x + b * *y = gcd(a, b) = *d */

{

	LINT q, r, x1, x2, y1, y2;
if (b == 0) 
	{
		*d = a, *x = 1, *y = 0;
		return;
	}
	x2 = 1, x1 = 0, y2 = 0, y1 = 1;

while (b > 0) 
	{
		q = a / b, r = a - q * b;
		*x = x2 - q * x1, *y = y2 - q * y1;
		a = b, b = r;
		x2 = x1, x1 = *x, y2 = y1, y1 = *y;
	}
	*d = a, *x = x2, *y = y2;
}

void Vinere(
	__in LINT E, 
	__in LINT N,
	__out LINT *D)
{
	LINT potential_D[n];
	LINT limitD = root(root(N/3)) - 1;

	podh_drob(E, N, potential_D);

	for (LINT M = 2; M < N; M++)
	{
		LINT LC = mexp(M, E, N);
		for (int i = 1; i < n; i++)
		{
			if (limitD < potential_D[i]) break;
			LINT M2 = mexp(LC, potential_D[i], N);

			if (M == M2)
			{
				*D = potential_D[i];
				return;
			}
		}
	}
	return;
}

void Vulnerable_Generator(
	__in vector <LINT> primes_vector, 
	__out LINT *E, 
	__out LINT *N)
{
	LINT p;
	LINT q;
	LINT NOD;
	LINT koef;

	for (int i = 0; i < (primes_vector.size() - 1); i++)
	{
		p = primes_vector[i];
		for (int j = i + 1; j < primes_vector.size(); j++)
		{
			q = primes_vector[j];

		    *N = mul(p,q);
			LINT eiler_func = mul((p - 1),(q - 1));
			
			LINT limitD = root(root((*N) / 3)) - 1;

			
			for (;limitD > 0; limitD--)
			{
				extended_euclid(limitD, eiler_func, E, &koef, &NOD); // E - ответ

				if ((*E > E_LOWER_VALUE) && (*E < E_UPPER_VALUE) && (NOD == 1) && (gcd(*E,eiler_func) == 1))
				{
					cout << "Eiler: " << eiler_func.decstr() << endl;
					cout << "Finished generating E, D and N." << endl;
					cout << "E is " << (*E).decstr() << endl;
					cout << "D is " << limitD.decstr() << endl;
					cout << "N is " << (*N).decstr() << endl;
					cout << "p is " << p.decstr() << endl;
					cout << "q is " << q.decstr() << endl;
					return;
				}
			}
		}
	}
}

int _tmain(int argc, _TCHAR* argv[])
{
	vector <LINT> primes_vector;
	ifstream FILE_read;
	FILE_read.open("11-20.txt");

	LINT E;
	LINT N;
	LINT D;

	int reader;

	while (FILE_read >> reader)
		primes_vector.push_back(reader);

	Vulnerable_Generator(primes_vector, &E, &N);

	cout << "Starting Vinere attack with E: " << E.decstr();
	cout << " and N: " << N.decstr() << endl;

	Vinere(E, N, &D);
	cout << "The key is: " << D.decstr() << endl;

	return 0;
}

