// LongArithmetics.cpp: определяет точку входа для консольного приложения.
//

#pragma once
#include "stdafx.h"

#define E_LOWER_VALUE "200000000"
#define E_UPPER_VALUE "1000000000000000000000000000000"
const int n = 50;

using namespace std;

/* 
Function convergent

This function calculates continued fraction's 
denominator for specified E and N.

The answer is stored two LINT array Q[n]

Arguments:

E - the first part of public key		
N - the second part of public key
*Q - pointer to denominator array
		
Return value: none
 */




void convergent(LINT E, LINT N, LINT* Q)
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


/*
Function extended_euclid

This function calculates coefficients a, b and GCD(a,b) 
in comparison a*x+b*y=GCD(a,b)

The answer is stored three LINT: x,y,d

Arguments:

a - the first coefficient in comparison
b - the second coefficient in comparison
x - the first variable in comparison
y - the second variable in comparison
d - greatest common divisor of a and b

Return value: none
*/

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




/*
Function Vinere

This function calculates the private part of key 
with using open part of key (E,N)

The answer is the private part of key D

Arguments:

E - the first part of public key
N - the second part of public key
D - the part of private key

Return value: none
*/

void Vinere(
	__in LINT E, 
	__in LINT N,
	__out LINT *D)
{
	LINT potential_D[n];
	LINT limitD = root(root(N/3)) - 1;

	convergent(E, N, potential_D);

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




/*
Function Vulnerable_Generator

This function recieves public part of key from
couple p,q

The answer is the public part of key E,N

Arguments:

primes_vector - vector primes
*E - the first part of public key
*N - the second part of public key

Return value: none
*/

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

