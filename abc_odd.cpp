//-----------------------------------------------------------------------------
//     Project:     Implementing Simple Matrix Encryption
//-----------------------------------------------------------------------------

//******************************************************************************
//      Copyright (C) 2024   Tao Chengdong
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
//(at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file LICENCE) for more details.
//*******************************************************************************

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define WORD int
#define FIELD 128

#define AROW 2
#define ACOL 2
#define BROW 2
#define BCOL 2
#define CROW 2
#define CCOL 2
#define VARIABLE 4
#define EQUATION 8
#define CENTRAL_MAP_SIZE 10
#define PUBLIC_KEY_SIZE 80
#define SECRET_KEY_SIZE 192

WORD add(WORD a, WORD b)
{
	return ((a + b) % FIELD);
}

WORD sub(WORD a, WORD b)
{
	return ((a - b) % FIELD);
}

WORD mul(WORD a, WORD b)
{
	/* multiply two elements of GF(2^m)*/
	return ((a * b) % FIELD);
}

WORD inv(WORD a, WORD b)
{
	int t, nt, r, nr, q, tmp;
	if (b < 0)
		b = -b;
	if (a < 0)
		a = b - (-a % b);
	t = 0;
	nt = 1;
	r = b;
	nr = a % b;
	while (nr != 0)
	{
		q = r / nr;
		tmp = nt;
		nt = t - q * nt;
		t = tmp;
		tmp = nr;
		nr = r - q * nr;
		r = tmp;
	}
	if (r > 1)
		return -1; /* No inverse */
	if (t < 0)
		t += b;
	return t;
}

// compute  A^-1
int matrixinv(WORD a[], int n)
{
	int *is, *js, i, j, k, l, u, v;
	WORD d, p;
	is = (int *)malloc(n * sizeof(int));
	js = (int *)malloc(n * sizeof(int));
	for (k = 0; k <= n - 1; k++)
	{
		d = 0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j;
				p = (a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d == 0)
		{
			free(is);
			free(js);
			printf("err**not inv\n");
			return (0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		l = k * n + k;
		a[l] = inv(a[l], FIELD);
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = mul(a[u], a[l]);
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = add(a[u], (mul(a[i * n + k], a[k * n + j])));
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k;
				a[u] = mul(a[u], a[l]);
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}
	free(is);
	free(js);
	return (1);
}

// A is m*n matrix, B is n*k matrix, C=A*B
void matrixmul(WORD a[], WORD b[], int m, int n, int k, WORD c[])
{
	int i, j, l, u;
	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= k - 1; j++)
		{
			u = i * k + j;
			c[u] = 0;
			for (l = 0; l <= n - 1; l++)
				c[u] = add(c[u], (mul(a[i * n + l], b[l * k + j])));
		}
	return;
}

// compute transpose(A)
int matrixtranspose(WORD *a, int n, WORD *b)
{
	int i, j;
	for (i = 0; i <= n - 1; i++)
	{
		for (j = 0; j <= n - 1; j++)
			b[i * n + j] = a[j * n + i];
	}
	return (1);
}

// c=(a_1*x_1+...+a_n*x_n)*(b_1*x_1+...+b_n*x_n)
// input: a=[a_1,..,a_n] and b=[b_1,...,b_n]
// output: the coefficents of c=[c_1,...]
int tensorproduct(WORD a[], WORD b[], WORD c[], int n, int r)
{
	int i, j, k, temp = 0;
	WORD t = 0;
	for (i = 0; i < n; i++)
	{
		for (j = i; j < n; j++)
		{
			if (i == j)
			{
				for (k = 0; k < r; k++)
				{
					t = add(t, mul(a[k * n + i], b[k * n + i]));
				}
				c[temp++] = t;
				t = 0;
			}
			if (i != j)
			{
				for (k = 0; k < r; k++)
				{
					t = add(t, add(mul(a[k * n + i], b[k * n + j]), mul(a[k * n + j], b[k * n + i])));
				}
				c[temp++] = t;
				t = 0;
			}
		}
	}
	return (1);
}

// compute the reduce echelon form of a matrix
int echelonform(WORD M[], int nrows, int ncols)
{
	int lead = 0, rix = 0, iix = 0, j, k;
	WORD temp, de, sb;
	for (rix = 0; rix < nrows; rix++)
	{
		if (lead >= ncols)
		{
			return 0;
		}
		iix = rix;
		while (M[iix * ncols + lead] == 0)
		{
			iix++;

			if (iix == nrows)
			{
				iix = rix;
				lead++;
				if (ncols == lead)
				{
					return 0;
				}
			}
		}
		for (j = 0; j < ncols; j++)
		{
			temp = M[rix * ncols + j];
			M[rix * ncols + j] = M[iix * ncols + j];
			M[iix * ncols + j] = temp;
		}
		de = M[rix * ncols + lead];
		if (de != 0)
		{
			for (j = 0; j < ncols; j++)
			{
				M[rix * ncols + j] = mul(M[rix * ncols + j], inv(de, FIELD));
			}
		}
		for (j = 0; j < nrows; j++)
		{
			if (j != rix)
			{
				sb = M[j * ncols + lead];
				for (k = 0; k < ncols; k++)
				{
					M[j * ncols + k] = sub(M[j * ncols + k], mul(sb, M[rix * ncols + k]));
				}
			}
		}
		lead++;
	}

	return 1;
}

// key pair generation

int keypair(WORD *sk, WORD *pk)
{
	WORD S[VARIABLE * VARIABLE], INVS[VARIABLE * VARIABLE], ST[VARIABLE * VARIABLE];
	WORD M[VARIABLE * VARIABLE], U[VARIABLE * VARIABLE], V[VARIABLE * VARIABLE];
	WORD T[EQUATION * EQUATION], INVT[EQUATION * EQUATION];
	WORD B[BROW * BCOL * VARIABLE], C[CROW * CCOL * VARIABLE], A[AROW * ACOL * VARIABLE];
	WORD TEMPA[ACOL * VARIABLE], TEMPB[BROW * VARIABLE], TEMPC[CROW * VARIABLE], TEMPF[CENTRAL_MAP_SIZE];
	WORD *F = new WORD[(AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE];
	WORD *FB = new WORD[(AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE];
	int i, j, s, t, ij = 0, eof, count, flag = 0; //
	WORD W = 0;

	printf("************************************************* \n");
	printf("Start ABC key pair generation \n");
	printf("************************************************* \n");

	// create L1 and check that L1 is not singular
	// initialize linear transformation  L1 in S
	eof = 0;
	while (eof == 0)
	{
		for (i = 0; i < VARIABLE; i++)
		{
			for (j = 0; j < VARIABLE; j++)
			{
				INVS[i * VARIABLE + j] = S[i * VARIABLE + j] = rand() % FIELD;
			}
		}
		eof = matrixinv(INVS, VARIABLE);
		//  if ( eof ) cout << "L1 was singular \n" ;
	}

	// create L2 and check that L2 is not singular
	// initialize linear transformation  L2 in T
	// and also put T^-1 on sk
	eof = 0;
	while (eof == 0)
	{
		for (i = 0; i < EQUATION; i++)
		{
			for (j = 0; j < EQUATION; j++)
			{
				INVT[i * EQUATION + j] = T[i * EQUATION + j] = rand() % FIELD;
			}
		}
		eof = matrixinv(INVT, EQUATION);
		//  if ( eof ) cout << "L1 was singular \n" ;
	}

	//  store T^(-1) in secrete key sk
	for (i = 0; i < EQUATION; i++)
	{
		for (j = 0; j < EQUATION; j++)
		{
			sk[ij++] = INVT[i * EQUATION + j];
		}
	}

	// generate the central map
	// generate A
	for (i = 0; i < AROW * ACOL * VARIABLE; i++)
	{
		A[i] = rand() % FIELD;
	}

	// generate B and stor B to sk
	for (i = 0; i < BROW * BCOL * VARIABLE; i++)
	{
		sk[ij++] = B[i] = rand() % FIELD;
	}

	// generate C and stor C to sk
	for (i = 0; i < CROW * CCOL * VARIABLE; i++)
	{
		sk[ij++] = C[i] = rand() % FIELD;
	}

	// the central map F
	// compute A*B
	count = 0;
	flag = 0;
	for (i = 0; i < AROW; i++)
	{
		for (j = 0; j < ACOL * VARIABLE; j++)
		{
			TEMPA[j] = A[i * ACOL * VARIABLE + j];
		}
		for (j = 0; j < BCOL; j++)
		{
			for (s = 0; s < BROW; s++)
			{
				for (t = 0; t < VARIABLE; t++)
				{
					TEMPB[s * VARIABLE + t] = B[j * VARIABLE + s * BCOL * VARIABLE + t];
				}
			}
			tensorproduct(TEMPA, TEMPB, TEMPF, VARIABLE, BROW);
			for (s = 0; s < VARIABLE; s++)
			{
				for (t = s; t < VARIABLE; t++)
				{
					sk[ij++] = F[flag++] = TEMPF[count++];
				}
			}
			count = 0;
		}
	}

	// compute A*C
	for (i = 0; i < AROW; i++)
	{
		for (s = 0; s < ACOL; s++)
		{
			for (t = 0; t < VARIABLE; t++)
			{
				TEMPA[s * VARIABLE + t] = A[i * ACOL * VARIABLE + s * VARIABLE + t];
			}
		}
		for (j = 0; j < CCOL; j++)
		{
			for (s = 0; s < CROW; s++)
			{
				for (t = 0; t < VARIABLE; t++)
				{
					TEMPC[s * VARIABLE + t] = C[j * VARIABLE + s * CCOL * VARIABLE + t];
				}
			}
			tensorproduct(TEMPA, TEMPC, TEMPF, VARIABLE, CROW);
			for (s = 0; s < VARIABLE; s++)
			{
				for (t = s; t < VARIABLE; t++)
				{
					sk[ij++] = F[flag++] = TEMPF[count++];
				}
			}
			count = 0;
		}
	}

	// compute FoS=Transpose(S)*F_i*S, i in 0...i<AROW*(BCOL+CCOL)-1
	matrixtranspose(S, VARIABLE, ST);
	count = 0;
	flag = 0;
	int num = 0;
	for (i = 0; i < AROW * (BCOL + CCOL); i++)
	{
		for (s = 0; s < VARIABLE * VARIABLE; s++)
		{
			M[s] = WORD(0);
		}
		for (s = 0; s < VARIABLE; s++)
		{
			for (t = s; t < VARIABLE; t++)
			{
				M[s * VARIABLE + t] = F[count++];
			}
		}
		matrixmul(ST, M, VARIABLE, VARIABLE, VARIABLE, U);
		matrixmul(U, S, VARIABLE, VARIABLE, VARIABLE, V);

		for (s = 0; s < VARIABLE; s++)
		{
			for (t = s; t < VARIABLE; t++)
			{
				if (t == s)
				{
					FB[flag++] = V[s * VARIABLE + t];
				}
				if (t != s)
				{
					FB[flag++] = add((WORD)V[s * VARIABLE + t], (WORD)V[t * VARIABLE + s]);
				}
			}
		}
	}

	// compute ToFoS
	matrixmul(T, FB, EQUATION, EQUATION, CENTRAL_MAP_SIZE, pk);

	// store S^-1 to the secret key sk
	for (i = 0; i < VARIABLE; i++)
	{
		for (j = 0; j < VARIABLE; j++)
		{
			sk[ij++] = INVS[i * VARIABLE + j];
		}
	}
	delete F;
	delete FB;
	return (1);
}

int encryption(WORD *pk, WORD *plaintext, WORD *ciphertext, int plaintextlen)
{
	int i, j, k, ij = 0, eof = 0;
	WORD temp, x[VARIABLE];
	if (plaintextlen != VARIABLE)
		return (0);

	printf("************************************************* \n");
	printf("Start ABC encryption \n");
	printf("************************************************* \n");

	for (k = 0; k < VARIABLE; k++)
		x[k] = WORD(plaintext[k]);

	for (k = 0; k < EQUATION; k++)
	{
		temp = 0;
		for (i = 0; i < VARIABLE; i++)
		{
			for (j = i; j < VARIABLE; j++)
			{
				temp = add(temp, mul(WORD(pk[ij++]), mul(x[i], x[j])));
			}
		}
		ciphertext[k] = temp;
	}
	return (1);
}

int decryption(WORD *sk, WORD *decrypttext, WORD *ciphertext, int ciphertextlen)
{
	int i, j, k, ij = 0, eof = 0;
	WORD x[AROW * ACOL + VARIABLE], y[EQUATION], y1[EQUATION], y2[EQUATION], X[VARIABLE];
	WORD coeff[ACOL * (BCOL + CCOL) * (ACOL * AROW + VARIABLE)];
	WORD TP, xn;
	if (ciphertextlen != EQUATION)
		return (0);

	printf("************************************************* \n");
	printf("Start ABC decryption \n");
	printf("************************************************* \n");

	for (k = 0; k < EQUATION; k++)
		y[k] = WORD(ciphertext[k]);
	for (k = 0; k < AROW * ACOL + VARIABLE; k++)
		x[k] = WORD(0);
	for (i = 0; i < ACOL * (BCOL + CCOL) * (ACOL * AROW + VARIABLE); i++)
	{
		coeff[i] = WORD(0);
	}

	WORD Sqrttable[FIELD];
	for (i = 0; i < FIELD; i++)
	{
		Sqrttable[i] = mul(i, i);
	}

	// compute first T^{-1}
	k = 0;
	for (i = 0; i < EQUATION; i++)
	{
		y1[i] = 0;
		for (j = 0; j < EQUATION; j++)
		{
			y1[i] = add(y1[i], mul(WORD(sk[ij++]), y[j]));
		}
	}

	// next find F^{-1}
	for (k = 0; k < ACOL; k++)
	{
		for (i = 0; i < BCOL; i++)
		{
			for (j = 0; j < AROW; j++)
			{
				coeff[(k * BCOL + i) * (ACOL * AROW + VARIABLE) + k * AROW + j] = y1[j * BCOL + i];
			}
		}
		for (i = 0; i < CCOL; i++)
		{
			for (j = 0; j < AROW; j++)
			{
				coeff[ACOL * BCOL * (ACOL * AROW + VARIABLE) + (k * CCOL + i) * (ACOL * AROW + VARIABLE) + k * AROW + j] = y1[AROW * BCOL + j * CCOL + i];
			}
		}
	}
	for (i = 0; i < ACOL * BCOL; i++)
	{
		for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++)
		{
			coeff[i * (AROW * ACOL + VARIABLE) + j] = WORD(sk[ij++]);
		}
	}

	for (i = 0; i < ACOL * CCOL; i++)
	{
		for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++)
		{
			coeff[ACOL * BCOL * (AROW * ACOL + VARIABLE) + i * (AROW * ACOL + VARIABLE) + j] = WORD(sk[ij++]);
		}
	}

	echelonform(coeff, ACOL * (BCOL + CCOL), AROW * ACOL + VARIABLE);

	/*
		printf( "\n coeff: \n");
		for (i=0; i<ACOL*(BCOL+CCOL)*(ACOL*AROW+VARIABLE); i++){
			printf( "%d ", coeff[i]) ;
			if ((i+1)%((ACOL*AROW+VARIABLE))==0){	printf("\n");printf("\n");}
		}
	*/

	WORD WM[VARIABLE];
	for (i = 0; i < VARIABLE - 1; i++)
	{
		WM[i] = coeff[(i + AROW * ACOL) * (ACOL * AROW + VARIABLE) + (ACOL * AROW + VARIABLE - 1)];
	}
	WM[VARIABLE - 1] = 1;
	for (k = 0; k < EQUATION; k++)
	{
		TP = 0;
		for (i = 0; i < VARIABLE; i++)
		{
			for (j = i; j < VARIABLE; j++)
			{
				TP = add(TP, mul(sk[ij++], mul(WM[i], WM[j])));
			}
		}
		y2[k] = TP;
	}

	for (i = 0; i < EQUATION; i++)
	{
		if (y2[i] != 0)
		{
			break;
		}
	}

	TP = mul(inv(y2[i], FIELD), y1[i]);

	// printf("\n\ni:= %d ", i);

	for (i = 0; i < FIELD; i++)
	{
		if (Sqrttable[i] == TP)
		{
			xn = (WORD)(i);
		}
	}
	//	printf("\n\ni:= %d ", i);

	X[VARIABLE - 1] = xn;
	for (i = 0; i < VARIABLE - 1; i++)
	{
		X[i] = mul(WM[i], X[VARIABLE - 1]);
	}

	// compute S^-1
	for (i = 0; i < VARIABLE; i++)
	{
		decrypttext[i] = 0;
		for (j = 0; j < VARIABLE; j++)
		{
			decrypttext[i] = add(decrypttext[i], mul(WORD(sk[ij++]), X[j]));
			//	if (decrypttext[i]<0) {decrypttext[i] = decrypttext[i] +FIELD;}
		}
	}
	return (1);
}

int main(int argc, int argv[])
{
	int plaintextlen = VARIABLE;
	int i;
	WORD *sk = new WORD[SECRET_KEY_SIZE];
	WORD *pk = new WORD[PUBLIC_KEY_SIZE];
	WORD plaintext[VARIABLE];
	WORD ciphertext[EQUATION];
	WORD decrypttext[EQUATION];

	srand((int)time(0));

	// test
	WORD t1 = 2, t2 = 3, t3;
	printf("\n t3: ");
	t3 = inv(t1, FIELD);
	printf("%d \n", t3);

	WORD ta[9], tb[9], tc[9];
	for (i = 0; i < 9; i++)
	{
		tb[i] = ta[i] = rand() % FIELD;
	}
	printf("\n ta: \n");
	for (i = 0; i < 9; i++)
	{
		printf("%d ", ta[i]);
		if ((i + 1) % 3 == 0)
		{
			printf("\n");
		}
	}
	printf("\n");

	matrixinv(tb, 9);
	printf("\n tb: \n");
	for (i = 0; i < 9; i++)
	{
		printf("%d ", tb[i]);
		if ((i + 1) % 3 == 0)
		{
			printf("\n");
		}
	}
	printf("\n");

	matrixmul(ta, tb, 9, 9, 9, tc);

	printf("\n tc: \n");
	for (i = 0; i < 9; i++)
	{
		printf("%d ", tc[i]);
		if ((i + 1) % 3 == 0)
		{
			printf("\n");
		}
	}
	printf("\n");

	keypair(sk, pk);
	printf("\n plaintext: \n");

	for (i = 0; i < VARIABLE; i++)
	{
		plaintext[i] = rand() % FIELD;
		printf("%d ", plaintext[i]);
	}
	printf("\n");
	encryption(pk, plaintext, ciphertext, plaintextlen);
	printf("\n ciphertext: ");
	for (i = 0; i < EQUATION; i++)
	{
		printf("%d ", ciphertext[i]);
	}
	printf("\n");
	decryption(sk, decrypttext, ciphertext, EQUATION);
	printf("\n decrypttext: \n");
	for (i = 0; i < VARIABLE; i++)
	{
		printf("%d ", decrypttext[i]);
	}
	printf("\n");
	delete sk;
	delete pk;
	return (1);
}