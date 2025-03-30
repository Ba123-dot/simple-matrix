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

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef unsigned char WORD;
#define FIELD 256

#define AROW 8
#define ACOL 8
#define BROW 8
#define BCOL 8
#define CROW 8
#define CCOL 8
#define VARIABLE 64
#define EQUATION 128
#define CENTRAL_MAP_SIZE 2080
#define PUBLIC_KEY_SIZE 266240
#define SECRET_KEY_SIZE 294912

// #define SECRET_KEY_SIZE 282624
/*
这里如果假设所有的row的变量为s的话那么
#define VARIABLE =s*s
#define EQUATION =2*s*s
#define CENTRAL_MAP_SIZE= s^2*(s^2+1)/2
#define PUBLIC_KEY_SIZE =2*s^2*s^2*(s^2+1)/2
#define SECRET_KEY_SIZE =6*(s^4)+2*s^2*s^2*(s^2+1)/2 // 当我们不需要S时，可以变成6*(s^4), 包括B、C的2s^2, T 的 (2s^2)^2=4s^4, 故而是6s^4.
*/

WORD Logtable[FIELD] = {
    0, 0, 25, 1, 50, 2, 26, 198, 75, 199, 27, 104, 51, 238, 223,
    3, 100, 4, 224, 14, 52, 141, 129, 239, 76, 113, 8, 200, 248, 105,
    28, 193, 125, 194, 29, 181, 249, 185, 39, 106, 77, 228, 166, 114, 154,
    201, 9, 120, 101, 47, 138, 5, 33, 15, 225, 36, 18, 240, 130, 69,
    53, 147, 218, 142, 150, 143, 219, 189, 54, 208, 206, 148, 19, 92, 210,
    241, 64, 70, 131, 56, 102, 221, 253, 48, 191, 6, 139, 98, 179, 37,
    226, 152, 34, 136, 145, 16, 126, 110, 72, 195, 163, 182, 30, 66, 58,
    107, 40, 84, 250, 133, 61, 186, 43, 121, 10, 21, 155, 159, 94, 202,
    78, 212, 172, 229, 243, 115, 167, 87, 175, 88, 168, 80, 244, 234, 214,
    116, 79, 174, 233, 213, 231, 230, 173, 232, 44, 215, 117, 122, 235, 22,
    11, 245, 89, 203, 95, 176, 156, 169, 81, 160, 127, 12, 246, 111, 23,
    196, 73, 236, 216, 67, 31, 45, 164, 118, 123, 183, 204, 187, 62, 90,
    251, 96, 177, 134, 59, 82, 161, 108, 170, 85, 41, 157, 151, 178, 135,
    144, 97, 190, 220, 252, 188, 149, 207, 205, 55, 63, 91, 209, 83, 57,
    132, 60, 65, 162, 109, 71, 20, 42, 158, 93, 86, 242, 211, 171, 68,
    17, 146, 217, 35, 32, 46, 137, 180, 124, 184, 38, 119, 153, 227, 165,
    103, 74, 237, 222, 197, 49, 254, 24, 13, 99, 140, 128, 192, 247, 112,
    7};

WORD Alogtable[FIELD] = {
    1, 3, 5, 15, 17, 51, 85, 255, 26, 46, 114, 150, 161, 248, 19,
    53, 95, 225, 56, 72, 216, 115, 149, 164, 247, 2, 6, 10, 30, 34,
    102, 170, 229, 52, 92, 228, 55, 89, 235, 38, 106, 190, 217, 112, 144,
    171, 230, 49, 83, 245, 4, 12, 20, 60, 68, 204, 79, 209, 104, 184,
    211, 110, 178, 205, 76, 212, 103, 169, 224, 59, 77, 215, 98, 166, 241,
    8, 24, 40, 120, 136, 131, 158, 185, 208, 107, 189, 220, 127, 129, 152,
    179, 206, 73, 219, 118, 154, 181, 196, 87, 249, 16, 48, 80, 240, 11,
    29, 39, 105, 187, 214, 97, 163, 254, 25, 43, 125, 135, 146, 173, 236,
    47, 113, 147, 174, 233, 32, 96, 160, 251, 22, 58, 78, 210, 109, 183,
    194, 93, 231, 50, 86, 250, 21, 63, 65, 195, 94, 226, 61, 71, 201,
    64, 192, 91, 237, 44, 116, 156, 191, 218, 117, 159, 186, 213, 100, 172,
    239, 42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 155, 182, 193, 88,
    232, 35, 101, 175, 234, 37, 111, 177, 200, 67, 197, 84, 252, 31, 33,
    99, 165, 244, 7, 9, 27, 45, 119, 153, 176, 203, 70, 202, 69, 207,
    74, 222, 121, 139, 134, 145, 168, 227, 62, 66, 198, 81, 243, 14, 18,
    54, 90, 238, 41, 123, 141, 140, 143, 138, 133, 148, 167, 242, 13, 23,
    57, 75, 221, 124, 132, 151, 162, 253, 28, 36, 108, 180, 199, 82, 246,
    1};

const WORD R1[VARIABLE] = {0};
const WORD R2[VARIABLE] = {1};
const WORD R3[VARIABLE] = {3};

WORD add(WORD a, WORD b) { return a ^ b; }

WORD sub(WORD a, WORD b) { return a ^ b; }

WORD mul(WORD a, WORD b)
{
    /* multiply two elements of GF(2^m)*/
    if (a && b)
        return Alogtable[(Logtable[a] + Logtable[b]) % (FIELD - 1)];
    else
        return 0;
}

WORD div(WORD a, WORD b)
{
    int j;
    if (b == 0)
    {
        printf("Division by zero\n");
        abort();
    }
    if (a == 0)
        return (0);

    if ((j = Logtable[a] - Logtable[b]) < 0)
        j += (FIELD - 1);

    return (Alogtable[j]);
}

WORD inv(WORD in)
{
    /* 0 is self inverting */
    if (in == 0)
        return 0;
    else
        return Alogtable[((FIELD - 1) - Logtable[in])];
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
            // printf("err**not inv\n");
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
        a[l] = inv(a[l]);
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
                        a[u] = a[u] ^ (mul(a[i * n + k], a[k * n + j]));
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

// A is m*n matrix, B is n*k matrix, C=A*B in this code A[i][j]=A[i*n+j]
void matrixmul(WORD a[], WORD b[], int m, int n, int k, WORD c[])
{
    int i, j, l, u;
    for (i = 0; i <= m - 1; i++)
        for (j = 0; j <= k - 1; j++)
        {
            u = i * k + j;
            c[u] = 0;
            for (l = 0; l <= n - 1; l++)
                c[u] = c[u] ^ (mul(a[i * n + l], b[l * k + j]));
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
                for (k = 0; k < r; k++) // 这个r 是怎么来的？
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
                    t = add(t, add(mul(a[k * n + i], b[k * n + j]),
                                   mul(a[k * n + j], b[k * n + i])));
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
                M[rix * ncols + j] = div(M[rix * ncols + j], de);
            }
        }
        for (j = 0; j < nrows; j++)
        {
            if (j != rix)
            {
                sb = M[j * ncols + lead];
                for (k = 0; k < ncols; k++)
                {
                    M[j * ncols + k] =
                        sub(M[j * ncols + k], mul(sb, M[rix * ncols + k]));
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
    WORD S[VARIABLE * VARIABLE], INVS[VARIABLE * VARIABLE],
        ST[VARIABLE * VARIABLE];
    WORD M[VARIABLE * VARIABLE], U[VARIABLE * VARIABLE], V[VARIABLE * VARIABLE];
    WORD T[EQUATION * EQUATION], INVT[EQUATION * EQUATION];
    WORD B[BROW * BCOL * VARIABLE], C[CROW * CCOL * VARIABLE],
        A[AROW * ACOL * VARIABLE];
    WORD TEMPA[ACOL * VARIABLE], TEMPB[BROW * VARIABLE], TEMPC[CROW * VARIABLE],
        TEMPF[CENTRAL_MAP_SIZE];
    WORD *F = new WORD[(AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE];
    // WORD *FB = new WORD[(AROW * BCOL + AROW * CCOL) * CENTRAL_MAP_SIZE];
    int i, j, s, t, ij = 0, eof, count, flag = 0; //
    WORD W = 0;

    // printf("************************************************* \n");
    // printf("Start ABC key pair generation \n");
    // printf("************************************************* \n");

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

    //  store T^(-1) in secrete key sk //sk 前 equation*equation 个位置存储 T^{-1} 的数据
    for (i = 0; i < EQUATION; i++)
    {
        for (j = 0; j < EQUATION; j++)
        {
            sk[ij++] = INVT[i * EQUATION + j];
        }
    }

    // generate the central map // 为什么这里的A要生成？
    // generate A
    for (i = 0; i < AROW * ACOL * VARIABLE; i++)
    {
        A[i] = rand() % FIELD;
    }

    // generate B and store B to sk // sk 从equation*equation+1 开始放 B 的 数据
    for (i = 0; i < BROW * BCOL * VARIABLE; i++)
    {
        sk[ij++] = B[i] = rand() % FIELD;
    }

    // generate C and store C to sk // sk 从equation*equation+variable*variable 开始放 C 的数据
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
                    TEMPB[s * VARIABLE + t] =
                        B[j * VARIABLE + s * BCOL * VARIABLE + t];
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
                TEMPA[s * VARIABLE + t] =
                    A[i * ACOL * VARIABLE + s * VARIABLE + t];
            }
        }
        for (j = 0; j < CCOL; j++)
        {
            for (s = 0; s < CROW; s++)
            {
                for (t = 0; t < VARIABLE; t++)
                {
                    TEMPC[s * VARIABLE + t] =
                        C[j * VARIABLE + s * CCOL * VARIABLE + t];
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

    // compute ToF
    matrixmul(T, F, EQUATION, EQUATION, CENTRAL_MAP_SIZE, pk);
    delete F;
    return (1);
}

int encryption(WORD *pk, WORD *plaintext, WORD *ciphertext, int plaintextlen)
{
    int i, j, k, ij = 0, eof = 0;
    WORD temp, x[VARIABLE];
    if (plaintextlen != VARIABLE)
        return (0);

    // printf("************************************************* \n");
    // printf("Start ABC encryption \n");
    // printf("************************************************* \n");

    for (k = 0; k < VARIABLE; k++)
        x[k] = WORD(plaintext[k]);

    for (k = 0; k < EQUATION; k++)
    {
        temp = 0;
        for (i = 0; i < VARIABLE; i++)
        {
            for (j = i; j < VARIABLE; j++)
            {
                temp = temp ^ mul(WORD(pk[ij++]), mul(x[i], x[j]));
            }
        }
        ciphertext[k] = temp;
    }
    return (1);
}

int decryption(WORD *sk, WORD *decrypttext, WORD *ciphertext,
               int ciphertextlen)
{
    int i, j, k, ij = 0, eof = 0;
    WORD x[AROW * ACOL + VARIABLE], y[EQUATION], y1[EQUATION], y2[EQUATION],
        X[VARIABLE];                                             // y 用来记录ciphtertext; y1=T^{-1}(y); X=F^{-1}(y1) 问题是y2 是什么？
    WORD coeff[ACOL * (BCOL + CCOL) * (ACOL * AROW + VARIABLE)]; // 存储E1,E2,B,C
    WORD TP, xn;                                                 // 这两个是什么？
    if (ciphertextlen != EQUATION)
        return (0);

    // printf("************************************************* \n");
    // printf("Start ABC decryption \n");
    // printf("************************************************* \n");

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

    // next find F^{-1}// 为什么看上去这里的 coeff 要记录更多的data。 这里的 coeff 是什么？？？
    for (k = 0; k < ACOL; k++)
    {
        // 计算E1,coeff 在这里放了 s*s 个数
        for (i = 0; i < BCOL; i++)
        {
            for (j = 0; j < AROW; j++)
            {
                coeff[(k * BCOL + i) * (ACOL * AROW + VARIABLE) + k * AROW +
                      j] = y1[j * BCOL + i];
            }
        }
        // 计算E2，coeff 在这里放了 s*s 个数
        for (i = 0; i < CCOL; i++)
        {
            for (j = 0; j < AROW; j++)
            {
                coeff[ACOL * BCOL * (ACOL * AROW + VARIABLE) +
                      (k * CCOL + i) * (ACOL * AROW + VARIABLE) + k * AROW +
                      j] = y1[AROW * BCOL + j * CCOL + i];
            }
        }
    } // 总共放了 s*s*s 个数

    for (i = 0; i < ACOL * BCOL; i++)
    {
        for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++)
        {
            coeff[i * (AROW * ACOL + VARIABLE) + j] = WORD(sk[ij++]);
        } // 这里的 sk 开始放 B 的 data
    }

    for (i = 0; i < ACOL * CCOL; i++)
    {
        for (j = AROW * ACOL; j < (ACOL * AROW + VARIABLE); j++)
        {
            coeff[ACOL * BCOL * (AROW * ACOL + VARIABLE) +
                  i * (AROW * ACOL + VARIABLE) + j] = WORD(sk[ij++]);
        } // 这里的sk 开始放 C 的data
    }

    echelonform(coeff, ACOL * (BCOL + CCOL), AROW * ACOL + VARIABLE);
    // coeff 现在是 row=2s^2=m, col= 2s^2=m
    /*
            printf( "\n coeff: \n");
            for (i=0; i<ACOL*(BCOL+CCOL)*(ACOL*AROW+VARIABLE); i++){
                    printf( "%d ", coeff[i]) ;
                    if ((i+1)%((ACOL*AROW+VARIABLE))==0){
       printf("\n");printf("\n");}
            }
    */

    WORD WM[VARIABLE];
    for (i = 0; i < VARIABLE - 1; i++)
    {
        WM[i] = coeff[(i + AROW * ACOL) * (ACOL * AROW + VARIABLE) +
                      (ACOL * AROW + VARIABLE - 1)];
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
    TP = mul(inv(y2[i]), y1[i]);
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
    for (i = 0; i < VARIABLE; i++)
    {
        X[i] = mul(WM[i], X[VARIABLE - 1]);
    }

    for (i = 0; i < VARIABLE; i++)
    {
        decrypttext[i] = X[i];
    }
    return (1);
}

int encryption_new(WORD *pk, WORD *plaintext, WORD *ciphertext, int plaintextlen)
{
    WORD ciphertext1[EQUATION], ciphertext2[EQUATION], ciphertext3[EQUATION];
    WORD plaintext1[VARIABLE], plaintext2[VARIABLE], plaintext3[VARIABLE];
    int k;
    for (k = 0; k < VARIABLE; k++)
        plaintext1[k] = WORD(plaintext[k] ^ R1[k]);
    plaintext2[k] = WORD(plaintext[k] ^ R2[k]);
    plaintext3[k] = WORD(plaintext[k] ^ R3[k]);

    encryption(pk, plaintext1, ciphertext1, plaintextlen);
    encryption(pk, plaintext2, ciphertext2, plaintextlen);
    encryption(pk, plaintext3, ciphertext3, plaintextlen);

    for (int i = 0; i < EQUATION; i++)
    {
        ciphertext[i] = ciphertext1[i];
        ciphertext[EQUATION + i] = ciphertext2[i];
        ciphertext[2 * EQUATION + i] = ciphertext3[i];
    }
    return (1);
}

bool is_equal(WORD *a, WORD *b)
{
    for (int i = 0; i < VARIABLE; i++)
    {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

int decryption_new(WORD *sk, WORD *decrypttext, WORD *ciphertext, int ciphertextlen)
{
    WORD ciphertext1[EQUATION], ciphertext2[EQUATION], ciphertext3[EQUATION];
    WORD decrypttext1[VARIABLE], decrypttext2[VARIABLE], decrypttext3[VARIABLE];
    for (int i = 0; i < EQUATION; i++)
    {
        ciphertext1[i] = ciphertext[i];
        ciphertext2[i] = ciphertext[EQUATION + i];
        ciphertext3[i] = ciphertext[2 * EQUATION + i];
    }

    decryption(sk, decrypttext1, ciphertext1, EQUATION);
    decryption(sk, decrypttext2, ciphertext2, EQUATION);
    decryption(sk, decrypttext3, ciphertext3, EQUATION);

    for (int i = 0; i < VARIABLE; i++)
    {
        decrypttext1[i] = decrypttext1[i];
        decrypttext2[i] = decrypttext2[i] ^ R2[i];
        decrypttext3[i] = decrypttext3[i] ^ R3[i];
    }
    if (is_equal(decrypttext2, decrypttext3))
    {
        for (int i = 0; i < VARIABLE; i++)
        {
            decrypttext[i] = decrypttext2[i];
        }
    }
    else
    {
        for (int i = 0; i < VARIABLE; i++)
        {
            decrypttext[i] = decrypttext1[i];
        }
    }
    return (1);
}

int main(int argc, int argv[])
{
    // Define test time
    unsigned int times = 1000;
    std::cout << "Test times: " << times << std::endl;

    int plaintextlen = VARIABLE;
    int i;
    WORD *sk = new WORD[SECRET_KEY_SIZE];
    WORD *pk = new WORD[PUBLIC_KEY_SIZE];
    WORD plaintext[VARIABLE];
    // WORD ciphertext[EQUATION];
    WORD ciphertext[3 * EQUATION];
    // WORD decrypttext[EQUATION];
    WORD decrypttext[VARIABLE];
    srand((int)time(0));
    // keypair(sk, pk);

    for (int t = 0; t < times; t++)
    {
        keypair(sk, pk);
        for (int j = 0; j < VARIABLE; j++)
        {
            plaintext[j] = rand() % FIELD;
        }

        encryption_new(pk, plaintext, ciphertext, plaintextlen);
        decryption_new(sk, decrypttext, ciphertext, EQUATION);

        for (int k = 0; k < VARIABLE; k++)
        {
            if (decrypttext[k] != plaintext[k])
            {
                // std::cout << decrypttext[k] << "  " << plaintext[k]
                //           << std::endl;
                std::cout << "Decrypt Error!" << std::endl;
                break;
            }
        }
    }

    delete[] sk;
    delete[] pk;
    return (1);
}