#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define PI 3.14159265359


double **CGUPDATE(double **U, double **R, double **P, double **AP, int n, int m) {
    int i, j;
    double nalpha, dalpha, alpha, nbeta, beta;
    double **Rk = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        Rk[i] = (double *) malloc(n * sizeof(double));

    double **Pk = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        Pk[i] = (double *) malloc(n * sizeof(double));

    nalpha = 0;
    dalpha = 0;
    nbeta = 0;

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            nalpha = nalpha + R[i][j] * R[i][j];

        }
    }
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            dalpha = dalpha + P[i][j] * AP[i][j];

        }
    }
    alpha = nalpha / dalpha;

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            U[i][j] = U[i][j] + (alpha * P[i][j]);

        }
    }
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            Rk[i][j] = R[i][j] - (alpha * AP[i][j]);

        }
    }

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            nbeta = nbeta + (Rk[i][j] * Rk[i][j]);

        }
    }

    beta = nbeta / nalpha;

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            Pk[i][j] = Rk[i][j] + (beta * P[i][j]);

        }
    }
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= m; j++) {
            P[i][j] = Pk[i][j];
            R[i][j] = Rk[i][j];

        }
    }
    return U, R, P;

}


double Fxy1(double w) {
    double value;

    value = exp(w * PI);
    return value;
}

double Fxy2(double w) {
    double value;

    value = -exp(w * PI);
    return value;
}

double Gxy1(double w) {
    double value;

    value = cos(w * PI);
    return value;
}

double Gxy2(double w) {
    double value;

    value = -cos(w * PI);
    return value;
}


double BDYVAL(int option, double w) {
    int i, j;
    double value;

    if (option == 1) {
        value = Fxy2(w);
    } else if (option == 2) {
        value = Fxy1(w);
    } else if (option == 3) {
        value = Gxy1(w);
    } else if (option == 4) {
        value = Gxy2(w);
    }
    return value;
}


double ERROR_METRIC(double **input_array, int N, int option) {
    int i;

    int n = sizeof(input_array) / sizeof(double);
    int m = N / n;

    double *v = (double *) malloc(N * sizeof(double));
    for (i = 0; i < n; i++) {
        for (unsigned j = 0; j < m; j++) {
            v[i + n * j] = input_array[i][j];
        }
    }

// Variables used to store the sum of absolute values, sum of squared values
// the maximum absolute value, the current absolute value
    double xSumAbsValues, xSumSquared, xMaxAbsValue, currentAbsValue;

// Variables used to store the results
    double L1norm, L2norm, Linfnorm;

// Initialize variables
    xSumAbsValues = 0.0;
    xSumSquared = 0.0;
    xMaxAbsValue = -1.0;

// Loop used to calculate xSumAbsValues, xSumSquared and find xMaxAbsValue
    for (i = 0; i < N; i = i + 1) {

// Calculate the absolute value of x[i]
        currentAbsValue = (v[i] < 0) ? (-v[i]) : (v[i]);

// Check if the current absolute value is the maximum
        if (currentAbsValue > xMaxAbsValue)
            xMaxAbsValue = currentAbsValue;

// Update the sum of absolute values
        xSumAbsValues = xSumAbsValues + currentAbsValue;

// Update the sum of squared values
        xSumSquared = xSumSquared + v[i] * v[i];

    }

    if (option == 1) {
// Find the L1
        L1norm = xSumAbsValues;
        printf("The L1-norm of x is %15.4f \n", L1norm);
        return L1norm;
    } else if (option == 2) {
// Find the L2 norm
        L2norm = sqrt(xSumSquared);
        printf("The L2-norm of x is %15.4f \n", L2norm);
        return L2norm;
    } else if (option == 3) {
        Linfnorm = xMaxAbsValue;
        printf("The Linf-norm of x is %15.4f \n", Linfnorm);
        return Linfnorm;
    }

}

double **LAPLACEWCG() {
    int i, j, m, n, cnt;
    double err, rx, ry, ave, a, b, h, hx, hy, tol, max1;
    


    hx = 0.1;
    hy = 0.1;
    h = 0.1;
    printf("Enter the number of rows: ");
    scanf("%lf", &a);

    printf("Enter the number of columns: ");
    scanf("%lf", &b);

    n = (a / hy) + 1;
    m = (b / hx) + 1;
	
	double *X = (double *) malloc(m * sizeof(double));
    double *Y = (double *) malloc(n * sizeof(double));

    double **R = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        R[i] = (double *) malloc(n * sizeof(double));

    double **P = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        P[i] = (double *) malloc(n * sizeof(double));

    double **AP = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        AP[i] = (double *) malloc(n * sizeof(double));

    double **U = (double **) malloc(m * sizeof(double));
    for (i = 0; i < n; i++)
        U[i] = (double *) malloc(n * sizeof(double));

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            U[i][j] = 1;
        }
    }

    for (j = 0; j < m; j++) {
        X[j] = j * hx;
    }
    X[m] = a;

    for (j = 0; j < n; j++) {
        Y[j] = (b - (j * hy));
    }
    Y[n] = 0;

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            R[i][j] = 0.0;
            P[i][j] = 0.0;
            AP[i][j] = 0.0;
        }
    }

    rx = (1 / (hx * hx));
    ry = (1 / (hy * hy));

    ave = (a * (BDYVAL(1, 0) + BDYVAL(2, 0)) + b * (BDYVAL(3, 0) + BDYVAL(4, 0))) / (2 * a + 2 * b);

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            U[i][j] = ave * U[i][j];
        }
    }

    for (i = 0; i <= n; i++) {
        U[i][0] = BDYVAL(3, Y[i]);
        U[i][m] = BDYVAL(4, Y[i]);
    }

    for (j = 0; j <= m; j++) {
        U[0][j] = BDYVAL(1, X[j]);
        U[n][j] = BDYVAL(2, X[j]);
    }

    U[0][0] = (U[0][1] + U[1][0]) / 2;
    U[0][m] = (U[0][m - 1] + U[1][m]) / 2;
    U[n][0] = (U[n - 1][0] + U[n][1]) / 2;
    U[n][m] = (U[n - 1][m] + U[n][m - 1]) / 2;

    for (j = 1; j < m; j++) {
        for (i = 1; i < n; i++) {
            R[i][j] = (rx * U[i][j + 1] + rx * U[i][j - 1] + ry * U[i + 1][j] + ry * U[i - 1][j]
                       - 2 * (rx + ry) * U[i][j]);
        }
    }

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            P[i][j] = R[i][j];
        }
    }

    err = ERROR_METRIC(R, m * n, 1);

    while ((err > tol) && (cnt <= max1)) {
        for (j = 2; j < m; j++) {
            for (i = 2; i < n; i++) {
                if (j == 2) {
                    if (i == 2) {
                        AP[i][j] = -rx * P[i][j + 1] - ry * P[i + 1][j] + 2 * (rx + ry) * P[i][j];

                    } else if (i == n - 1) {
                        AP[i][j] = -rx * P[i][j + 1] - ry * P[i - 1][j] + 2 * (rx + ry) * P[i][j];
                    } else {
                        AP[i][j] = -rx * P[i][j + 1] - ry * P[i + 1][j] - ry * P[i - 1][j] + 2 * (rx + ry) * P[i][j];
                    }

                } else if (j == m - 1) {
                    if (i == 2) {
                        AP[i][j] = -rx * P[i][j - 1] - ry * P[i + 1][j] + 2 * (rx + ry) * P[i][j];
                    } else if (i == n - 1) {
                        AP[i][j] = -rx * P[i][j - 1] - ry * P[i - 1][j] + 2 * (rx + ry) * P[i][j];
                    } else {
                        AP[i][j] = -rx * P[i][j - 1] - ry * P[i + 1][j] - ry * P[i - 1][j] + 2 * (rx + ry) * P[i][j];
                    }
                } else if (i == n - 1) {
                    AP[i][j] = -rx * P[i][j + 1] - ry * P[i][j - 1] - ry * P[i - 1][j] + 2 * (rx + ry) * P[i][j];
                } else if (i == 2) {
                    AP[i][j] = -rx * P[i][j + 1] - ry * P[i][j - 1] - ry * P[i + 1][j] + 2 * (rx + ry) * P[i][j];
                } else {
                    AP[i][j] = -rx * P[i][j + 1] - rx * P[i][j - 1] - ry * P[i + 1][j] - ry * P[i - 1][j]
                               + 2 * (rx + ry) * P[i][j];
                }
            }
        }
        R = CGUPDATE(U, R, P, AP, n, m);
        P = CGUPDATE(U, R, P, AP, n, m);
        U = CGUPDATE(U, R, P, AP, n, m);
        err = ERROR_METRIC(R, m * n, 1);
        cnt = cnt + 1;
    }

    if (cnt >= max1) {
        printf("Maximum number of iterations exceeded");
    }
	 for (i = 0; i <= n; i++) {
        for (j = 0; j <= m; j++) {
            printf("U: \n %lf", U[i][j]);
        }
    }
    return U;
}

int main() {
	
    LAPLACEWCG();
	
    return 0;
}




