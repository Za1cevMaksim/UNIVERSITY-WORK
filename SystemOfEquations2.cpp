#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#define n 5
#define dt pow(10, -2)
#define eps pow(10, -10)

double* F(double X[n])
{
	double f[n] = { 0 };
	//f[0] = (pow(X[0], 3) + X[1] - 1) * dt;
	//f[1] = (pow(X[1], 3) - X[0] + 1) * dt;
	f[0] = (X[0] + X[2] * cos(3 * M_PI / 2 - X[3]) + 0.25) * dt;
	f[1] = (X[2] + X[2] * sin(3 * M_PI / 2 - X[3]) - 0.3) * dt;
	f[2] = (X[1] + X[2] * cos(3 * M_PI / 2 + X[4]) - 0.25) * (-dt);
	f[3] = (X[2] + X[2] * sin(3 * M_PI / 2 + X[4]) - 0.3) * dt;
	f[4] = ((X[3] + X[4]) * X[2] + (X[1] - X[0]) - M_PI / 4) * dt;
	return f;
}
bool Condition(double X[n])
{
	double f[n] = { 0 };
	//f[0] = pow(X[0], 3) + X[1] - 1;
	//f[1] = pow(X[1], 3) - X[0] + 1;
	f[0] = X[0] + X[2] * cos(3 * M_PI / 2 - X[3]) + 0.25;
	f[1] = X[2] + X[2] * sin(3 * M_PI / 2 - X[3]) - 0.3;
	f[2] = X[1] + X[2] * cos(3 * M_PI / 2 + X[4]) - 0.25;
	f[3] = X[2] + X[2] * sin(3 * M_PI / 2 + X[4]) - 0.3;
	f[4] = (X[3] + X[4]) * X[2] + (X[1] - X[0]) - M_PI / 4;
	for (int i = 0; i < n; ++i)
	{
		if (fabs(f[i]) > eps)
			return true;
	}
	return false;
}

int main(int argc, char* argv[])
{
	double X[n] = { 0 };
	for (int i = 0; i < n; ++i)
		X[i] = 0.5;

	while (Condition(X))
	{
		for (int i = 0; i < n; ++i)
			X[i] -= F(X)[i];
	}

	for (int i = 0; i < n; ++i)
		printf("X[%d] = %f\n", i, X[i]);
}
