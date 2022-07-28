#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#define SIZE 4
#define SIZE2 2
struct node
{
	float a;
	float b;
	float c;
	float d;
	node* next;
};
float d(float x1, float x2, node* splain1, node* splain2)
{
	float ans = (splain1->a * x1 * x1 * x1 + splain1->b * x1 * x1 + splain1->c * x1 + splain1->d - (splain2->a * x2 * x2 * x2 + splain2->b * x2 * x2 + splain2->c * x2 + splain2->d));
	ans *= ans;
	ans += (x1 - x2) * (x1 - x2);
	ans = sqrt(ans);
	printf("\nx1=%.2f, x2=%.2f, d=%.2f\n", x1, x2, ans);
	return ans;
}
float y(float x, node* splain)
{
	return splain->a * x * x * x + splain->b * x * x + splain->c * x + splain->d;
}
bool spusk(node* splain1, node* splain2, float x1, float x2, float L)
{
	float X[SIZE2] = {};
	X[0] = (x2 - x1) / 2;
	X[1] = (x2 - x1) / 2;
	float eps = 0.0001;
	float x0 = X[0] + 1, x = X[1] + 1;
	int i = 0;
	float first = 2 * (X[0] - X[1]) + 2 * (y(X[0], splain1) - y(X[1], splain2)) * (3 * splain1->a * X[0] * X[0] + 2 * splain1->b * X[0] + splain1->c);
	float second = -2 * (X[0] - X[1]) - 2 * (y(X[0], splain1) - y(X[1], splain2)) * (3 * splain2->a * X[1] * X[1] + 2 * splain2->b * X[1] + splain2->c);
	while (abs(first) > eps && abs(second) > eps)
	{
		i += 1;
		x0 = X[0];
		x = X[1];
		X[0] -= first * L;
		X[1] -= second * L;
		if (X[0] == x0&& X[1] == x)
		{
			printf("stopped inf cycle");
			break;
		}
		first = 2 * (X[0] - X[1]) + 2 * (y(X[0], splain1) - y(X[1], splain2)) * (3 * splain1->a * X[0] * X[0] + 2 * splain1->b * X[0] + splain1->c);
		second = -2 * (X[0] - X[1]) - 2 * (y(X[0], splain1) - y(X[1], splain2)) * (3 * splain2->a * X[1] * X[1] + 2 * splain2->b * X[1] + splain2->c);
	}
	d(X[0], X[1], splain1, splain2);
	return 1;
}
float gam(int n, float** f, float* gam)
{
	gam[0] = 0, gam[n - 1] = 0;
	float** k = (float**)malloc(sizeof(float*) * (n-2));
	for (int i = 0; i < n-2; i++)
	{
		k[i] = (float*)malloc(sizeof(float) * (n));
	}
	for (int i = 0; i < n - 2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j >= i && j <= i + 2)
			{
				if ((j == i || j == i + 2)&&j!=0)
				{
					k[i][j] = (f[j][0]-f[j-1][0]) / 6;
				}
				if (j == i + 1 && j!=n-1)
				{
					k[i][j] = (f[j + 1][0] - f[j - 1][0]) / 3;
				}
			}
			else
			{
				k[i][j] = 0;
			}
		}
	}
	k[0][0] = 0, k[n-3][n-1] = 0;
	float* b = (float*)malloc(sizeof(float) * (n - 2));
	for (int i = 0; i < n - 2; i++)
	{
		b[i] = (f[i + 2][1] - f[i + 1][1]) / (f[i + 2][0] - f[i + 1][0]) - (f[i + 1][1] - f[i][1]) / (f[i + 1][0] - f[i][0]);
	}
	for (int i = 1; i < n - 2; i++)
	{
		k[i][i + 1] +=k[i][i]*(-k[i-1][i+1]/k[i-1][i]);
		b[i] -=b[i-1]/k[i-1][i];
		k[i][i] = 0;
	}
	for (int i = n - 4; i >= 0; i--)
	{
		b[i] -= b[i + 1] * (k[i][i+2]/k[i+1][i+2]);
		k[i][i + 2] = 0;
	}
	for (int i = 1; i < n - 1; i++)
	{
		gam[i] = b[i-1] / k[i-1][i];
	}
	return 0;
}
void splain(float** f, int i, float* gamma,node* splain)
{
	float h = f[i + 1][0] - f[i][0];
	splain->a=(-gamma[i] + gamma[i + 1]) / (6 * h);
	splain->b = (gamma[i] * f[i + 1][0] - gamma[i + 1] * f[i][0]) / (2 * h);
	splain->c = (gamma[i] * (h * h - 3 * f[i + 1][0] * f[i + 1][0]) + gamma[i + 1] * (3 * f[i][0] * f[i][0] - h * h)) / (6 * h) + (f[i + 1][1] - f[i][1]) / h;
	splain->d=(gamma[i]*(f[i+1][0]* f[i + 1][0]* f[i + 1][0]- f[i + 1][0]*h*h)+gamma[i+1]*(h*h* f[i][0]- f[i][0]* f[i][0]* f[i][0]))/(6*h)+(f[i][1]*f[i+1][0]-f[i+1][1]*f[i][0])/h;
	printf("(%.2f)x^3+(%.2f)x^2+(%.2f)x+(%.2f)", splain->a, splain->b,splain->c, splain->d);
}
float min(float** f1, float** f2, node* splain1,node* splain2,float minr, int n)
{
	int i = 0, j = 0;
	float x1 = f1[i][0], x2=f2[j][0];
	float eps = 0.01;
	while (x1 < f1[i + 1][0])
	{
		while (x2 < f1[i + 1][0] && x2<f2[j+1][0])
		{
			if (abs(2 * (x1 - x2 + (splain1->a * x1 * x1 * x1 + splain1->b * x1 * x1 + splain1->c * x1 + splain1->d-(splain2->a * x2 * x2 * x2 + splain2->b *x2 * x2 + splain2->c * x2 + splain2->d )*(3*splain1->a * x1 * x1 + 2*splain1->b * x1  + splain1->c))))<eps && abs(2 * (x1 - x2 + (splain1->a * x1 * x1 * x1 + splain1->b * x1 * x1 + splain1->c * x1 + splain1->d - (splain2->a * x2 * x2 * x2 + splain2->b * x2 * x2 + splain2->c * x2 + splain2->d) * (3 * splain2->a * x2 * x2 + 2 * splain2->b * x2 + splain2->c))))<eps)
			{
				if (d(x1, x2, splain1, splain2) < eps)
				{
					printf("\ncross: x=%.2f, y=%.2f", x1, splain1->a * x1 * x1 * x1 + splain1->b * x1 * x1 + splain1->c * x1 + splain1->d);
					minr = 0;
					return minr;
				}
				if (d(x1, x2, splain1, splain2) < minr)
				{
					minr = d(x1, x2, splain1, splain2);
				}
			}
			x2 += 0.0001;
		}
		x1 += 0.0001;
		printf("%.4f", x1);
if (abs(x1 - f1[i + 1][0]) < 0.0001)
{
	i += 1;
	if (i == n - 1)
	{
		i -= 1;
	}
	if (i == n)
	{
		return minr;
	}
	splain1 = splain1->next;
	printf("\nipassed%d\n", i);
}
if (abs(x2 - f2[j + 1][0]) < 0.0001 && abs(x1 - f2[i + 1][0]) < 0.0001)
{
	j += 1;
	if (j == n - 1)
	{
		j -= 1;
	}
	if (j == n)
	{
		return minr;
	}
	splain2 = splain2->next;
	printf("\njpassed%d\n", j);
}
x2 = f2[j][0];
	}
	return minr;
}
float F(node* splain1, node* splain2, float x)
{
	return (splain1->a - splain2->a) * x * x * x + (splain1->b - splain2->b) * x * x + (splain1->c - splain2->c) * x + splain1->d - splain2->d;
}
bool dih(float x1, float x2, node* splain1, node* splain2)
{
	float eps = 0.0001;
	float x = (x2 + x1) / 2;
	if (F(splain1, splain2, x1) < 0 && F(splain1, splain2, x2) > 0)
	{
		while (abs(F(splain1, splain2, x)) > eps)
		{
			if (F(splain1, splain2, x) > 0)
			{
				x2 = x;
			}
			else
			{
				x1 = x;
			}
			x = (x2 + x1) / 2;
		}
		printf("\n(%.2f)x^3+(%.2f)x^2+(%.2f)x+(%.2f)\n", splain1->a, splain1->b, splain1->c, splain1->d);
		printf("(%.2f)x^3+(%.2f)x^2+(%.2f)x+(%.2f)\n", splain2->a, splain2->b, splain2->c, splain2->d);
		printf("cross: x=%.2f, y=%.2f\n", x, splain1->a * x * x * x + splain1->b * x * x + splain1->c * x + splain1->d);
		return 1;
	}
	else if (F(splain1, splain2, x1) > 0 && F(splain1, splain2, x2) < 0)
	{
		while (abs(F(splain1, splain2, x)) > eps)
		{
			if (F(splain1, splain2, x) > 0)
			{
				x1 = x;
			}
			else
			{
				x2 = x;
			}
			x = (x2 + x1) / 2;
		}
		printf("\n(%.2f)x^3+(%.2f)x^2+(%.2f)x+(%.2f\n)", splain1->a, splain1->b, splain1->c, splain1->d);
		printf("(%.2f)x^3+(%.2f)x^2+(%.2f)x+(%.2f)\n", splain2->a, splain2->b, splain2->c, splain2->d);
		printf("cross: x=%.2f, y=%.2f\n", x, splain1->a * x * x * x + splain1->b * x * x + splain1->c * x + splain1->d);
		return 1;
	}
	else
	{
		return 0;
	}
}
bool newton(node* splain1, node* splain2, float x1,float x2)
{
	float H[SIZE2][SIZE2] = {};
	float detH;
	float time;

	float X[SIZE2] = {};
	X[0] = (x1+x2)/2;
	X[1] = (x1+x2)/2;
	float eps = 0.01;
	float x0, x;
	float first = 2 * (X[0] - X[1]) + 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (3 * splain1->a * X[0] * X[0] + 2 * splain1->b * X[0] + splain1->c);
	float second = -2 * (X[0] - X[1]) - 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (3 * splain2->a * X[1] * X[1] + 2 * splain2->b * X[1] + splain2->c);
	while (second> eps || first>eps)
	{
		H[0][0] = 2 + 2 * (splain1->c + 2 * splain1->b * X[0] + 3 * splain1->a * X[0] * X[0]) * (splain1->c + 2 * splain1->b * X[0] + 3 * splain1->a * X[0] * X[0]) + 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (6 * splain1->a * X[0] + 2 * splain1->b);
		H[1][1] = 2 +2 * (splain2->c + 2 * splain2->b * X[1] + 3 * splain2->a * X[1] * X[1]) * (splain2->c + 2 * splain2->b * X[1] + 3 * splain2->a * X[1] * X[1]) - 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (6 * splain2->a * X[1] + 2 * splain2->b);
		H[0][1] = -2 - 2 * (3 * splain2->a * X[1] * X[1] + 2 * splain2->b * X[1] + splain2->c) * (3 * splain1->a * X[0] * X[0] + 2 * splain1->b * X[0] + splain1->c);
		H[1][0] = H[0][1];
		detH = H[0][0] * H[1][1] - H[1][0] * H[0][1];
		time = H[0][1];
		H[0][1] = -H[1][0] / detH;
		H[1][0] = -time / detH;
		time = H[1][1];
		H[1][1] =H[0][0]/ detH;
		H[0][0] = H[1][1]/detH;
		x0 = X[0];
		x = X[1];
		X[0] -=(H[0][0] * first + H[0][1] * second);
		X[1] -= (H[1][0] * first + H[1][1]*second);
		if (x0 == X[0] && x == X[1])
		{
			printf("\nStopped inf cycle: %.5f, %.5f", X[0], X[1]);
			break;
		}
		first = 2 * (X[0] - X[1]) + 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (3 * splain1->a * X[0] * X[0] + 2 * splain1->b * X[0] + splain1->c);
		second = -2 * (X[0] - X[1]) - 2 * (splain1->a * X[0] * X[0] * X[0] + splain1->b * X[0] * X[0] + splain1->c * X[0] + splain1->d - (splain2->a * X[1] * X[1] * X[1] + splain2->b * X[1] * X[1] + splain2->c * X[1] + splain2->d)) * (3 * splain2->a * X[1] * X[1] + 2 * splain2->b * X[1] + splain2->c);
	}
	//if ((X[0] >= x2 || X[0] <= x1) || (X[1] >= x2 || X[1] <= x1))
	//{
	//	return 0;
	//}
	if (abs(d(X[0], X[1], splain1, splain2)) < eps)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
bool dihper(node* splain1, node* splain2, float** f1, float** f2, int n)
{
	int j = 0;
	float x1, x2;
	for (int i = 0; i < n-1; i++)
	{
		while (j<n-1)
		{
			if (f1[i][0] == f2[j][0])
			{
				x1 = f1[i][0];
				if (f1[i + 1][0] == f2[j + 1][0])
				{
					x2 = f1[i + 1][0];
				}
				if (f1[i + 1][0] > f2[j + 1][0])
				{
					x2 = f2[j + 1][0];
				}
				else
				{
					x2 = f1[i + 1][0];
				}
				j += 1;
				if(dih(x1,x2,splain1,splain2) == 0)
				{
					printf("spusk with: 0.001, 0.0005");
					spusk(splain1, splain2, 0, 1, 0.001);
					spusk(splain1, splain2, 0, 1, 0.0005);
					printf("newton:");
					newton(splain1, splain2, x1, x2);
				}
				else
				{
					printf("We have a cross");
					return 0;
				}
				splain2 = splain2->next;
				break;
			}
			if (f1[i][0] > f2[j][0])
			{
				x1 = f2[j][0];
				x2 = f1[i][0];
				if(dih(x1, x2, splain1, splain2) == 0)
				{
					printf("spusk with: 0.001, 0.0005");
					spusk(splain1, splain2, 0, 1, 0.001);
					spusk(splain1, splain2, 0, 1, 0.0005);
					printf("newton:");
					newton(splain1, splain2, x1, x2);
				}
				else
				{
					printf("We have a cross");
					return 0;
				}
				j +=1;
				splain2 = splain2->next;
			}
			else
			{
				x1 = f1[i][0];
				x2 = f2[j][0];
				if(dih(x1, x2, splain1, splain2) == 0)
				{
					printf("spusk with: 0.001, 0.0005");
					spusk(splain1, splain2, 0, 1, 0.001);
					spusk(splain1, splain2, 0, 1, 0.0005);
					printf("newton:");
					newton(splain1, splain2, x1, x2);
				}
				else
				{
					printf("We have a cross");
					return 0;
				}
				break;
			}
		}
		splain1 = splain1->next;
	}
	return 1;
}
int main()
{
	int  n = 4; float xx;
	float f3[SIZE][SIZE2] = { {1,3}, { 5, 4},{6,-3},{8,5} };
	float f4[SIZE][SIZE2] = { {1,7},{5, -5},{6,-5},{8,6} };
	//float tim;
	float** f1 = (float**)malloc(sizeof(float*) * n);
	for (int i = 0; i < n; i++)
	{
		f1[i] = (float*)malloc(sizeof(float) * 2);
	}
	srand((unsigned)time(0));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (j == 0)
			{
				/*scanf_s("%f", &xx);*/
				f1[i][j] = (float)i;
			}
			else
			{
				/*scanf_s("%f", &xx)*/;
				f1[i][j] = (float)(rand() % 21 - 10);
			}
			printf("%8.2f", f1[i][j]);
		}
		printf("\n");
	}
	printf("\n\n\n");
	float** f2 = (float**)malloc(sizeof(float*) * n);
	for (int i = 0; i < n; i++)
	{
		f2[i] = (float*)malloc(sizeof(float) * 2);
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			if (j == 0)
			{
				/*scanf_s("%f", &xx);*/
				f2[i][j] = (float)i;
			}
			else
			{
				/*scanf_s("%f", &xx);*/
				f2[i][j] = (float)(rand() % 21 - 10);
			}
			printf("%8.2f", f2[i][j]);
		}
		printf("\n");
	}
	float x;
	float* gamma1 = (float*)malloc(sizeof(float) * n);
	float* gamma2 = (float*)malloc(sizeof(float) * n);
	gam(n, f1, gamma1);
	gam(n, f2, gamma2);
	node* splain1= (node*)malloc(sizeof(node));
	node* splain2=(node*)malloc(sizeof(node));
	node* splaintime = splain1;
	for (int j = 0; j < n - 1; j++)
	{
		splain(f1, j, gamma1, splaintime);
		splaintime->next = (node*)malloc(sizeof(node));
		splaintime = splaintime->next;
		splaintime->next = nullptr;
		printf("\n");
	}
	splaintime = nullptr;
	printf("\n\n\n");
	splaintime = splain2;
	for (int j = 0; j < n - 1; j++)
	{
		splain(f2, j, gamma2, splaintime);
		splaintime->next = (node*)malloc(sizeof(node));
		splaintime = splaintime->next;
		splaintime->next = nullptr;
		printf("\n");
	}
	splaintime = nullptr;
	float minr = 100000;
	///dihper(splain1,splain2,f1,f2,n);
	if (dih(0, 1, splain1, splain2) ==0)
	{
		spusk(splain1, splain2,0, 1, 0.001);
		spusk(splain1, splain2,0, 1, 0.0005);
		spusk(splain1, splain2,0, 1, 0.0001);
		printf("newton");
		newton(splain1, splain2, 0, 1);
	}
}
