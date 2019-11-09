#include<iostream>
#include<fstream>
#include"test_function.h"

using namespace std;

//定义x的n次方
double power(double x, int n) {
	double res = 1;
	for (int i = 0; i < n; ++i) {
		res = res * x;
	}
	return res;
}



#define N 2  //对于第二题 N=2,可修改得到更高次的高斯积分

//消元法求解线性方程组,并将结果保留在b[]中
void Gaussian_elimination(double* A[], double* b,int n) {
	//对角化
	for (int i = 0; i < n; ++i) {
		if (A[i][i] == 0) { cout << "error"; break; }
		for (int j = i + 1; j < n; ++j) {
			double c = -1*A[j][i] / A[i][i];
			for (int k = i + 1; k < n; ++k) {
				A[j][k] = A[j][k] + c * A[i][k];
			}
			b[j] = b[j] + b[i]*c;
		}
	}
	//回代求解x
	for (int i = (n - 1); i >= 0; --i) {
		for (int j = i + 1; j < n; ++j) {
			b[i] = b[i] - A[i][j] * b[j];
		}
		b[i] = b[i] / A[i][i];
	}
}

//定义多项式类,其中最高次系数为1
class polynomial {
public:
	//定义n次多项式
	polynomial(int m) {
		coefficient = new double[m + 1];
		coefficient[m] = 1;
		n = m;
	}
	~polynomial() { delete[] coefficient; }
	//赋值
	void setvalue(double* a) {
		for (int i = 0; i < n; ++i) {
			coefficient[i] = a[i];
		}
	}
	//求值f(x)
	double getvalue(double x) {
		double res = 1;
		for (int i = n - 1; i >= 0; --i) {
			res = res * x + coefficient[i];
		}
		return res;
	}
	//在区间[a,b]用二分法求零点
	double zero(double a, double b) {
		double fa = getvalue(a);
		while ((b - a) / 2 > (1e-10)) {
			double c = a + (b - a) / 2;
			double temp = getvalue(c);
			if (temp * fa < 0) {
				b = c;
			}
			else {
				a = c;
				fa = temp;
			}
		}
		return (b + a) / 2;
	}

	//在[0,1]区间内计算所有零点
	double* zero_point() {
		double* res = new double[n];
		//找可能存在零点的区间
		double** section = new double*[n];//存储可能零点区间
		for (int i = 0; i < n; ++i) section[i] = new double[2];
		
		int point_num = 0;
		for (int i = 1; point_num < (n-1); i*=2) {
			point_num = 0;
			double step= 1.0/ (i*(n + 1));
			for (int j = 0; j < (i*(n + 1));++j) {
				if (getvalue(j*step)*getvalue((j + 1)*step) < 0) {
					section[point_num][0] = j * step;
					section[point_num][1] = (j+1) * step;
					point_num++;
				}
			}
		}
		for (int i = 0; i < n; ++i) {
			res[i] = zero(section[i][0], section[i][1]);
		}
		return res;
	}


private:
	double* coefficient;
	int n;
};

//测试函数部分
void Test(double f_x[], double f_A[], double s_x[], double s_A[], double t_x[], double t_A[], double(*f)(double)) {
	double res = f_A[0] * f(f_x[0]) + f_A[1] * f(f_x[1]);
	cout << "第一题结果：" << res << endl;
	res = 0;
	for (int i = 0; i <= N; ++i) {
		res += s_A[i] * f(s_x[i]);
	}
	cout << "第二题结果：" << res << endl;
	res = t_A[0] * f(t_x[0]) + t_A[1] * f(t_x[1]);
	cout << "第三题结果：" << res << endl;
}

//主函数
int main() {
//第一题部分：
	double f_x[2];
	double f_A[2];
	f_x[0] = (35 - 2 * sqrt(70)) / 63.0;
	f_x[1] = (35 + 2 * sqrt(70)) / 63.0;
	f_A[0] = (1 / 3.0) - (sqrt(70) / 150);
	f_A[1] = (1 / 3.0) + (sqrt(70) / 150);
	cout << "第一题：" << endl << "节点为：x0=" << f_x[0] << ", x1=" << f_x[1] << endl;
	cout << "系数为：A0=" << f_A[0] << ", A1=" << f_A[1] << endl;


//第二题部分：
	double* b = new double[N+1];
	double** A=new double*[N+1];
	//初始化方程
	for (int i = 0; i <= N; ++i) {
		A[i] = new double[N];
		for (int j = 0; j <= N; ++j) {
			A[i][j] = 1 / (i + j + 1.5);
		}
		b[i] = -1 / (N + i + 2.5);
	}
	//计算系数an
	Gaussian_elimination(A, b, N+1);
	//得到高斯节点，保存在zero数组中
	polynomial p{ N + 1 };
	p.setvalue(b);
	double* zero = p.zero_point();
    //计算节点系数An
	for (int i = 0; i <= N; ++i) {
		for (int j = 0; j <= N; ++j) {
			A[i][j] = power(zero[j]*10, i);
		}
		b[i] = (power(10, i)) / (1.5 + i);
	}
	Gaussian_elimination(A, b, N + 1);
	cout << "第二题：" << endl << "节点为：";
	for (int i = 0; i <= N; ++i) {
		cout << "x" << i << "=" << zero[i];
		if (i != N) cout << ", ";
	}
	cout << endl;
	cout << "系数为：";
	for (int i = 0; i <= N; ++i) {
		cout << "A" << i << "=" << b[i];
		if (i != N) cout << ", ";
	}
	cout << endl;
//第三题部分：
	double t_x[2];
	double t_A[2];
	t_x[0] = sqrt(10) / 5;
	t_x[1] = -sqrt(10) / 5;
	t_A[0] = 4 / 3.0;
	t_A[1] = 4 / 3.0;
	cout << "第三题：" << endl << "节点为：x0=" << t_x[0] << ", x1=" << t_x[1] << endl;
	cout << "系数为：A0=" << t_A[0] << ", A1=" << t_A[1] << endl;

//Test部分(自由发挥部分)
	cout << endl << "Test部分：" << endl;

	cout << "f(x)=x^5 积分：" << endl;
	Test(f_x, f_A, zero, b, t_x, t_A, power_5);
	cout << "f(x)=x^6 积分：" << endl;
	Test(f_x, f_A, zero, b, t_x, t_A, power_6);
	cout << "f(x)=Sin(x) 积分：" << endl;
	Test(f_x, f_A, zero, b, t_x, t_A, sin);
	cout << "f(x)=exp(x) 积分：" << endl;
	Test(f_x, f_A, zero, b, t_x, t_A, exp);
}

