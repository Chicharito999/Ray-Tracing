#ifndef POINT_H
#define POINT_H

using namespace std;

class point {
private:
	// 成员变量
	double u;
	double v;
	double w;

public:
	// 构造、析构函数 
	point(double x, double y, double z) { u = x; v = y; w = z; }
	~point() {};

	// 方法
	void sety(double y) { v = y; }
	double x() { return u; }
	double y() { return v; }
	double z() { return w; }
};

#endif 
