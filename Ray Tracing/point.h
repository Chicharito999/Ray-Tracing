#ifndef POINT_H
#define POINT_H

using namespace std;

class point {
private:
	// ��Ա����
	double u;
	double v;
	double w;

public:
	// ���졢�������� 
	point(double x, double y, double z) { u = x; v = y; w = z; }
	~point() {};

	// ����
	void sety(double y) { v = y; }
	double x() { return u; }
	double y() { return v; }
	double z() { return w; }
};

#endif 
