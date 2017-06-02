#ifndef SPHERE_H
#define SPHERE_H


#include "point.h"

using namespace std;


class sphere {
private:
	// ��Ա����
	point *cent;
	double rad;
	double Ka[3];       // ���ʷ��价����
	double Kd[3];       // ���ʷ����������
	double Ks[3];       // ���ʷ��侵�淴���
	double n_exp;       // ���淴��߹�ϵ��
	double refl;		// �����ϵ��
	double tran;		// �����ϵ��
	double idx_ref;		// ������

public:
	// ���졢��������
	sphere(point *c, double r) { cent = c; rad = r; }
	~sphere() {};

	// ���� 
	void setIa(double r, double g, double b) { Ka[0] = r; Ka[1] = g; Ka[2] = b; }
	void setId(double r, double g, double b) { Kd[0] = r; Kd[1] = g; Kd[2] = b; }
	void setIr(double r, double g, double b) { Ks[0] = r; Ks[1] = g; Ks[2] = b; }
	void setn(double n) { n_exp = n; }
	void setReflection(double reflection) { refl = reflection; }
	void setTransmission(double transmission) { tran = transmission; }
	void setIndexOfRefraction(double idxRefraction) { idx_ref = idxRefraction; }
	point* center() { return cent; }
	double radius() { return rad; }
	double getIa(int color) { return Ka[color]; }
	double getId(int color) { return Kd[color]; }
	double getIr(int color) { return Ks[color]; }
	double getn() { return n_exp; }
	double getReflection() { return refl; }
	double getTransmission() { return tran; }
	double getIndexOfRefraction() { return idx_ref; }
};
#endif 
