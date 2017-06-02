#pragma once
#ifndef BOARD_H
#define BOARD_H


#include "point.h"

using namespace std;


class board {
private:
	// ��Ա����
	point *cent;
	double side;
	point *pa;
	point *pb;
	point *pc;
	double Ka[3];       // ���ʷ��价����
	double Kd[3];       // ���ʷ����������
	double Ks[3];       // ���ʷ��侵�淴���
	double n_exp;       // ���淴��߹�ϵ��
	double refl;		// �����ϵ��
	double tran;		// �����ϵ��
	double idx_ref;		// ������

public:
	// ���졢��������
	board(point *pt1, point *pt2, point *pt3) { pa = pt1; pb = pt2; pc = pt3; }
	~board() {};

	// ����
	void setIa(double r, double g, double b) { Ka[0] = r; Ka[1] = g; Ka[2] = b; }
	void setId(double r, double g, double b) { Kd[0] = r; Kd[1] = g; Kd[2] = b; }
	void setIr(double r, double g, double b) { Ks[0] = r; Ks[1] = g; Ks[2] = b; }
	void setn(double n) { n_exp = n; }
	void setReflection(double reflection) { refl = reflection; }
	void setTransmission(double transmission) { tran = transmission; }
	void setIndexOfRefraction(double idxRefraction) { idx_ref = idxRefraction; }
	point* center() { return cent; }
	point* point1() { return pa; }
	point* point2() { return pb; }
	point* point3() { return pc; } 
	double sides() { return side; }
	double getIa(int color) { return Ka[color]; }
	double getId(int color) { return Kd[color]; }
	double getIr(int color) { return Ks[color]; }
	double getn() { return n_exp; }
	double getReflection() { return refl; }
	double getTransmission() { return tran; }
	double getIndexOfRefraction() { return idx_ref; }
};
#endif 