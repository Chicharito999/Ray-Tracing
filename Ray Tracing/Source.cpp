#include <windows.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "point.h"
#include "sphere.h"
#include "light.h"
#include "board.h"
#include "cylinder.h"
//���ڳߴ�
#define WIDTH 800
#define HEIGHT 800
//������RGB��ȫ�ֻ�����RGB
double background[3];
double ambientlight[3];
//��Դ��λ�ú�RGB
point *p1 = new point(20, 20, -5);
light *lightSource = new light(p1, 1.0, 1.0, 1.0);
//��������ĵ�Ͱ뾶
point *p2 = new point(15.0, -5.0, -80.0);
sphere *sph1 = new sphere(p2, 12.0);
//����������ĵ�Ͱ뾶
point *p3 = new point(-8.0, -12.0, -70.0);
sphere *sph2 = new sphere(p3, 8.0);
//Բ�������ĵ㡢�뾶�͸߶�
point *pc = new point(-20.0, -20.0, -50.0);
cylinder *cyl = new cylinder(pc, 6.0, 6.0);
//���֧��
point *p5, *p6, *p7, *p8;
board *bd1, *bd2;

sphere *spheres[2];
cylinder *olist[1];
board *blist[2];

int i;
double tf, tg;
//��������
double* raytrace(point*, point*, int);
void setvalues();
bool visible(point*, point*, point*, point*);

// ��ʼ������������ʼ��Opengl windows
void init() {
	glMatrixMode(GL_PROJECTION); // �����ͶӰ�任�������
	glLoadIdentity();
	gluOrtho2D(0, WIDTH, 0, HEIGHT); // ��������ͶӰ����

	glClearColor(1.0, 1.0, 1.0, 1.0); // ���ñ�����ɫ
}
// display function - ����Ļ���ػ���
void display(void) {
	/* ������Ļ��ɫ*/
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POINTS);
	//	�����ػ���
	point *p0 = new point(0, 0, 0);
	for (int y = 0; y<HEIGHT; y++) {
		for (int x = 0; x<WIDTH; x++) {
			double vx = ((double)x + 0.5) / ((double)WIDTH / 2.0) - 1.0;
			double vy = (((double)y + 0.5) / ((double)HEIGHT / 2.0) - 1.0)*.899;
			//vx,vy ����������Ļ�Ľ���
			point *p1 = new point(vx, vy, -1);
			double *pixel = new double[3];
			pixel = raytrace(p0, p1, 3);//����raytrace����
			glColor3f(pixel[0], pixel[1], pixel[2]);
			glVertex2i(x, y);//�����������

		}
	}
	glEnd();

	glutSwapBuffers();
}
//���̺���
void keyboard(unsigned char key, int x, int y) {}
//������̺���
void specialKey(int key, int x, int y) {}
//���������
void mouseButton(int button, int state, int x, int y) {}
// idle function
void idleFunc() { glutPostRedisplay(); } // �ػ溯��

int main(int argc, char** argv) {

	glutInit(&argc, argv);

	// ���ó�ʼ��ʾģʽ
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// ���ô��ڳߴ�
	glutInitWindowSize(WIDTH, HEIGHT);
	// ����һ�����ڣ�����ΪProject 3
	glutCreateWindow("Project 3");

	/*����setvalue������ʼ������*/
	setvalues();
	/* �����߼�ͼ�οռ� */
	init();

	/* ���ռ������� */
	glutKeyboardFunc(keyboard);
	/* ���ռ������� */
	glutSpecialFunc(specialKey);
	/* �����������*/
	glutMouseFunc(mouseButton);
	/* ����display function */
	glutDisplayFunc(display);
	/* ����idle function */
	glutIdleFunc(idleFunc);

	glutMainLoop();
	return (0);
}

void setvalues()
{

	//������ɫ
	background[0] = 1.0;
	background[1] = 0.5;
	background[2] = 1.0;
	//������	
	ambientlight[0] = 1.0;
	ambientlight[1] = 1.0;
	ambientlight[2] = 1.0;
	// SPHERE -- �������������
	//���ò��ʷ��价����Ka
	sph1->setIa(0.231, 0.231, 0.231);
	//���ò��ʷ����������Kd
	sph1->setId(0.278, 0.278, 0.278);
	//���ò��ʷ��侵�淴���Ks
	sph1->setIr(0.774, 0.774, 0.774);
	//���ø߹�ϵ��n
	sph1->setn(89.6);
	//���÷����ϵ��reflection
	sph1->setReflection(0.65);
	//���������ϵ��transmission
	sph1->setTransmission(0.0);
	//����������refractive index
	sph1->setIndexOfRefraction(1.0);

	// SPHERE -- �ƹ�������������
	//���ò��ʷ��价����Ka
	sph2->setIa(0.3, 0.16, 0.12);
	//���ò��ʷ����������Kd
	sph2->setId(0.8, 0.4, 0.35);
	//���ò��ʷ��侵�淴���Ks
	sph2->setIr(0.1, 0.1, 0.1);
	//���ø߹�ϵ��n
	sph2->setn(5.0);
	//���÷����ϵ��reflection
	sph2->setReflection(0.65);
	//���������ϵ��transmission
	sph2->setTransmission(0.0);
	//����������refractive index
	sph2->setIndexOfRefraction(1.0);

	// CYLINDER -- ��ɫ�ƹ�Բ������������
	//���ò��ʷ��价����Ka
	cyl->setIa(1.0, 0.0, 0.0);
	//���ò��ʷ����������Kd
	cyl->setId(1.0, 0.0, 0.1);
	//���ò��ʷ��侵�淴���Ks
	cyl->setIr(0.1, 0.1, 0.1);
	//���ø߹�ϵ��n
	cyl->setn(5.0);
	//���÷����ϵ��reflection
	cyl->setReflection(0.65);
	//���������ϵ��transmission
	cyl->setTransmission(0.0);
	//����������refractive index
	cyl->setIndexOfRefraction(1.0);

	// Board -- �ƹ��ɫ  ����������
	//���֧��
	p5 = new point(-40.0, -20.0, -20.0);
	p6 = new point(40.0, -20.0, -20.0);
	p7 = new point(40.0, -20.0, -100.0);
	p8 = new point(-40.0, -20.0, -100.0);
	bd1 = new board(p5, p6, p7);
	//���ò��ʷ��价����Ka
	bd1->setIa(0.3, 0.3, 0.02);
	//���ò��ʷ����������Kd
	bd1->setId(0.8, 0.8, 0.1);
	//���ò��ʷ��侵�淴���Ks
	bd1->setIr(0.1, 0.1, 0.1);
	//���ø߹�ϵ��n
	bd1->setn(5.0);
	//���÷����ϵ��reflection
	bd1->setReflection(0.65);
	//���������ϵ��transmission
	bd1->setTransmission(0.0);
	//����������refractive index
	bd1->setIndexOfRefraction(1.0);

	bd2 = new board(p5, p7, p8);
	//���ò��ʷ��价����Ka
	bd2->setIa(0.3, 0.3, 0.02);
	//���ò��ʷ����������Kd
	bd2->setId(0.8, 0.8, 0.1);
	//���ò��ʷ��侵�淴���Ks
	bd2->setIr(0.1, 0.1, 0.1);
	//���ø߹�ϵ��n
	bd2->setn(5.0);
	//���÷����ϵ��reflection
	bd2->setReflection(0.65);
	//���������ϵ��transmission
	bd2->setTransmission(0.0);
	//����������refractive index
	bd2->setIndexOfRefraction(1.0);

}

/*ray tracing ����*/
double* raytrace(point* p0, point* p1, int count) {

	if (count>3) {	//��ʵ��󷴵�����
		double* actual = new double[3];
		actual[0] = 0; actual[1] = 0; actual[2] = 0;
		return actual;
	}

	/*������*/
	spheres[0] = sph1;
	spheres[1] = sph2;
	/*����Բ���Ͱ�*/
	olist[0] = cyl;	//for cylinder
	blist[0] = bd1;	//for board
	blist[1] = bd2;
	/*������������*/
	double dx = p1->x() - p0->x();
	double dy = p1->y() - p0->y();
	double dz = p1->z() - p0->z();
	double t[50];
	double red, green, blue;

	//for sphere��Ѱ��ray-sphere���㣬���ù�ʽ���
	for (i = 0; i<2; i++) {
		double a, b, c, t1, t2;
		point *center = new point(spheres[i]->center()->x(), spheres[i]->center()->y(), spheres[i]->center()->z());
		a = dx*dx + dy*dy + dz*dz;
		b = 2.0*(dx*(p0->x() - center->x()) + dy*(p0->y() - center->y()) + dz*(p0->z() - center->z()));
		c = (p0->x() - center->x())*(p0->x() - center->x()) + (p0->y() - center->y())*(p0->y() - center->y()) + (p0->z() - center->z())*(p0->z() - center->z()) - (spheres[i]->radius())*(spheres[i]->radius());
		t[i] = 1000;

		delete center;

		if (((b*b) - (4.0*a*c)) >= 0) {
			t1 = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
			t2 = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
			double x, y, z;
			if (t1 <= t2) {
				t[i] = t1;
			}
			else {
				t[i] = t2;
			}
		}
	}

	//for board��Ѱ��ray-board���㣬���ù�ʽ���
	for (int i = 0; i<2; i++) {
		double a, b, c, d, intersect;
	
		point *ab = new point(blist[i]->point2()->x() - blist[i]->point1()->x(), blist[i]->point2()->y() - blist[i]->point1()->y(), blist[i]->point2()->z() - blist[i]->point1()->z());
		point *ad = new point(blist[i]->point3()->x() - blist[i]->point2()->x(), blist[i]->point3()->y() - blist[i]->point2()->y(), blist[i]->point3()->z() - blist[i]->point2()->z());
		a = ab->y()*ad->z() - ab->z()*ad->y();
		b = ab->z()*ad->x() - ab->x()*ad->z();
		c = ab->x()*ad->y() - ab->y()*ad->x();
		d = -(blist[i]->point3()->x()*a + blist[i]->point3()->y()*b + blist[i]->point3()->z()*c);

		
		intersect = -((a*p0->x() + b*p0->y() + c*p0->z() + d) / (a*dx + b*dy + c*dz));

		point *p = new point(p0->x() + intersect*dx, p0->y() + intersect*dy, p0->z() + intersect*dz);

		delete ab, ad;
		t[i + 2] = 1000;
		//��ʵ�ɼ��Բ��Ҹ��½�������
		if (visible(p, blist[i]->point1(), blist[i]->point2(), blist[i]->point3())) {
			if (visible(p, blist[i]->point2(), blist[i]->point1(), blist[i]->point3())) {
				//if (visible(p, blist[i]->point3(), blist[i]->point1(), blist[i]->point2())) {
					t[i + 2] = intersect;
				//}
			}
		}
		delete p;
	}

	//for cylinder��Ѱ��ray-cylinder���㣬���ù�ʽ���
	for (i = 0; i<1; i++) {

		double a, b, c, t1, t2;

		point *center = new point(olist[i]->center()->x(), olist[i]->center()->y(), olist[i]->center()->z());
		a = dx*dx + dz*dz;
		b = 2.0*(dx*(p0->x() - center->x()) + dz*(p0->z() - center->z()));
		c = (p0->x() - center->x())*(p0->x() - center->x()) + (p0->z() - center->z())*(p0->z() - center->z()) - (olist[i]->radius())*(olist[i]->radius());

		delete center;

		if (((b*b) - (4.0*a*c)) >= 0) {
			t1 = (-b + sqrt(b*b - 4.0*a*c)) / (2.0*a);
			t2 = (-b - sqrt(b*b - 4.0*a*c)) / (2.0*a);
			double x, y, z;
			if (t1 <= t2) {
				t[2 + i] = t1;
				tf = t1;
				tg = t2;
			}
			else {
				t[2 + i] = t2;
				tf = t2;
				tg = t1;
			}

			//���ִ���Բ������ϱ�����±���
			if (((p0->y() + t[2 + i] * p1->y()) > (cyl->center()->y() + cyl->height())) && (((p0->y() + tg*p1->y()) > (cyl->center()->y() + cyl->height()))))
				t[2 + i] = -1;
			else if (((p0->y() + t[2 + i] * p1->y()) >= (cyl->center()->y() + cyl->height())) && (((p0->y() + tg*p1->y()) <= (cyl->center()->y() + cyl->height()))))
				t[2 + i] = ((cyl->center()->y() + cyl->height()) - (p0->y())) / (p1->y());

			if (((p0->y() + t[2 + i] * p1->y()) < (cyl->center()->y())) && (((p0->y() + tg*p1->y()) < (cyl->center()->y()))))
				t[2 + i] = -1;
			else if (((p0->y() + t[2 + i] * p1->y()) <= (cyl->center()->y())) && (((p0->y() + tg*p1->y()) >= (cyl->center()->y()))))
				t[2 + i] = ((cyl->center()->y()) - (p0->y())) / (p1->y());
		}
	}



	//Ѱ��t��С�Ľ��㣬���ɼ������t
	double min = t[0];
	int mindex = 0;
	for (int i = 1; i <= 3; i++) {
		if (t[i]<min && t[i]>0.02) {
			min = t[i];
			mindex = i;
		}
	}

	//���û�н��㣬��Ϊ������ɫ
	if (min == 1000) {
		double* res = new double[3];
		res[0] = background[0];
		res[1] = background[1];
		res[2] = background[2];
		return res;
	}

	//ͨ��ray��ʽ���ҵ���������
	double x = p0->x() + t[mindex] * dx;
	double y = p0->y() + t[mindex] * dy;
	double z = p0->z() + t[mindex] * dz;
	point* newp0 = new point(x, y, z);
	double vmag = sqrt((p0->x() - x)*(p0->x() - x) + (p0->y() - y)*(p0->y() - y) + (p0->z() - z)*(p0->z() - z));
	point* view = new point((p0->x() - x) / vmag, (p0->y() - y) / vmag, (p0->z() - z) / vmag);//V���ߵĵ�λ������
	double cy = y;
	if (mindex<2) {
		sphere* sphere = spheres[mindex];
		//����ķ���
		point* normal = new point((x - sphere->center()->x()) / sphere->radius(), (y - sphere->center()->y()) / sphere->radius(), (z - sphere->center()->z()) / sphere->radius());//N����ĵ�λ������
		double Lmod = sqrt((lightSource->location()->x() - x)*(lightSource->location()->x() - x) + (lightSource->location()->y() - y)*(lightSource->location()->y() - y) + (lightSource->location()->z() - z)*(lightSource->location()->z() - z));
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L������ߵĵ�λ������
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos��N*L����������������ߵļн�cosֵ
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos��R*V��������������ߵļн�cosֵ


		double specular = pow(rv, sphere->getn());//�߹�ϵ����cos��R*V����n���� 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//�������ߵķ����ߺ����ߵ�������
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos��N*V��:�����������ߵļн�cosֵ
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//���ߵķ����

		double Kn;
		if (count % 2 == 1) { Kn = sphere->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / sphere->getIndexOfRefraction(); }
		double anglei = acos(nv);//�����������ߵļн�
		double angler = asin(Kn*sin(anglei));//�����
		double tx, ty, tz;

		if (nv>0) {
			tx = -Kn*view->x() + (Kn*cos(anglei) - cos(angler))*normal->x();
			ty = -Kn*view->y() + (Kn*cos(anglei) - cos(angler))*normal->y();
			tz = -Kn*view->z() + (Kn*cos(anglei) - cos(angler))*normal->z();
		}
		else {
			tx = -Kn*view->x() - (Kn*cos(anglei) - cos(angler))*normal->x();
			ty = -Kn*view->y() - (Kn*cos(anglei) - cos(angler))*normal->y();
			tz = -Kn*view->z() - (Kn*cos(anglei) - cos(angler))*normal->z();
		}
		point* T = new point(tx + x, ty + y, tz + z);//���ߵ������

		refl = raytrace(newp0, R, count + 1);
		refr = raytrace(newp0, T, count + 1);//�ݹ����

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		if (mindex == 0) {

		}

		//��������ص��RGB��ɫ
		red = ambientlight[0] * sphere->getIa(0) +
			sphere->getId(0)*lightSource->red()*nl +
			sphere->getIr(0)*lightSource->red()*specular +
			sphere->getReflection()*refl[0] +
			sphere->getTransmission()*refr[0];

		green = ambientlight[1] * sphere->getIa(1) +
			sphere->getId(1)*lightSource->green()*nl +
			sphere->getIr(1)*lightSource->green()*specular +
			sphere->getReflection()*refl[1] +
			sphere->getTransmission()*refr[1];

		blue = ambientlight[2] * sphere->getIa(2) +
			sphere->getId(2)*lightSource->blue()*nl +
			sphere->getIr(2)*lightSource->blue()*specular +
			sphere->getReflection()*refl[2] +
			sphere->getTransmission()*refr[2];


		delete normal, L, R, T;

	}
	else if ((mindex >= 2) && (mindex<4)) {
		board* board = blist[mindex - 2];
		// ��ķ�����		
		point* normal = new point(0.0, 1.0, 0.0);//N����ĵ�λ������
		double Lmod = sqrt((lightSource->location()->x() - x)*(lightSource->location()->x() - x) + (lightSource->location()->y() - y)*(lightSource->location()->y() - y) + (lightSource->location()->z() - z)*(lightSource->location()->z() - z));
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L������ߵĵ�λ������
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos��N*L����������������ߵļн�cosֵ
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos��R*V��������������ߵļн�cosֵ

		double specular = pow(rv, board->getn());//�߹�ϵ����cos��R*V����n���� 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//�������ߵķ����ߺ����ߵ�������
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos��N*V��:�����������ߵļн�cosֵ
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//���ߵķ����

		double Kn;
		if (count % 2 == 1) { Kn = board->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / board->getIndexOfRefraction(); }
		double anglei = acos(nv);//�����������ߵļн�
		double angler = asin(Kn*sin(anglei));//�����
		double tx, ty, tz;
		if (nv>0) {
			tx = -Kn*view->x() + (Kn*cos(anglei) - cos(angler))*normal->x();
			ty = -Kn*view->y() + (Kn*cos(anglei) - cos(angler))*normal->y();
			tz = -Kn*view->z() + (Kn*cos(anglei) - cos(angler))*normal->z();
		}
		else {
			tx = -Kn*view->x() + (Kn*cos(anglei) + cos(angler))*normal->x();
			ty = -Kn*view->y() + (Kn*cos(anglei) + cos(angler))*normal->y();
			tz = -Kn*view->z() + (Kn*cos(anglei) + cos(angler))*normal->z();
		}
		point* T = new point(tx + x, ty + y, tz + z);//���ߵ������


		refl = raytrace(newp0, R, count + 1);//�ݹ����
		refr = raytrace(newp0, T, count + 1);

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		//��������ص��RGB��ɫ
		red = ambientlight[0] * board->getIa(0) +
			board->getId(0)*lightSource->red()*nl +
			board->getIr(0)*lightSource->red()*specular +
			board->getReflection()*refl[0] +
			board->getTransmission()*refr[0];

		green = ambientlight[1] * board->getIa(1) +
			board->getId(1)*lightSource->green()*nl +
			board->getIr(1)*lightSource->green()*specular +
			board->getReflection()*refl[1] +
			board->getTransmission()*refr[1];

		blue = ambientlight[2] * board->getIa(2) +
			board->getId(2)*lightSource->blue()*nl +
			board->getIr(2)*lightSource->blue()*specular +
			board->getReflection()*refl[2] +
			board->getTransmission()*refr[2];


		delete normal, L, R, T;

	}

	else if (mindex == 4) {
		//����Բ���ķ�����
		cylinder* cylinder = cyl;
		point* normal;
		if ((p0->y() == (cylinder->center()->y() + cylinder->height())) || (p0->y() == cylinder->center()->y()))
		{
			normal = new point(0.0, 1.0, 0.0);
		}
		else
		{
			cylinder->center()->sety(cy);
			normal = new point((x - cylinder->center()->x()) / cylinder->radius(), (y - cylinder->center()->y()) / cylinder->radius(), (z - cylinder->center()->z()) / cylinder->radius());
		}

		double Lmod = sqrt((lightSource->location()->x() - x)*(lightSource->location()->x() - x) + (lightSource->location()->y() - y)*(lightSource->location()->y() - y) + (lightSource->location()->z() - z)*(lightSource->location()->z() - z));
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L������ߵĵ�λ������
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos��N*L����������������ߵļн�cosֵ
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos��R*V��������������ߵļн�cosֵ


		double specular = pow(rv, cylinder->getn());//�߹�ϵ����cos��R*V����n���� 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//�������ߵķ����ߺ����ߵ�������
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos��N*V��:�����������ߵļн�cosֵ
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//���ߵķ����

		double Kn;
		if (count % 2 == 1) { Kn = cylinder->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / cylinder->getIndexOfRefraction(); }
		double anglei = acos(nv);//�����������ߵļн�
		double angler = asin(Kn*sin(anglei));//�����
		double tx, ty, tz;
		if (nv>0) {
			tx = -Kn*view->x() + (Kn*cos(anglei) - cos(angler))*normal->x();
			ty = -Kn*view->y() + (Kn*cos(anglei) - cos(angler))*normal->y();
			tz = -Kn*view->z() + (Kn*cos(anglei) - cos(angler))*normal->z();
		}
		else {
			tx = -Kn*view->x() - (Kn*cos(anglei) - cos(angler))*normal->x();
			ty = -Kn*view->y() - (Kn*cos(anglei) - cos(angler))*normal->y();
			tz = -Kn*view->z() - (Kn*cos(anglei) - cos(angler))*normal->z();
		}
		point* T = new point(tx + x, ty + y, tz + z);//���ߵ������


		refl = raytrace(newp0, R, count + 1);//�ݹ����
		refr = raytrace(newp0, T, count + 1);

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		if (mindex == 0) {

		}

		//��������ص��RGB��ɫ
		red = ambientlight[0] * cylinder->getIa(0) +
			cylinder->getId(0)*lightSource->red()*nl +
			cylinder->getIr(0)*lightSource->red()*specular +
			cylinder->getReflection()*refl[0] +
			cylinder->getTransmission()*refr[0];

		green = ambientlight[1] * cylinder->getIa(1) +
			cylinder->getId(1)*lightSource->green()*nl +
			cylinder->getIr(1)*lightSource->green()*specular +
			cylinder->getReflection()*refl[1] +
			cylinder->getTransmission()*refr[1];

		blue = ambientlight[2] * cylinder->getIa(2) +
			cylinder->getId(2)*lightSource->blue()*nl +
			cylinder->getIr(2)*lightSource->blue()*specular +
			cylinder->getReflection()*refl[2] +
			cylinder->getTransmission()*refr[2];


		delete normal, L, R, T;

	}

	delete newp0, view;

	double *result = new double[3];
	result[0] = red;
	result[1] = green;
	result[2] = blue;

	return result;
}

/*�ɼ��Ժ���*/
bool visible(point* p1, point* p2, point* a, point* b) {
	double u, v, w;
	u = b->x() - a->x();
	v = b->y() - a->y();
	w = b->z() - a->z();
	double u1, v1, w1;
	u1 = p1->x() - a->x();
	v1 = p1->y() - a->y();
	w1 = p1->z() - a->z();

	point *cp1 = new point(v*w1 - w*v1, w*u1 - u*w1, u*v1 - v*u1);

	u1 = p2->x() - a->x();
	v1 = p2->y() - a->y();
	w1 = p2->z() - a->z();

	point *cp2 = new point(v*w1 - w*v1, w*u1 - u*w1, u*v1 - v*u1);

	double dotproduct = cp1->x()*cp2->x() + cp1->y()*cp2->y() + cp1->z()*cp2->z();

	delete cp1, cp2;
	if (dotproduct >= 0) { return true; }
	else { return false; }
}