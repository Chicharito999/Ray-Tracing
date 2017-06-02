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
//窗口尺寸
#define WIDTH 800
#define HEIGHT 800
//背景光RGB和全局环境光RGB
double background[3];
double ambientlight[3];
//光源的位置和RGB
point *p1 = new point(20, 20, -5);
light *lightSource = new light(p1, 1.0, 1.0, 1.0);
//银球的中心点和半径
point *p2 = new point(15.0, -5.0, -80.0);
sphere *sph1 = new sphere(p2, 12.0);
//其他球的中心点和半径
point *p3 = new point(-8.0, -12.0, -70.0);
sphere *sph2 = new sphere(p3, 8.0);
//圆柱的中心点、半径和高度
point *pc = new point(-20.0, -20.0, -50.0);
cylinder *cyl = new cylinder(pc, 6.0, 6.0);
//板的支点
point *p5, *p6, *p7, *p8;
board *bd1, *bd2;

sphere *spheres[2];
cylinder *olist[1];
board *blist[2];

int i;
double tf, tg;
//声明函数
double* raytrace(point*, point*, int);
void setvalues();
bool visible(point*, point*, point*, point*);

// 初始化函数——初始化Opengl windows
void init() {
	glMatrixMode(GL_PROJECTION); // 下面对投影变换矩阵操作
	glLoadIdentity();
	gluOrtho2D(0, WIDTH, 0, HEIGHT); // 定义正交投影矩阵

	glClearColor(1.0, 1.0, 1.0, 1.0); // 设置背景颜色
}
// display function - 逐屏幕像素绘制
void display(void) {
	/* 清理屏幕颜色*/
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POINTS);
	//	逐像素绘制
	point *p0 = new point(0, 0, 0);
	for (int y = 0; y<HEIGHT; y++) {
		for (int x = 0; x<WIDTH; x++) {
			double vx = ((double)x + 0.5) / ((double)WIDTH / 2.0) - 1.0;
			double vy = (((double)y + 0.5) / ((double)HEIGHT / 2.0) - 1.0)*.899;
			//vx,vy 是视线与屏幕的交点
			point *p1 = new point(vx, vy, -1);
			double *pixel = new double[3];
			pixel = raytrace(p0, p1, 3);//调用raytrace函数
			glColor3f(pixel[0], pixel[1], pixel[2]);
			glVertex2i(x, y);//绘制这个像素

		}
	}
	glEnd();

	glutSwapBuffers();
}
//键盘函数
void keyboard(unsigned char key, int x, int y) {}
//特殊键盘函数
void specialKey(int key, int x, int y) {}
//鼠标点击函数
void mouseButton(int button, int state, int x, int y) {}
// idle function
void idleFunc() { glutPostRedisplay(); } // 重绘函数

int main(int argc, char** argv) {

	glutInit(&argc, argv);

	// 设置初始显示模式
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// 设置窗口尺寸
	glutInitWindowSize(WIDTH, HEIGHT);
	// 创建一个窗口，命名为Project 3
	glutCreateWindow("Project 3");

	/*调用setvalue函数初始化数据*/
	setvalues();
	/* 设置逻辑图形空间 */
	init();

	/* 接收键盘输入 */
	glutKeyboardFunc(keyboard);
	/* 接收键盘输入 */
	glutSpecialFunc(specialKey);
	/* 接受鼠标输入*/
	glutMouseFunc(mouseButton);
	/* 分配display function */
	glutDisplayFunc(display);
	/* 分配idle function */
	glutIdleFunc(idleFunc);

	glutMainLoop();
	return (0);
}

void setvalues()
{

	//背景颜色
	background[0] = 1.0;
	background[1] = 0.5;
	background[2] = 1.0;
	//环境光	
	ambientlight[0] = 1.0;
	ambientlight[1] = 1.0;
	ambientlight[2] = 1.0;
	// SPHERE -- 银球的其余属性
	//设置材质反射环境光Ka
	sph1->setIa(0.231, 0.231, 0.231);
	//设置材质反射漫反射光Kd
	sph1->setId(0.278, 0.278, 0.278);
	//设置材质反射镜面反射光Ks
	sph1->setIr(0.774, 0.774, 0.774);
	//设置高光系数n
	sph1->setn(89.6);
	//设置反射光系数reflection
	sph1->setReflection(0.65);
	//设置折射光系数transmission
	sph1->setTransmission(0.0);
	//设置折射率refractive index
	sph1->setIndexOfRefraction(1.0);

	// SPHERE -- 哑光粉球的其余属性
	//设置材质反射环境光Ka
	sph2->setIa(0.3, 0.16, 0.12);
	//设置材质反射漫反射光Kd
	sph2->setId(0.8, 0.4, 0.35);
	//设置材质反射镜面反射光Ks
	sph2->setIr(0.1, 0.1, 0.1);
	//设置高光系数n
	sph2->setn(5.0);
	//设置反射光系数reflection
	sph2->setReflection(0.65);
	//设置折射光系数transmission
	sph2->setTransmission(0.0);
	//设置折射率refractive index
	sph2->setIndexOfRefraction(1.0);

	// CYLINDER -- 红色哑光圆柱的其余属性
	//设置材质反射环境光Ka
	cyl->setIa(1.0, 0.0, 0.0);
	//设置材质反射漫反射光Kd
	cyl->setId(1.0, 0.0, 0.1);
	//设置材质反射镜面反射光Ks
	cyl->setIr(0.1, 0.1, 0.1);
	//设置高光系数n
	cyl->setn(5.0);
	//设置反射光系数reflection
	cyl->setReflection(0.65);
	//设置折射光系数transmission
	cyl->setTransmission(0.0);
	//设置折射率refractive index
	cyl->setIndexOfRefraction(1.0);

	// Board -- 哑光灰色  的其余属性
	//板的支点
	p5 = new point(-40.0, -20.0, -20.0);
	p6 = new point(40.0, -20.0, -20.0);
	p7 = new point(40.0, -20.0, -100.0);
	p8 = new point(-40.0, -20.0, -100.0);
	bd1 = new board(p5, p6, p7);
	//设置材质反射环境光Ka
	bd1->setIa(0.3, 0.3, 0.02);
	//设置材质反射漫反射光Kd
	bd1->setId(0.8, 0.8, 0.1);
	//设置材质反射镜面反射光Ks
	bd1->setIr(0.1, 0.1, 0.1);
	//设置高光系数n
	bd1->setn(5.0);
	//设置反射光系数reflection
	bd1->setReflection(0.65);
	//设置折射光系数transmission
	bd1->setTransmission(0.0);
	//设置折射率refractive index
	bd1->setIndexOfRefraction(1.0);

	bd2 = new board(p5, p7, p8);
	//设置材质反射环境光Ka
	bd2->setIa(0.3, 0.3, 0.02);
	//设置材质反射漫反射光Kd
	bd2->setId(0.8, 0.8, 0.1);
	//设置材质反射镜面反射光Ks
	bd2->setIr(0.1, 0.1, 0.1);
	//设置高光系数n
	bd2->setn(5.0);
	//设置反射光系数reflection
	bd2->setReflection(0.65);
	//设置折射光系数transmission
	bd2->setTransmission(0.0);
	//设置折射率refractive index
	bd2->setIndexOfRefraction(1.0);

}

/*ray tracing 函数*/
double* raytrace(point* p0, point* p1, int count) {

	if (count>3) {	//核实最大反弹次数
		double* actual = new double[3];
		actual[0] = 0; actual[1] = 0; actual[2] = 0;
		return actual;
	}

	/*对于球*/
	spheres[0] = sph1;
	spheres[1] = sph2;
	/*对于圆柱和板*/
	olist[0] = cyl;	//for cylinder
	blist[0] = bd1;	//for board
	blist[1] = bd2;
	/*计算视线向量*/
	double dx = p1->x() - p0->x();
	double dy = p1->y() - p0->y();
	double dz = p1->z() - p0->z();
	double t[50];
	double red, green, blue;

	//for sphere：寻找ray-sphere交点，利用公式求解
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

	//for board：寻找ray-board交点，利用公式求解
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
		//核实可见性并且更新交点坐标
		if (visible(p, blist[i]->point1(), blist[i]->point2(), blist[i]->point3())) {
			if (visible(p, blist[i]->point2(), blist[i]->point1(), blist[i]->point3())) {
				//if (visible(p, blist[i]->point3(), blist[i]->point1(), blist[i]->point2())) {
					t[i + 2] = intersect;
				//}
			}
		}
		delete p;
	}

	//for cylinder：寻找ray-cylinder交点，利用公式求解
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

			//区分处理圆柱体的上表面和下表面
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



	//寻找t最小的交点，即可见交点的t
	double min = t[0];
	int mindex = 0;
	for (int i = 1; i <= 3; i++) {
		if (t[i]<min && t[i]>0.02) {
			min = t[i];
			mindex = i;
		}
	}

	//如果没有交点，设为背景颜色
	if (min == 1000) {
		double* res = new double[3];
		res[0] = background[0];
		res[1] = background[1];
		res[2] = background[2];
		return res;
	}

	//通过ray的式子找到交点坐标
	double x = p0->x() + t[mindex] * dx;
	double y = p0->y() + t[mindex] * dy;
	double z = p0->z() + t[mindex] * dz;
	point* newp0 = new point(x, y, z);
	double vmag = sqrt((p0->x() - x)*(p0->x() - x) + (p0->y() - y)*(p0->y() - y) + (p0->z() - z)*(p0->z() - z));
	point* view = new point((p0->x() - x) / vmag, (p0->y() - y) / vmag, (p0->z() - z) / vmag);//V视线的单位法向量
	double cy = y;
	if (mindex<2) {
		sphere* sphere = spheres[mindex];
		//球体的法线
		point* normal = new point((x - sphere->center()->x()) / sphere->radius(), (y - sphere->center()->y()) / sphere->radius(), (z - sphere->center()->z()) / sphere->radius());//N交点的单位法向量
		double Lmod = sqrt((lightSource->location()->x() - x)*(lightSource->location()->x() - x) + (lightSource->location()->y() - y)*(lightSource->location()->y() - y) + (lightSource->location()->z() - z)*(lightSource->location()->z() - z));
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L入射光线的单位法向量
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos（N*L）法向量与入射光线的夹角cos值
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos（R*V）反射光线与视线的夹角cos值


		double specular = pow(rv, sphere->getn());//高光系数：cos（R*V）的n次幂 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//计算视线的反射线和视线的折射线
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos（N*V）:法向量与视线的夹角cos值
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//视线的反射点

		double Kn;
		if (count % 2 == 1) { Kn = sphere->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / sphere->getIndexOfRefraction(); }
		double anglei = acos(nv);//法向量与视线的夹角
		double angler = asin(Kn*sin(anglei));//折射角
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
		point* T = new point(tx + x, ty + y, tz + z);//视线的折射点

		refl = raytrace(newp0, R, count + 1);
		refr = raytrace(newp0, T, count + 1);//递归计算

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		if (mindex == 0) {

		}

		//计算该像素点的RGB颜色
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
		// 板的法向量		
		point* normal = new point(0.0, 1.0, 0.0);//N交点的单位法向量
		double Lmod = sqrt((lightSource->location()->x() - x)*(lightSource->location()->x() - x) + (lightSource->location()->y() - y)*(lightSource->location()->y() - y) + (lightSource->location()->z() - z)*(lightSource->location()->z() - z));
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L入射光线的单位法向量
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos（N*L）法向量与入射光线的夹角cos值
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos（R*V）反射光线与视线的夹角cos值

		double specular = pow(rv, board->getn());//高光系数：cos（R*V）的n次幂 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//计算视线的反射线和视线的折射线
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos（N*V）:法向量与视线的夹角cos值
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//视线的反射点

		double Kn;
		if (count % 2 == 1) { Kn = board->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / board->getIndexOfRefraction(); }
		double anglei = acos(nv);//法向量与视线的夹角
		double angler = asin(Kn*sin(anglei));//折射角
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
		point* T = new point(tx + x, ty + y, tz + z);//视线的折射点


		refl = raytrace(newp0, R, count + 1);//递归计算
		refr = raytrace(newp0, T, count + 1);

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		//计算该像素点的RGB颜色
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
		//计算圆柱的法向量
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
		point* L = new point((lightSource->location()->x() - x) / Lmod, (lightSource->location()->y() - y) / Lmod, (lightSource->location()->z() - z) / Lmod);//L入射光线的单位法向量
		double nl = L->x()*normal->x() + L->y()*normal->y() + L->z()*normal->z();//cos（N*L）法向量与入射光线的夹角cos值
		double rv = (2 * normal->x()*nl - L->x())*view->x() + (2 * normal->y()*nl - L->y())*view->y() + (2 * normal->z()*nl - L->z())*view->z();//cos（R*V）反射光线与视线的夹角cos值


		double specular = pow(rv, cylinder->getn());//高光系数：cos（R*V）的n次幂 
		if (!(specular > 0.0)) specular = 0.0;

		if (nl < 0.0) nl = 0.0;

		//计算视线的反射线和视线的折射线
		double *refl = new double[3];
		double *refr = new double[3];

		double nv = normal->x()*view->x() + normal->y()*view->y() + normal->z()*view->z();//cos（N*V）:法向量与视线的夹角cos值
		double viewreflx = (2 * normal->x()*nv - view->x());
		double viewrefly = (2 * normal->y()*nv - view->y());
		double viewreflz = (2 * normal->z()*nv - view->z());

		point* R = new point(x + viewreflx, y + viewrefly, z + viewreflz);//视线的反射点

		double Kn;
		if (count % 2 == 1) { Kn = cylinder->getIndexOfRefraction() / 1.0; }
		else { Kn = 1.0 / cylinder->getIndexOfRefraction(); }
		double anglei = acos(nv);//法向量与视线的夹角
		double angler = asin(Kn*sin(anglei));//折射角
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
		point* T = new point(tx + x, ty + y, tz + z);//视线的折射点


		refl = raytrace(newp0, R, count + 1);//递归计算
		refr = raytrace(newp0, T, count + 1);

		if (!(refl[0] >= 0.0 && refl[0] <= 100.0)) refl[0] = 0.0;
		if (!(refl[1] >= 0.0 && refl[1] <= 100.0)) refl[1] = 0.0;
		if (!(refl[2] >= 0.0 && refl[2] <= 100.0)) refl[2] = 0.0;

		if (!(refr[0] >= 0.0 && refr[0] <= 100.0)) refr[0] = 0.0;
		if (!(refr[1] >= 0.0 && refr[1] <= 100.0)) refr[1] = 0.0;
		if (!(refr[2] >= 0.0 && refr[2] <= 100.0)) refr[2] = 0.0;

		if (mindex == 0) {

		}

		//计算该像素点的RGB颜色
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

/*可见性函数*/
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