# Ray-Tracing
:sunny:Realize at least two different method such as Phong model, Gouraud model<br>
__________________________________________________________________________________________
Author:赵明福                                        Student ID：201400301087                            E-mail:1109646702@qq.com<br>
## Assignment
Realize your ray-tracing algorithm and compare with the effect of opengl.<br>
　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片24.png) 
　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片25.png)
## Analysis
### 反向光线追踪：
　　现实生活中，我们看到的物体的色彩是由太阳发出的光线，经过物体的反射等最终到达人的眼睛。<br>
　　　　　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片57.png)<br>
　　在上图中，没有到达摄像机的射线都没有画出来，在现实生活中，没有进入我们眼睛的光线也不会被看到。所以我们不用以光源为起点追踪每一条光线，因为大多数光线最终都不会到达观察点，相反我们可以从摄像机的每个像素发出一条射线，来追踪他们到达的地方。<br>
　　　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片58.png)<br> 
　　在上图中，起点代表摄像机，矩形平面代表屏幕，我们以摄像机为起点向屏幕上的每个像素点发射射线，找到射线与空间中物体的第一个相交点，计算该点的颜色值，最后将该颜色值渲染到屏幕上。<br>
```cpp
逐像素渲染伪代码:
for (x,y) in screen
  {
     建立一个由摄像机穿过这个像素点的射线;
     找到这条射线碰到的第一个物体;
     测定在射线与物体的交点的颜色;
     将颜色画到像素上;
  }
```  
```cpp
逐像素渲染代码:
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
```    
### 光线类型：
* 环境光、漫反射光、镜面反射光：<br>
　　　　　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片59.png)<br>
* 反射光：如果物体表面具有反射性质，则部分光将会被反射出去，继续在场景中前进；同时根据物体的反射特性会显示部分反射光颜色<br>
　　　　　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片60.png)<br>
* 折射光：当物体表面具有折射性质并且部分透明，部分光线将会进入物体继续传播；同时根据物体的透明度会显示部分折射光颜色<br>
　　　　　　　　![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片61.png)<br>
```cpp
颜色值计算代码:根据表面性质(反射率、折射率)，和不同类型光线计算得出的颜色值，来确定交点的颜色值，即当前像素点的颜色值。
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
```    
### 交点计算：
* Ray-sphere：寻找ray-sphere交点，利用公式求解<br> 
```cpp
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
``` 
* Ray-cylinder：寻找ray-cylinder交点，利用公式求解<br> 
```cpp
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
``` 
* Ray-board：寻找ray-board交点，利用公式求解<br> 
```cpp
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
``` 
### Ray-tracing 算法：
```cpp
ray-tracing伪代码（最多递归3次）:
for each ray
  {
     compute ray-object intersection；
     if not exit{set up background color;}
     else{
       compute its reflection ray;
       ray-tracing(reflection ray);
       compute its refracted ray;
       ray-tracing(refracted ray);
       compute its phong model color(ambient、diffuse、specula)；
       pixel color = phong model color+ray-tracing(reflection ray)+ray-tracing(refracted ray);
     }    
  }

```  

## Display
![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片62.png)<br>
                                    　　　　　　　　　　　　　效果图<br>
　　可以在效果图里看到在黄色板子下面若隐若现有一个小圆柱体，这是因为下方确实放置了一个小圆柱体，由于上方的黄板具有一定的折射特性：透明度Transmission和折射率IndexOfRefraction，所以透过黄板可以大约看到圆柱的轮廓；除此之外，该算法还结合了反射光，也能略微看到一定的反射效果。无论是折射还是反射，这里设置的最大反弹次数为3次，即最多可以递归求解反射光或折射光3次，这也是为了提高速度，如果反弹次数过多的话，时间复杂度会非常大，渲染时间会非常长。<br>
　　
