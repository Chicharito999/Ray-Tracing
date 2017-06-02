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
### 古罗着色 冯氏着色 Blinn冯氏着色：
* 古罗着色：又叫做逐顶点着色，故名思意跟顶点有关，也就是在我们的顶点着色器中根据每个顶点上的入射向量L、法向量N、观察向量V等直接计算出每个顶点该有的颜色，然后传递给片段着色器进行插值着色，可想而知由于顶点是离散的，片段是连续的，所以引起着色效果的不光滑很容易理解。<br> 
* 冯氏着色：与古罗着色对应即是在片段着色器中，对法向量与坐标进行插值，然后再通过冯氏反射模型计算出每个像素点的颜色值，从而使离散的顶点计算出来的离散的颜色变得连续而光滑。我们直接把环境光、漫反射光、镜面反射光的计算拿到片段着色器中计算即可完成修改，那么法向量、观察向量、入射向量同理需要传递给片段着色器，而不再是直接传递一个颜色。<br> 
* Blinn冯氏着色：类似于与冯氏着色模型，只是Blinn-Phong模型镜面光的计算，采用了半角向量，这个向量是光照向量L和观察向量V的取中向量H，通过计算这个取中向量H与法向量N的夹角来得到镜面发射光强度，如下图所示：<br> 
![](https://github.com/Chicharito999/ImageCache/raw/master/image/图片28.png)<br>
## Code
* 光照计算
```cg
      output OUT;
      float3 N = normalize(IN.normal);//计算法向量
      float3 P = IN.objectPos;
      float3 L = normalize(LightPosition - P);//计算入射向量
      float NdotL = max(dot(N,L),0);//入射向量与法向量夹角
      float3 ambient = Ka * I;//环境光
      float3 diffuse = Kd * I * NdotL;//漫反射光
      float3 V = normalize(eyePosition - P);//眼睛-物体连线向量  

 //Phong光照模型
      float3 R=reflect(-L,N);//计算反射光线向量
      R=normalize(R);
      float NdotH = pow(max(dot(V,R), 0), shininess);//反射光与视角的夹角
 //BlinnPhong光照模型
      //  float3 H = normalize(L+V);
      //  float NdotH = pow(max(dot(N,H), 0), shininess);

      if(NdotL<=0)
           NdotH = 0.0;
      float3 specular = Ks*I*NdotH;//镜面发射光
      float3 color = ambient + diffuse + specular;//所有成分相加
      OUT.color.xyz= color;
      OUT.color.w = 1.0;
      return OUT;
```
* 古罗着色<br>
01vs.cg:<br>
```cg
struct output
{
      float4 position : POSITION; 
      float4 color     : COLOR; 
};
 
output vs_main( float4 position : POSITION,
                     float3 normal   : NORMAL,
                     uniform float4x4 MV, // 在相机坐标系中计算，所以要用到ModelView变换矩阵
                     uniform float4x4 MVP // ModelViewProjection变换矩阵
                    )
{

   //光照计算

}
```
01fs.cg:<br>
```cg
float4 fs_main( float4 color    : COLOR ) : COLOR
{
      return color;
}  
```
　　在顶点着色器中对所有顶点进行光照计算，得到每个顶点的颜色，在将输出的所有顶点颜色作为参数输入片段着色器，片段着色器根据每个三角形片的顶点颜色差值得到所有像素的颜色，从而实现了古罗着色模型。<br>
* 冯氏着色<br>
02vs.cg:<br>
```cg
struct output
{
      float4 position  : POSITION;    
      float3 objectPos : TEXCOORD0;  //顶点坐标传入fragment shader ，然后插值
      float3 normal     : TEXCOORD1;//法向量传入fragment shader，然后插值
};
 
output vs_main( float4 position : POSITION,
                  float3 normal   : NORMAL,
                  uniform float4x4 MV,
                  uniform float4x4 MVP
                    )
{
      output OUT;
      OUT.position = mul(MVP, position);//顶点位置转换到裁剪坐标
      OUT.objectPos = mul(MV, position).xyz;//转换到相机坐标
      OUT.normal = mul(MV, float4(normal,0.0)).xyz;//计算相机坐标下的法向量
 
      return OUT;
}
```
02fs.cg:<br>
```cg
struct input{//传入每个顶点的在相机坐标的位置信息和法向量，插值得到每个像素点的
      float3 objectPos: TEXCOORD0;   
      float3 normal   : TEXCOORD1;
};
 
struct output{
      float4 color     : COLOR;
};
 
output fs_main( in input IN )//为每个像素点计算color
{
     //光照计算

}
```
　　在顶点着色器中将所有顶点的坐标转换到相机坐标、计算出每个顶点的法向量并将它们传入片段着色器，片段着色器自动插值得到所有像素点的相机坐标和法向量，最后通过光照计算得到每个像素点的颜色<br>
* Blinn冯氏着色<br>
　　实现方式类似于冯氏着色，只是在计算镜面反射时用到的夹角不同<br>

## Display
视频链接：http://www.bilibili.com/video/av9903444/<br>
截图：<br>
![](https://github.com/Chicharito999/ImageCache/raw/master/image/Gouraud.png)<br>
                                    　　　　　　　　　Gouraud Shading<br>
![](https://github.com/Chicharito999/ImageCache/raw/master/image/Phong.png)<br>
                                　  　　　　　　　　  Phong Shading<br>
![](https://github.com/Chicharito999/ImageCache/raw/master/image/BlinnPhong.png) <br>
                                   　　　　　　　　　 BlinnPhong Shading<br>
　　通过Gouraud Shading和Phong Shading的对比可以看出，由于Gouraud Shading是基于顶点计算的光照，其余光照的元素由顶点插值得到，这样插值后的光照相比于基于片元的Phong Shading显得不是很真实，Phong Shading能够获得更为平滑的光照效果。<br>
　　通过Phong Shading和BlinnPhong Shading的对比可以看出，BlinnPhong Shading在明暗交接处的变化是渐变的不像Phong Shading那样变化得那么突然，这是因为在计算镜面反射光时用到的角度不同。　
