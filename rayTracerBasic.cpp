#include <fstream>
#include <cmath>
#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>
using namespace std;
#define PI 3.14159265358979323846

void generateBMP(int,int,unsigned char*);
class Vector3 {
    public:
  double x,y,z;
  Vector3(){}
  Vector3(double x_, double y_, double z_) : x(x_),y(y_),z(z_)
  {

    }
  Vector3 operator + (const Vector3& v) const { return Vector3(x+v.x, y+v.y, z+v.z); }
  Vector3 operator - (const Vector3& v) const { return Vector3(x-v.x, y-v.y, z-v.z); }
  Vector3 operator * (double d) const { return Vector3(x*d, y*d, z*d); }
  Vector3 operator / (double d) const { return Vector3(x/d, y/d, z/d); }
  Vector3 normalize() const {
    double mg = sqrt(x*x + y*y + z*z);
    return Vector3(x/mg,y/mg,z/mg);
  }
  double Magnitude() const
  {
      return sqrt(x*x + y*y + z*z);
  }

  Vector3 reverseVector() const
  {
      return Vector3(0.0-x, 0.0-y, 0.0-z);
  }
};

inline double dot(const Vector3& a, const Vector3& b) {

  return (a.x*b.x + a.y*b.y + a.z*b.z);
}
bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
float discr = b * b - 4 * a * c;
if (discr < 0) return false;
else if (discr == 0) x0 = x1 = - 0.5 * b / a;
else {
float q = (b > 0) ?
-0.5 * (b + sqrt(discr)) :
-0.5 * (b - sqrt(discr));
x0 = q / a;
x1 = c / q;
}
if (x0 > x1) std::swap(x0, x1);

return true;
}

class Light
{
    public:
        Vector3 position;
        Light(Vector3& pos) : position(pos)
        {

        }
};

class Ray
{
    public:
        Vector3 po;
        Vector3 d;
        Ray(Vector3& point, Vector3& _direction) : po(point), d(_direction)
        {

        }
};

class Material
{

};

struct Color
{
    int red;
    int green;
    int blue;
};

class Shapes
{

public:
    virtual Vector3 getNormal(const Vector3& point) {cout<<"Base"<<endl;}
    virtual bool intersect(const Ray &ray,double& t3,double& t4) {cout<<"Base"<<endl;}
    virtual Vector3 getPosition(){}
    virtual Color getColor(){}
    virtual void setColor(int,int,int){}
    virtual void setReflectivity(float _reflectivity){}
    virtual void setTransparency(float){}
    virtual float getReflectivity(){}
    virtual float getTranparency(float){}
};

class Sphere : public Shapes
{
    private:
        Vector3 position;
        double r;
        Color color;
        float reflectivity;
        float transparency;
    public:

    Sphere(double x,double y,double z,int _r)
    {
        position.x = x;
        position.y = y;
        position.z = z;

        r=_r;
        color.red = 255;
        color.green = 255;
        color.blue = 255;
    }

    void setReflectivity(float _reflectivity)
    {
        reflectivity = _reflectivity;
    }

    float getReflectivity()
    {

        return reflectivity;
    }

    void setTransparency(float trans)
    {
        transparency = trans;
    }

    float getTranparency()
    {

        return transparency;
    }
    void setColor(int r,int g,int b)
    {
        color.red = r;
        color.green = g;
        color.blue = b;
    }

    Color getColor()
    {
        return color;
    }

    Vector3 getPosition()
    {
        return position;
    }
    Vector3 getNormal(const Vector3& point)
    {
        return (point - position);
    }


    bool intersect(const Ray& ray, double& t0,double &t1)
    {
        double radius2 = r*r;
         Vector3 L = position - ray.po;
        double tc = dot(L, ray.d);

        if ( tc < 0.0 ) return false;
        double d2 = (tc*tc) - dot(L,L);


        if ( d2 > radius2) return false;

        //solve for t1c
        float t1c = sqrt( radius2 - d2 );

        //solve for intersection points
        t1 = tc - t1c;
        t0 = tc + t1c;

        return true;
    }
};




void clamp(int x, int y,int z) {
    /*Clamping function to clamp values from 0 to 255*/
  x = (x > 255) ? 255 : (x < 0) ? 0 : x;
  y = (y > 255) ? 255 : (y < 0) ? 0 : y;
  z = (z > 255) ? 255 : (z < 0) ? 0 : z;
}


void render(vector<Shapes*> shapeList)
{
    const int Height = 500;
    const int Width = 500;
    vector<Light> lights;



    Vector3 eye(Width/2,Height/2,-40);
    Vector3 lightPosition(Width/2,Height/2,-40);

     Light point(lightPosition);
    lights.push_back(point);


    unsigned char *img = NULL;
    img = (unsigned char *)malloc(3*Width*Height);
    memset(img,0,3*Width*Height);

    vector<Shapes*>::iterator it;
    vector<Shapes*>::iterator it2;
    for(int y=0;y<Height;y++)
    {
        for(int x=0;x<Width;x++)
        {
            double middleX = x+0.5;
            double middleY = y+0.5;

            Vector3 point(x,y,0);
            Vector3 direction(middleX-eye.x,middleY-eye.y,0-eye.z);

            Vector3 normalizedDirection = direction.normalize();


            Ray ray(point,normalizedDirection);



            for(it=shapeList.begin();it!=shapeList.end();it++)
            {

                Shapes* intersectingObject =NULL;
                double intersection1 = 999999;
                double intersection2 = 999999;
                bool isIntersecting = (*it)->intersect(ray,intersection1,intersection2);
                float minimumDistance = 9999999;
                if (isIntersecting)
                {



                    if(intersection1<0)
                        intersection1 = intersection2;
        if(intersection1<minimumDistance)
                    {
                        minimumDistance = intersection1;
                        intersectingObject = (*it);

                    }

                }
                bool shadow = false;

                if(intersectingObject!=NULL)
                {

                    Vector3 pointOfIntersection = ray.po + (ray.d*minimumDistance);

                    Vector3 normalOfIntersection = intersectingObject->getNormal(pointOfIntersection).reverseVector();

                    Vector3 directionOfShadowRay = lightPosition - pointOfIntersection ;

                    Vector3 normalizedDirectionShadow = directionOfShadowRay.normalize();
                    Vector3 pointOfIntersectionRev = pointOfIntersection.reverseVector();

                     Vector3 normalizedIntersection = normalOfIntersection.normalize();

                    Ray secondRay(pointOfIntersection,normalizedDirectionShadow);


                    double temp1,temp2;
                    for(it2=shapeList.begin();it2!=shapeList.end();it2++)
                    {
                        if(it!=it2)
                        {
                            bool isIntersecting = (*it2)->intersect(secondRay,temp1,temp2);
                            if(isIntersecting)
                            {

                                break;
                            }
                        }
                    }


                    int i=x;
                    int j=(Height-1)-y;
                    Color newColor = (*it)->getColor();
                    double alpha = max(0.0,abs(dot(normalizedIntersection,normalizedDirectionShadow)));
                   // double colorScale = alpha * 1 + (1 - alpha) * 0;
                    int r =  (newColor.red * alpha);
                    int g = (newColor.green * alpha);
                    int b = (newColor.blue * alpha);
                    clamp(r,g,b);
                    if(alpha==0 && shadow)
                    {


                    }
                    else if(alpha==0)
                    {
                        continue;

                    }




                    if(shadow)
                    {
                        r = r * 0.9;
                        g = g * 0.9;
                       b = b * 0.9;


                        img[(i+j*Width)*3+2] = (unsigned char)(255);
                        img[(i+j*Width)*3+1] = (unsigned char)(255);
                        img[(i+j*Width)*3+0] = (unsigned char)(255);

                    }
                    else
                    {

                        img[(i+j*Width)*3+2] = (unsigned char)(r);
                        img[(i+j*Width)*3+1] = (unsigned char)(g);
                        img[(i+j*Width)*3+0] = (unsigned char)(b);
                    }


                }
                else
                {

                }

            }

        }
    }
    generateBMP(Width,Height,img);

}


float mix(const float &a, const float &b, const float &mix)
{
return b * mix + a * (1 - mix);
}

void generateBMP(int w,int h,unsigned char *img)
{

    FILE *f;

    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int


    unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
    unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
    unsigned char bmppad[3] = {0,0,0};

    bmpfileheader[ 2] = (unsigned char)(filesize    );
    bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
    bmpfileheader[ 4] = (unsigned char)(filesize>>16);
    bmpfileheader[ 5] = (unsigned char)(filesize>>24);

    bmpinfoheader[ 4] = (unsigned char)(       w    );
    bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
    bmpinfoheader[ 6] = (unsigned char)(       w>>16);
    bmpinfoheader[ 7] = (unsigned char)(       w>>24);
    bmpinfoheader[ 8] = (unsigned char)(       h    );
    bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
    bmpinfoheader[10] = (unsigned char)(       h>>16);
    bmpinfoheader[11] = (unsigned char)(       h>>24);

    f = fopen("imgBasic.bmp","wb");
    fwrite(bmpfileheader,1,14,f);
    fwrite(bmpinfoheader,1,40,f);
    for(int i=0; i<h; i++)
    {
        fwrite(img+(w*(h-i-1)*3),3,w,f);
        fwrite(bmppad,1,(4-(w*3)%4)%4,f);
    }

    free(img);
    fclose(f);
}

int main()
{



    vector<Shapes*> shapeList;
    Shapes* sphere = new Sphere(125,125,100,50);
    sphere->setColor(255,125,0);
    sphere->setTransparency(1);
    sphere->setReflectivity(0.5);
    Shapes* sphere2 = new Sphere(375,375,80,50);
    sphere2->setTransparency(1);
    sphere2->setReflectivity(0);
    sphere2->setColor(0,255,255);
        shapeList.push_back(sphere);
    shapeList.push_back(sphere2);
    render(shapeList);


}
