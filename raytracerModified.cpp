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

float tweak(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

double dot(const Vector3& a, const Vector3& b) {

  return (a.x*b.x + a.y*b.y + a.z*b.z);
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
    double red;
    double green;
    double blue;
    Color(){}
    Color(double x_, double y_, double z_) : red(x_),green(y_),blue(z_)
  {

    }
      Color operator + (const Color& v) const { return Color(red+v.red, green+v.green, blue+v.blue); }
  Color operator - (const Color& v) const { return Color(red-v.red, green-v.green, blue-v.blue); }
    Color operator * (Color v) const { return Color(red*v.red, green*v.green, blue*v.blue); }
  Color operator * (double d) const { return Color(red*d, green*d, blue*d); }

  Color operator / (double d) const { return Color(red/d, green/d, blue/d); }
};

class Box
{
    public:
    float left;
    float right;
    float top;
    float bottom;
    float near;
    float far;
};
class Shapes
{

public:
    virtual Vector3 getNormal(const Vector3& point) {cout<<"Base"<<endl;}
    virtual bool intersect(const Ray &ray,double& t3,double& t4) {cout<<"Base"<<endl;}
    virtual Vector3 getPosition(){}
    virtual Color getColor(){}
    virtual void setColor(double,double,double){}
    virtual void setReflectivity(float _reflectivity){}
    virtual void setTransparency(float){}
    virtual float getReflectivity(){cout<<"Base"<<endl;}
    virtual float getTranparency(){cout<<"Base"<<endl;}
    virtual Box getBoundingBox() { cout<<"Base"<<endl;}
    virtual void setBoundingBox(Box b){cout<<"base"<<endl;}
};

class Sphere : public Shapes
{
    private:
        Vector3 position;
        double r;
        Color color;
        float reflectivity;
        float transparency;
        Box boundingBox;
    public:

    Sphere(double x,double y,double z,int _r)
    {
        position.x = x;
        position.y = y;
        position.z = z;

        r=_r;
        color.red = 1;
        color.green = 1;
        color.blue = 1;

        boundingBox.left = x - r;
        boundingBox.right = x + r;
        boundingBox.top = y + r;
        boundingBox.bottom = y-r;
        boundingBox.near = z - r;
        boundingBox.far = z+r;
    }

    Box getBoundingBox()
    {
        return boundingBox;
    }
    void setBoundingBox(Box b)
    {
        boundingBox = b;
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
    void setColor(double r,double g,double b)
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


    bool intersect(const Ray& ray, double& intersection1,double &intersection2)
    {

        //This intersection is based on the formula
        double radius2 = r*r;
         Vector3 Lvector = position - ray.po;
        double temp1 = dot(Lvector, ray.d);

        if ( temp1 < 0.0 ) return false;
        double temp2 =   dot(Lvector,Lvector)- (temp1*temp1) ;


        if ( temp2 > radius2) return false;


        float quad = sqrt( radius2 - temp2 );

        //solve for intersection points
        intersection1 = temp1 - quad;
        intersection2 = temp1 + quad;

        return true;
    }
};

class Rectangle : public Shapes
{
    private:

        Color color;
        Vector3 position;
        float reflectivity;
        float transparency;
        Box boundingBox;
    public:

    Rectangle(double left,double right,double top,int bottom,double near,double far)
    {


        position.x = (left + right) / 2;
        position.y = (top + bottom) / 2;
        position.z= (near + far) / 2;
        color.red = 1;
        color.green = 1;
        color.blue = 1;

        boundingBox.left = left;
        boundingBox.right = right;
        boundingBox.top = top;
        boundingBox.bottom = bottom;
        boundingBox.near = near;
        boundingBox.far = far;
    }

    Box getBoundingBox()
    {
        return boundingBox;
    }
    void setBoundingBox(Box b)
    {
        boundingBox = b;
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
    void setColor(double r,double g,double b)
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
        return position;
    }


    bool intersect(const Ray& ray, double& intersection1,double &intersection2)
    {

        float tmin = (boundingBox.left - ray.po.x) / ray.d.x;
        float tmax = (boundingBox.right - ray.po.x) / ray.d.x;

        if (tmin > tmax) swap(tmin, tmax);

        float tymin = (boundingBox.bottom - ray.po.y) / ray.d.y;
        float tymax = (boundingBox.top - ray.po.y) / ray.d.y;

        if (tymin > tymax) swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax))
        return false;

        if (tymin > tmin)
        tmin = tymin;

        if (tymax < tmax)
        tmax = tymax;

        float tzmin = (boundingBox.near - ray.po.z) / ray.d.z;
        float tzmax = (boundingBox.far - ray.po.z) / ray.d.z;

        if (tzmin > tzmax) swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax))
        return false;

        if (tzmin > tmin)
        tmin = tzmin;

        if (tzmax < tmax)
        tmax = tzmax;

        return true;
    }
};

void clamp(int x, int y,int z) {
    /*Clamping function to clamp values from 0 to 255*/
  x = (x > 255) ? 255 : (x < 0) ? 0 : x;
  y = (y > 255) ? 255 : (y < 0) ? 0 : y;
  z = (z > 255) ? 255 : (z < 0) ? 0 : z;
}

Color trace(Ray& ray, vector<Shapes*>& shapeList, Vector3& lightPosition, int depth,int xpixel,int ypixel)
{
    vector<Shapes*>::iterator it;
    vector<Shapes*>::iterator it2;
    Shapes* intersectingObject =NULL;
    float minimumDistance = 99999999;
    for(it=shapeList.begin();it!=shapeList.end();it++)
    {

        //Acceleration structure - Bounding volume
        Box shapeVolumeBox = (*it)->getBoundingBox();
        Vector3 rayDir = ray.d;
        if(xpixel < (shapeVolumeBox.left - shapeVolumeBox.far))
            continue;
        if(ypixel<(shapeVolumeBox.top - shapeVolumeBox.far))
            continue;


        double intersection1 = 999999;
        double intersection2 = 999999;
        bool isIntersecting = (*it)->intersect(ray,intersection1,intersection2);

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
    }


    if(intersectingObject!=NULL)
    {
        Color newColor(0,0,0);
        Vector3 pointOfIntersection = ray.po + (ray.d*minimumDistance);

        Vector3 normalOfIntersection = intersectingObject->getNormal(pointOfIntersection).reverseVector();

        Vector3 normalizedIntersection = normalOfIntersection.normalize();




        float deviate = 1e-2;

        bool isInside = false;

        if (dot(ray.d,normalizedIntersection) > 0)
        {
            normalizedIntersection = normalizedIntersection.reverseVector();
            isInside = true;
        }



        if(depth>0 && (intersectingObject->getTranparency()>0 || intersectingObject->getReflectivity()>0) )
        {
            float facing = -dot(ray.d,normalizedIntersection);

            double reflectiveWaves = tweak(pow(1 - facing, 3), 1, 0.1);

            Vector3 reflectedRayDir = ray.d - normalizedIntersection * 2 * dot(ray.d,normalizedIntersection);

            Vector3 reflectedRayDirNormalize = reflectedRayDir.normalize();

            Vector3 reflectedRayPoint = pointOfIntersection + normalizedIntersection * deviate;

            Ray reflectionRay(reflectedRayPoint,reflectedRayDirNormalize);

            Color reflectionColor = trace(reflectionRay,shapeList,lightPosition,depth-1,xpixel,ypixel);

            Color refraction(0,0,0);

            if(intersectingObject->getTranparency()>0)
            {
                float InsideOrOutside = 1.1;
                float deviateRayC = (isInside) ? InsideOrOutside : 1 / InsideOrOutside;
                float refrectiveInd = -dot(normalizedIntersection,ray.d);
                float k = 1 - deviateRayC * deviateRayC * (1 - refrectiveInd * refrectiveInd);
                Vector3 refrdir = ray.d * deviateRayC + normalizedIntersection * (deviateRayC * refrectiveInd - sqrt(k));
                Vector3 refNormalized = refrdir.normalize();
                Vector3 point = pointOfIntersection - normalizedIntersection * deviate;
                Ray insideRay(point,refNormalized);
                refraction = trace(insideRay, shapeList,lightPosition, depth - 1,xpixel,ypixel);
            }

            newColor = (reflectionColor * reflectiveWaves + refraction * (1 - reflectiveWaves) * intersectingObject->getTranparency() ) * intersectingObject->getColor();

        }
        else
        {

            float shadow = 1;

            Vector3 directionOfShadowRay = lightPosition - pointOfIntersection ;

            Vector3 normalizedDirectionShadow = directionOfShadowRay.normalize();



            for(it2=shapeList.begin();it2!=shapeList.end();it2++)
            {

                if(intersectingObject!=(*it2))
                {

                    double temp1,temp2;
                    Vector3 diffusePoint = pointOfIntersection + normalizedIntersection * deviate;
                    Vector3 reversing = normalizedDirectionShadow.reverseVector();
                    Ray secondRay(diffusePoint,normalizedDirectionShadow);
                    bool isIntersecting = (*it2)->intersect(secondRay,temp1,temp2);
                    if(isIntersecting)
                    {


                        shadow = 0;
                        break;
                    }
                }
            }

            Color lightColor(1,1,1);
            newColor = intersectingObject->getColor() * shadow * std::max(double(0), dot(normalizedDirectionShadow,normalizedIntersection)) * lightColor;


        }

        return newColor;
        }
        else
        {
            Color background(1,1,1);
            return background;
        }


}

void render(vector<Shapes*> shapeList)
{

    //Initializing variables
    const int Height = 768;
    const int Width = 1024;

    int depth = 1;


    Vector3 eye(0,0,0);
    Vector3 lightPosition(0,20,-30);

    //Setting image (bitmap) properties
    unsigned char *img = NULL;
    img = (unsigned char *)malloc(3*Width*Height);
    memset(img,0,3*Width*Height);


    //calculating angle from field of view = 30
    float angle = tan(PI * 0.5 * 30 / 180.);
    for(int y=0;y<Height;y++)
    {
        for(int x=0;x<Width;x++)
        {

            //Calculating direction of the ray based on the x and y
            float middleX = (2 * ((x + 0.5) * (1/float(Width))) - 1) * angle * (Width / float(Height));
            float middleY = (1 - 2 * ((y + 0.5) * (1 / float(Height)))) * angle;


            Vector3 direction(middleX,middleY,-1);
            Vector3 normalizedDirection = direction.normalize();


            Ray ray(eye,normalizedDirection);

            Color thisRayColor = trace(ray,shapeList,lightPosition,depth,x,y);


            //Writing color to file
            int i=x;
            int j=(Height-1)-y;

            //Converting from float to int coordinates
            int r = max(0, min(255, (int)floor(thisRayColor.red * 256.0)));
            int g = max(0, min(255, (int)floor(thisRayColor.green * 256.0)));
            int b = max(0, min(255, (int)floor(thisRayColor.blue * 256.0)));


            img[(i+j*Width)*3+2] = (unsigned char)(r);
            img[(i+j*Width)*3+1] = (unsigned char)(g);
            img[(i+j*Width)*3+0] = (unsigned char)(b);



        }
    }
    generateBMP(Width,Height,img);

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

    f = fopen("img.bmp","wb");
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

    Shapes* sphere = new Sphere(0,-10000,-100,10000);
    sphere->setColor(0.5,0.5,0.5);
    sphere->setTransparency(0);
    sphere->setReflectivity(2);

    Shapes* sphere2 = new Sphere( 0, 0, -20, 4);
    sphere2->setTransparency(0.5);
    sphere2->setReflectivity(1);
    sphere2->setColor(0.1,1,0.1);


    Shapes* sphere3 = new Sphere( 3, 1, -10, 2);
    sphere3->setTransparency(0);
    sphere3->setReflectivity(1);
    sphere3->setColor(0.8,0.4,0.3);

    Shapes* sphere4 = new Sphere( -2, 1, -10, 1);
    sphere4->setTransparency(0);
    sphere4->setReflectivity(1);
    sphere4->setColor(0.2,0.3,1);


    Rectangle* rect1 = new Rectangle(0,3,0,3,-30,-40);
    rect1->setTransparency(0);
    rect1->setReflectivity(0);
    rect1->setColor(1,1,1);

    shapeList.push_back(sphere);
    shapeList.push_back(sphere2);
    shapeList.push_back(sphere3);
    shapeList.push_back(sphere4);
    shapeList.push_back(rect1);



    render(shapeList);


}
