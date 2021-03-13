#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>

#include <cmath>

#define M_PI 3.14159265359

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <random>
static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0,1);

 
 class Vector {
 public:
     explicit Vector(double x=0, double y=0, double z=0) {
         coord[0] = x;
         coord[1] = y;
         coord[2] = z;
     };
     double operator[](int i) const { return coord[i];};
     double &operator[](int i) { return coord[i];};
     double sqrNorm() const { 
         return coord[0]*coord[0] +coord[1]*coord[1] +coord[2]*coord[2];}
     Vector get_normalized() {
         double n = sqrt(sqrNorm());
         return Vector(coord[0]/n,coord[1]/n,coord[2]/n );
     }
     Vector get_normalized1() {
         Vector result(*this);
         result.get_normalized();
         return result;
     }
     
     Vector& operator+=(const Vector& a) {
     coord[0]+=a[0];
     coord[1]+=a[1];
     coord[2]+=a[2];

     return *this;
     }

private:
     double coord[3];
 };
 Vector operator+(const Vector& a, const Vector& b ) {
     return Vector(a[0]+b[0],a[1]+b[1],a[2]+b[2]);
 }
 Vector operator-(const Vector& a, const Vector& b ) {
     return Vector(a[0]-b[0],a[1]-b[1],a[2]-b[2]);
 }
 Vector operator*(double a, const Vector& b) {
     return Vector(a*b[0], a*b[1], a*b[2]);
 }

 Vector operator*(const Vector& a, const Vector& b) {
     return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
 }
 
 Vector operator*(const Vector& a, double b) {
     return Vector(a[0]*b, a[1]*b, a[2]*b);
 }

  Vector operator-(const Vector& a) {
     return Vector(-a[0], -a[1], -a[2]);
 }

 Vector operator/(const Vector& a, double b) {
     return Vector(a[0]/b, a[1]/b, a[2]/b);
 }
 double dot(const Vector& a, const Vector& b) {
     return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
 }
 
  Vector cross(const Vector& a, const Vector& b) {
     return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
 }

 double sqr(const float a){
    return a*a;
}

Vector random_cos(const Vector &N){
	double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2*M_PI*r1)*sqrt(1-r2);
    double y = sin(2*M_PI*r1)*sqrt(1-r2);
    double z = sqrt(r2);
    Vector T1;
    if (N[0] < N[1] && N[0] < N[2]){
        T1 = Vector(0, N[2], -N[1]);
    }
    else{
        if (N[1] < N[2] && N[1] < N[0]){
            T1 = Vector(N[2], 0, -N[0]);
        }
        else{
            T1 = Vector(N[1], -N[0], 0);
        }

    }
    T1 = T1.get_normalized();
    Vector T2 = cross(N, T1);
    return z*N + x*T1 + y*T2;

}
class Ray{
     public:
     Ray(const Vector& C, const Vector& u): C(C), u(u) {}
     Vector C, u;

 };
 class Sphere{
 public:
     Sphere(const Vector& O, double R, Vector couleur, bool isMirror = false, bool transp = false, bool creux = false, double n = 1.4): O(O), R(R), albedo(couleur), isMirror(isMirror), transp(transp), creux(creux), n(n) {
     }

     bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
         double a = 1;
         double b = 2*dot(r.u, r.C - O);
         double c = (r.C-O).sqrNorm() - R*R;

         double delta = b*b - 4*a*c;
         if (delta<0) return false;
         double t1=(-b-sqrt(delta))/(2*a);
         double t2=(-b+sqrt(delta))/(2*a);

         if (t2<0) return false;

         if (t1>0) 
            t = t1;
         else 
            t = t2;

         P = r.C + t*r.u; // la direction P est définié par l'origine du rayon et t* la direction du rayon
         N = (P - O).get_normalized();
         return true;
     }

     Vector O;
     double R;
     Vector albedo;
     bool isMirror;
     bool transp, creux;
	 double n;

 };

class Scene{
public:
    Scene() {};
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &n, double &t, int &objectId){
        t = 1E10;
        bool has_inter = false;
        for(int i = 0; i < objects.size(); i++){
            Vector localP, localN;
            double localt;
            if (objects[i].intersect(r, localP, localN, localt) && localt < t){
                has_inter = true;
                t = localt;
                P = localP;
                N = localN;
                if (objects[i].creux) {
                    N = -N;
                }
                albedo = objects[i].albedo;
                mirror = objects[i].isMirror;
                transp = objects[i].transp;
                n = objects[i].n;
                objectId = i;
            }
        }
        return has_inter;
     }

     Vector getColor(const Ray& r, int rebond, bool last_diffu) {
            Vector P, N, albedo;
            double t, n2;
            bool mirror, transp;
            int objectId;
            bool inter=intersect(r, P, N, albedo, mirror, transp, n2, t, objectId);
            Vector color(0, 0, 0);
            if (rebond > 10) return color;

            if (inter) {
                if (objectId == 0){
                    if (rebond == 0 || !last_diffu)
                        return Vector(intensite_lum, intensite_lum, intensite_lum)/(4*M_PI*M_PI*objects[0].R*objects[0].R);
                    else
                        return Vector(0., 0., 0.);
                        }
                if (mirror) {
                    Vector reflectedDir = r.u - 2*dot(r.u, N)*N;
                    Ray reflectedRay(P + 1E-5*N, reflectedDir);
                    return getColor(reflectedRay, rebond + 1, false);
                } else {
                    if (transp) {
                        double n1 = 1;
                        Vector N2 = N;
                        if (dot(r.u, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        Vector Tt = n1 / n2 * (r.u - dot(r.u, N2)*N2);
                        double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));
                        if (rad < 0) {
                            Vector reflectedDir = r.u - 2*dot(r.u, N)*N;
                            Ray reflectedRay(P + 1E-5*N, reflectedDir);
                            return getColor(reflectedRay, rebond + 1, false);

                        }
                        Vector Tn = -sqrt(rad)*N2;
                        Vector refractedDir = Tt + Tn;
                        Ray refractedRay(P - 1E-5*N2, refractedDir);
                        return getColor(refractedRay, rebond + 1, false);
                    }
                    else {
                        /*Vector PL = L - P;
                        double d = sqrt(PL.sqrNorm());
                        Vector shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+1E-5*N, PL/d);
                        bool shadowMirror, shadowTransp;
                        int shadowId;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt, shadowId);
                        if (shadowInter && shadowt < d) {
                        } else {
                            color = I/(4*M_PI*d*d)*albedo/M_PI*std::max(0., dot(N, PL/d));
                        }*/
                        // Eclairage direct
                        Vector PL = lumiere - P;
                        PL = PL.get_normalized();
                        Vector random_w = random_cos(-PL);
                        Vector xprime = random_w*objects[0].R + objects[0].O;
                        Vector Pxprime = xprime - P;
                        double distance = sqrt(Pxprime.sqrNorm());
                        Pxprime = Pxprime/distance;

                        Vector shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+1E-4*N, Pxprime);
                        bool shadowMirror, shadowTransp;
                        int shadowId;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt, shadowId);
                        if (shadowInter && shadowt < distance - 2E-4) {
                            color = Vector(0., 0., 0.);
                        } else {
                            double Jacobien = std::max(1E-8, dot(random_w, -Pxprime))/(distance*distance);
                            double probabilite = std::max(1E-8, dot(-PL, random_w))/(M_PI*objects[0].R*objects[0].R);
                            color = intensite_lum/(4*M_PI*M_PI*objects[0].R*objects[0].R)*albedo/M_PI*std::max(0., dot(N, Pxprime))*Jacobien/probabilite;
                        }

                        // Eclairage indirect
                        Vector wi = random_cos(N);
                        Ray Rayw_inc(P + 1E-5*N, wi);
                        color += albedo*getColor(Rayw_inc, rebond + 1, true);
                    }
                }
            };
            return color;
        }

    std::vector<Sphere> objects;
    Vector lumiere;
    double intensite_lum;
};
 

 //Test de Résolution Monte Carlo 
 /*
void intergreCos() {
    int N=1000;
    double sigma = 0.25;
    double s = 0;
    for (int i = 0; i < N; i++) {
        double u1 = uniform(engine);
        double u2 = uniform(engine);
        double xi = sigma * cos(2*M_PI*u1)*sqrt(-2* log(u2));
        double p = 1/(sigma*sqrt(2*M_PI))*exp(-xi*xi/ (2*sigma*sigma));
        s+=pow(cos(xi), 10) /p/N;
    }
    std::cout << s << std::endl;


}
*/

int main() {

    int W = 300;
    int H = 300;
    double fov = 60 * M_PI/180;
    const int nrays = 200 ;

    Sphere Slum(Vector(15,70,-30),30, Vector(1.,1.,1.));

    Sphere s1(Vector(0,0,-55),15, Vector(1,1,1));
    Sphere s2(Vector(0,-2000-20,0),2000, Vector(1,1,1)); //sol
    Sphere s3(Vector(0,-2000+100,0),2000, Vector(1,1,1)); //plafond
    Sphere s4(Vector(-2000-50,0,0),2000, Vector(0,1,0)); //mur gauche
    Sphere s5(Vector(2000+50,0,0),2000, Vector(0,0,1)); //mur droit
    Sphere s6(Vector(0,0,-2000-100),2000, Vector(0,1,1)); //fond

    Scene s;
    s.objects.push_back(Slum);
    s.objects.push_back(s1);
    s.objects.push_back(s2);
    s.objects.push_back(s3);
    s.objects.push_back(s4);
    s.objects.push_back(s5);
    s.objects.push_back(s6);
    
    s.lumiere = Vector(15,70,-30);
    s.intensite_lum = 5E9;
    
    
    std::vector<unsigned char> image(W*H * 3,0);

    #pragma omp parallel for 
    for (int i = 0; i < H; i++) {
        std::cout << i << '\n';
        for (int j = 0; j < W; j++) {

           
            
            Vector color(0., 0., 0.);
            for (int k=0; k<nrays; k++) {

                double r1 = uniform(engine);
                double r2 =uniform(engine);

                double dx = 0.25*sqrt(-2*log(r1))*cos(2*M_PI*r2);
                double dy = 0.25*sqrt(-2*log(r1))*(2*M_PI*r2);
                
                Vector u(j - W/2 +0.5 +dx, i-H/2+0.5+dy, -W/(2*tan(fov/2)));
                u = u.get_normalized();
                Ray r(Vector(0,0,0),u);

                color += s.getColor(r,0, false);
            }
			color = color/nrays;
                
            image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0],1/2.2)));
            image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1],1/2.2)));
            image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2],1/2.2)));
            

            //std::min(255., std::max(0., intensite_pixel));
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    
    return 0;
};