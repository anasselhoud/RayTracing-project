#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.14159265359

#include <algorithm>
#include <complex>
#include <random>
#include <string>
#include <iostream>
#include <stdio.h>
#include <list>


static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);


class Vector{
    public:
        explicit Vector(double x=0, double y=0, double z=0) {
            coords[0] = x;
            coords[1] = y;
            coords[2] = z;
        };
        double operator[](int i) const {return coords[i]; };
        double &operator[](int i) {return coords[i]; };
        double sqrNorm() {
            return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
        }
        Vector get_normalized(){
            double n = sqrt(sqrNorm());
            return Vector(coords[0]/n, coords[1]/n, coords[2]/n);
        }
        Vector& operator+=(const Vector& a){
            coords[0] += a[0];
            coords[1] += a[1];
            coords[2] += a[2];
            return *this;
        }

    private:
        double coords[3];
    };

Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b){
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector& a){
    return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(double a, const Vector& b){
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, double b){
    return Vector(b*a[0], b*a[1], b*a[2]);
}
Vector operator*(const Vector& a, const Vector& b){
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, double b){
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}
Vector cross(const Vector& a, const Vector& b){
    return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}
double dot(const Vector& a, const Vector& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double sqr(const float a){
    return a*a;
}

Vector random_cos(const Vector &N){
    double u1 = uniform(engine);
    double u2 = uniform(engine);
    double x = cos(2*M_PI*u1)*sqrt(1-u2);
    double y = sin(2*M_PI*u1)*sqrt(1-u2);
    double z = sqrt(u2);
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


class Ray
{
    public:
        explicit Ray(const Vector& C, const Vector& u): C(C), u(u) {
        };
        Vector C, u;
};

class Sphere
{
    public:
        explicit Sphere(const Vector& O, double R, const Vector& albedo, bool isMirror = false, bool transp = false, double n = 1.4, bool isCreux = false): O(O), R(R), albedo(albedo), isMirror(isMirror), transp(transp), n(n), isCreux(isCreux){
        };
        bool intersect(const Ray& r, Vector& P, Vector& N, double& t){
            double a = 1;
            double b = 2*dot(r.u, r.C-O);
            double c = (r.C-O).sqrNorm() - R*R;
            double delta = b*b - 4*a*c;
            if (delta < 0) return false;

            double sqDelta = sqrt(delta);
            double t2 = (-b + sqDelta)/ (2*a);
            if (t2 < 0) return false;

            double t1 = (-b - sqDelta) / (2*a);
            if (t1 > 0)
                t = t1;
            else
                t = t2;

            P = r.C + t*r.u;
            N = (P - O).get_normalized();
            return true;
        };
        bool isMirror, transp, isCreux;
        Vector albedo;
        double n;
    private:
        Vector O;
        double R;
};

class Scene{
public:
    Scene() {};
    bool intersect(const Ray& r, Vector& P, Vector& N, Vector &albedo, bool &mirror, bool &transp, double &n, double &t){
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
                if (objects[i].isCreux) {
                    N = -N;
                }
                albedo = objects[i].albedo;
                mirror = objects[i].isMirror;
                transp = objects[i].transp;
                n = objects[i].n;
            }
        }
        return has_inter;
     }

     Vector getColor(const Ray& r, int rebond) {
            Vector P, N;
            double t, n2;
            bool mirror, transp;
            Vector albedo;
            int object_id;

            bool has_inter=intersect(r, P, N, albedo, mirror, transp, n2, t);

            Vector color(0, 0, 0);
            if (rebond > 10) return color;

            if (has_inter) {
                if (mirror) {
                    Vector reflect_direct = r.u - 2*dot(r.u, N)*N;
                    Ray reflect_ray(P + .0001*N, reflect_direct);
                    return getColor(reflect_ray, rebond + 1);
                } else {
                    if (transp) {

                        double n1 = 1;
                        Vector N2 = N;
                        if (dot(r.u, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        Vector tang = n1 / n2 * (r.u - dot(r.u, N2)*N2);
                        double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(r.u, N2)));
                        if (rad < 0) {
                            Vector reflect_direct = r.u - 2*dot(r.u, N)*N;
                            Ray reflect_ray(P + 0.0001*N, reflect_direct);
                            return getColor(reflect_ray, rebond + 1);

                        }
                        Vector nor = -sqrt(rad)*N2;
                        Vector reflect_direct = tang + nor;
                        Ray reflect_ray(P - 0.0001*N2, reflect_direct);
                        return getColor(reflect_ray, rebond + 1);
                    }
                    else {
                         Vector PL = lumiere - P;
                        double d = sqrt(PL.sqrNorm());
                        Vector shadowP, shadowN, shadowAlbedo;
                        double shadowt, shadown2;
                        Ray shadowRay(P+0.0001*N, PL/d);
                        bool shadowMirror, shadowTransp;
                        bool shadowInter = intersect(shadowRay, shadowP, shadowN, shadowAlbedo, shadowMirror, shadowTransp, shadown2, shadowt);
                        if (shadowInter && shadowt < d) {
                        } else {
                            color = intensite_lum/(4*M_PI*d*d)*albedo/M_PI*std::max(0., dot(N, PL/d));
                        }
                        // Eclairage indirect
                        Vector direction_random = random_cos(N);
                        Ray rayon_random(P + 1E-3*N, direction_random);
                        color += albedo*getColor(rayon_random, rebond + 1);
                    }
                }
            };
            return color;
        }

    std::vector<Sphere> objects;
    Vector lumiere;
    double intensite_lum;
};

int main() {
    int W = 500;
    int H = 500;
    int nrays = 50;
    double fov = 60 * M_PI / 180;

    Scene s;
    Vector C(0, 0, 55);
   

     Sphere S1(Vector(0,0,0),10, Vector(1,1,0));
    Sphere S1plus(Vector(-10, -5, 15), 2, Vector(1., 0., 0.), true, true);
    Sphere S2(Vector(-10, -5, 30), 2, Vector(1., 0., 0.), false, false);
    Sphere S3(Vector(10, -7, 10), 3, Vector(1., 0., 0.), true, false);

  Sphere s2(Vector(0,-1000,0),990, Vector(1,1,1)); //sol
    Sphere s3(Vector(0,-2000+100,0),2000, Vector(1,1,1)); //plafonds
    Sphere s4(Vector(-2000-50,0,0),2000, Vector(0,1,0)); //mur gauche
    Sphere s5(Vector(2000+50,0,0),2000, Vector(0,0,1)); //mur droit
    Sphere s6(Vector(0,0,-2000-100),2000, Vector(0,1,1)); //fond
    s.objects.push_back(S1);
    s.objects.push_back(S2);
    s.objects.push_back(S3);
    s.objects.push_back(S1plus);
    s.objects.push_back(s2);
    s.objects.push_back(S3);
    s.objects.push_back(s4);
    s.objects.push_back(s5);
    s.objects.push_back(s6);
    
    s.intensite_lum = 5E9;
    s.lumiere = Vector(-10, 20, 40);

    std::vector<unsigned char> image(W*H * 3, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        std::cout << i << '\n';
        for (int j = 0; j < W; j++) {

            Vector color(0,0,0);
            for (int k=0; k<nrays; k++){
                double r1 = uniform(engine);
                double r2 =uniform(engine);

                double dx = 0.25*sqrt(-2*log(r1))*cos(2*M_PI*r2);
                double dy = 0.25*sqrt(-2*log(r1))*(2*M_PI*r2);

                Vector u(j - W/2 +0.5 +dx, i-H/2+0.5+dy, -W/(2*tan(fov/2)));
                u = u.get_normalized();
                Ray r(C, u);
                color += s.getColor(r, 0);
            }
            color = color / nrays;
            image[((H - i - 1)* W + j) * 3 + 0] = std::min(255., std::pow(color[0], 1/2.2));
            image[((H - i - 1)* W + j) * 3 + 1] = std::min(255., std::pow(color[1], 1/2.2));
            image[((H - i - 1)* W + j) * 3 + 2] = std::min(255., std::pow(color[2], 1/2.2));

             //std::min(255., std::max(0., intensite_pixel));
            }
        }
    cout << "Saving the image" << endl;
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    cout << "Finished saving. Check it!" << endl;

    return 0;
}
