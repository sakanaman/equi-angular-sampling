#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#define EQUI_ANG

struct Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    double length() const {return sqrt(this->dot(*this));}
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};
struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()
struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    {                     // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
Sphere spheres[] = {
    //Scene: radius, position, emission, color, material
    // Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   //Left
    // Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), //Rght
    // Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         //Back
    // Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               //Frnt
    // Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         //Botm
    // Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Top
    Sphere(16.5, Vec(27, 30.5, 78), Vec(), Vec(1, 1, 1) * .999, SPEC),        //Mirr
    Sphere(16.5, Vec(73, 30.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR)        //Glas
    //Sphere(8, Vec(50, 52, 81.6), Vec(12, 12, 12), Vec(), DIFF)     //Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id))
        return Vec();                // if miss, return black
    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl
    if (++depth > 5)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.e; //R.R.
    if (obj.refl == DIFF)
    { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    else if (obj.refl == SPEC) // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.e + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                           radiance(reflRay, depth, Xi) * RP
                                                       : radiance(Ray(x, tdir), depth, Xi) * TP)
                                    : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}

//param
constexpr double kEPS = 1e-6;
constexpr double Max_dist = 200.0;
constexpr double sigma_s = 0.0015;
constexpr double sigma_a = 0.0;
constexpr double sigma_t = sigma_a + sigma_s;
Vec light_power = Vec(500.0, 500.0, 500.0);
Vec point_light = Vec(50, 52, 91.6);
Vec Albedo(1.0,1.0,1.0);

double equianglar_sample(const Ray& ray, const Vec& light_pos, const double u,
                         const double tNear,const double tFar, double* pdf)
{
    const double delta = ray.d.dot(light_pos - ray.o);
    const double D = ((ray.o + ray.d * delta) - light_pos).length();

    double a = tNear - delta;
    double b = tFar - delta;

    //if(a > b) printf("それはヤバイで\n");
    //printf("%f\n", u);
    double t;
    if(D > kEPS)
    {
        double thetaA = atan(a/D);
        double thetaB = atan(b/D);

        t = D * tan((1 - u)*thetaA + u*thetaB);
        *pdf = D/abs(thetaA - thetaB)/(D*D + t*t);
    }
    else
    {
        //Dがほぼ0の時
        t = a*b / (b + (a - b)*u);
        *pdf = a*b / (b - a) / (t*t);
    }
    
    //*contribution = (*contribution) + (light_power * (1.0/(D*D + t*t)) * sigma_s * exp(-sigma_t*(t + delta + sqrt(D*D + t*t))) * (1.0/(*pdf)));
    //printf("%f\n",delta + t); -->値は正常
    return delta + t;
}

double distance_sample(const Ray& ray, const Vec& light_pos, const double u,
                        const double tNear,const double tFar, double* pdf)
{
    double t = tNear - (1.0/sigma_t) * log(1.0 - u * (1.0 - exp(-sigma_t * (tFar - tNear))));
    *pdf = sigma_t/(exp(sigma_t*(t - tNear)) - exp(sigma_t*(t - tFar)));

    return t;
}

Vec radiance_density(const Ray &r, int depth, unsigned short *Xi)
{
    
}

Vec radiance_equiangular(const Ray &r, int depth, unsigned short *Xi)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object

    Vec contribute = Vec(0.0, 0.0, 0.0);

    double tNear, tFar;
    if (!intersect(r, t, id)) //何にも当たらなかった
    {
        tNear = 0.0;
        tFar = Max_dist;
    }
    else
    {
        tNear = 0.0;
        tFar = t;
        // printf("in"); -->入ってる

        //点光源と衝突点でNEE
        const Sphere &obj = spheres[id];
        Vec hitPos = r.o + r.d * t;
        Vec h2l = (point_light-hitPos).norm();
        Vec n = (hitPos - obj.p).norm();
        double l = (point_light - hitPos).length();
        double costerm = h2l.dot(n);

        double shadow_t;
        int shadow_id = 0;
        Ray shadow_r(hitPos, h2l);
        bool is_shadowhit = intersect(shadow_r, shadow_t, shadow_id);
        if(is_shadowhit == false)
        {
            contribute = contribute + (Albedo * (1.0/M_PI)).mult(light_power * costerm * (1.0/(l*l))) * exp(-sigma_t*(l + t));
        }
        // contribute = Vec(1.0, 0, 0);
        // return contribute;-->ちゃんと衝突してる。
    }

    /// Equi-Angular Sampleの開始だ!!

    double pdf, scatter_dist;
    scatter_dist = equianglar_sample(r, point_light, erand48(Xi), 
                                    tNear, tFar, &pdf);
    //printf("%f\n",pdf);

    Vec scatter_point = r.o + r.d * scatter_dist;
    Vec scatter_dir = (point_light - scatter_point).norm();

    double shadow_t;
    int shadow_id = 0;
    Ray shadow_r(scatter_point, scatter_dir);
    bool is_shadowhit = intersect(shadow_r, shadow_t, shadow_id);
    if(!is_shadowhit || (is_shadowhit && shadow_t > (point_light - scatter_point).length()))
    {
        // printf("in"); --> 入ってる
        double l2s = (point_light - scatter_point).length();
        //printf("%f\n", l2s);

        Vec plus = light_power * (1.0/(l2s * l2s)) * sigma_s * exp(-sigma_t*(scatter_dist + l2s)) * (1.0/pdf);
        //printf("%f\n", -sigma_t*(scatter_dist + l2s));
        //printf("(%f, %f, %f)\n", plus.x, plus.y, plus.z);

        contribute = contribute + light_power * (1.0/(l2s * l2s)) * sigma_s * exp(-sigma_t*(scatter_dist + l2s)) * (1.0/pdf);
    }
    //printf("(%f, %f, %f)\n",contribute.x,contribute.y,contribute.z);
    return contribute;
}

int main(int argc, char *argv[])
{
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
    Ray cam(Vec(50, 52, 245.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)       // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance_equiangular(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}