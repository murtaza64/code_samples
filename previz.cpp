#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <float.h>
#include "SETTINGS.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"

using namespace std;

#define PLASTIC 0
#define MIRROR 1
#define GLASS 2

#define SS_GRID_SIZE 0.005

typedef struct cylinderInfo {
  Matrix3d R;
  VEC3 translation;
} cylinderInfo;

typedef struct cylinder {
  VEC3 leftVertex;
  VEC3 rightVertex;
  float radius;
  VEC3 color;
  cylinderInfo info;
} cylinder;

typedef struct {
  VEC3 origin;
  VEC3 dir;
  int in_glass;
} ray;

typedef struct {
  VEC3 center;
  VEC3 color;
  float radius;
  int material;
} sphere;

typedef struct {
  VEC3 pos;
  VEC3 color;
} light;

typedef struct {
  VEC3 a;
  VEC3 b;
  VEC3 c;
  int material;
  VEC3 color;
} triangle;

typedef enum {NONE, SPHERE, CYLINDER, TRIANGLE} intersectType;

// Stick-man classes
DisplaySkeleton displayer;    
Skeleton* skeleton;
Motion* motion;

int windowWidth = 640;
int windowHeight = 480;

VEC3 eye(-6, 0.5, 1);
VEC3 lookingAt(5, 0.5, 1);
VEC3 up(0,1,0);

// scene geometry
vector<VEC3> sphereCenters;
vector<float> sphereRadii;
vector<VEC3> sphereColors;
vector<cylinder> cylinders;
vector<sphere> spheres;
vector<triangle> tris;
light lights[2];
int nlights;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void calculateCylinderInfo(){
  for (int i = 0; i < cylinders.size(); i++) {
    cylinder cyl = cylinders[i];
    VEC3 translation = (cyl.leftVertex + cyl.rightVertex) / 2;
    VEC3 z = {0,0,1};
    VEC3 rotAxisDir = -(cyl.leftVertex - cyl.rightVertex).cross(z);
    rotAxisDir = rotAxisDir/rotAxisDir.norm();
    // std::cout << "rotAxisDir:\n" << rotAxisDir << endl;
    float a = rotAxisDir[0];
    float b = rotAxisDir[1];
    float c = rotAxisDir[2];
    float d = sqrt(b*b + c*c);
    Matrix3d Rx;
    Matrix3d Rx_inv;
    if (d == 0) {
      Rx << 1, 0, 0, 0, 1, 0, 0, 0, 1;
      Rx_inv << 1, 0, 0, 0, 1, 0, 0, 0, 1; 
    }
    else {
      Rx << 1, 0, 0, 0, c/d, -b/d, 0, b/d, c/d;
      Rx_inv << 1, 0, 0, 0, c/d, b/d, 0, -b/d, c/d;
    }
    Matrix3d Ry;
    Ry << d, 0, -a, 0, 1, 0, a, 0, d;
    Matrix3d Ry_inv;
    Ry_inv << d, 0, a, 0, 1, 0, -a, 0, d;
    float acos_arg = (cyl.leftVertex - cyl.rightVertex).dot(z)/(cyl.leftVertex - cyl.rightVertex).norm();
    float theta = acos((cyl.leftVertex - cyl.rightVertex).dot(z)/(cyl.leftVertex - cyl.rightVertex).norm());
    // printf("theta = %f\n", theta);
    // cout << "acos arg = " << acos_arg << endl;
    Matrix3d Rz;
    Rz << cos(theta), -sin(theta), 0, sin(theta), cos(theta), 0, 0, 0, 1;
    Matrix3d R = Rx_inv * Ry_inv * Rz * Ry * Rx;
    cylinders[i].info.R = R;
    cylinders[i].info.translation = translation;
  }
}

// bool raySphereIntersect(const VEC3& center, 
//                         const float radius, 
//                         const VEC3& rayPos, 
//                         const VEC3& rayDir,
//                         float& t)
// {
//   const VEC3 op = center - rayPos;
//   const float eps = 1e-8;
//   const float b = op.dot(rayDir);
//   float det = b * b - op.dot(op) + radius * radius;

//   // determinant check
//   if (det < 0) 
//     return false; 
  
//   det = sqrt(det);
//   t = b - det;
//   if (t <= eps)
//   {
//     t = b + det;
//     if (t <= eps)
//       t = -1;
//   }

//   if (t < 0) return false;
//   return true;
// }

bool rayCylinderIntersect(ray r, cylinder cyl, float& t, VEC3& p, VEC3& n) {
  VEC3 translation = cyl.info.translation;
  Matrix3d R = cyl.info.R;

  VEC3 o_new = R.transpose() * (r.origin - translation);
  VEC3 dir_new = R.transpose() * r.dir;
  // printf("GOT TO HERE!\n");
  // std::cout << "o_new:\n" << o_new << endl;
  // std::cout << "dir_new:\n" << dir_new << endl;
  float qa = dir_new[0]*dir_new[0] + dir_new[1]*dir_new[1];
  float qb = 2*(dir_new[0]*o_new[0] + dir_new[1]*o_new[1]);
  float qc = o_new[0]*o_new[0] + o_new[1]*o_new[1] - cyl.radius*cyl.radius;
  float det = qb*qb - 4*qa*qc;
  // printf("GOT TO HERE 2\n");
  if (det < 0) {
    // printf("no intersect\n");
    return false;
  }
  else {
    // printf("det > 0\n");
    float t1 = (-qb - sqrt(qb*qb - 4*qa*qc)) / (2 * qa);
    float t2 = (-qb + sqrt(qb*qb - 4*qa*qc)) / (2 * qa);
    // std::cout << "o_new:\n" << o_new << endl;
    // std::cout << "dir_new:\n" << dir_new << endl;
    VEC3 p1 = o_new + t1*dir_new;
    VEC3 p2 = o_new + t2*dir_new;
    VEC3 lv_new = (R.transpose() * (cyl.leftVertex - translation));
    VEC3 rv_new = (R.transpose() * (cyl.rightVertex - translation));
    float z_min = lv_new[2];
    float z_max = rv_new[2];
    // cout << "R:\n" << R << endl;
    // printf("z_min = %f, z_max = %f\n", z_min, z_max);
    // cout << "p1 = " << p1 << endl;
    if (z_min > z_max) {
      float temp = z_min;
      z_min = z_max;
      z_max = temp;
    }
    if (p1[2] < z_max && p1[2] > z_min && t1 > 0.001) {
      t = t1;
      p = R*p1 + translation;
      n = (R*(p1 - VEC3(0, 0, p1[2]))).normalized();
      // printf("true\n");
      return true;
    }
    else if (p2[2] < z_max && p2[2] > z_min && t2 > 0.001) {
      t = t2;
      p = R*p2 + translation;
      n = (R*(p2 - VEC3(0, 0, p2[2]))).normalized();
      // printf("true\n");
      return true;
    }
    else {
      return false;
    }
  }
}

bool rayTriangleIntersect(ray r, triangle tri, float& t_out) {
    float x_a = tri.a[0];
    float y_a = tri.a[1];
    float z_a = tri.a[2];
    float x_b = tri.b[0];
    float y_b = tri.b[1];
    float z_b = tri.b[2];
    float x_c = tri.c[0];
    float y_c = tri.c[1];
    float z_c = tri.c[2];
    float x_e = r.origin[0];
    float y_e = r.origin[1];
    float z_e = r.origin[2];
    float x_d = r.dir[0];
    float y_d = r.dir[1];
    float z_d = r.dir[2];
    float M = (x_a-x_b)*((y_a-y_c)*z_d - y_d*(z_a-z_c))
            + (y_a-y_b)*(x_d*(z_a-z_c) - (x_a-x_c)*z_d)
            + (z_a-z_b)*((x_a-x_c)*y_d - (y_a-y_c)*x_d);
    
    float t = -((z_a-z_c)*((x_a-x_b)*(y_a-y_e) - (x_a-x_e)*(y_a-y_b))
            +  (y_a-y_c)*((x_a-x_e)*(z_a-z_b) - (x_a-x_b)*(z_a-z_e))
            +  (x_a-x_c)*((y_a-y_b)*(z_a-z_e) - (y_a-y_e)*(z_a-z_b)))/M;
    // printf("t for triangle %i = %f", i, t);
    if (t <= 0) {
      return false;
    }
    float gamma =(z_d*((x_a-x_b)*(y_a-y_e) - (x_a-x_e)*(y_a-y_b))
                + y_d*((x_a-x_e)*(z_a-z_b) - (x_a-x_b)*(z_a-z_e))
                + x_d*((y_a-y_b)*(z_a-z_e) - (y_a-y_e)*(z_a-z_b)))/M;
    if (gamma < 0 || gamma > 1) {
      return false;
    }
    float beta =((x_a-x_e)*((y_a-y_c)*z_d - y_d*(z_a-z_c))
               + (y_a-y_e)*(x_d*(z_a-z_c) - (x_a-x_c)*z_d)
               + (z_a-z_e)*((x_a-x_c)*y_d - (y_a-y_c)*x_d))/M;
    if (beta < 0 || beta > 1 - gamma) {
      return false;
    }
    if (t > 0.00001) {
      t_out = t;
      return true;

    }
}

bool raySphereIntersect(ray r, sphere s, float& t_out) {
  float a = r.dir.dot(r.dir);
  float b = 2 * r.dir.dot(r.origin - s.center);
  float c = (r.origin - s.center).dot(r.origin - s.center) - s.radius*s.radius;
  float t, t1, t2;
  if (b*b >= 4*a*c) {
    t1 = (-b - sqrt(b*b - 4*a*c)) / (2 * a);
    t2 = (-b + sqrt(b*b - 4*a*c)) / (2 * a);
    if (t1 <= 0.00001) {
      t = t2;
    }
    else {
      t = fmin(t1, t2);
    }
    if (t > 0.0001) {
        t_out = t;
        return true;
      }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// void rayColor(ray r, VEC3& pixelColor) 
// {
//   pixelColor = VEC3(1,1,1);

//   // look for intersections
//   int hitID = -1;
//   intersectType hitType;
//   float tMinFound = FLT_MAX;
//   for (int y = 0; y < sphereCenters.size(); y++) {
//     float tMin = FLT_MAX;
//     if (raySphereIntersect(sphereCenters[y], sphereRadii[y], r.origin, r.dir, tMin)) { 
//       // is the closest so far?
//       if (tMin < tMinFound) {
//         tMinFound = tMin;
//         hitID = y;
//         hitType = SPHERE;
//       }
//     }
//   }
  
//   for (int i = 0; i < cylinders.size(); i++) {
//     float tMin = FLT_MAX;
//     if (rayCylinderIntersect(r, cylinders[i], tMin)) {
//       if (tMin < tMinFound) {
//         tMinFound = tMin;
//         hitID = i;
//         hitType = CYLINDER;
//       }
//     }
//   }
  
  
  
//   // No intersection, return white
//   if (hitID == -1)
//     return;

//   // set to the sphere color
//   pixelColor = VEC3(1, 0, 0);
// }

bool shadow_ray_intersect_tri(ray r, light li) {
  //o + td = li.pos
  //o[0] + td[0] = li.pos[0]
  float t_min = (li.pos[0] - r.origin[0])/r.dir[0];
    int hitID = -1;
  intersectType hitType;
  float tMinFound = FLT_MAX;
  for (int i = 0; i < spheres.size(); i++) {
    float tMin = FLT_MAX;
    if (raySphereIntersect(r, spheres[i], tMin)) { 
      // is the closest so far?
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = SPHERE;
      }
    }
  }
  
  for (int i = 0; i < cylinders.size(); i++) {
    float tMin = FLT_MAX;
    VEC3 p;
    VEC3 n;
    if (rayCylinderIntersect(r, cylinders[i], tMin, p, n)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = CYLINDER;

      }
    }
  }

  for (int i = 0; i < tris.size(); i++) {
    float tMin = FLT_MAX;
    if (rayTriangleIntersect(r, tris[i], tMin)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = TRIANGLE;
      }
    }
  }
  if (hitID != -1) {
    return true;
  }
  return false;
}

float fresnel(float costheta, float cosphi, int into_glass) {
  float alpha_1 = 1.5;
  float alpha_2 = 1;
  if (into_glass) {
    alpha_1 = 1;
    alpha_2 = 1.5;
  }
  // float p_par = (alpha_2*costheta - alpha_1*cosphi)/(alpha_2*costheta + alpha_1*cosphi);
  // float p_perp = (alpha_1*costheta - alpha_2*cosphi)/(alpha_2*costheta + alpha_1*cosphi);
  // return 0.5 * (p_par*p_par + p_perp*p_perp);
  float R_0 = pow(((alpha_2 -1)/(alpha_2 + 1)), 2);
  float R = R_0 + (1-R_0)*pow((1-costheta), 5);
  return R;
}

VEC3 rayColor(ray r, int nlights, light* lights, int depth) {
  // printf("r\n");
  float t_min = INFINITY;
  VEC3 color = {0,0,0};
  if (depth == 10) {
    return color;
  }

  int hitID = -1;
  intersectType hitType;
  float tMinFound = FLT_MAX;
  for (int i = 0; i < spheres.size(); i++) {
    float tMin = FLT_MAX;
    if (raySphereIntersect(r, spheres[i], tMin)) { 
      // is the closest so far?
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = SPHERE;
      }
    }
  }
  VEC3 p_found;
  VEC3 n_found;
  VEC3 p_cyl;
  VEC3 n_cyl;
  for (int i = 0; i < cylinders.size(); i++) {
    float tMin = FLT_MAX;
    if (rayCylinderIntersect(r, cylinders[i], tMin, p_cyl, n_cyl)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = CYLINDER;
        p_found = p_cyl;
        n_found = n_cyl;
      }
    }
  }

  for (int i = 0; i < tris.size(); i++) {
    float tMin = FLT_MAX;
    if (rayTriangleIntersect(r, tris[i], tMin)) {
      if (tMin < tMinFound) {
        tMinFound = tMin;
        hitID = i;
        hitType = TRIANGLE;
      }
    }
  }
  // printf("done intersections\n");
  // printf("hitID = %i\n", hitID);

  if (hitID == -1) {
    return color;
  }

  float t = tMinFound;

  if (hitType == CYLINDER) {
    cylinder cyl = cylinders[hitID];
    VEC3 p = p_found;
    VEC3 n = n_found;
    color = {0,0,0};
    VEC3 e = (r.origin-p)/(r.origin-p).norm();
    for (int j = 0; j < nlights; j++) {
      light li = lights[j];
      ray shadowray;
      // shadowray.origin = p;
      // shadowray.dir = (li.pos - p)/(li.pos - p).norm();
      VEC3 surface1 = n.cross(e);
      VEC3 surface2 = surface1.cross(n);
      int shadowrays_blocked = 0;
      for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
          shadowray.origin = p + x*SS_GRID_SIZE*surface1 + y*SS_GRID_SIZE*surface2
            + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface1 
            + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface2;
          shadowray.dir = (li.pos - shadowray.origin).normalized();
          if (shadow_ray_intersect_tri(shadowray, li)) {
            shadowrays_blocked++;
          }
        }
      }
      // printf("shadowrays blocked: %d\n", shadowrays_blocked);
      float shadow_ratio = ((float) (25 - shadowrays_blocked)) / ((float) 25);
      // shadow_ratio = 1;
      VEC3 l = (li.pos - p)/(li.pos - p).norm();
      color += shadow_ratio * cyl.color.cwiseProduct(li.color * fmax(0, n.dot(l)));
      VEC3 refl = -l + 2*(n.dot(l))*n;
      color += shadow_ratio * cyl.color.cwiseProduct(li.color * pow(fmax(0, refl.dot(e)), 10));
    }
  }

  if (hitType == SPHERE) {
    sphere s = spheres[hitID];
    VEC3 p = r.origin + t * r.dir;
    VEC3 n = (p - s.center)/(p - s.center).norm();
    if (s.material == PLASTIC) {
      color = {0,0,0};
      VEC3 e = (r.origin-p)/(r.origin-p).norm();
      for (int j = 0; j < nlights; j++) {
        light li = lights[j];
        ray shadowray;
        // shadowray.origin = p;
        // shadowray.dir = (li.pos - p)/(li.pos - p).norm();
        VEC3 surface1 = n.cross(e);
        VEC3 surface2 = surface1.cross(n);
        int shadowrays_blocked = 0;
        for (int x = -2; x <= 2; x++) {
          for (int y = -2; y <= 2; y++) {
            shadowray.origin = p + x*SS_GRID_SIZE*surface1 + y*SS_GRID_SIZE*surface2
              + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface1 
              + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface2;
            shadowray.dir = (li.pos - shadowray.origin).normalized();
            if (shadow_ray_intersect_tri(shadowray, li)) {
              shadowrays_blocked++;
            }
          }
        }
        // if (shadowrays_blocked > 0 && shadowrays_blocked < 25) {
        //   printf("shadowrays blocked: %d\n", shadowrays_blocked);
        // }
        float shadow_ratio = ((float) (25 - shadowrays_blocked)) / ((float) 25);
        VEC3 l = (li.pos - p)/(li.pos - p).norm();
        color += shadow_ratio * s.color.cwiseProduct(li.color * fmax(0, n.dot(l)));
        VEC3 refl = -l + 2*(n.dot(l))*n;
        color += shadow_ratio * s.color.cwiseProduct(li.color * pow(fmax(0, refl.dot(e)), 10));
      }
    }
    else if (s.material == MIRROR) {
      ray reflect;
      VEC3 e = (r.origin-p)/(r.origin-p).norm();
      reflect.origin = p;
      reflect.dir = (-e + 2*(n.dot(e))*n)/(-e + 2*(n.dot(e))*n).norm();
      color = rayColor(reflect, nlights, lights, depth+1);
    }
    else if (s.material == GLASS) {
      ray refract;
      refract.origin = p;
      VEC3 in = (p - r.origin)/(p - r.origin).norm();
      float costheta;
      float cosphi;
      float k_reflect;
      if (r.in_glass) {
        // printf("in glass\n");
        n = -n;
        // printf("%f\n", n[2]);
        VEC3 dir1 = (1.5/1) * (in - n * in.dot(n));
        VEC3 dir2 = sqrt(1 - pow(1.5/1, 2) * (1 - pow(in.dot(n), 2))) * n;
        refract.dir = (dir1 - dir2)/(dir1 - dir2).norm();
        refract.in_glass = 0;
        costheta = -in.dot(n);
        cosphi = sqrt(1 - pow((1.5/1 * sqrt(1 - costheta*costheta)), 2)); 
        k_reflect = fresnel(costheta, cosphi, 0);
      }
      else {
        // printf("hit glass\n");
        VEC3 dir1 = (1/1.5) * (in - n * in.dot(n));
        VEC3 dir2 = sqrt(1 - pow(1/1.5, 2) * (1 - pow(in.dot(n), 2))) * n;
        refract.dir = (dir1 - dir2)/(dir1 - dir2).norm();
        refract.in_glass = 1;
        costheta = -in.dot(n);
        cosphi = sqrt(1 - pow((1/1.5 * sqrt(1 - costheta*costheta)), 2));
        k_reflect = fresnel(costheta, cosphi, 1);
      }

      ray reflect;
      VEC3 e = (r.origin-p)/(r.origin-p).norm();
      reflect.origin = p;
      reflect.dir = (-e + 2*(n.dot(e))*n)/(-e + 2*(n.dot(e))*n).norm();


      VEC3 refract_color = rayColor(refract, nlights, lights, depth+1);
      VEC3 reflect_color = rayColor(reflect, nlights, lights, depth+1);
      color = k_reflect*reflect_color + (1-k_reflect)*refract_color;
    }
  }

  // printf("t_min = %f\n", t_min);
  if (hitType == TRIANGLE) {
    triangle tri = tris[hitID];

    VEC3 p = r.origin + t * r.dir;
    VEC3 n = -(tri.c-tri.a).cross(tri.b-tri.a)/((tri.c-tri.a).cross(tri.b-tri.a)).norm();
    if (tri.material == PLASTIC) {
      color = {0,0,0};
      VEC3 e = (r.origin-p)/(r.origin-p).norm();
      for (int j = 0; j < nlights; j++) {
        light li = lights[j];
        ray shadowray;
        // shadowray.origin = p;
        // shadowray.dir = (li.pos - p)/(li.pos - p).norm();
        VEC3 surface1 = n.cross(e);
        VEC3 surface2 = surface1.cross(n);
        int shadowrays_blocked = 0;
        for (int x = -2; x <= 2; x++) {
          for (int y = -2; y <= 2; y++) {
            shadowray.origin = p + x*SS_GRID_SIZE*surface1 + y*SS_GRID_SIZE*surface2
              + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface1 
              + (((float) rand()/(float) RAND_MAX)*SS_GRID_SIZE - SS_GRID_SIZE/2)*surface2;
            shadowray.dir = (li.pos - shadowray.origin).normalized();
            if (shadow_ray_intersect_tri(shadowray, li)) {
              shadowrays_blocked++;
            }
          }
        }
        float shadow_ratio = ((float) (25 - shadowrays_blocked)) / ((float) 25);
        VEC3 l = (li.pos - p)/(li.pos - p).norm();
        color += shadow_ratio * tri.color.cwiseProduct(li.color * fmax(0, n.dot(l)));
        VEC3 refl = -l + 2*(n.dot(l))*n;
        color += shadow_ratio * tri.color.cwiseProduct(li.color * pow(fmax(0, refl.dot(e)), 10));
      }
    }
    else if (tri.material == MIRROR) {
        ray reflect;
        VEC3 e = (r.origin-p)/(r.origin-p).norm();
        reflect.origin = p;
        reflect.dir = (-e + 2*(n.dot(e))*n)/(-e + 2*(n.dot(e))*n).norm();
        color = rayColor(reflect, nlights, lights, depth+1);
      }
      // color = tri.color;
  }

  color[0] = fmin(color[0], 1.0);
  color[1] = fmin(color[1], 1.0);
  color[2] = fmin(color[2], 1.0);
  return color;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
float clamp(float value)
{
  if (value < 0.0)      return 0.0;
  else if (value > 1.0) return 1.0;
  return value;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void renderImage(int& xRes, int& yRes, const string& filename, VEC3** prev_frames) 
{
  // allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];

  // compute image plane
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 4.0f / 3.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();
  VEC3* curr_frame = (VEC3*) malloc(sizeof(VEC3) * xRes * yRes);
  // int n_prev_frames = 0;
  // n_prev_frames += prev_frames[0] != NULL;
  // n_prev_frames += prev_frames[1] != NULL;
  // n_prev_frames += prev_frames[2] != NULL;
  // n_prev_frames += prev_frames[3] != NULL;
  // printf("n_prev_frames = %d\n", n_prev_frames);
  for (int y = 0; y < yRes; y++) 
    for (int x = 0; x < xRes; x++) 
    {
      // generate the ray, making x-axis go left to right
      const float ratioX = 1.0f - ((xRes - 1) - x) / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt + 
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      const VEC3 rayDir = (rayHitImage - eye).normalized();

      // get the color
      ray r;
      r.origin = eye;
      r.dir = rayDir;
      VEC3 color = {0,0,0};
      if (prev_frames[3] == NULL) {
        color = rayColor(r, nlights, lights, 0);
      }
      else {
        color += 0.5 * rayColor(r, nlights, lights, 0);
        color += 0.25 * prev_frames[0][y * xRes + x];
        color += 0.125 * prev_frames[1][y * xRes + x];
        color += 0.0625 * prev_frames[2][y * xRes + x];
        color += 0.0625 * prev_frames[3][y * xRes + x];
      }
      // VEC3 color = rayColor(r, nlights, lights, 0);
      curr_frame[y * xRes + x] = color;

      // set, in final image
      ppmOut[3 * (y * xRes + x)] = clamp(color[0]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 1] = clamp(color[1]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 2] = clamp(color[2]) * 255.0f;
    }
  free(prev_frames[3]);
  prev_frames[3] = prev_frames[2];
  prev_frames[2] = prev_frames[1];
  prev_frames[1] = prev_frames[0];
  prev_frames[0] = curr_frame;
  writePPM(filename, xRes, yRes, ppmOut);

  delete[] ppmOut;
}

//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void setSkeletonsToSpecifiedFrame(int frameIndex)
{
  if (frameIndex < 0)
  {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL)
  {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames())
    {
      cout << " We hit the last frame! You might want to pick a different sequence. " << endl;
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    }
    else 
      postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Build a list of spheres in the scene
//////////////////////////////////////////////////////////////////////////////////
void buildScene()
{
  sphereCenters.clear();
  sphereRadii.clear();
  sphereColors.clear();
  // spheres.clear();
  displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);
  cylinders.clear();

  // lights[0].color = VEC3(1, 0.3, 1);
  // lights[0].pos = VEC3(0, 12, -3);
  lights[0].color = VEC3(1, 1, 1);
  lights[0].pos = VEC3(-10, 10, 0);
  nlights = 1;

  // 1.4, 0, 0.7
  cylinder c0;
  c0.leftVertex = VEC3(1.3, 0, 1.8);
  c0.rightVertex = VEC3(1.4, 0, 0.7);
  c0.radius = 0.35;
  c0.color = VEC3(0.7,0.25,0.5);
  cylinders.push_back(c0);
  // sphereCenters.push_back(VEC3(5, 0.5, 1));
  // sphereRadii.push_back(2);

  sphere s0;
  s0.center = VEC3(5, -1, 2);
  s0.radius = 4;
  s0.material = PLASTIC;
  s0.color = VEC3{0, 1, 0};
  // spheres.push_back(s0);

  

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
  vector<float>& lengths     = displayer.lengths();

  // build a sphere list, but skip the first bone, 
  // it's just the origin
  int totalBones = rotations.size();
  for (int x = 1; x < totalBones; x++)
  {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    // get the direction vector
    VEC3 direction = (rightVertex - leftVertex).head<3>();
    const float magnitude = direction.norm();
    direction *= 1.0 / magnitude;

    // how many spheres?
    const float sphereRadius = 0.05;
    const int totalSpheres = magnitude / (2.0 * sphereRadius);
    const float rayIncrement = magnitude / (float)totalSpheres;

    // store the spheres
    // sphereCenters.push_back(leftVertex.head<3>());
    // sphereRadii.push_back(0.05);
    // sphereColors.push_back(VEC3(1,0,0));
    
    // sphereCenters.push_back(rightVertex.head<3>());
    // sphereRadii.push_back(0.05);
    // sphereColors.push_back(VEC3(1,0,0));
    // for (int y = 0; y < totalSpheres; y++)
    // {
    //   VEC3 center = ((float)y + 0.5) * rayIncrement * direction + leftVertex.head<3>();
    //   sphereCenters.push_back(center);
    //   sphereRadii.push_back(0.05);
    //   sphereColors.push_back(VEC3(1,0,0));
    // } 
    cylinder c;
    c.leftVertex = leftVertex.head<3>();
    c.rightVertex = rightVertex.head<3>();
    c.color = VEC3(0.9, 0.7, 0.5);
    c.radius = 0.05;
    cylinders.push_back(c);
  }
  calculateCylinderInfo();
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  string skeletonFilename("36.asf");
  string motionFilename("36_24.amc");
  
  // load up skeleton stuff
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);

  // load up the motion
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));

  // Note we're going 8 frames at a time, otherwise the animation
  // is really slow.
  // for (int x = -6; x < 8; x++) {
  //   for (int z = -4; z < 4; z++) {
  //     triangle t0;
  //     t0.a = VEC3(2*x, 0, 2*z);
  //     t0.b = VEC3(2*x, 0, 2*(z+1));
  //     t0.c = VEC3(2*(x+1), 0, 2*z);
  //     t0.color = VEC3(0.25, 0.25, 0.4);
  //     t0.material = PLASTIC;
  //     tris.push_back(t0);
  //     triangle t1;
  //     t1.a = VEC3(2*(x+1), 0, 2*(z+1));
  //     t1.c = VEC3(2*x, 0, 2*(z+1));
  //     t1.b = VEC3(2*(x+1), 0, 2*z);
  //     t1.color = VEC3(0.1, 0.1, 0.1);
  //     t1.material = PLASTIC;
  //     tris.push_back(t1);
  //   }
  // }
  triangle t0;
  t0.a = VEC3(-12, 0, -8);
  t0.b = VEC3(-12, 0, 8);
  t0.c = VEC3(16, 0, -8);
  t0.color = VEC3(0.25, 0.25, 0.4);
  t0.material = PLASTIC;
  tris.push_back(t0);
  triangle t1;
  t1.a = VEC3(16, 0, 8);
  t1.c = VEC3(-12, 0, 8);
  t1.b = VEC3(16, 0, -8);
  t1.color = VEC3(0.25, 0.25, 0.4);
  t1.material = PLASTIC;
  tris.push_back(t1);
  

  //LEFT Z IS POSITIVE
  sphere s0;
  s0.center = VEC3(1.4, 0, 0.7);
  s0.radius = 0.45;
  s0.material = MIRROR;
  s0.color = VEC3(0, 1, 0);
  spheres.push_back(s0);
  sphere s1;
  s1.center = VEC3(0.5, 0, 0.7);
  s1.radius = 0.25;
  s1.material = PLASTIC;
  s1.color = VEC3(0, 0.75, 0.25);
  spheres.push_back(s1);
  sphere s2;
  s2.center = VEC3(1, 0, -0.17);
  s2.radius = 0.33;
  s2.material = PLASTIC;
  s2.color = VEC3(0, 0.75, 0.75);
  spheres.push_back(s2);
  sphere s3;
  s3.center = VEC3(0.8, -0.1, -1);
  s3.radius = 0.28;
  s3.material = MIRROR;
  s3.color = VEC3(0, 1, 0);
  spheres.push_back(s3);
  s3.center = VEC3(1.4, -0.2, -0.8);
  s3.radius = 0.35;
  s3.material = PLASTIC;
  s3.color = VEC3(0.8, 0.2, 0);
  spheres.push_back(s3);

  VEC3** prev_frames = (VEC3**) malloc(4 * sizeof(VEC3*));
  prev_frames[0] = NULL;
  prev_frames[1] = NULL;
  prev_frames[2] = NULL;
  prev_frames[3] = NULL; 
  for (int x = 0; x < 2400; x += 8)
  {
    setSkeletonsToSpecifiedFrame(x);
    buildScene();
    // VEC3 lookingAt(5, 0.5, 1);
    if (x > 400 && x <= 750) {
      lookingAt = VEC3(5, 0.5, 1-(((float) x-400)/350)*2.3);
    }
    if (x > 1200 && x < 1800) {
      lookingAt = VEC3(5, 0.5, -1.3+(((float) x-1200)/600)*2.3);
    }

    char buffer[256];
    sprintf(buffer, "./frames/frame.%04i.ppm", x / 8);
    renderImage(windowWidth, windowHeight, buffer, prev_frames);
    cout << "Rendered " + to_string(x / 8) + " frames" << endl;
  }
  free(prev_frames[0]);
  free(prev_frames[1]);
  free(prev_frames[2]);
  free(prev_frames[3]);
  free(prev_frames);

  return 0;
}
