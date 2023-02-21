#ifndef HITTABLE_H
#define HITTABLE_H

#include "vec3.h"
#include "ray.h"

struct hit_record {
  point3 point;
  vec3 norm;
  double t;
  bool front_face;

  inline void set_face_normal(const ray& r, const vec3& out_norm) {
    front_face = dot(r.direction(), out_norm) < 0;
    norm = front_face ? out_norm : -out_norm;
  }
};

class hittable {
  public:
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

#endif