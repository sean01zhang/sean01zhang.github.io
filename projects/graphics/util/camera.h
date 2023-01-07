#ifndef CAMERA_H
#define CAMERA_H

#include "rtweekend.h"

class camera {
public:
  camera() {
    double aspect = 16.0 / 9.0;
    double vp_height = 2.0;
    double vp_width = vp_height * aspect;
    double focal_length = 1.0;

    origin = vec3(0, 0, 0);
    horizontal = vec3(vp_width, 0, 0);
    vertical = vec3(0, vp_height, 0);
    lower_left =
        origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
  }

  ray get_ray(double u, double v) const {
    return ray(origin, lower_left + horizontal * u + vertical * v - origin);
  }

public:
  vec3 origin;
  vec3 horizontal;
  vec3 vertical;
  vec3 lower_left;
};

#endif
