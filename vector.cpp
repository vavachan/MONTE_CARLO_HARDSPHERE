#include "vector.hpp"
#include<math.h>
 Vector v_add(Vector v1,Vector v2) {
    Vector v_sum;
    v_sum.x=v1.x+v2.x;
    v_sum.y=v1.y+v2.y;
    v_sum.z=v1.z+v2.z;
    return v_sum;
}
 Vector v_sub(Vector v1,Vector v2) {
    Vector v_diff;
    v_diff.x=v1.x-v2.x;
    v_diff.y=v1.y-v2.y;
    v_diff.z=v1.z-v2.z;
    return v_diff;
}

 long double v_dot(Vector v1,Vector v2) {
    return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}
 Vector VWrap(Vector v,Vector box) {
    long double x,y,z;
    long double twob=2*box.x;
    v.x=fmod(v.x,twob);
    v.y=fmod(v.y,twob);
    v.z=fmod(v.z,twob);
    if (v.x >= box.x) {
        v.x -=  twob;
    }
    else if (v.x < -box.x) {
        v.x += twob;
    }

    if (v.y >= box.y) {
        v.y -= twob;
    }
    else if (v.y < -box.y) {
        v.y += twob;
    }

    if (v.z >= box.z) {
        v.z -=twob;
    }
    else if (v.z < -box.z) {
        v.z += twob;
    }

    return v;
}

 Vector v_Sadd(Vector v1,Vector v2,long double factor) {
    Vector v_sum;
    v_sum.x=v1.x+factor*v2.x;
    v_sum.y=v1.y+factor*v2.y;
    v_sum.z=v1.z+factor*v2.z;
    return v_sum;
}

 Vector v_scale(Vector v,long double s) {
    v.x=v.x/s;
    v.y=v.y/s;
    v.z=v.z/s;
    return v;
}
