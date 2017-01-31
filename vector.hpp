#ifndef VECTOR_INCLUDED_H
#define VECTOR_INCLUDED_H
class Vector {
public:
  long double x=0,y=0,z=0;
};
#endif
Vector v_add(Vector v1,Vector v2);
Vector v_sub(Vector v1,Vector v2); 
long double v_dot(Vector v1,Vector v2);
Vector VWrap(Vector v,Vector box); 
Vector v_Sadd(Vector v1,Vector v2,long double factor); 
Vector v_scale(Vector v,long double s); 
