#ifndef VECTOR_INCLUDED_H
#define VECTOR_INCLUDED_H
class Vector {
public:
   double x=0,y=0,z=0;
};
#endif
Vector v_add(Vector v1,Vector v2);
Vector v_sub(Vector v1,Vector v2); 
double v_dot(Vector v1,Vector v2);
Vector VWrap(Vector v,Vector box); 
Vector v_Sadd(Vector v1,Vector v2,double factor); 
Vector v_scale(Vector v,double s); 
