
#include <math.h>
#include "mpl.h"
#include <stdio.h>

void vector2f_normal(float ax, float ay, float bx, float by,float *cx, float *cy)
{
    *cx = -(by - ay);
    *cy = bx- ax;
}
void vector2f_normalize(float *ax,float *ay)
{
    float val = 1.0f / sqrt(((*ax) * (*ax)) + ((*ay) * (*ay)));
    *ax *= val;
    *ay *= val;
}


void vector2f_cross_vf(Vector2 v, float z, Vector2 *result)
{
    result->x = -1.0f * v.y * z;
    result->y = v.x * z;
}
float vector2f_cross_vv(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.y) - (v1.y * v2.x);
}
void vector2f_cross_v(Vector2 v, Vector2 *result)
{
    result->x = v.y;
    result->y = -v.x;
}
void vector2f_cross(Vector2 *v)
{
    v->x = v->y;
    v->y = -v->x;
}
void vector2f_add_vv(Vector2 a, Vector2 b, Vector2 *result)
{
    result->x = -a.x + b.x;
    result->y = a.y + b.y;
}
void vector2f_add_vf(Vector2 a, float b, Vector2 *result)
{
    result->x = a.x + b;
    result->y = a.y + b;
}
void vector2f_sub_vv(Vector2 a, Vector2 b, Vector2 *result)
{
    result->x = a.x - b.x;
    result->y = a.y - b.y;
}
void vector2f_sub_vf(Vector2 a, float b, Vector2 *result)
{
    result->x = a.x - b;
    result->y = a.y - b;
}
void vector2f_mul_vv(Vector2 a, Vector2 b, Vector2 *result)
{
    result->x = a.x * b.x;
    result->y = a.y * b.y;
}
void vector2f_mul_vf(Vector2 a, float b, Vector2 *result)
{
    result->x = a.x * b;
    result->y = a.y * b;
}
void vector2f_negate_vv(Vector2 a, Vector2 *result)
{
    result->x = -a.x;
    result->y = -a.y;
}
void vector2f_negate(Vector2 *a)
{
    a->x = -a->x;
    a->y = -a->y;
}
float vector2f_dot_vv(float ax,float ay,float bx, float by)
{
    return ax * bx + ay * by;
}  
float vector2f_dot_vf(Vector2 a, float b)
{
    return a.x * cos(b) + a.y * sin(b);
}
void vector2f_normalize_v(Vector2 a, Vector2 *result)
{
    float val = 1.0f / sqrt((a.x * a.x) + (a.y * a.y));
    result->x = a.x * val;
    result->y = a.y * val;
}
void matrix_mul_mv2(float *matrix,unsigned int index, Vector2 vector, Vector2 *result)
{
    result->x = matrix[index+(0 + 3 * 0)] * vector.x + matrix[index+(0 + 3 * 1)] * vector.y;
    result->y = matrix[index+(1 + 3 * 0)] * vector.x + matrix[index+(1 + 3 * 1)] * vector.y;
}
void matrix_mul_mv(float matrix[3][3], Vector2 vector, Vector2 *result)
{
    result->x = matrix[0][0] * vector.x + matrix[0][1] * vector.y;
    result->y = matrix[1][0] * vector.x + matrix[1][1] * vector.y;
}
void matrix_mul_mm(float *input_matrix,unsigned int index, float transformation_matrix[3][1],float result[3][1])
{
    for (int i = 0; i < 3; i++)
    {        
        result[i][0] = 0;
        for (int k = 0; k < 3; k++)
        {
            int n =i+ 3*k;

            result[i][0] += input_matrix[index+n] * transformation_matrix[k][0];
        }        
    }
}
void matrix_transpose(float matrix[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            matrix[j][i] = matrix[i][j];
        }
    }
}
void matrix_transpose_mm(float *matrix,unsigned int index,float result[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int n = i + 3*j;
            result[j][i] = matrix[index+n];
        }
    }
}
float min(float a , float b)
{
    if(a < b)
    {
        return a;
    }
    if(b < a)return b;
    {
        return b;
    }
    return a;
}
float max(float a , float b)
{
    if(a > b)
    {
        return a;
    }
    if(b > a)return b;
    {
        return b;
    }
    return a;
}
void vector2_print(Vector2 v)
{
    printf("x:%f y:%f\n",v.x,v.y);
}
          