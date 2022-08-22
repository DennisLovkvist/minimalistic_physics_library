#pragma once
#include "mpl.h"
#define FLT_MAX 3.4028235E+38F /* max value */
#define FLT_MIN -3.4028235E+38F /* min positive value */
#define PI 3.14159265358979323846

void vector2f_normal(float ax, float ay, float bx, float by,float *cx, float *cy);
void vector2f_normalize(float *ax,float *ay);
float vector2f_dot_vv(float ax,float ay,float bx, float by);
void matrix_mul_mv(float matrix[3][3], Vector2 vector, Vector2 *result);
void vector2f_cross_vf(Vector2 v, float z, Vector2 *result);
float vector2f_cross_vv(Vector2 v1, Vector2 v2);
void vector2f_cross_v(Vector2 v, Vector2 *result);
float min(float a , float b);
float max(float a , float b);
void matrix_mul_mm(float *input_matrix,unsigned int index, float transformation_matrix[3][1],float result[3][1]);
void matrix_mul_mv2(float *matrix,unsigned int index, Vector2 vector, Vector2 *result);
void matrix_transpose_mm(float *matrix,unsigned int index,float result[3][3]);