#pragma once
#include "mpl.h"
#define FLT_MAX 3.4028235E+38F /* max value */
#define FLT_MIN -3.4028235E+38F /* min positive value */
#define PI 3.14159265358979323846

void vector2f_normal(Vector2 a, Vector2 b, Vector2 *result);
void vector2f_cross_vf(Vector2 v, float z, Vector2 *result);
float vector2f_cross_vv(Vector2 v1, Vector2 v2);
void vector2f_cross_v(Vector2 v, Vector2 *result);
void vector2f_cross(Vector2 *v);
void vector2f_add_vv(Vector2 a, Vector2 b, Vector2 *result);
void vector2f_add_vf(Vector2 a, float b, Vector2 *result);
void vector2f_sub_vv(Vector2 a, Vector2 b, Vector2 *result);
void vector2f_sub_vf(Vector2 a, float b, Vector2 *result);
void vector2f_mul_vv(Vector2 a, Vector2 b, Vector2 *result);
void vector2f_mul_vf(Vector2 a, float b, Vector2 *result);
void vector2f_negate_vv(Vector2 a, Vector2 *result);
void vector2f_negate(Vector2 *a);
float vector2f_dot_vv(Vector2 a, Vector2 b);
float vector2f_dot_vf(Vector2 a, float b);
void vector2f_normalize_v(Vector2 a, Vector2 *result);
void vector2f_normalize(Vector2 *a);
void matrix_mul_mv(float matrix[3][3], Vector2 vector, Vector2 *result);
void matrix_mul_mm(float input_matrix[3][3], float transformation_matrix[3][1],float result[3][1]);
void matrix_transpose(float matrix[3][3]);
void matrix_transpose_mm(float matrix[3][3],float result[3][3]);
float min(float a , float b);
float max(float a , float b);
void vector2_print(Vector2 v);