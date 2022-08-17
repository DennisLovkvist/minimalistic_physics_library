
#include <math.h>
#include "mpl.h"

void vector2f_normal(Vector2 a, Vector2 b, Vector2 *result)
{
    result->x = -(b.y - a.y);
    result->y = b.x- a.x;
}
void vector2f_cross_vf(Vector2 v, float z, Vector2 *result)
{
    result->x = -1.0f * v.y * z;
    result->y = v.x * z;
}
void vector2f_cross_vv(Vector2 v1, Vector2 v2, float *result)
{
    *result = (v1.x * v2.y) - (v1.y * v2.x);
}
void vector2f_cross_v(Vector2 v, Vector2 *result)
{
    result->x = v.y;
    result->y = -v.x;
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
    result->x = a.x + b.x;
    result->y = a.y + b.y;
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
void vector2f_negate(Vector2 a, Vector2 *result)
{
    result->x = -a.x;
    result->y = -a.y;
}
void vector2f_dot_vv(Vector2 a, Vector2 b, float *result)
{
    *result = a.x * b.x + a.y * b.y;
}
void vector2f_dot_vf(Vector2 a, float b, float *result)
{
    *result = a.x * cos(b) + a.y * sin(b);
}
void vector2f_normalize_v(Vector2 a, Vector2 *result)
{
    float val = 1.0f / sqrt((a.x * a.x) + (a.y * a.y));
    result->x = a.x * val;
    result->y = a.y * val;
}
void vector2f_normalize(Vector2 *a)
{
    float val = 1.0f / sqrt((a->x * a->x) + (a->y * a->y));
    a->x *= val;
    a->y *= val;
}
void matrix_mul_mv(float matrix[3][3], Vector2 vector, Vector2 *result)
{
    result->x = matrix[0][0] * vector.x + matrix[0][1] * vector.y;
    result->y = matrix[1][0] * vector.x + matrix[1][1] * vector.y;
}
void matrix_mul_mm(float input_matrix[3][3], float transformation_matrix[3][1],float result[3][1])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < 3; k++)
            {
                result[i][j] += input_matrix[i][k] * transformation_matrix[k][j];
            }
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
void matrix_transpose_mm(float matrix[3][3],float result[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[j][i] = matrix[i][j];
        }
    }
}
        