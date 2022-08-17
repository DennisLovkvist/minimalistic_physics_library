
#include "mpl.h"
#include "delo_math.h"
#include <math.h>
void mpl_polygon_get_centroid(Vector2 *vertices, unsigned int vertex_count, Vector2 *centroid)
{
    centroid->x = vertices[0].x;
    centroid->y = vertices[0].y;
    
    for (int i = 1; i < vertex_count; i++)
    {
        centroid->x += vertices[i].x;
        centroid->y += vertices[i].y;
    }
    centroid->x /= vertex_count;
    centroid->y /= vertex_count;
}
void mpl_matrix_set(RotationMatrix *rotation_matrix,float tx, float ty, float r)
{    
    rotation_matrix->matrix[0][0] = cos(r);
    rotation_matrix->matrix[0][1] = sin(r);
    rotation_matrix->matrix[0][2] = tx;

    rotation_matrix->matrix[1][0] = -sin(r);
    rotation_matrix->matrix[1][1] = cos(r);
    rotation_matrix->matrix[1][2] = ty;

    rotation_matrix->matrix[2][0] = 0;
    rotation_matrix->matrix[2][1] = 0;
    rotation_matrix->matrix[2][2] = 1;
}
void mpl_polygon_transform(Polygon *polygon)
{
    float tx = polygon->position.x - polygon->position_last.x;
    float ty = polygon->position.y - polygon->position_last.y;
    float r = polygon->orientation - polygon->orientation_last;    
    
    for (int i = 0; i < polygon->vertex_count; i++)
    {
        mpl_matrix_set(&polygon->matrix,-polygon->position_last.x, -polygon->position_last.y, 0);
        polygon->matrix.input[0][0] = polygon->m_vertices[i].x;
        polygon->matrix.input[1][0] = polygon->m_vertices[i].y;
        polygon->matrix.input[2][0] = 1;
        matrix_mul_mm(polygon->matrix.matrix,polygon->matrix.input,polygon->matrix.output);
        polygon->m_vertices[i].x = polygon->matrix.output[0][0];
        polygon->m_vertices[i].y = polygon->matrix.output[1][0];

        mpl_matrix_set(&polygon->matrix,0, 0, r);
        polygon->matrix.input[0][0] = polygon->m_vertices[i].x;
        polygon->matrix.input[1][0] = polygon->m_vertices[i].y;
        polygon->matrix.input[2][0] = 1;
        matrix_mul_mm(polygon->matrix.matrix,polygon->matrix.input,polygon->matrix.output);
        polygon->m_vertices[i].x = polygon->matrix.output[0][0];
        polygon->m_vertices[i].y = polygon->matrix.output[1][0];

        polygon->matrix.input[0][0] = polygon->normals[i].x;
        polygon->matrix.input[1][0] = polygon->normals[i].y;
        polygon->matrix.input[2][0] = 1;
        matrix_mul_mm(polygon->matrix.matrix,polygon->matrix.input,polygon->matrix.output);
        polygon->normals[i].x = polygon->matrix.output[0][0];
        polygon->normals[i].y = polygon->matrix.output[1][0];

        mpl_matrix_set(&polygon->matrix,polygon->position_last.x + tx, polygon->position_last.y + ty, 0);
        polygon->matrix.input[0][0] = polygon->m_vertices[i].x;
        polygon->matrix.input[1][0] = polygon->m_vertices[i].y;
        polygon->matrix.input[2][0] = 1;
        matrix_mul_mm(polygon->matrix.matrix,polygon->matrix.input,polygon->matrix.output);
        polygon->m_vertices[i].x = polygon->matrix.output[0][0];
        polygon->m_vertices[i].y = polygon->matrix.output[1][0];
    }
    polygon->position_last.x = polygon->position.x;
    polygon->position_last.y = polygon->position.y;
    polygon->orientation_last = polygon->orientation;
}
void mpl_polygon_init(Polygon *polygon, Vector2 *vertices, unsigned int vertex_count)
{
    polygon->orientation = polygon->orientation_last = 0.0;
    polygon->radius = 0;
    vertex_count = (vertex_count > 4) ? 4:vertex_count;
    polygon->vertex_count = vertex_count;

    for (int i = 0; i < vertex_count; i++)
    {
        polygon->m_vertices[i].x = vertices[i].x;
        polygon->m_vertices[i].y = vertices[i].y;
    }    
    for (int i = 0; i < 2; i++)
    {
        polygon->m_edges[i].m_index_a = i;
        polygon->m_edges[i].m_index_b = (i + 1 < vertex_count) ? i + 1 : 0;    
        Vector2 v0 = polygon->m_vertices[polygon->m_edges[i].m_index_b];
        Vector2 v1 = polygon->m_vertices[polygon->m_edges[i].m_index_a];
        vector2f_normal(v0,v1,&polygon->normals[i]);
        vector2f_normalize(&polygon->normals[i]);
    }
    mpl_polygon_get_centroid(polygon->m_vertices,vertex_count, &polygon->position);
    
    polygon->position_last.x = polygon->position.x;
    polygon->position_last.y = polygon->position.y;

    mpl_polygon_transform(polygon);
    
}