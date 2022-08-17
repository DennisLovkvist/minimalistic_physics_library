
#include <stdlib.h>
#include <stdio.h>
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
void mpl_polygon_get_support_point(Vector2 *vertices,unsigned int vertex_count, Vector2 dir,Vector2 *best_vertex)
{
    float bestProjection = FLT_MIN;   
    best_vertex->x = 0;
    best_vertex->y = 0;

    for (int i = 0; i < vertex_count; ++i)
    {
        Vector2 v = vertices[i];
        float projection;
        vector2f_dot_vv(v, dir, &projection);

        if (projection > bestProjection)
        {
            best_vertex->x = v.x;
            best_vertex->y = v.y;
            bestProjection = projection;
        }
    }
}
void mpl_polygon_get_reverse_support_point(Vector2 *vertices,unsigned int vertex_count, Vector2 dir,Vector2 *best_vertex)
{
    float bestProjection = FLT_MAX;
    best_vertex->x = 0;
    best_vertex->y = 0;
    for (int i = 0; i < vertex_count; ++i)
    {
        Vector2 v = vertices[i];
        float projection;
        vector2f_dot_vv(v, dir, &projection);

        if (projection < bestProjection)
        {            
            best_vertex->x = v.x;
            best_vertex->y = v.y;
            bestProjection = projection;
        }
    }
}
void mpl_polygon_get_maximum_vertex(Vector2 observer, Vector2 *vertices,unsigned int vertex_count,Vector2 *max_vertex)
{    
    mpl_polygon_get_centroid(vertices,vertex_count,max_vertex);
    Vector2 v = *max_vertex;
    float dx = observer.x - v.x;
    float dy = observer.y - v.y;
    float center_angle = atan2(dy, dx) * 180 / (float)PI;
    float delta_max = FLT_MIN;

    for (int a = 0; a < vertex_count; a++)
    {
        v = vertices[a];
        dx = observer.x - v.x;
        dy = observer.y - v.y;
        float angle = atan2(dy, dx) * 180 / (float)PI;
        float delta = abs(angle-center_angle);

        if (delta > delta_max)
        {
            delta_max = delta;
            *max_vertex = v;
        }
    }
}
int mpl_polygon_get_next_vertex_index(int index,int capacity)
{
    int next;

    if (index + 1 >= capacity)
    {
        next = 0;
    }
    else
    {
        next = index + 1;
    }
    return next;
}
int mpl_polygon_get_prev_vertex_index(int index, int capacity)
{
    int next;

    if (index - 1 < 0)
    {
        next = capacity-1;
    }
    else
    {
        next = 0;
    }
    return next;
}
float mpl_polygon_get_average_radius(Vector2* vertices,unsigned int vertex_count)
{
    Vector2 centroid;
    mpl_polygon_get_centroid(vertices,vertex_count,&centroid);
    float radius = 0;

    for (int a = 0; a < vertex_count; a++)
    {
        Vector2 v = vertices[a];
        float dx = centroid.x - v.x;
        float dy = centroid.y - v.y;
        float dist = sqrt(dx * dx + dy * dy);
        radius += dist;
    }
    radius /= vertex_count;

    return radius;
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
void mpl_rigid_body_compute_mass(RigidBody *rigid_body,float density)
{
    float area = 0;
    unsigned int vertex_count = rigid_body->polygon.vertex_count;
    Vector2 c;
    mpl_polygon_get_centroid(rigid_body->polygon.m_vertices,vertex_count,&c);
    float dx = 0;
    float dy = 0;
    float longest_radius = 0;
    float radius = 0;

    for (int i = 0; i < vertex_count; i++)
    {
        Vector2 v = rigid_body->polygon.m_vertices[i];
        int next = mpl_polygon_get_next_vertex_index(i,vertex_count);
        dx = c.x - v.x;
        dy = c.y - v.y;
        float d1 = sqrt(dx * dx + dy * dy);
        radius += d1;
        dx = rigid_body->polygon.m_vertices[next].x - v.x;
        dy = rigid_body->polygon.m_vertices[next].y - v.y;
        float d2 = sqrt(dx * dx + dy * dy);
        area += ((d1 * d2) * 0.5f);

        if (radius > longest_radius)
        {
            longest_radius = radius;
        }
    }
    radius /= vertex_count;

    rigid_body->polygon.radius = radius;
    rigid_body->longest_radius = longest_radius;

    rigid_body->mass = density * area;
    rigid_body->inertia = rigid_body->mass * radius * radius;
    rigid_body->invMass = 1.0f / rigid_body->mass;
    rigid_body->invInertia = 1.0f / rigid_body->inertia;
    
}
void mpl_rigid_body_init(RigidBody *rigid_body)
{
    Vector2 vertices[4];
    vertices[0].x = 100;vertices[0].y = 100;
    vertices[1].x = 200;vertices[1].y = 100;
    vertices[2].x = 200;vertices[2].y = 200;
    vertices[3].x = 100;vertices[3].y = 200;

    mpl_polygon_init(&rigid_body->polygon,vertices,4);

    rigid_body->active = 1;
    rigid_body->is_colliding = 0;
    rigid_body->is_static = 0;
    rigid_body->lock_orientation = 0;
    rigid_body->ignore_gravity = 0;
    rigid_body->restitution = 0.2f;
    rigid_body->static_frictionf =  0.5f;
    rigid_body->dynamic_frictionf = 0.3f;
    rigid_body->angularVelocity = 0;
    rigid_body->torque = 0;
    rigid_body->orient = 0;
    rigid_body->mass = 0;
    rigid_body->invMass = 0;
    rigid_body->inertia = 0;
    rigid_body->invInertia = 0;
    rigid_body->aprox_point_area = 1;
    rigid_body->longest_radius = 0;
    rigid_body->damping = 0.7f;
    rigid_body->velocity_epsilon = 0.001f;
    rigid_body->velocity.x = rigid_body->velocity.y =  0;
    rigid_body->force.x = rigid_body->force.y = 0;

    mpl_rigid_body_compute_mass(rigid_body,10);

}
void mpl_rigid_body_set_position(RigidBody *rigid_body,float x, float y)
{
    rigid_body->polygon.position.x = x;
    rigid_body->polygon.position.y = y;
    mpl_polygon_transform(&rigid_body->polygon);
}
void mpl_rigid_body_apply_force(RigidBody *rigid_body,Vector2 f)
{
    vector2f_add_vv(rigid_body->force, f, &rigid_body->force);
}

void mpl_rigid_body_apply_impulse(RigidBody *rigid_body,Vector2 impulse, Vector2 contactVector)
{    
    rigid_body->velocity.x += rigid_body->invMass * impulse.x;
    rigid_body->velocity.y += rigid_body->invMass * impulse.y;
    float n;
    vector2f_cross_vv(contactVector, impulse,&n);
    if (!rigid_body->lock_orientation) 
    {
        rigid_body->angularVelocity += rigid_body->invInertia * n;
    }
}
void mpl_rigid_body_integrate_velocity(RigidBody *rigid_body,float dt)
{
    rigid_body->polygon.position.x += rigid_body->velocity.x * rigid_body->damping * dt;
    rigid_body->polygon.position.y += rigid_body->velocity.y * rigid_body->damping * dt;
    if (!rigid_body->lock_orientation) 
    {
        rigid_body->polygon.orientation += rigid_body->angularVelocity * rigid_body->damping * dt;
    }
}
void mpl_rigid_body_integrate_forces(RigidBody *rigid_body,float dt)
{
    rigid_body->velocity.x += (rigid_body->force.x * rigid_body->invMass) * (dt);
    rigid_body->velocity.y += (rigid_body->force.y * rigid_body->invMass) * (dt);
    if(!rigid_body->lock_orientation)
    {
        rigid_body->angularVelocity += rigid_body->torque * rigid_body->invInertia * (dt);
    }
}
void mpl_rigid_body_clear_forces(RigidBody *rigid_body)
{
    rigid_body->force.x = 0;
    rigid_body->force.y = 0;
    
}
void update(RigidBody *bodies, unsigned int body_count, int iterations ,float dt, float G)
{
    bodies[0].force.y = 100000;
    for (int i = 0; i < body_count; ++i)
    {
        if (!bodies[i].is_static && bodies[i].active)
        {
            mpl_rigid_body_integrate_forces(&bodies[i],dt);
        }
    }  


        
    /*
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].EMPTY)
        {
            manifolds[i].Init(dt);
        }
    }
            
            for (int j = 0; j < iterations; ++j)
            {
                for (int i = 0; i < manifold_count; ++i)
                {
                    if (!manifolds[i].EMPTY)
                    {
                        manifolds[i].ApplyImpulse();
                    }
                }
            }
*/
            
    for (int i = 0; i < body_count; ++i)
    {
        if (!bodies[i].is_static && bodies[i].active)
        {
            mpl_rigid_body_integrate_velocity(&bodies[i],dt);
        }
    }
         /*   
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].EMPTY)
        {
            manifolds[i].PositionalCorrection();
        }
    }
    */
   printf("%f\n",bodies[0].force.y);
   printf("%f\n",bodies[0].velocity.y);

    for (int i = 0; i < body_count; ++i)
    {
        mpl_rigid_body_clear_forces(&bodies[i]);
    }
    for (int i = 0; i < body_count; ++i)
    {
        mpl_polygon_transform(&bodies[i].polygon);
    }
}