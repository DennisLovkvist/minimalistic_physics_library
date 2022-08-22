
#include <stdlib.h>
#include <stdio.h>
#include "mpl.h"
#include "delo_math.h"
#include <math.h>
#include "immintrin.h"
void mpl_polygon_get_centroid(float *vertices_x,float *vertices_y,unsigned int index, unsigned int vertex_count, float *centroid_x,float *centroid_y)
{
    *centroid_x = vertices_x[index];
    *centroid_y = vertices_y[index];
    
    for (int i = 1; i < vertex_count; i++)
    {
        *centroid_x += vertices_x[index+i];
        *centroid_y += vertices_y[index+i];
    }
    *centroid_x /= vertex_count;
    *centroid_y /= vertex_count;
}
void mpl_polygon_get_support_point(float *vertices_x,float *vertices_y,unsigned int index,unsigned int vertex_count, Vector2 dir,Vector2 *best_vertex)
{
    float best_projection = FLT_MIN;   
    best_vertex->x = 0;
    best_vertex->y = 0;

    for (int i = 0; i < vertex_count; ++i)
    {
        Vector2 v;
        v.x = vertices_x[index+i];
        v.y = vertices_y[index+i];
        float projection = vector2f_dot_vv(v.x,v.y, dir.x,dir.y);

        if (projection > best_projection)
        {
            best_vertex->x = v.x;
            best_vertex->y = v.y;
            best_projection = projection;
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
        float projection = vector2f_dot_vv(v.x,v.y, dir.x,dir.y);

        if (projection < bestProjection)
        {            
            best_vertex->x = v.x;
            best_vertex->y = v.y;
            bestProjection = projection;
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
void mpl_matrix_set(float *matrix,unsigned int index,float tx, float ty, float r)
{    
    
    matrix[index+(0 + 3*0)] = cos(r);
    matrix[index+(0 + 3*1)] = sin(r);
    matrix[index+(0 + 3*2)] = tx;

    matrix[index+(1 + 3*0)] = -sin(r);
    matrix[index+(1 + 3*1)] = cos(r);
    matrix[index+(1 + 3*2)] = ty;

    matrix[index+(2 + 3*0)] = 0;
    matrix[index+(2 + 3*1)] = 0;
    matrix[index+(2 + 3*2)] = 1;
}
void mpl_polygon_transform(RigidBody*rigid_body, unsigned int index)
{
    float tx = rigid_body->position_x[index] - rigid_body->position_last_x[index];
    float ty = rigid_body->position_y[index] - rigid_body->position_last_y[index];
    float r = rigid_body->orientation[index] - rigid_body->orientation_last[index];    
    
    for (int i = 0; i < rigid_body->vertex_count[index]; i++)
    {
        int v_index = (index * 4)+i;

        float normal_x = rigid_body->normals_x[v_index];
        float normal_y = rigid_body->normals_y[v_index];
        float position_x = rigid_body->position_x[index];
        float position_y = rigid_body->position_y[index];
        float position_last_x = rigid_body->position_last_x[index];
        float position_last_y = rigid_body->position_last_y[index];

        float input[3][1];
        float output[3][1];

        int m_index = index*(3*3);
        mpl_matrix_set(rigid_body->matrix,m_index,-position_last_x, -position_last_y, 0);
        input[0][0] = rigid_body->vertices_x[v_index];
        input[1][0] = rigid_body->vertices_y[v_index];
        input[2][0] = 1;
        matrix_mul_mm(rigid_body->matrix,m_index,input,output);

        rigid_body->vertices_x[v_index] = output[0][0];
        rigid_body->vertices_y[v_index] = output[1][0];
        
        mpl_matrix_set(rigid_body->matrix,m_index,0, 0, r);
        input[0][0] = rigid_body->vertices_x[v_index];
        input[1][0] = rigid_body->vertices_y[v_index];
        input[2][0] = 1;
        matrix_mul_mm(rigid_body->matrix,m_index,input,output);
        rigid_body->vertices_x[v_index] = output[0][0];
        rigid_body->vertices_y[v_index] = output[1][0];

        input[0][0] = normal_x;
        input[1][0] = normal_y;
        input[2][0] = 1;
        matrix_mul_mm(rigid_body->matrix,m_index,input,output);
        rigid_body->normals_x[v_index] = output[0][0];
        rigid_body->normals_y[v_index] = output[1][0];

        mpl_matrix_set(rigid_body->matrix,m_index,rigid_body->position_last_x[index] + tx, rigid_body->position_last_y[index] + ty, 0);
        input[0][0] = rigid_body->vertices_x[v_index];
        input[1][0] = rigid_body->vertices_y[v_index];
        input[2][0] = 1;
        matrix_mul_mm(rigid_body->matrix,m_index,input,output);
        rigid_body->vertices_x[v_index] = output[0][0];
        rigid_body->vertices_y[v_index] = output[1][0];
    }
    rigid_body->position_last_x[index] = rigid_body->position_x[index];
    rigid_body->position_last_y[index] = rigid_body->position_y[index];
    rigid_body->orientation_last[index] = rigid_body->orientation[index];
}
void mpl_polygon_init(RigidBody *rigid_body, unsigned int index,Vector2 *vertices, unsigned int vertex_count)
{
    rigid_body->orientation[index] = rigid_body->orientation_last[index] = 0.0;
    rigid_body->radius[index] = 0;
    vertex_count = (vertex_count > 4) ? 4:vertex_count;
    rigid_body->vertex_count[index] = vertex_count;

    for (int i = 0; i < vertex_count; i++)
    {
        int v_index = (index * 4)+i;
        rigid_body->vertices_x[v_index] = vertices[i].x;
        rigid_body->vertices_y[v_index] = vertices[i].y;
    }    
    for (int i = 0; i < vertex_count; i++)
    {
        int v_index = (index * 4)+i;
        rigid_body->edge_index_a[v_index] = i;
        rigid_body->edge_index_b[v_index] = (i + 1 < vertex_count) ? i + 1 : 0;  

        int edge_index = rigid_body->edge_index_b[v_index];
        int v_edge_index = (index * 4)+edge_index;
        Vector2 v0 = 
        {
            rigid_body->vertices_x[v_edge_index],
            rigid_body->vertices_y[v_edge_index]
        };
        edge_index = rigid_body->edge_index_a[v_index];
        v_edge_index = (index * 4)+edge_index;
        Vector2 v1 = 
        {
            rigid_body->vertices_x[v_edge_index],
            rigid_body->vertices_y[v_edge_index]
        };

        vector2f_normal(v0.x,v0.y,v1.x,v1.y,&rigid_body->normals_x[v_index],&rigid_body->normals_y[v_index]);

        vector2f_normalize(&rigid_body->normals_x[v_index],&rigid_body->normals_y[v_index]);
    }
    mpl_polygon_get_centroid(rigid_body->vertices_x,rigid_body->vertices_y,index * 4,vertex_count, &rigid_body->position_x[index],&rigid_body->position_y[index]);
    
    rigid_body->position_last_x[index] = rigid_body->position_x[index];
    rigid_body->position_last_y[index] = rigid_body->position_y[index];

    mpl_polygon_transform(rigid_body,index);
    
}
void mpl_rigid_body_compute_mass(RigidBody *rigid_body,unsigned int index, float density)
{
    float area = 0;
    unsigned int vertex_count = rigid_body->vertex_count[index];
    Vector2 c;
    mpl_polygon_get_centroid(rigid_body->vertices_x,rigid_body->vertices_y,index*4,vertex_count,&c.x,&c.y);
    float dx = 0;
    float dy = 0;
    float longest_radius = 0;
    float radius = 0;

    for (int i = 0; i < vertex_count; i++)
    {
        int v_index = (index * 4)+i;
        Vector2 v = 
        {
            rigid_body->vertices_x[v_index],
            rigid_body->vertices_y[v_index],
        };
        int next = mpl_polygon_get_next_vertex_index(i,vertex_count);
        dx = c.x - v.x;
        dy = c.y - v.y;
        float d1 = sqrt(dx * dx + dy * dy);
        radius += d1;
        dx = rigid_body->vertices_x[v_index+next] - v.x;
        dy = rigid_body->vertices_y[v_index+next] - v.y;
        float d2 = sqrt(dx * dx + dy * dy);
        area += ((d1 * d2) * 0.5f);

        if (radius > longest_radius)
        {
            longest_radius = radius;
        }
    }
    radius /= vertex_count;

    rigid_body->radius[index] = radius;
    rigid_body->longest_radius[index] = longest_radius;

    rigid_body->mass[index] = density * area;
    rigid_body->inertia[index] = rigid_body->mass[index] * radius * radius;
    rigid_body->inverse_mass[index] = 1.0f / rigid_body->mass[index];
    rigid_body->inverse_inertia[index] = 1.0f / rigid_body->inertia[index];
    
}
void mpl_rigid_body_compute_mass_circle(RigidBody *rigid_body,unsigned int index,float density,float radius)
{
    float area = PI * radius*radius;

    rigid_body->radius[index] = radius;
    rigid_body->longest_radius[index] = radius;

    rigid_body->mass[index] = density * area;
    rigid_body->inertia[index] = rigid_body->mass[index] * radius * radius;
    rigid_body->inverse_mass[index] = 1.0f / rigid_body->mass[index];
    rigid_body->inverse_inertia[index] = 1.0f / rigid_body->inertia[index];
}
void mpl_polygon_make_square(RigidBody *rigid_body,unsigned int index,unsigned int width,unsigned int height)
{
    rigid_body->vertex_count[index] = 4;
    Vector2 vertices[4];
    float half_width = (float)width*0.5;
    float half_height = (float)height*0.5;

    vertices[0].x = -half_width;vertices[0].y = -half_height;
    vertices[1].x = half_width;vertices[1].y = -half_height;
    vertices[2].x = half_width;vertices[2].y = half_height;
    vertices[3].x = -half_width;vertices[3].y = half_height;
    mpl_polygon_init(rigid_body,index,vertices,4);
}
void mpl_polygon_make_triangle(RigidBody *rigid_body,unsigned int index, unsigned int width,unsigned int height)
{
    rigid_body->vertex_count[index] = 3;
    Vector2 vertices[3];
    float half_width = (float)width*0.5;
    float half_height = (float)height*0.5;

    vertices[0].x = 0;vertices[0].y = -half_height;
    vertices[1].x = half_width;vertices[1].y = half_height;
    vertices[2].x = -half_width;vertices[2].y = half_height;

    mpl_polygon_init(rigid_body,index,vertices,3);
}
void mpl_rigid_bodies_init(RigidBody *rigid_body, unsigned int capacity)
{
    rigid_body->is_polygon = malloc(sizeof(unsigned int)*capacity);
    rigid_body->is_colliding = malloc(sizeof(unsigned int)*capacity);
    rigid_body->is_static = malloc(sizeof(unsigned int)*capacity);
    rigid_body->active = malloc(sizeof(unsigned int)*capacity);
    rigid_body->lock_orientation = malloc(sizeof(unsigned int)*capacity);
    rigid_body->ignore_gravity = malloc(sizeof(unsigned int)*capacity);
    rigid_body->restitution = malloc(sizeof(float)*capacity);
    rigid_body->static_friction = malloc(sizeof(float)*capacity);
    rigid_body->dynamic_friction = malloc(sizeof(float)*capacity);
    rigid_body->velocity_x = malloc(sizeof(float)*capacity);
    rigid_body->velocity_y = malloc(sizeof(float)*capacity);
    rigid_body->force_x = malloc(sizeof(float)*capacity);
    rigid_body->force_y = malloc(sizeof(float)*capacity);
    rigid_body->angular_velocity = malloc(sizeof(float)*capacity);
    rigid_body->torque = malloc(sizeof(float)*capacity);
    rigid_body->orient = malloc(sizeof(float)*capacity);
    rigid_body->mass = malloc(sizeof(float)*capacity);
    rigid_body->inverse_mass = malloc(sizeof(float)*capacity);
    rigid_body->inertia = malloc(sizeof(float)*capacity);
    rigid_body->inverse_inertia = malloc(sizeof(float)*capacity);
    rigid_body->aprox_point_area = malloc(sizeof(float)*capacity);
    rigid_body->longest_radius = malloc(sizeof(float)*capacity);
    rigid_body->damping = malloc(sizeof(float)*capacity);
    rigid_body->velocity_epsilon = malloc(sizeof(float)*capacity);

    rigid_body->matrix = malloc(sizeof(float)*capacity*3*3);
    rigid_body->vertex_count = malloc(sizeof(unsigned int)*capacity);
    rigid_body->vertices_x = malloc(sizeof(float)*capacity*4);
    rigid_body->vertices_y = malloc(sizeof(float)*capacity*4);

    rigid_body->edge_vertex_a_x = malloc(sizeof(float)*capacity*4);
    rigid_body->edge_vertex_b_y = malloc(sizeof(float)*capacity*4);
    rigid_body->edge_index_a = malloc(sizeof(int)*capacity*4);
    rigid_body->edge_index_b = malloc(sizeof(int)*capacity*4);

    rigid_body->normals_x = malloc(sizeof(float)*capacity*4);
    rigid_body->normals_y = malloc(sizeof(float)*capacity*4);

    rigid_body->position_x = malloc(sizeof(float)*capacity);
    rigid_body->position_y = malloc(sizeof(float)*capacity);
    rigid_body->position_last_x = malloc(sizeof(float)*capacity);
    rigid_body->position_last_y = malloc(sizeof(float)*capacity);

    rigid_body->orientation = malloc(sizeof(float)*capacity);
    rigid_body->orientation_last = malloc(sizeof(float)*capacity);
    rigid_body->radius = malloc(sizeof(float)*capacity);

    rigid_body->contact_points_x = malloc(sizeof(float)*capacity*4);
    rigid_body->contact_points_y = malloc(sizeof(float)*capacity*4);
}
void mpl_rigid_body_init(RigidBody *rigid_body, unsigned int index,unsigned int shape,unsigned int width,unsigned int height)
{    
    if(shape == 0)
    {
        rigid_body->is_polygon[index] = 0;
        rigid_body->vertex_count[index] = 0;
    }
    else if(shape == 1)
    {
        rigid_body->is_polygon[index] = 1;
        mpl_polygon_make_triangle(rigid_body,index,width, height);        
    }
    else if(shape == 2)
    {
        rigid_body->is_polygon[index] = 1;
        mpl_polygon_make_square(rigid_body,index,width, height);
    }
    else
    {
        rigid_body->is_polygon[index] = 1;
        mpl_polygon_make_square(rigid_body,index,width, height);
    }
    rigid_body->active[index] = 1;
    rigid_body->is_colliding[index] = 0;
    rigid_body->is_static[index] = 0;
    rigid_body->lock_orientation[index] = 0;
    rigid_body->ignore_gravity[index] = 0;
    rigid_body->restitution[index] = 0.2f;
    rigid_body->static_friction[index] =  0.4f;
    rigid_body->dynamic_friction[index] = 0.2f;
    rigid_body->angular_velocity[index] = 0;
    rigid_body->torque[index] = 0;
    rigid_body->orient[index] = 0;
    rigid_body->mass[index] = 0;
    rigid_body->inverse_mass[index] = 0;
    rigid_body->inertia[index] = 0;
    rigid_body->inverse_inertia[index] = 0;
    rigid_body->aprox_point_area[index] = 1;
    rigid_body->longest_radius[index] = 0;
    rigid_body->damping[index] = 0.7f;
    rigid_body->velocity_epsilon[index] = 0.001f;
    rigid_body->velocity_x[index] = rigid_body->velocity_y[index] =  0;
    rigid_body->force_x[index] = rigid_body->force_y[index] =  0;
    rigid_body->contact_points_x[index*4+0] = rigid_body->contact_points_y[index*4+0] =  0;
    rigid_body->contact_points_x[index*4+1] = rigid_body->contact_points_y[index*4+1] =  0;
    rigid_body->contact_points_x[index*4+2] = rigid_body->contact_points_y[index*4+2] =  0;
    rigid_body->contact_points_x[index*4+3] = rigid_body->contact_points_y[index*4+3] =  0;


    if(rigid_body->is_polygon[index])
    {
        mpl_rigid_body_compute_mass(rigid_body,index,0.001);
    }
    else
    {
        mpl_rigid_body_compute_mass_circle(rigid_body,index,0.001,width/2);
    }

}
void mpl_rigid_body_set_position(RigidBody *rigid_body,unsigned int index,float x, float y)
{
    rigid_body->position_x[index] = x;
    rigid_body->position_y[index] = y;
    mpl_polygon_transform(rigid_body,index);
}
void mpl_rigid_body_apply_force(RigidBody *rigid_body, unsigned int index,Vector2 f)
{
    rigid_body->force_x[index] += f.x;
    rigid_body->force_y[index] += f.y;
}

void mpl_rigid_body_apply_impulse(RigidBody *rigid_body,unsigned int index,Vector2 impulse, Vector2 contactVector)
{    
    if(isnan(impulse.x) || isnan(impulse.y))
    {
        return;
    }
    rigid_body->velocity_x[index] += rigid_body->inverse_mass[index] * impulse.x;
    rigid_body->velocity_y[index] += rigid_body->inverse_mass[index] * impulse.y;
    
    float n = vector2f_cross_vv(contactVector, impulse);
    if (!rigid_body->lock_orientation[index]) 
    {
        rigid_body->angular_velocity[index] += rigid_body->inverse_inertia[index] * n;
    }
}
void mpl_rigid_body_integrate_velocity(RigidBody *rigid_body,unsigned int index,float dt)
{
    rigid_body->position_x[index] += rigid_body->velocity_x[index] * rigid_body->damping[index] * dt;
    rigid_body->position_y[index] += rigid_body->velocity_y[index] * rigid_body->damping[index] * dt;

    if (!rigid_body->lock_orientation) 
    {
        rigid_body->orientation[index] += rigid_body->angular_velocity[index] * rigid_body->damping[index] * dt;
    }
}
void mpl_rigid_body_integrate_forces(RigidBody *rigid_body,unsigned int index,float dt)
{
    rigid_body->velocity_x[index] += (rigid_body->force_x[index] * rigid_body->inverse_mass[index]) * (dt);
    rigid_body->velocity_y[index] += (rigid_body->force_y[index] * rigid_body->inverse_mass[index]) * (dt);

    if(!rigid_body->lock_orientation[index])
    {
        rigid_body->angular_velocity[index] += rigid_body->torque[index] * rigid_body->inverse_inertia[index] * (dt);
    }
}
void mpl_rigid_body_clear_forces(RigidBody *rigid_body,unsigned int index)
{
    rigid_body->force_x[index] = 0;
    rigid_body->force_y[index] = 0;    
}
void mpl_manifold_init(Manifold *manifold,RigidBody *rigid_body,unsigned int index,float dt)
{
    
    manifold->e = min(rigid_body->restitution[manifold->rigid_body_a], rigid_body->restitution[manifold->rigid_body_b]);
    manifold->static_friction = sqrt(rigid_body->static_friction[manifold->rigid_body_a] * rigid_body->static_friction[manifold->rigid_body_b]);
    manifold->dynamic_friction = sqrt(rigid_body->dynamic_friction[manifold->rigid_body_a] * rigid_body->dynamic_friction[manifold->rigid_body_b]);
    
}
void mpl_manifold_apply_impulse(Manifold *manifold,RigidBody *rigid_body)
{

    unsigned int rigid_body_a = manifold->rigid_body_a;
    unsigned int rigid_body_b = manifold->rigid_body_b;

    for (int i = 0; i < manifold->contact_count; ++i)
    {
        Vector2 ra = 
        {
            rigid_body->position_x[rigid_body_a]-manifold->contacts[i].x,
            rigid_body->position_y[rigid_body_a]-manifold->contacts[i].y
        };

        Vector2 rb = 
        {
            rb.x = rigid_body->position_x[rigid_body_b]-manifold->contacts[i].x,
            rb.y = rigid_body->position_y[rigid_body_b]-manifold->contacts[i].y
        };        
        
        Vector2 cross_rb_av;
        vector2f_cross_vf(rb,rigid_body->angular_velocity[rigid_body_b],&cross_rb_av);
        Vector2 cross_ra_av;
        vector2f_cross_vf(ra,rigid_body->angular_velocity[rigid_body_a],&cross_ra_av);
        
        Vector2 rv;
        rv.x = rigid_body->velocity_x[rigid_body_b] + cross_rb_av.x - rigid_body->velocity_x[rigid_body_a] - cross_ra_av.x;
        rv.y = rigid_body->velocity_y[rigid_body_b] + cross_rb_av.y - rigid_body->velocity_y[rigid_body_a] - cross_ra_av.y;

        float contact_velocity = vector2f_dot_vv(rv.x,rv.y, manifold->normal.x,manifold->normal.y);
        
        if (contact_velocity > 0)
        {
            return;
        }

        float raCrossN = vector2f_cross_vv(ra, manifold->normal);
        float rbCrossN = vector2f_cross_vv(rb, manifold->normal);

        float invMassSum = rigid_body->inverse_mass[rigid_body_a] + rigid_body->inverse_mass[rigid_body_b] + (raCrossN * raCrossN) * rigid_body->inverse_inertia[rigid_body_a] + (rbCrossN * rbCrossN) * rigid_body->inverse_inertia[rigid_body_b];
        
        float j = -(1.0f + manifold->e) * contact_velocity;                
        j /= invMassSum;
        j /= manifold->contact_count;


        Vector2 impulse;
        impulse.x = manifold->normal.x * j;
        impulse.y = manifold->normal.y * j;

        mpl_rigid_body_apply_impulse(rigid_body,rigid_body_b,impulse,rb);
        impulse.x = -impulse.x;
        impulse.y = -impulse.y;
        mpl_rigid_body_apply_impulse(rigid_body,rigid_body_a,impulse,ra);

        float dot = vector2f_dot_vv(rv.x,rv.y, manifold->normal.x,manifold->normal.y);
        Vector2 t;
        t.x = rv.x - (manifold->normal.x * dot);
        t.y = rv.y - (manifold->normal.y * dot);
        vector2f_normalize(&t.x,&t.y);

        float jt = -vector2f_dot_vv(rv.x,rv.y, t.x,t.y);
        jt /= invMassSum;
        jt /= manifold->contact_count;
        
        if (jt == 0.0f)
        {
            return;
        }

        Vector2 tangent_impulse;

        if (abs(jt) < j * manifold->static_friction)
        {
            tangent_impulse.x = t.x * jt;
            tangent_impulse.y = t.y * jt;
        }
        else
        {
            tangent_impulse.x = t.x * -j * manifold->dynamic_friction;
            tangent_impulse.y = t.y * -j * manifold->dynamic_friction;
        }
        mpl_rigid_body_apply_impulse(rigid_body,rigid_body_b,tangent_impulse,rb);
        tangent_impulse.x = -tangent_impulse.x;
        tangent_impulse.y = -tangent_impulse.y;
        mpl_rigid_body_apply_impulse(rigid_body,rigid_body_a,tangent_impulse,ra);    
    }
}
void mpl_manifold_positional_correction(Manifold *manifold,RigidBody *rigid_body)
{
    unsigned int rigid_body_a = manifold->rigid_body_a;
    unsigned int rigid_body_b = manifold->rigid_body_b;

    const float k_slop = 0.05f;
    const float percent = 0.4f;

    float f = max(manifold->penetration - k_slop, 0.0f)/(rigid_body->inverse_mass[rigid_body_a] + rigid_body->inverse_mass[rigid_body_b]);
    Vector2 correction;
    correction.x = f * manifold->normal.x * percent;
    correction.y = f * manifold->normal.y * percent;
    
    rigid_body->position_x[rigid_body_a] -= correction.x * rigid_body->inverse_mass[rigid_body_a];
    rigid_body->position_y[rigid_body_a] -= correction.y * rigid_body->inverse_mass[rigid_body_a];

    rigid_body->position_x[rigid_body_b] += correction.x * rigid_body->inverse_mass[rigid_body_b];
    rigid_body->position_y[rigid_body_b] += correction.y * rigid_body->inverse_mass[rigid_body_b];
}
void mpl_manifold_reset(Manifold *manifold)
{
    manifold->contact_count = 0;
    manifold->penetration = 0;
}
float mpl_find_axis_least_penetration(RigidBody *rigid_bodies, unsigned int a, unsigned int b,int *face)
{
    float best_distance = FLT_MIN;
    int best_index = 0;
    for (int i = 0; i < rigid_bodies->vertex_count[a]; ++i)
    {
        int v_index_a = a*4;
        int v_index_b = b*4;

        Vector2 normal = {rigid_bodies->normals_x[v_index_a+i],rigid_bodies->normals_y[v_index_a+i]};
        Vector2 normal_flipped = {-normal.x,-normal.y};

        Vector2 support_point = {0,0};
        mpl_polygon_get_support_point(rigid_bodies->vertices_x,rigid_bodies->vertices_y,v_index_b,rigid_bodies->vertex_count[b],normal_flipped,&support_point);

    
        Vector2 vertex = {rigid_bodies->vertices_x[v_index_a+i],rigid_bodies->vertices_y[v_index_a+i]};

        Vector2 diff = {support_point.x - vertex.x,support_point.y - vertex.y};

        float dot = vector2f_dot_vv(normal.x,normal.y, diff.x,diff.y);

        if (dot > best_distance)
        {
            best_distance = dot;
            best_index = i;
        }
    }
    *face = best_index;
    return best_distance;
}

unsigned int mpl_gt(float a, float b)
{
    return a >= b * 0.95f + a * 0.0f;
}
int mpl_find_incident_face(RigidBody *rigid_bodies, unsigned int reference, unsigned int incident, int reference_index)
{
    Vector2 reference_normal;

    int v_reference = reference*4;
    int v_incident = incident*4;
    
    reference_normal.x = rigid_bodies->normals_x[v_reference+reference_index];
    reference_normal.y = rigid_bodies->normals_y[v_reference+reference_index];      

    matrix_mul_mv2(rigid_bodies->matrix,reference*3*3,reference_normal,&reference_normal);

    float transposed[3][3];
    matrix_transpose_mm(rigid_bodies->matrix,reference*3*3,transposed);
    matrix_mul_mv(transposed,reference_normal,&reference_normal);

    int incidentFace = 0;
    float minDot = FLT_MAX;
    for (int i = 0; i < rigid_bodies->vertex_count[incident]; ++i)
    {        
        float dot = vector2f_dot_vv(reference_normal.x,reference_normal.y,rigid_bodies->normals_x[v_incident+i],rigid_bodies->normals_y[v_incident+i]);

        if (dot < minDot)
        {
            minDot = dot;
            incidentFace = i;
        }
    }
    return incidentFace;
}
int mpl_clip(Vector2 n, float c, Vector2 *face)
{
    int sp = 0;
    Vector2 points[2];
    points[0].x = face[0].x;
    points[0].y = face[0].y;
    points[1].x = face[1].x;
    points[1].y = face[1].y;

    float d1 = vector2f_dot_vv(n.x,n.y, face[0].x,face[0].y)-c;
    float d2 = vector2f_dot_vv(n.x,n.y, face[1].x,face[1].y)-c;

    if (d1 <= 0.0f) points[sp++] = face[0];
    if (d2 <= 0.0f) points[sp++] = face[1];

    // If the points are on different sides of the plane
    if (d1 * d2 < 0.0f)
    {
        float alpha = d1 / (d1 - d2);
        points[sp].x = face[0].x  + alpha * (face[1].x  - face[0].x );
        points[sp].y = face[0].y  + alpha * (face[1].y  - face[0].y );
        ++sp;
    }
    face[0] = points[0];
    face[1] = points[1];

    return sp;
}
unsigned int mpl_broad_phase(RigidBody *rigid_bodies,unsigned int a, unsigned int b)
{
    float dx = rigid_bodies->position_x[a] - rigid_bodies->position_x[b];
    float dy = rigid_bodies->position_y[a] - rigid_bodies->position_y[b];
    float dist = sqrt(dx * dx + dy * dy);
    return (dist < rigid_bodies->longest_radius[a] + rigid_bodies->longest_radius[b]);
}
float mpl_dist_sqr(float ax,float ay,float bx,float by)
{
    float cx = ax-bx;
    float cy = ay-by;
    return vector2f_dot_vv(cx,cy,cx,cy);
}
unsigned int mpl_narrow_phase_c2c(Manifold *manifold, RigidBody *rigid_bodies,unsigned int a, unsigned int b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = a;
    manifold->rigid_body_b = b;
    Vector2 normal = {rigid_bodies->position_x[b]-rigid_bodies->position_x[a],rigid_bodies->position_y[b]-rigid_bodies->position_y[a]};

    float dist_sqr = normal.x * normal.x + normal.y * normal.y;
    float radius = rigid_bodies->longest_radius[a] + rigid_bodies->longest_radius[b];

    // Not in contact
    if (dist_sqr >= radius * radius)
    {
        manifold->contact_count = 0;
        manifold->is_empty = 1;
        return 0;
    }

    float distance = sqrt(dist_sqr);

    manifold->contact_count = 1;

    if (distance == 0.0f)
    {
        manifold->penetration = rigid_bodies->longest_radius[a];
        manifold->normal.x = 1;
        manifold->normal.y = 0;
        manifold->contacts[0].x = rigid_bodies->position_x[a];
        manifold->contacts[0].y = rigid_bodies->position_y[a];
    }
    else
    {
        manifold->penetration = radius - distance;
        manifold->normal.x = normal.x / distance; // Faster than using Normalized since we already performed sqrt
        manifold->normal.y = normal.y / distance;
        manifold->contacts[0].x = manifold->normal.x * rigid_bodies->longest_radius[a] + rigid_bodies->position_x[a];
        manifold->contacts[0].y = manifold->normal.y * rigid_bodies->longest_radius[a] + rigid_bodies->position_y[a];
    }

    return 1;
}
unsigned int mpl_narrow_phase_p2c(Manifold *manifold, RigidBody *rigid_bodies, unsigned int a,unsigned int b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = a;
    manifold->rigid_body_b = b;

   Vector2 normal;
   vector2f_normal(rigid_bodies->position_x[b],rigid_bodies->position_y[b], rigid_bodies->position_x[a],rigid_bodies->position_y[a], &normal.x,&normal.y);

    vector2f_normalize(&normal.x,&normal.y);
    vector2f_cross_v(normal,&normal);

    Vector2 normal_flipped = {-normal.x,-normal.y};

    float separation = FLT_MIN;
    int face_normal = 0;

    int v_index_a = a*4;

    for (int i = 0; i < rigid_bodies->vertex_count[a]; i++)
    {
        Vector2 delta = {delta.x = rigid_bodies->position_x[b] - rigid_bodies->position_x[b],delta.y = rigid_bodies->position_y[a] - rigid_bodies->position_y[a]};
        
        float s = vector2f_dot_vv(rigid_bodies->normals_x[v_index_a+i],rigid_bodies->normals_x[v_index_a+i],delta.x,delta.y);

        if (s > rigid_bodies->radius[b])
        {
            manifold->is_empty = 1;
            return 0;
        }
        if (s > separation)
        {
            separation = s;
            face_normal = i;
        }

    }
    int edge_index_a = rigid_bodies->edge_index_a[v_index_a+face_normal];
    int edge_index_b = rigid_bodies->edge_index_b[v_index_a+face_normal];

    Vector2 v1 = {rigid_bodies->vertices_x[v_index_a+edge_index_a],rigid_bodies->vertices_y[v_index_a+edge_index_a]};
    Vector2 v2 = {rigid_bodies->vertices_x[v_index_a+edge_index_b],rigid_bodies->vertices_y[v_index_a+edge_index_b]};

    Vector2 delta1 = {delta1.x = rigid_bodies->position_x[b] - v1.x,delta1.y = rigid_bodies->position_y[b] - v1.y};
    
    Vector2 delta2 = {v2.x - v1.x,v2.y - v1.y};

    float dot1 = vector2f_dot_vv(delta1.x,delta1.y, delta2.x,delta2.y);

    delta1.x = rigid_bodies->position_x[b] - v2.x;
    delta1.y = rigid_bodies->position_y[b] - v2.y;

    delta2.x = v1.x - v2.x;
    delta2.y = v1.y - v2.y;    
    float dot2 = vector2f_dot_vv(delta1.x,delta1.y, delta2.x,delta2.y);

    manifold->contact_count = 1;
    manifold->penetration = rigid_bodies->radius[b] - separation;

    if (separation < 0.0001f)
    {
        manifold->contact_count = 1;
        manifold->normal.x = rigid_bodies->normals_x[v_index_a+face_normal];
        manifold->normal.y = rigid_bodies->normals_y[v_index_a+face_normal];
        manifold->contacts[0].x = manifold->normal.x * rigid_bodies->radius[b] + rigid_bodies->position_x[b];
        manifold->contacts[0].y = manifold->normal.y * rigid_bodies->radius[b] + rigid_bodies->position_y[b];

        manifold->penetration = rigid_bodies->radius[b];
        return 0;
    }
    if (dot1 <= 0.0f)
    {
        if (mpl_dist_sqr(rigid_bodies->position_x[b],rigid_bodies->position_y[b], v1.x,v1.y) > rigid_bodies->radius[b] * rigid_bodies->radius[b])
        {
            manifold->is_empty = 1;
            return 0;
        }
        manifold->contacts[0] = v1;
        float r = atan2(rigid_bodies->position_y[b] - v1.y, rigid_bodies->position_x[b] - v1.x);
        Vector2 n;
        n.x = cos(r);
        n.y = sin(r);
        vector2f_normalize(&n.x,&n.y);

        manifold->normal = n;
        Vector2 surface;
        surface.x = rigid_bodies->position_x[b] + cos(r - PI) * rigid_bodies->radius[b];
        surface.y = rigid_bodies->position_y[b] + sin(r - PI) * rigid_bodies->radius[b];

        float dx = surface.x - v1.x;
        float dy = surface.y - v1.y;
        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else if (dot2 <= 0.0f)
    {
        if (mpl_dist_sqr(rigid_bodies->position_x[b],rigid_bodies->position_y[b], v2.x,v2.y) > rigid_bodies->radius[b] * rigid_bodies->radius[b])
        {
            manifold->is_empty = 1;
            return 0;
        }
        manifold->contacts[0] = v2;
        float r = atan2(rigid_bodies->position_y[b] - v2.y, rigid_bodies->position_x[b] - v2.x);
        Vector2 n = {cos(r),sin(r)};
        vector2f_normalize(&n.x,&n.y);
        manifold->normal = n;
        Vector2 surface = {rigid_bodies->position_x[b] + cos(r - PI) * rigid_bodies->radius[b],rigid_bodies->position_y[b] + sin(r - PI) * rigid_bodies->radius[b]};

        float dx = surface.x - v2.x;
        float dy = surface.y - v2.y;
        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else
    {
        manifold->normal.x = rigid_bodies->normals_x[v_index_a+face_normal];
        manifold->normal.y = rigid_bodies->normals_y[v_index_a+face_normal];
        manifold->contacts[0].x = rigid_bodies->position_x[b] - rigid_bodies->normals_x[v_index_a+face_normal] * (rigid_bodies->radius[b] - manifold->penetration);
        manifold->contacts[0].y = rigid_bodies->position_y[b] - rigid_bodies->normals_y[v_index_a+face_normal] * (rigid_bodies->radius[b] - manifold->penetration);
        
    }

    int v_index_b = manifold->rigid_body_b * 4;

    rigid_bodies->contact_points_x[v_index_b+0] = rigid_bodies->position_x[b];
    rigid_bodies->contact_points_y[v_index_b+0] = rigid_bodies->position_y[b];
    rigid_bodies->contact_points_x[v_index_b+1] = manifold->contacts[0].x;
    rigid_bodies->contact_points_y[v_index_b+1] = manifold->contacts[0].y;

    rigid_bodies->contact_points_x[v_index_b+2] = v1.x;
    rigid_bodies->contact_points_y[v_index_b+2] = v1.y;
    rigid_bodies->contact_points_x[v_index_b+3] = v2.x;
    rigid_bodies->contact_points_y[v_index_b+3] = v2.y;

    rigid_bodies->is_colliding[a] = 1;
    rigid_bodies->is_colliding[b] = 1;
    return 1;
}
unsigned int mpl_narrow_phase_p2p(Manifold *manifold, RigidBody *rigid_bodies, unsigned int a,unsigned int b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = a;
    manifold->rigid_body_b = b;

    manifold->contact_count = 0;

    int face_a = 0;
    int face_b = 0;    
    
    float penetration_a = mpl_find_axis_least_penetration(rigid_bodies, a,b,&face_a);
    float penetration_b = mpl_find_axis_least_penetration(rigid_bodies, b,a,&face_b);
  

    if (penetration_a >= 0.0f)
    {
        manifold->is_empty = 1;
        return 0;
    }
    if (penetration_b >= 0.0f)
    {
        manifold->is_empty = 1;
        return 0;
    }

    int reference_face_index = 0;
    unsigned int flip = 0; // Always point from a to b

    unsigned int *reference_polygon = &manifold->reference_polygon;
    unsigned int *incident_polygon = &manifold->incident_polygon;

    // Determine which shape contains reference face
    if (mpl_gt(penetration_a, penetration_b))
    {
        *reference_polygon = a;
        *incident_polygon = b;
        reference_face_index = face_a;
        flip = 0;
    }
    else
    {
        *reference_polygon = b;
        *incident_polygon = a;
        reference_face_index = face_b;
        flip = 1;
    }                           

    Vector2 reference_face[2];
    int v_reference = (*reference_polygon)*4;
    int index_a = rigid_bodies->edge_index_a[v_reference+reference_face_index];
    int index_b = rigid_bodies->edge_index_b[v_reference+reference_face_index]; 
              
    reference_face[0].x = rigid_bodies->vertices_x[v_reference+index_a];
    reference_face[0].y = rigid_bodies->vertices_y[v_reference+index_a];

    reference_face[1].x = rigid_bodies->vertices_x[v_reference+index_b];
    reference_face[1].y = rigid_bodies->vertices_y[v_reference+index_b];

    int incident_face_index = mpl_find_incident_face(rigid_bodies,*reference_polygon, *incident_polygon, reference_face_index);
    
    Vector2 incident_face[2];
    int v_incident = *incident_polygon*4;
    index_a = rigid_bodies->edge_index_a[v_incident + incident_face_index];
    index_b = rigid_bodies->edge_index_b[v_incident + incident_face_index];  
    incident_face[0].x = rigid_bodies->vertices_x[v_incident + index_a];
    incident_face[0].y = rigid_bodies->vertices_y[v_incident + index_a];
    incident_face[1].x = rigid_bodies->vertices_x[v_incident + index_b];
    incident_face[1].y = rigid_bodies->vertices_y[v_incident + index_b];


    // Calculate normals
    Vector2 side_plane_normal;
    side_plane_normal.x = reference_face[1].x - reference_face[0].x;
    side_plane_normal.y = reference_face[1].y - reference_face[0].y;
    vector2f_normalize(&side_plane_normal.x,&side_plane_normal.y);

    Vector2 reference_face_normal;
    reference_face_normal.x = side_plane_normal.y;
    reference_face_normal.y = -side_plane_normal.x;       
    
    float dist_from_origin = vector2f_dot_vv(reference_face_normal.x,reference_face_normal.y, reference_face[0].x,reference_face[0].y);
    float negative_side = -vector2f_dot_vv(side_plane_normal.x,side_plane_normal.y, reference_face[0].x,reference_face[0].y);
    float positive_side = vector2f_dot_vv(side_plane_normal.x,side_plane_normal.y, reference_face[1].x,reference_face[1].y);

    manifold->normal.x = reference_face_normal.x;
    manifold->normal.y = reference_face_normal.y;
    // Clip incident face to reference face side planes
    Vector2 side_plane_normal_negated = {-side_plane_normal.x,-side_plane_normal.y};

    if (mpl_clip(side_plane_normal_negated, negative_side, &incident_face[0]) < 2)
    {
        //return 0;
    }
    if (mpl_clip(side_plane_normal, positive_side, &incident_face[0]) < 2)
    {
        //return 0;
    }
    if (flip)
    {
        manifold->normal.x = -manifold->normal.x;
        manifold->normal.y = -manifold->normal.y;
    }


    int v_index_a = a*4;   
    rigid_bodies->contact_points_x[v_index_a+0] = incident_face[0].x;
    rigid_bodies->contact_points_y[v_index_a+0] = incident_face[0].y;
    rigid_bodies->contact_points_x[v_index_a+1] = incident_face[1].x;
    rigid_bodies->contact_points_y[v_index_a+1] = incident_face[1].y;

    rigid_bodies->contact_points_x[v_index_a+2] = reference_face[0].x;
    rigid_bodies->contact_points_y[v_index_a+2] = reference_face[0].y;
    rigid_bodies->contact_points_x[v_index_a+3] = reference_face[1].x;
    rigid_bodies->contact_points_y[v_index_a+3] = reference_face[1].y;

    // Keep points behind reference face
    int contact_count = 0;
    float separation = vector2f_dot_vv(reference_face_normal.x,reference_face_normal.y, incident_face[0].x,incident_face[0].y) - dist_from_origin;

    if (separation <= 0.0f)
    {
        manifold->contacts[contact_count] = incident_face[0];
        manifold->penetration = -separation;
        ++contact_count;
    }
    else
    {
        manifold->penetration = 0;
    }
    separation = vector2f_dot_vv(reference_face_normal.x,reference_face_normal.y, incident_face[1].x,incident_face[1].y) - dist_from_origin;

    if (separation <= 0.0f)
    {
        manifold->contacts[contact_count] = incident_face[1];
        manifold->penetration += -separation;
        ++contact_count;
        manifold->penetration /= contact_count;
    }
    manifold->contact_count = contact_count;
    manifold->reference_face_index = reference_face_index;
    manifold->incident_face_index = incident_face_index;

    rigid_bodies->is_colliding[a] = 1;            
    rigid_bodies->is_colliding[b] = 1;

    return 1;
}
void mpl_update(RigidBody *rigid_bodies, unsigned int body_count,Manifold *manifolds, unsigned int manifold_capacity, int iterations ,float dt, float G)
{
    unsigned int manifold_count = 0;
    Vector2 gravity;
    gravity.x = 0;
    gravity.y = G;

    for (int i = 0; i < body_count; i++)
    {
        if (!rigid_bodies->is_static[i] && rigid_bodies->active[i] && !rigid_bodies->ignore_gravity[i])
        {                    
            gravity.y = round(G * rigid_bodies->mass[i]);
            mpl_rigid_body_apply_force(rigid_bodies,i,gravity);
        }
    }  
    
    for (int i = 0; i < body_count; i++)
    {
        if (rigid_bodies->active[i])
        {
            for (int j = i; j < body_count; ++j)
            {
                if (rigid_bodies->active[j])
                {
                    if (i != j)
                    {
                        if (!(rigid_bodies->is_static[i] && rigid_bodies->is_static[j]))
                        {
                            if(manifold_count < manifold_capacity)
                            {
                                
                                if (mpl_broad_phase(rigid_bodies,i,j))
                                {
                                    Manifold m = manifolds[manifold_count];
                                    mpl_manifold_reset(&m);
                                    unsigned int flag = 0;
                                    if (!rigid_bodies->is_polygon[i])
                                    {
                                        if (!rigid_bodies->is_polygon[j])
                                        {
                                            flag =mpl_narrow_phase_c2c(&manifolds[manifold_count],rigid_bodies,i,j);
                                        }
                                        else
                                        {
                                           flag = mpl_narrow_phase_p2c(&manifolds[manifold_count],rigid_bodies,j,i);
                                        }
                                    }
                                    else
                                    {
                                        if (!rigid_bodies->is_polygon[j])
                                        {
                                            flag = mpl_narrow_phase_p2c(&manifolds[manifold_count],rigid_bodies,i,j);
                                        }
                                        else
                                        {
                                            flag = mpl_narrow_phase_p2p(&manifolds[manifold_count],rigid_bodies,i,j);
                                        }
                                    }

                                    if (flag) manifold_count++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


/*
    for (int i = 0; i < body_count; ++i)
    {
        __m256 delta_time = _mm256_set_ps(dt, dt, dt, dt, dt, dt, dt, dt);
        __m256 rb_torque = _mm256_set_ps(rigid_bodies[i].torque, rigid_bodies[i+1].torque, rigid_bodies[i+2].torque, rigid_bodies[i+3].torque, rigid_bodies[i+4].torque, rigid_bodies[i+5].torque, rigid_bodies[i+6].torque, rigid_bodies[i+7].torque);
        __m256 rb_inverse_inertia = _mm256_set_ps(rigid_bodies[i].inverse_inertia, rigid_bodies[i+1].inverse_inertia, rigid_bodies[i+2].inverse_inertia, rigid_bodies[i+3].inverse_inertia, rigid_bodies[i+4].inverse_inertia, rigid_bodies[i+5].inverse_inertia, rigid_bodies[i+6].inverse_inertia, rigid_bodies[i+7].inverse_inertia);
        
        __m256 rb_force_x = _mm256_set_ps(rigid_bodies[i].force.x, rigid_bodies[i+1].force.x, rigid_bodies[i+2].force.x, rigid_bodies[i+3].force.x, rigid_bodies[i+4].force.x, rigid_bodies[i+5].force.x, rigid_bodies[i+6].force.x, rigid_bodies[i+7].force.x);
        __m256 rb_inverse_mass = _mm256_set_ps(rigid_bodies[i].inverse_mass, rigid_bodies[i+1].inverse_mass, rigid_bodies[i+2].inverse_mass, rigid_bodies[i+3].inverse_mass, rigid_bodies[i+4].inverse_mass, rigid_bodies[i+5].inverse_mass, rigid_bodies[i+6].inverse_mass, rigid_bodies[i+7].inverse_mass);
        __m256 vel = _mm256_mul_ps(rb_force_x, rb_inverse_mass);
        vel = _mm256_mul_ps(vel, delta_time);
        __m256 rb_vel_x = _mm256_add_ps(rb_vel_x, vel);

         __m256 rb_force_y = _mm256_set_ps(rigid_bodies[i].force.y, rigid_bodies[i+1].force.y, rigid_bodies[i+2].force.y, rigid_bodies[i+3].force.y, rigid_bodies[i+4].force.y, rigid_bodies[i+5].force.y, rigid_bodies[i+6].force.y, rigid_bodies[i+7].force.y);
        vel = _mm256_mul_ps(rb_force_x, rb_inverse_mass);
        vel = _mm256_mul_ps(vel, delta_time);
        __m256 rb_vel_y = _mm256_add_ps(rb_vel_y, vel);

        
        

        __m256 rb_lock_orientation = _mm256_set_ps(rigid_bodies[i].lock_orientation, rigid_bodies[i+1].lock_orientation, rigid_bodies[i+2].lock_orientation, rigid_bodies[i+3].lock_orientation, rigid_bodies[i+4].lock_orientation, rigid_bodies[i+5].lock_orientation, rigid_bodies[i+6].lock_orientation, rigid_bodies[i+7].lock_orientation);
        __m256 rb_angular_velocity_x = _mm256_set_ps(rigid_bodies[i].velocity.x, rigid_bodies[i+1].velocity.x, rigid_bodies[i+2].velocity.x, rigid_bodies[i+3].velocity.x, rigid_bodies[i+4].velocity.x, rigid_bodies[i+5].velocity.x, rigid_bodies[i+6].velocity.x, rigid_bodies[i+7].velocity.x);
        __m256 rb_angular_velocity_y = _mm256_set_ps(rigid_bodies[i].velocity.y, rigid_bodies[i+1].velocity.y, rigid_bodies[i+2].velocity.y, rigid_bodies[i+3].velocity.y, rigid_bodies[i+4].velocity.y, rigid_bodies[i+5].velocity.y, rigid_bodies[i+6].velocity.y, rigid_bodies[i+7].velocity.y);

        __m256 angular_velocity = _mm256_mul_ps(rb_torque, rb_inverse_inertia);
        __m256 angular_velocity = _mm256_mul_ps(angular_velocity, delta_time);
        __m256 angular_velocity = _mm256_mul_ps(angular_velocity, rb_lock_orientation);
        rb_angular_velocity_x = _mm256_add_ps(rb_angular_velocity_x, angular_velocity);
        rb_angular_velocity_y = _mm256_add_ps(rb_angular_velocity_y, angular_velocity);


        
        if (!bodies[i].is_static && bodies[i].active)
        {
            mpl_rigid_body_integrate_forces(&bodies[i],dt);
        }
    }  */
     for (int i = 0; i < body_count; ++i)
    {
        if (!rigid_bodies->is_static[i] && rigid_bodies->active[i])
        {
            mpl_rigid_body_integrate_forces(rigid_bodies,i,dt);
        }
    } 
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].is_empty)
        {
            mpl_manifold_init(&manifolds[i],rigid_bodies,i,dt);
        }
    }
            
    for (int j = 0; j < iterations; ++j)
    {
        for (int i = 0; i < manifold_count; ++i)
        {
            if (!manifolds[i].is_empty)
            {
               mpl_manifold_apply_impulse(&manifolds[i],rigid_bodies);
            }
        }
    }
            
    for (int i = 0; i < body_count; ++i)
    {
        if (!rigid_bodies->is_static[i] && rigid_bodies->active[i])
        {
            mpl_rigid_body_integrate_velocity(rigid_bodies,i,dt);
        }
    }
      
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].is_empty)
        {
            mpl_manifold_positional_correction(&manifolds[i],rigid_bodies);
        }
    }
    

    for (int i = 0; i < body_count; ++i)
    {
        mpl_rigid_body_clear_forces(rigid_bodies,i);
    }
    for (int i = 0; i < body_count; ++i)
    {
        mpl_polygon_transform(rigid_bodies,i);
    }
}
void mpl_rigid_body_set_static(RigidBody *rb, unsigned int index)
{
    rb->inertia[index] = 0.0f;
    rb->inverse_inertia[index] = 0.0f;
    rb->mass[index] = 0.0f;
    rb->inverse_mass[index] = 0.0f;
    rb->lock_orientation[index] = 1;
}
