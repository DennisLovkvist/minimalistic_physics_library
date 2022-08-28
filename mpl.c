#include <stdlib.h>
#include <stdio.h>
#include "mpl.h"
#include <math.h>

//math functions begin
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
float vector2f_cross_vv(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.y) - (v1.y * v2.x);
}
void vector2f_cross_v(Vector2 v, Vector2 *result)
{
    result->x = v.y;
    result->y = -v.x;
}
float vector2f_dot_vv(Vector2 a, Vector2 b)
{
    return a.x * b.x + a.y * b.y;
}  
float vector2f_dot_vv2(float ax,float ay,float bx, float by)
{
    return ax * bx + ay * by;
}  
float vector2f_dot_vf(Vector2 a, float b)
{
    return a.x * cos(b) + a.y * sin(b);
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
        result[i][0] = 0;
        for (int k = 0; k < 3; k++)
        {
            result[i][0] += input_matrix[i][k] * transformation_matrix[k][0];
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
float min(float a , float b)
{
    return (a < b) ? a:b;
}
float max(float a , float b)
{
    return (a > b) ? a:b;
}
//math functions end
//helper functions begin
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
    float best_projection = FLT_MIN;   
    best_vertex->x = 0;
    best_vertex->y = 0;

    for (int i = 0; i < vertex_count; ++i)
    {
        Vector2 v;
        v.x = vertices[i].x;
        v.y = vertices[i].y;
        float projection = vector2f_dot_vv(v, dir);

        if (projection > best_projection)
        {
            best_vertex->x = v.x;
            best_vertex->y = v.y;
            best_projection = projection;
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
void mpl_matrix_set(float matrix[3][3],float tx, float ty, float r)
{    
    matrix[0][0] = cos(r);
    matrix[0][1] = sin(r);
    matrix[0][2] = tx;

    matrix[1][0] = -sin(r);
    matrix[1][1] = cos(r);
    matrix[1][2] = ty;

    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = 1;
}
void mpl_polygon_transform(Polygon *polygon)
{    
    Vector2 *vertices = polygon->vertices;
    Vector2 *normals = polygon->normals;
    Vector2 position = polygon->position;
    Vector2 position_last = polygon->position_last;
    float orientation = polygon->orientation;

    float tx = position.x - position_last.x;
    float ty = position.y - position_last.y;
    float r = orientation - polygon->orientation_last;    
    float (*matrix)[3] = polygon->matrix;
    float input[3][1];
    float output[3][1];
    for (int i = 0; i < polygon->vertex_count; i++)
    {        
        mpl_matrix_set(matrix,-position_last.x, -position_last.y, 0);
        input[0][0] = vertices[i].x;
        input[1][0] = vertices[i].y;
        input[2][0] = 1;
        matrix_mul_mm(matrix,input,output);
        vertices[i].x = output[0][0];
        vertices[i].y = output[1][0];

        mpl_matrix_set(matrix,0, 0, r);
        input[0][0] = vertices[i].x;
        input[1][0] = vertices[i].y;
        input[2][0] = 1;
        matrix_mul_mm(matrix,input,output);
        vertices[i].x = output[0][0];
        vertices[i].y = output[1][0];

        input[0][0] = normals[i].x;
        input[1][0] = normals[i].y;
        input[2][0] = 1;
        matrix_mul_mm(matrix,input,output);
        normals[i].x = output[0][0];
        normals[i].y = output[1][0];

        mpl_matrix_set(matrix,position_last.x + tx, position_last.y + ty, 0);
        input[0][0] = vertices[i].x;
        input[1][0] = vertices[i].y;
        input[2][0] = 1;
        matrix_mul_mm(matrix,input,output);
        vertices[i].x = output[0][0];
        vertices[i].y = output[1][0];
    }
    polygon->position_last.x = position.x;
    polygon->position_last.y = position.y;
    polygon->orientation_last = orientation;
}
int mpl_find_incident_face(Polygon *polygon_reference, Polygon *polygon_incident, int reference_index)
{
    Vector2 reference_normal;
    
    reference_normal.x = polygon_reference->normals[reference_index].x;
    reference_normal.y = polygon_reference->normals[reference_index].y;

    matrix_mul_mv(polygon_reference->matrix,reference_normal,&reference_normal);

    float transposed[3][3];
    matrix_transpose_mm(polygon_reference->matrix,transposed);
    matrix_mul_mv(transposed,reference_normal,&reference_normal);

    int incident_face = 0;
    float min_dot = FLT_MAX;
    for (int i = 0; i < polygon_incident->vertex_count; ++i)
    {        
        float dot = vector2f_dot_vv(reference_normal,polygon_incident->normals[i]);

        if (dot < min_dot)
        {
            min_dot = dot;
            incident_face = i;
        }
    }
    return incident_face;
}
int mpl_clip(Vector2 n, float c, Vector2 *face)
{
    int sp = 0;
    Vector2 points[2];
    points[0].x = face[0].x;
    points[0].y = face[0].y;
    points[1].x = face[1].x;
    points[1].y = face[1].y;

    float d1 = vector2f_dot_vv(n, face[0])-c;
    float d2 = vector2f_dot_vv(n, face[1])-c;

    if (d1 <= 0.0f) points[sp++] = face[0];
    if (d2 <= 0.0f) points[sp++] = face[1];

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
float mpl_dist_sqr(Vector2 a, Vector2 b)
{
    Vector2 c = {a.x-b.x,a.y-b.y};
    return vector2f_dot_vv(c,c);
}
//helper functions end
//physics init functions begin
void mpl_polygon_init(Polygon *polygon, Vector2 *vertices, unsigned int vertex_count)
{
    polygon->orientation = polygon->orientation_last = 0.0;
    polygon->radius = 0;
    vertex_count = (vertex_count > 4) ? 4:vertex_count;
    polygon->vertex_count = vertex_count;

    for (int i = 0; i < vertex_count; i++)
    {
        polygon->vertices[i].x = vertices[i].x;
        polygon->vertices[i].y = vertices[i].y;
    }    
    for (int i = 0; i < vertex_count; i++)
    {
        polygon->edges[i].m_index_a = i;
        polygon->edges[i].m_index_b = (i + 1 < vertex_count) ? i + 1 : 0;  

        Vector2 v0 = polygon->vertices[polygon->edges[i].m_index_b];
        Vector2 v1 = polygon->vertices[polygon->edges[i].m_index_a];

        vector2f_normal(v0,v1,&polygon->normals[i]);
        vector2f_normalize(&polygon->normals[i]);
    }
    mpl_polygon_get_centroid(polygon->vertices,vertex_count, &polygon->position);
    
    polygon->position_last.x = polygon->position.x;
    polygon->position_last.y = polygon->position.y;

    mpl_polygon_transform(polygon);
    
}
void mpl_rigid_body_compute_mass(RigidBody *rigid_body,float density)
{
    float area = 0;
    unsigned int vertex_count = rigid_body->polygon.vertex_count;
    Vector2 c;
    mpl_polygon_get_centroid(rigid_body->polygon.vertices,vertex_count,&c);
    float dx = 0;
    float dy = 0;
    float longest_radius = 0;
    float radius = 0;

    for (int i = 0; i < vertex_count; i++)
    {
        Vector2 v = rigid_body->polygon.vertices[i];
        int next = mpl_polygon_get_next_vertex_index(i,vertex_count);
        dx = c.x - v.x;
        dy = c.y - v.y;
        float d1 = sqrt(dx * dx + dy * dy);
        radius += d1;
        dx = rigid_body->polygon.vertices[next].x - v.x;
        dy = rigid_body->polygon.vertices[next].y - v.y;
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
    rigid_body->inverse_mass = 1.0f / rigid_body->mass;
    rigid_body->inverse_inertia = 1.0f / rigid_body->inertia;
    
}
void mpl_rigid_body_compute_mass_circle(RigidBody *rigid_body,float density,float radius)
{
    float area = PI * radius*radius;

    rigid_body->polygon.radius = radius;
    rigid_body->longest_radius = radius;

    rigid_body->mass = density * area;
    rigid_body->inertia = rigid_body->mass * radius * radius;
    rigid_body->inverse_mass = 1.0f / rigid_body->mass;
    rigid_body->inverse_inertia = 1.0f / rigid_body->inertia;
}
void mpl_polygon_make_square(Polygon *polygon, unsigned int width,unsigned int height)
{
    polygon->vertex_count = 4;
    Vector2 vertices[4];
    float half_width = (float)width*0.5;
    float half_height = (float)height*0.5;

    vertices[0].x = -half_width;vertices[0].y = -half_height;
    vertices[1].x = half_width;vertices[1].y = -half_height;
    vertices[2].x = half_width;vertices[2].y = half_height;
    vertices[3].x = -half_width;vertices[3].y = half_height;
    mpl_polygon_init(polygon,vertices,4);
}
void mpl_polygon_make_triangle(Polygon *polygon, unsigned int width,unsigned int height)
{
    polygon->vertex_count = 3;
    Vector2 vertices[3];
    float half_width = (float)width*0.5;
    float half_height = (float)height*0.5;

    vertices[0].x = 0;vertices[0].y = -half_height;
    vertices[1].x = half_width;vertices[1].y = half_height;
    vertices[2].x = -half_width;vertices[2].y = half_height;

    mpl_polygon_init(polygon,vertices,3);
}
void mpl_rigid_body_init(RigidBody *rigid_body,RigidBodyMeta *rigid_body_meta, unsigned int shape,unsigned int width,unsigned int height,unsigned int max_collisions)
{    
    if(shape == 0)
    {
        rigid_body_meta->is_polygon = 0;
        rigid_body->polygon.vertex_count = 0;
    }
    else if(shape == 1)
    {
        rigid_body_meta->is_polygon = 1;
        mpl_polygon_make_triangle(&rigid_body->polygon,width, height);        
    }
    else if(shape == 2)
    {
        rigid_body_meta->is_polygon = 1;
        mpl_polygon_make_square(&rigid_body->polygon,width, height);
    }
    else
    {
        rigid_body_meta->is_polygon = 1;
        mpl_polygon_make_square(&rigid_body->polygon,width, height);
    }
    rigid_body_meta->active = 1;
    rigid_body_meta->collisions = 0;
    rigid_body_meta->max_collisions = max_collisions;
    rigid_body_meta->is_static = 0;
    rigid_body_meta->ignore_gravity = 0;
    rigid_body->lock_orientation = 0;
    rigid_body->restitution = 0.6f;
    rigid_body->static_friction =  0.2f;
    rigid_body->dynamic_friction = 0.4f;
    rigid_body->angular_velocity = 0;
    rigid_body->torque = 0;
    rigid_body->mass = 0;
    rigid_body->inverse_mass = 0;
    rigid_body->inertia = 0;
    rigid_body->inverse_inertia = 0;
    rigid_body->longest_radius = 0;
    rigid_body->damping = 0.95f;
    rigid_body->velocity_epsilon = 0.001f;
    rigid_body->velocity.x = rigid_body->velocity.y =  0;
    rigid_body->force.x = rigid_body->force.y = 0;
    rigid_body->polygon.orientation = rigid_body->polygon.orientation_last = 0;

    if(rigid_body_meta->is_polygon)
    {
        mpl_rigid_body_compute_mass(rigid_body,0.001);
    }
    else
    {
        mpl_rigid_body_compute_mass_circle(rigid_body,0.0001,width/2);
    }

}
void mpl_rigid_body_set_static(RigidBody *rb,RigidBodyMeta *rbm)
{
    rb->inertia = 0.0f;
    rb->inverse_inertia = 0.0f;
    rb->mass = 0.0f;
    rb->inverse_mass = 0.0f;
    rb->lock_orientation = 1;
    rbm->is_static = 1;
}
//physics init functions end
//integration begin
void mpl_rigid_body_set_position(RigidBody *rigid_body,float x, float y)
{
    rigid_body->polygon.position.x = x;
    rigid_body->polygon.position.y = y;
    mpl_polygon_transform(&rigid_body->polygon);
}
void mpl_rigid_body_apply_force(RigidBody *rigid_body,Vector2 f)
{
    rigid_body->force.x += f.x;
    rigid_body->force.y += f.y;
}

void mpl_rigid_body_apply_impulse(RigidBody *rigid_body,Vector2 impulse, Vector2 contactVector)
{    
    if(isnan(impulse.x) || isnan(impulse.y))
    {
        return;
    }
    rigid_body->velocity.x += rigid_body->inverse_mass * impulse.x;
    rigid_body->velocity.y += rigid_body->inverse_mass * impulse.y;
    
    float n = vector2f_cross_vv(contactVector, impulse);
    rigid_body->angular_velocity += rigid_body->inverse_inertia * n * (!rigid_body->lock_orientation);
    
}
void mpl_rigid_body_clear_forces(RigidBody *rigid_body)
{
    rigid_body->force.x = 0;
    rigid_body->force.y = 0;
    
}
//integration end
//collision begin
void mpl_manifold_init(Manifold *manifold,float dt)
{
    manifold->e = min(manifold->rigid_body_a->restitution, manifold->rigid_body_b->restitution);
    manifold->static_friction = sqrt(manifold->rigid_body_a->static_friction * manifold->rigid_body_b->static_friction);
    manifold->dynamic_friction = sqrt(manifold->rigid_body_a->dynamic_friction * manifold->rigid_body_b->dynamic_friction);
    
}
void mpl_manifold_apply_impulse(Manifold *manifold)
{
    Polygon *polygon_a = &manifold->rigid_body_a->polygon;
    Polygon *polygon_b = &manifold->rigid_body_b->polygon;

    RigidBody *rigid_body_a = manifold->rigid_body_a;
    RigidBody *rigid_body_b = manifold->rigid_body_b;

    for (int i = 0; i < manifold->contact_count; ++i)
    {
        Vector2 ra;
        ra.x = polygon_a->position.x-manifold->contacts[i].x;
        ra.y = polygon_a->position.y-manifold->contacts[i].y;

        Vector2 rb;
        rb.x = polygon_b->position.x-manifold->contacts[i].x;
        rb.y = polygon_b->position.y-manifold->contacts[i].y;         
        
        Vector2 cross_rb_av;
        vector2f_cross_vf(rb,rigid_body_b->angular_velocity,&cross_rb_av);
        Vector2 cross_ra_av;
        vector2f_cross_vf(ra,rigid_body_a->angular_velocity,&cross_ra_av);
        
        Vector2 rv;
        rv.x = rigid_body_b->velocity.x + cross_rb_av.x - rigid_body_a->velocity.x - cross_ra_av.x;
        rv.y = rigid_body_b->velocity.y + cross_rb_av.y - rigid_body_a->velocity.y - cross_ra_av.y;

        float contact_velocity = vector2f_dot_vv(rv, manifold->normal);
        
        if (contact_velocity > 0)
        {
            return;
        }

        float raCrossN = vector2f_cross_vv(ra, manifold->normal);
        float rbCrossN = vector2f_cross_vv(rb, manifold->normal);

        float inverse_mass_sum = rigid_body_a->inverse_mass + rigid_body_b->inverse_mass + (raCrossN * raCrossN) * rigid_body_a->inverse_inertia + (rbCrossN * rbCrossN) * rigid_body_b->inverse_inertia;
        
        float j = -(1.0f + manifold->e) * contact_velocity;                
        j /= inverse_mass_sum;
        j /= manifold->contact_count;


        Vector2 impulse;
        impulse.x = manifold->normal.x * j;
        impulse.y = manifold->normal.y * j;        

        float dot = vector2f_dot_vv(rv, manifold->normal);
        Vector2 t;
        t.x = rv.x - (manifold->normal.x * dot);
        t.y = rv.y - (manifold->normal.y * dot);
        vector2f_normalize(&t);

        float jt = -vector2f_dot_vv(rv, t);
        jt /= inverse_mass_sum;
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

        mpl_rigid_body_apply_impulse(rigid_body_b,impulse,rb);
        mpl_rigid_body_apply_impulse(rigid_body_b,tangent_impulse,rb);
        
        impulse.x = -impulse.x;
        impulse.y = -impulse.y;
        tangent_impulse.x = -tangent_impulse.x;
        tangent_impulse.y = -tangent_impulse.y;

        mpl_rigid_body_apply_impulse(rigid_body_a,impulse,ra);
        mpl_rigid_body_apply_impulse(rigid_body_a,tangent_impulse,ra);    
    }
}
void mpl_manifold_positional_correction(Manifold *manifold)
{
    RigidBody *rigid_body_a = manifold->rigid_body_a;
    RigidBody *rigid_body_b = manifold->rigid_body_b;

    const float k_slop = 0.01f;
    const float percent = 0.2f;

    float f = max(manifold->penetration - k_slop, 0.0f)/(rigid_body_a->inverse_mass + rigid_body_b->inverse_mass);
    Vector2 correction;
    correction.x = f * manifold->normal.x * percent;
    correction.y = f * manifold->normal.y * percent;
    
    rigid_body_a->polygon.position.x -= correction.x * rigid_body_a->inverse_mass;
    rigid_body_a->polygon.position.y -= correction.y * rigid_body_a->inverse_mass;

    rigid_body_b->polygon.position.x += correction.x * rigid_body_b->inverse_mass;
    rigid_body_b->polygon.position.y += correction.y * rigid_body_b->inverse_mass;
}
void mpl_manifold_reset(Manifold *manifold)
{
    manifold->contact_count = 0;
    manifold->penetration = 0;
}
float mpl_find_axis_least_penetration(Polygon *polygon_a, Polygon *polygon_b,int *face)
{
    float best_distance = FLT_MIN;
    int best_index = 0;
    for (int i = 0; i < polygon_a->vertex_count; ++i)
    {
        Vector2 normal;
        normal.x = polygon_a->normals[i].x;
        normal.y = polygon_a->normals[i].y;

        Vector2 normal_flipped;
        normal_flipped.x = -normal.x;
        normal_flipped.y = -normal.y;

        Vector2 support_point;
        support_point.x = 0;
        support_point.y = 0;
        mpl_polygon_get_support_point(polygon_b->vertices,polygon_b->vertex_count,normal_flipped,&support_point);

        Vector2 delta;
        delta.x = support_point.x - polygon_a->vertices[i].x;
        delta.y = support_point.y - polygon_a->vertices[i].y;

        float dot = vector2f_dot_vv(normal, delta);

        if (dot > best_distance)
        {
            best_distance = dot;
            best_index = i;
        }
    }
    *face = best_index;
    return best_distance;
}
//collision begin

unsigned int mpl_broad_phase(RigidBody *A, RigidBody *B)
{
    float dx = A->polygon.position.x - B->polygon.position.x;
    float dy = A->polygon.position.y - B->polygon.position.y;
    float dist = sqrt(dx * dx + dy * dy);

    return (dist < A->longest_radius + B->longest_radius);
}
unsigned int mpl_gt(float a, float b)
{
    return a >= b * 0.95f + a * 0.0f;
}

unsigned int mpl_narrow_phase_c2c(Manifold *manifold, RigidBody *a, RigidBody *b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = a;
    manifold->rigid_body_b = b;
    Vector2 normal;
    normal.x = b->polygon.position.x - a->polygon.position.x;
    normal.y = b->polygon.position.y - a->polygon.position.y;

    float dist_sqr = normal.x * normal.x + normal.y * normal.y;
    float radius = a->polygon.radius + b->polygon.radius;

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
        manifold->penetration = a->polygon.radius;
        manifold->normal.x = 1;
        manifold->normal.y = 0;
        manifold->contacts[0] = a->polygon.position;
    }
    else
    {
        manifold->penetration = radius - distance;
        manifold->normal.x = normal.x / distance;
        manifold->normal.y = normal.y / distance;
        manifold->contacts[0].x = manifold->normal.x * a->polygon.radius + a->polygon.position.x;
        manifold->contacts[0].y = manifold->normal.y * a->polygon.radius + a->polygon.position.y;
    }

    return 1;
}
unsigned int mpl_narrow_phase_p2c(Manifold *manifold, RigidBody *a, RigidBody *b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = a;
    manifold->rigid_body_b = b;

    Vector2 b_position = b->polygon.position;

    float b_radius = b->polygon.radius;

    Vector2 normal;
    vector2f_normal(b_position, b_position,&normal);
    vector2f_normalize(&normal);
    vector2f_cross_v(normal,&normal);

    Vector2 normal_flipped = {-normal.x,-normal.y};

    Vector2 support_point; 
    mpl_polygon_get_support_point(a->polygon.vertices,a->polygon.vertex_count, normal_flipped,&support_point);

    Vector2 v2;

    float separation = FLT_MIN;
    int face_normal = 0;

    int va_count = a->polygon.vertex_count;
    Vector2 *normals_a = a->polygon.normals;
    Vector2 *vertices_a = a->polygon.vertices;
    for (int i = 0; i < va_count; i++)
    {
        Vector2 temp = {b_position.x - vertices_a[i].x,b_position.y - vertices_a[i].y};
        float s = vector2f_dot_vv(normals_a[i],temp);
        if (s > b_radius)
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
    support_point = vertices_a[a->polygon.edges[face_normal].m_index_a];
    v2 = vertices_a[a->polygon.edges[face_normal].m_index_b];

    float dot1 = vector2f_dot_vv2(b_position.x- support_point.x,b_position.y - support_point.y,v2.x- support_point.x,v2.y - support_point.y);
    float dot2 = vector2f_dot_vv2(b_position.x- v2.x,b_position.y - v2.y,support_point.x-v2.x,support_point.y-v2.y);

    manifold->contact_count = 1;
    manifold->penetration = b->polygon.radius - separation;

    if (separation < 0.0001f)
    {
        manifold->contact_count = 1;
        manifold->normal = normals_a[face_normal];
        manifold->contacts[0].x = manifold->normal.x * b_radius + b_position.x;
        manifold->contacts[0].y = manifold->normal.y * b_radius + b_position.y;
        manifold->penetration = b_radius;
        return 1;
    }
    if (dot1 <= 0.0f)
    {
        if (mpl_dist_sqr(b_position,support_point) > b_radius * b_radius)
        {
            manifold->is_empty = 1;
            return 0;
        }
        
        manifold->contacts[0] = support_point;
        float r = atan2(b_position.y - support_point.y, b_position.x - support_point.x);
        Vector2 n = {cos(r), sin(r)};
        vector2f_normalize(&n);
        manifold->normal = n;
        Vector2 surface = {b_position.x + cos(r-PI) * b_radius,b_position.y + sin(r-PI) * b_radius};

        float dx = surface.x - support_point.x;
        float dy = surface.y - support_point.y;
        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else if (dot2 <= 0.0f)
    {
        if (mpl_dist_sqr(b_position, v2) > b_radius * b_radius)
        {
            manifold->is_empty = 1;
            return 0;
        }
        manifold->contacts[0] = v2;
        float r = atan2(b_position.y - v2.y, b_position.x - v2.x);
        Vector2 n = {cos(r),sin(r)};
        vector2f_normalize(&n);
        manifold->normal = n;
        Vector2 contact_point = {b_position.x + cos(r-PI) * b_radius,b_position.y + sin(r-PI) * b_radius};

        float dx = contact_point.x - v2.x;
        float dy = contact_point.y - v2.y;

        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else
    {
        manifold->normal = normals_a[face_normal];
        manifold->contacts[0].x = b_position.x - normals_a[face_normal].x * (b_radius - manifold->penetration);
        manifold->contacts[0].y = b_position.y - normals_a[face_normal].y * (b_radius - manifold->penetration);
    }
    return 1;
}
unsigned int mpl_narrow_phase_p2p(Manifold *manifold, RigidBody *rigid_body_a, RigidBody *rigid_body_b)
{
    manifold->is_empty = 0;
    manifold->rigid_body_a = rigid_body_a;
    manifold->rigid_body_b = rigid_body_b;

    manifold->contact_count = 0;

    int face_a = 0;
    int face_b = 0;                                            

    float penetration_a = mpl_find_axis_least_penetration(&rigid_body_a->polygon, &rigid_body_b->polygon,&face_a);
    float penetration_b = mpl_find_axis_least_penetration(&rigid_body_b->polygon, &rigid_body_a->polygon,&face_b);

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
    unsigned int flip = 0;

    Polygon *reference_polygon = manifold->reference_polygon;
    Polygon *incident_polygon = manifold->incident_polygon;

    if (mpl_gt(penetration_a, penetration_b))
    {
        reference_polygon = &rigid_body_a->polygon;
        incident_polygon = &rigid_body_b->polygon;
        reference_face_index = face_a;
        flip = 0;
    }
    else
    {
        reference_polygon = &rigid_body_b->polygon;
        incident_polygon = &rigid_body_a->polygon;
        reference_face_index = face_b;
        flip = 1;
    }

    Vector2 reference_face[2];
    reference_face[0] = reference_polygon->vertices[reference_polygon->edges[reference_face_index].m_index_a];
    reference_face[1] = reference_polygon->vertices[reference_polygon->edges[reference_face_index].m_index_b];

    int incident_face_index = mpl_find_incident_face(reference_polygon, incident_polygon, reference_face_index);
    
    Vector2 incident_face[2];
    incident_face[0] = incident_polygon->vertices[incident_polygon->edges[incident_face_index].m_index_a];
    incident_face[1] = incident_polygon->vertices[incident_polygon->edges[incident_face_index].m_index_b];

    Vector2 side_plane_normal;
    side_plane_normal.x = reference_face[1].x - reference_face[0].x;
    side_plane_normal.y = reference_face[1].y - reference_face[0].y;
    vector2f_normalize(&side_plane_normal);

    Vector2 reference_face_normal;
    reference_face_normal.x = side_plane_normal.y;
    reference_face_normal.y = -side_plane_normal.x;       
    
    float dist_from_origin = vector2f_dot_vv(reference_face_normal, reference_face[0]);
    float negative_side = -vector2f_dot_vv(side_plane_normal, reference_face[0]);
    float positive_side = vector2f_dot_vv(side_plane_normal, reference_face[1]);

    manifold->normal.x = reference_face_normal.x;
    manifold->normal.y = reference_face_normal.y;

    Vector2 side_plane_normal_negated = {-side_plane_normal.x,-side_plane_normal.y};

    if (mpl_clip(side_plane_normal_negated, negative_side, &incident_face[0]) < 2)
    {
        return 0;
    }
    if (mpl_clip(side_plane_normal, positive_side, &incident_face[0]) < 2)
    {
        return 0;
    }
    if (flip)
    {
        manifold->normal.x = -manifold->normal.x;
        manifold->normal.y = -manifold->normal.y;
    }

    int contact_count = 0;
    float separation = vector2f_dot_vv(reference_face_normal, incident_face[0]) - dist_from_origin;

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
    separation = vector2f_dot_vv(reference_face_normal, incident_face[1]) - dist_from_origin;

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
    return 1;
}
void mpl_update(RigidBody *rigid_bodies,RigidBodyMeta *rigid_bodies_meta, unsigned int body_count,Manifold *manifolds,unsigned int manifold_capacity, int iterations ,float dt, float G)
{
    unsigned int manifold_count = 0;
    Vector2 gravity;
    gravity.x = 0;
    gravity.y = G;
    for (int i = 0; i < body_count; i++)
    {
        RigidBody *rigid_body = &rigid_bodies[i];
        RigidBodyMeta *rigid_body_meta = &rigid_bodies_meta[i];

        if (!rigid_body_meta->is_static && rigid_body_meta->active && !rigid_body_meta->ignore_gravity)
        {                    
            gravity.y = round(G * rigid_body->mass)*dt;
            rigid_body->force.x += gravity.x;
            rigid_body->force.y += gravity.y;
        }
    }       
    for (size_t k = 0; k < iterations; k++)
    {        
        for (int i = 0; i < body_count; i++)
        {
            RigidBody *rigid_body_a = &rigid_bodies[i];
            RigidBodyMeta *rigid_body_meta_a = &rigid_bodies_meta[i];
            rigid_body_meta_a->collisions = 0;
            for (int j = i+1; j < body_count; ++j)
            {
                RigidBody *rigid_body_b = &rigid_bodies[j];
                RigidBodyMeta *rigid_body_meta_b = &rigid_bodies_meta[j];

                if (rigid_body_meta_a->is_static + rigid_body_meta_b->is_static != 2)
                {
                    if(manifold_count < manifold_capacity)
                    {                                
                        if (mpl_broad_phase(rigid_body_a, rigid_body_b))
                        {
                            Manifold *manifold = &manifolds[manifold_count];
                            mpl_manifold_reset(manifold);
                            int collision_flag = 0;

                            if (!rigid_body_meta_a->is_polygon)
                            {
                                if (!rigid_body_meta_b->is_polygon)
                                {
                                    collision_flag =mpl_narrow_phase_c2c(manifold,rigid_body_a, rigid_body_b);
                                }
                                else
                                {
                                    collision_flag = mpl_narrow_phase_p2c(manifold,rigid_body_b, rigid_body_a);
                                }
                            }
                            else
                            {
                                if (!rigid_body_meta_b->is_polygon)
                                {
                                    collision_flag = mpl_narrow_phase_p2c(manifold,rigid_body_a, rigid_body_b);
                                }
                                else
                                {
                                    collision_flag = mpl_narrow_phase_p2p(manifold,rigid_body_a, rigid_body_b);
                                }
                            }
                            rigid_body_meta_a->collisions+=collision_flag*!rigid_body_meta_a->is_static;                       
                            manifold_count+=collision_flag;

                            if(rigid_body_meta_a->collisions > rigid_body_meta_a->max_collisions)
                            {
                                break;
                            }
                        }
                    }
                }                
            }   
        }
    }
    for (int i = 0; i < body_count; ++i)
    {
        RigidBody *rigid_body = &rigid_bodies[i];
        RigidBodyMeta *rigid_body_meta = &rigid_bodies_meta[i];
        
        int flag = rigid_body_meta->active*(!rigid_body_meta->is_static);
        rigid_body->velocity.x += ((rigid_body->force.x * rigid_body->inverse_mass) * (dt)*flag);
        rigid_body->velocity.y += ((rigid_body->force.y * rigid_body->inverse_mass) * (dt)*flag);
        rigid_body->angular_velocity += (rigid_body->torque * rigid_body->inverse_inertia * (dt))*(!rigid_body->lock_orientation);
        
    }  
    for (int i = 0; i < manifold_count; ++i)
    {
        Manifold *manifold = &manifolds[i];
        if (!manifold->is_empty)
        {
            RigidBody *rigid_body_a = manifold->rigid_body_a;
            RigidBody *rigid_body_b = manifold->rigid_body_b;
            manifold->e = min(rigid_body_a->restitution, rigid_body_b->restitution);
            manifold->static_friction = sqrt(rigid_body_a->static_friction * rigid_body_b->static_friction);
            manifold->dynamic_friction = sqrt(rigid_body_a->dynamic_friction * rigid_body_b->dynamic_friction);
        }
    }  
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].is_empty)
        {
            mpl_manifold_apply_impulse(&manifolds[i]);
        }
    }    
    for (int i = 0; i < body_count; ++i)
    {
        RigidBodyMeta *rigid_body_meta = &rigid_bodies_meta[i];
        RigidBody *rigid_body = &rigid_bodies[i];
        if (!rigid_body_meta->is_static && rigid_body_meta->active)
        {            
            rigid_body->polygon.position.x += rigid_body->velocity.x * rigid_body->damping * dt;
            rigid_body->polygon.position.y += rigid_body->velocity.y * rigid_body->damping * dt;
            rigid_body->polygon.orientation += rigid_body->angular_velocity * rigid_body->damping * dt*(!rigid_body->lock_orientation);
        }
    }       
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].is_empty)
        {
            mpl_manifold_positional_correction(&manifolds[i]);
        }
    }   
    for (int i = 0; i < body_count; ++i)
    {
        mpl_rigid_body_clear_forces(&rigid_bodies[i]);
    }
    for (int i = 0; i < body_count; ++i)
    {
        mpl_polygon_transform(&rigid_bodies[i].polygon);        
    }
}
//collision end