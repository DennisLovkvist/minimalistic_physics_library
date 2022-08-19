
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
void mpl_polygon_get_reverse_support_point(Vector2 *vertices,unsigned int vertex_count, Vector2 dir,Vector2 *best_vertex)
{
    float bestProjection = FLT_MAX; 
    best_vertex->x = 0;
    best_vertex->y = 0;
    for (int i = 0; i < vertex_count; ++i)
    {
        Vector2 v = vertices[i];
        float projection = vector2f_dot_vv(v, dir);

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
    for (int i = 0; i < vertex_count; i++)
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
        radius = d1;
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
void mpl_rigid_body_compute_mass_circle(RigidBody *rigid_body,float density,float radius)
{
    float area = PI * radius*radius;

    rigid_body->polygon.radius = radius;
    rigid_body->longest_radius = radius;

    rigid_body->mass = density * area;
    rigid_body->inertia = rigid_body->mass * radius * radius;
    rigid_body->invMass = 1.0f / rigid_body->mass;
    rigid_body->invInertia = 1.0f / rigid_body->inertia;
    
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
void mpl_rigid_body_init(RigidBody *rigid_body, unsigned int shape,unsigned int width,unsigned int height)
{    
    if(shape == 0)
    {
        rigid_body->is_polygon = 0;
        rigid_body->polygon.vertex_count = 0;
    }
    else if(shape == 1)
    {
        rigid_body->is_polygon = 1;
        mpl_polygon_make_triangle(&rigid_body->polygon,width, height);        
    }
    else if(shape == 2)
    {
        rigid_body->is_polygon = 1;
        mpl_polygon_make_square(&rigid_body->polygon,width, height);
    }
    else
    {
        rigid_body->is_polygon = 1;
        mpl_polygon_make_square(&rigid_body->polygon,width, height);
    }
    rigid_body->active = 1;
    rigid_body->is_colliding = 0;
    rigid_body->is_static = 0;
    rigid_body->lock_orientation = 0;
    rigid_body->ignore_gravity = 0;
    rigid_body->restitution = 0.2f;
    rigid_body->static_frictionf =  0.4f;
    rigid_body->dynamic_frictionf = 0.2f;
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
    rigid_body->polygon.contact_points[0].x = rigid_body->polygon.contact_points[0].y =  0;
    rigid_body->polygon.contact_points[1].x = rigid_body->polygon.contact_points[1].y =  0;
    rigid_body->polygon.contact_points[2].x = rigid_body->polygon.contact_points[2].y =  0;
    rigid_body->polygon.contact_points[3].x = rigid_body->polygon.contact_points[3].y =  0;

    if(rigid_body->is_polygon)
    {
        mpl_rigid_body_compute_mass(rigid_body,10);
    }
    else
    {
        mpl_rigid_body_compute_mass_circle(rigid_body,10,width/2);
    }

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
    if(isnan(impulse.x) || isnan(impulse.y))
    {
        return;
    }
    rigid_body->velocity.x += rigid_body->invMass * impulse.x;
    rigid_body->velocity.y += rigid_body->invMass * impulse.y;
    
    float n = vector2f_cross_vv(contactVector, impulse);
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
void mpl_manifold_init(Manifold *manifold,float dt)
{
    manifold->e = min(manifold->A->restitution, manifold->B->restitution);
    manifold->static_friction = sqrt(manifold->A->static_frictionf * manifold->B->static_frictionf);
    manifold->dynamic_friction = sqrt(manifold->A->dynamic_frictionf * manifold->B->dynamic_frictionf);
    
}
void mpl_manifold_apply_impulse(Manifold *manifold)
{
    for (int i = 0; i < manifold->contact_count; ++i)
    {
        Vector2 ra;
        ra.x = manifold->A->polygon.position.x - manifold->contacts[i].x;
        ra.y = manifold->A->polygon.position.y - manifold->contacts[i].y;

        Vector2 rb;
        rb.x = manifold->contacts[i].x-manifold->B->polygon.position.x;
        rb.y = manifold->contacts[i].y-manifold->B->polygon.position.y;         
        
        Vector2 cross_rb_av;
        vector2f_cross_vf(rb,manifold->B->angularVelocity,&cross_rb_av);
        Vector2 cross_ra_av;
        vector2f_cross_vf(ra,manifold->A->angularVelocity,&cross_ra_av);
        Vector2 rv;
        rv.x = manifold->B->velocity.x + cross_rb_av.x - manifold->A->velocity.x - cross_ra_av.x;
        rv.y = manifold->B->velocity.y + cross_rb_av.y - manifold->A->velocity.y - cross_ra_av.y;

        float contactVel = vector2f_dot_vv(rv, manifold->normal);
        
        if (contactVel > 0)
        {
            return;
        }

        float raCrossN = vector2f_cross_vv(ra, manifold->normal);
        float rbCrossN = vector2f_cross_vv(rb, manifold->normal);

        float invMassSum = manifold->A->invMass + manifold->B->invMass + (raCrossN * raCrossN) * manifold->A->invInertia + (rbCrossN * rbCrossN) * manifold->B->invInertia;
        
        float j = -(1.0f + manifold->e) * contactVel;                
        j /= invMassSum;
        j /= manifold->contact_count;


        Vector2 impulse;
        impulse.x = manifold->normal.x * j;
        impulse.y = manifold->normal.y * j;

       mpl_rigid_body_apply_impulse(manifold->B,impulse,rb);
        impulse.x = -impulse.x;
        impulse.y = -impulse.y;
        mpl_rigid_body_apply_impulse(manifold->A,impulse,ra);

        float dot = vector2f_dot_vv(rv, manifold->normal);
        Vector2 t;
        t.x = rv.x - (manifold->normal.x * dot);
        t.y = rv.y - (manifold->normal.y * dot);
        vector2f_normalize(&t);

        float jt = -vector2f_dot_vv(rv, t);
        jt /= invMassSum;
        jt /= manifold->contact_count;
        
        if (jt == 0.0f)
        {
            return;
        }

        Vector2 tangentImpulse;

        if (abs(jt) < j * manifold->static_friction)
        {
            tangentImpulse.x = t.x * jt;
            tangentImpulse.y = t.y * jt;
        }
        else
        {
            tangentImpulse.x = t.x * -j * manifold->dynamic_friction;
            tangentImpulse.y = t.y * -j * manifold->dynamic_friction;
        }
        mpl_rigid_body_apply_impulse(manifold->B,tangentImpulse,rb);
        tangentImpulse.x = -tangentImpulse.x;
        tangentImpulse.y = -tangentImpulse.y;
        mpl_rigid_body_apply_impulse(manifold->A,tangentImpulse,ra);    
    }
}
void mpl_manifold_positional_correction(Manifold *manifold)
{
    const float k_slop = 0.05f;
    const float percent = 0.4f;

    float f = max(manifold->penetration - k_slop, 0.0f)/(manifold->A->invMass + manifold->B->invMass);
    Vector2 correction;
    correction.x = f*manifold->normal.x * percent;
    correction.y = f*manifold->normal.y * percent;
    
    manifold->A->polygon.position.x -= correction.x * manifold->A->invMass;
    manifold->A->polygon.position.y -= correction.y * manifold->A->invMass;

    manifold->B->polygon.position.x += correction.x * manifold->B->invMass;
    manifold->B->polygon.position.y += correction.y * manifold->B->invMass;
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
        mpl_polygon_get_support_point(polygon_b->m_vertices,polygon_b->vertex_count,normal_flipped,&support_point);

        Vector2 vertex;
        vertex.x = polygon_a->m_vertices[i].x;
        vertex.y = polygon_a->m_vertices[i].y;

        Vector2 diff;
        diff.x = support_point.x - vertex.x;
        diff.y = support_point.y - vertex.y;

        float dot = vector2f_dot_vv(normal, diff);

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
int mpl_find_incident_face(Vector2 v[2], Polygon *RefPoly, Polygon *IncPoly, int referenceIndex)
{
    Vector2 reference_normal;
    
    reference_normal.x = RefPoly->normals[referenceIndex].x;
    reference_normal.y = RefPoly->normals[referenceIndex].y;

    matrix_mul_mv(RefPoly->matrix.matrix,reference_normal,&reference_normal);

    float transposed[3][3];
    matrix_transpose_mm(RefPoly->matrix.matrix,transposed);
    matrix_mul_mv(transposed,reference_normal,&reference_normal);

    int incidentFace = 0;
    float minDot = FLT_MAX;
    for (int i = 0; i < IncPoly->vertex_count; ++i)
    {        
        float dot = vector2f_dot_vv(reference_normal,IncPoly->normals[i]);

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

    float d1 = vector2f_dot_vv(n, face[0])-c;
    float d2 = vector2f_dot_vv(n, face[1])-c;

    if (d1 <= 0.0f) points[sp++] = face[0];
    if (d2 <= 0.0f) points[sp++] = face[1];

    // If the points are on different sides of the plane
    if (d1 * d2 < 0.0f) // less than to ignore -0.0f
    {
        // Push intersection point
        float alpha = d1 / (d1 - d2);
        points[sp].x = face[0].x  + alpha * (face[1].x  - face[0].x );
        points[sp].y = face[0].y  + alpha * (face[1].y  - face[0].y );
        ++sp;
    }

    // Assign our new converted values
    face[0] = points[0];
    face[1] = points[1];

    return sp;
}
unsigned int mpl_broad_phase(RigidBody *A, RigidBody *B)
{
    float dx = A->polygon.position.x - B->polygon.position.x;
    float dy = A->polygon.position.y - B->polygon.position.y;
    float dist = sqrt(dx * dx + dy * dy);

    return (dist < A->longest_radius + B->longest_radius);
}
float mpl_dist_sqr(Vector2 a, Vector2 b)
{
    Vector2 c;
    vector2f_sub_vv(a,b,&c);
    return vector2f_dot_vv(c,c);
}
unsigned int mpl_narrow_phase_c2c(Manifold *manifold, RigidBody *a, RigidBody *b)
{
    manifold->EMPTY = 0;
    manifold->A = a;
    manifold->B = b;
    Vector2 normal;
    vector2f_sub_vv(b->polygon.position, a->polygon.position,&normal);

    float dist_sqr = normal.x * normal.x + normal.y * normal.y;
    float radius = a->polygon.radius + a->polygon.radius;

    // Not in contact
    if (dist_sqr >= radius * radius)
    {
        manifold->contact_count = 0;
        manifold->EMPTY = 1;
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
        manifold->normal.x = normal.x / distance; // Faster than using Normalized since we already performed sqrt
        manifold->normal.y = normal.y / distance;
        manifold->contacts[0].x = manifold->normal.x * a->polygon.radius + a->polygon.position.x;
        manifold->contacts[0].y = manifold->normal.y * a->polygon.radius + a->polygon.position.y;
    }

    return 1;
}





unsigned int mpl_narrow_phase_p2c(Manifold *manifold, RigidBody *a, RigidBody *b)
{
    manifold->EMPTY = 0;
    manifold->A = a;
    manifold->B = b;

   Vector2 normal;
   vector2f_normal(b->polygon.position, a->polygon.position, &normal);

    vector2f_normalize(&normal);
    vector2f_cross_v(normal,&normal);

    Vector2 normal_flipped;
    vector2f_negate_vv(normal,&normal_flipped);


    float separation = FLT_MIN;
    int face_normal = 0;

    for (int i = 0; i < a->polygon.vertex_count; i++)
    {
        Vector2 delta;
        delta.x = b->polygon.position.x - a->polygon.m_vertices[i].x;
        delta.y = b->polygon.position.y - a->polygon.m_vertices[i].y;
        
        float s = vector2f_dot_vv(a->polygon.normals[i],delta);

        if (s > b->polygon.radius)
        {
            manifold->EMPTY = 1;
            return 0;
        }
        if (s > separation)
        {
            separation = s;
            face_normal = i;
        }

    }
    Vector2 v1 = a->polygon.m_vertices[a->polygon.m_edges[face_normal].m_index_a];
    Vector2 v2 = a->polygon.m_vertices[a->polygon.m_edges[face_normal].m_index_b];

    Vector2 delta1;
    delta1.x = b->polygon.position.x - v1.x;
    delta1.y = b->polygon.position.y - v1.y;
    Vector2 delta2;
    delta2.x = v2.x - v1.x;
    delta2.y = v2.y - v1.y;
    float dot1 = vector2f_dot_vv(delta1, delta2);

    delta1.x = b->polygon.position.x - v2.x;
    delta1.y = b->polygon.position.y - v2.y;
    delta2.x = v1.x - v2.x;
    delta2.y = v1.y - v2.y;    
    float dot2 = vector2f_dot_vv(delta1, delta2);

    manifold->contact_count = 1;
    manifold->penetration = b->polygon.radius - separation;

    if (separation < 0.0001f)
    {
        manifold->contact_count = 1;
        manifold->normal = a->polygon.normals[face_normal];
        manifold->contacts[0].x = manifold->normal.x * b->polygon.radius + b->polygon.position.x;
        manifold->contacts[0].y = manifold->normal.y * b->polygon.radius + b->polygon.position.y;

        manifold->penetration = b->polygon.radius;
        return 0;
    }
    if (dot1 <= 0.0f)
    {
        if (mpl_dist_sqr(b->polygon.position, v1) > b->polygon.radius * b->polygon.radius)
        {
            manifold->EMPTY = 1;
            return 0;
        }
        manifold->contacts[0] = v1;
        float r = atan2(b->polygon.position.y - v1.y, b->polygon.position.x - v1.x);
        Vector2 n;
        n.x = cos(r);
        n.y = sin(r);
        vector2f_normalize(&n);

        manifold->normal = n;
        Vector2 surface;
        surface.x = b->polygon.position.x + cos(r - PI) * b->polygon.radius;
        surface.y = b->polygon.position.y + sin(r - PI) * b->polygon.radius;

        float dx = surface.x - v1.x;
        float dy = surface.y - v1.y;
        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else if (dot2 <= 0.0f)
    {
        if (mpl_dist_sqr(b->polygon.position, v2) > b->polygon.radius * b->polygon.radius)
        {
            manifold->EMPTY = 1;
            return 0;
        }
        manifold->contacts[0] = v2;
        float r = atan2(b->polygon.position.y - v2.y, b->polygon.position.x - v2.x);
        Vector2 n;
        n.x = cos(r);
        n.x = sin(r);
        vector2f_normalize(&n);
        manifold->normal = n;
        Vector2 surface;
        surface.x = b->polygon.position.x + cos(r - PI) * b->polygon.radius;
        surface.y = b->polygon.position.y + sin(r - PI) * b->polygon.radius;

        float dx = surface.x - v2.x;
        float dy = surface.y - v2.y;
        manifold->penetration = sqrt(dx * dx + dy * dy);
    }
    else
    {
        manifold->normal = a->polygon.normals[face_normal];

        manifold->contacts[0].x = b->polygon.position.x - a->polygon.normals[face_normal].x * (b->polygon.radius - manifold->penetration);
        manifold->contacts[0].y = b->polygon.position.y - a->polygon.normals[face_normal].y * (b->polygon.radius - manifold->penetration);
    }

    manifold->B->polygon.contact_points[0].x = b->polygon.position.x;
    manifold->B->polygon.contact_points[0].y = b->polygon.position.y;
    manifold->B->polygon.contact_points[1].x = manifold->contacts[0].x;
    manifold->B->polygon.contact_points[1].y = manifold->contacts[0].y;

    manifold->B->polygon.contact_points[2].x = v1.x;
    manifold->B->polygon.contact_points[2].y = v1.y;
    manifold->B->polygon.contact_points[3].x = v2.x;
    manifold->B->polygon.contact_points[3].y = v2.y;

    a->is_colliding = 1;            
    b->is_colliding = 1;
    return 1;
}









unsigned int mpl_narrow_phase_p2p(Manifold *manifold, RigidBody *A, RigidBody *B)
{
    manifold->EMPTY = 0;
    manifold->A = A;
    manifold->B = B;

    manifold->contact_count = 0;

    int faceA = 0;
    int faceB = 0;                                            

    float penetrationA = mpl_find_axis_least_penetration(&A->polygon, &B->polygon,&faceA);
    float penetrationB = mpl_find_axis_least_penetration(&B->polygon, &A->polygon,&faceB);


    if (penetrationA >= 0.0f)
    {
        manifold->EMPTY = 1;
        return 0;
    }
    if (penetrationB >= 0.0f)
    {
        manifold->EMPTY = 1;
        return 0;
    }


    int reference_face_index = 0;
    unsigned int flip = 0; // Always point from a to b

    Polygon *reference_polygon; // Reference
    Polygon *incident_polygon; // Incident

    // Determine which shape contains reference face
    if (mpl_gt(penetrationA, penetrationB))
    {
        reference_polygon = &A->polygon;
        incident_polygon = &B->polygon;
        reference_face_index = faceA;
        flip = 0;
    }
    else
    {
        reference_polygon = &B->polygon;
        incident_polygon = &A->polygon;
        reference_face_index = faceB;
        flip = 1;
    }

    Vector2 reference_face[2];
    reference_face[0] = reference_polygon->m_vertices[reference_polygon->m_edges[reference_face_index].m_index_a];
    reference_face[1] = reference_polygon->m_vertices[reference_polygon->m_edges[reference_face_index].m_index_b];

    Vector2 vertex1 = reference_face[0];
    Vector2 vertex2 = reference_face[1];

    int incident_face_index = mpl_find_incident_face(reference_face, reference_polygon, incident_polygon, reference_face_index);
    
    Vector2 incident_face[2];
    incident_face[0] = incident_polygon->m_vertices[incident_polygon->m_edges[incident_face_index].m_index_a];
    incident_face[1] = incident_polygon->m_vertices[incident_polygon->m_edges[incident_face_index].m_index_b];


    // Calculate normals
    Vector2 side_plane_normal;
    side_plane_normal.x = vertex2.x - vertex1.x;
    side_plane_normal.y = vertex2.y - vertex1.y;
    vector2f_normalize(&side_plane_normal);

    Vector2 reference_face_normal;
    reference_face_normal.x = side_plane_normal.y;
    reference_face_normal.y = -side_plane_normal.x;       
    
    float dist_from_origin = vector2f_dot_vv(reference_face_normal, vertex1);
    float negative_side = -vector2f_dot_vv(side_plane_normal, vertex1);
    float positive_side = vector2f_dot_vv(side_plane_normal, vertex2);

    manifold->normal.x = reference_face_normal.x;
    manifold->normal.y = reference_face_normal.y;
    // Clip incident face to reference face side planes
    Vector2 side_plane_normal_negated;
    vector2f_negate_vv(side_plane_normal,&side_plane_normal_negated);

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

    manifold->A->polygon.contact_points[0].x = incident_face[0].x;
    manifold->A->polygon.contact_points[0].y = incident_face[0].y;
    manifold->A->polygon.contact_points[1].x = incident_face[1].x;
    manifold->A->polygon.contact_points[1].y = incident_face[1].y;


    manifold->A->polygon.contact_points[2].x = reference_face[0].x;
    manifold->A->polygon.contact_points[2].y = reference_face[0].y;    
    manifold->A->polygon.contact_points[3].x = reference_face[1].x;
    manifold->A->polygon.contact_points[3].y = reference_face[1].y;


    // Keep points behind reference face
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
    manifold->reference_polygon = reference_polygon;
    manifold->incident_polygon = incident_polygon;
    manifold->reference_face_index = reference_face_index;
    manifold->incident_face_index = incident_face_index;

    A->is_colliding = 1;            
    B->is_colliding = 1;

    return 1;
}
void mpl_update(RigidBody *bodies, unsigned int body_count,Manifold *manifolds, int iterations ,float dt, float G)
{
    unsigned int manifold_count = 0;
    unsigned int manifold_capacity = 100;
    Vector2 gravity;
    gravity.x = 0;
    gravity.y = 0;

    for (int i = 0; i < body_count; i++)
    {
        if (!bodies[i].is_static && bodies[i].active && !bodies[i].ignore_gravity)
        {                    
            gravity.y = round(G * bodies[i].mass);
            mpl_rigid_body_apply_force(&bodies[i],gravity);
        }
    }    
    for (int i = 0; i < body_count; i++)
    {
        if (bodies[i].active)
        {
            for (int j = i; j < body_count; ++j)
            {
                if (bodies[j].active)
                {
                    if (i != j)
                    {
                        if (!(bodies[i].is_static && bodies[j].is_static))
                        {
                            if(manifold_count < manifold_capacity)
                            {
                                
                                if (mpl_broad_phase(&bodies[i], &bodies[j]))
                                {
                                    Manifold m = manifolds[manifold_count];
                                    mpl_manifold_reset(&m);
                                    unsigned int flag = 0;
                                    if (!bodies[i].is_polygon)
                                    {
                                        if (!bodies[j].is_polygon)
                                        {
                                            flag =mpl_narrow_phase_c2c(&manifolds[manifold_count],&bodies[i], &bodies[j]);
                                        }
                                        else
                                        {
                                           flag = mpl_narrow_phase_p2c(&manifolds[manifold_count],&bodies[j], &bodies[i]);
                                        }
                                    }
                                    else
                                    {
                                        if (!bodies[j].is_polygon)
                                        {
                                            flag = mpl_narrow_phase_p2c(&manifolds[manifold_count],&bodies[i], &bodies[j]);
                                        }
                                        else
                                        {
                                            flag = mpl_narrow_phase_p2p(&manifolds[manifold_count],&bodies[i], &bodies[j]);
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
    for (int i = 0; i < body_count; ++i)
    {
        if (!bodies[i].is_static && bodies[i].active)
        {
            mpl_rigid_body_integrate_forces(&bodies[i],dt);
        }
    }  
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].EMPTY)
        {
            mpl_manifold_init(&manifolds[i],dt);
        }
    }
            
    for (int j = 0; j < iterations; ++j)
    {
        for (int i = 0; i < manifold_count; ++i)
        {
            if (!manifolds[i].EMPTY)
            {
               mpl_manifold_apply_impulse(&manifolds[i]);
            }
        }
    }
            
    for (int i = 0; i < body_count; ++i)
    {
        if (!bodies[i].is_static && bodies[i].active)
        {
            mpl_rigid_body_integrate_velocity(&bodies[i],dt);
        }
    }
      
    for (int i = 0; i < manifold_count; ++i)
    {
        if (!manifolds[i].EMPTY)
        {
            mpl_manifold_positional_correction(&manifolds[i]);
        }
    }
    

    for (int i = 0; i < body_count; ++i)
    {
        mpl_rigid_body_clear_forces(&bodies[i]);
    }
    for (int i = 0; i < body_count; ++i)
    {
        mpl_polygon_transform(&bodies[i].polygon);
    }
}
void mpl_rigid_body_set_static(RigidBody *rb)
{
    rb->inertia = 0.0f;
    rb->invInertia = 0.0f;
    rb->mass = 0.0f;
    rb->invMass = 0.0f;
    rb->lock_orientation = 1;
}
