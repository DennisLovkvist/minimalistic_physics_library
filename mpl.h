#pragma once

#define FLT_MAX 3.4028235E+38F /* max value */
#define FLT_MIN -3.4028235E+38F /* min positive value */
#define PI 3.14159265358979323846

typedef struct Vector2 Vector2;
struct Vector2
{
    float x;
    float y;
};
typedef struct Edge Edge;
struct Edge
{
    Vector2 m_vertex_a,m_vertex_b;
    int m_index_a,m_index_b;
}; 
typedef struct Polygon Polygon;
struct Polygon
{
    unsigned int vertex_count;
    Vector2 vertices[4];
    Edge edges[4];
    Vector2 normals[4];
    float matrix[3][3];
    Vector2 position,position_last;
    float orientation,orientation_last,radius;
};
typedef struct RigidBody RigidBody;
struct RigidBody
{
    Polygon polygon;
    float restitution,static_friction,dynamic_friction;
    Vector2 velocity;
    Vector2 force;
    float angular_velocity,torque;
    float mass, inverse_mass, inertia, inverse_inertia;
    float longest_radius,damping,velocity_epsilon;
    unsigned int lock_orientation;
};
typedef struct RigidBodyMeta RigidBodyMeta;
struct RigidBodyMeta
{
    unsigned int is_polygon,is_static,active,ignore_gravity,max_collisions,collisions;
    //unsigned int is_colliding;
};
typedef struct Manifold Manifold;
struct Manifold
{
        float e,static_friction,dynamic_friction,penetration;
        RigidBody *rigid_body_a;
        RigidBody *rigid_body_b;
        int contact_count,reference_face_index,incident_face_index;
        Polygon *reference_polygon;
        Polygon *incident_polygon;
        Vector2 normal;
        Vector2 contacts[2];
        unsigned int is_empty;
};
void mpl_polygon_init(Polygon *polygon, Vector2 *vertices, unsigned int vertex_count);
void mpl_polygon_transform(Polygon *polygon);
void mpl_rigid_body_init(RigidBody *rigid_body,RigidBodyMeta *rigid_body_meta, unsigned int shape,unsigned int width,unsigned int height,unsigned int max_collisions);
void mpl_update(RigidBody *bodies,RigidBodyMeta *rbm, unsigned int body_count,Manifold *manifolds, unsigned int manifold_capacity, int iterations ,float dt, float G);
void mpl_rigid_body_set_position(RigidBody *rigid_body,float x, float y);
void mpl_rigid_body_set_static(RigidBody *rb,RigidBodyMeta *rbm);
void mpl_rigid_body_apply_force(RigidBody *rigid_body,Vector2 f);