#pragma once
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
typedef struct RotationMatrix RotationMatrix;
struct RotationMatrix
{
    float matrix[3][3];
    float input[3][1];
    float output[3][1];
};
typedef struct RigidBody RigidBody;
struct RigidBody
{
    unsigned int *is_polygon;
    unsigned int *is_colliding;
    unsigned int *is_static;
    unsigned int *active;
    unsigned int *lock_orientation;
    float *ignore_gravity;
    float *restitution;
    float *static_friction;
    float *dynamic_friction;
    float *velocity_x;
    float *velocity_y;
    float *force_x;
    float *force_y;
    float *angular_velocity;
    float *torque;
    float *orient;
    float *mass;
    float *inverse_mass;
    float *inertia;
    float *inverse_inertia;
    float *aprox_point_area;
    float *longest_radius;
    float *damping;
    float *velocity_epsilon;

    unsigned int *vertex_count;
    float *vertices_x;
    float *vertices_y;

    float *edge_vertex_a_x;
    float *edge_vertex_b_y;
    int *edge_index_a;
    int *edge_index_b;

    float *normals_x;
    float *normals_y;

    float *matrix;
    float *position_x;
    float *position_y;
    float *position_last_x;
    float *position_last_y;

    float *orientation;
    float *orientation_last;
    float *radius;

    float *contact_points_x;
    float *contact_points_y;
};
typedef struct Manifold Manifold;
struct Manifold
{
        float e,static_friction,dynamic_friction,penetration;
        unsigned int rigid_body_a;
        unsigned int rigid_body_b;
        int contact_count,reference_face_index,incident_face_index;
        unsigned int reference_polygon;
        unsigned int incident_polygon;
        Vector2 normal;
        Vector2 contacts[2];
        unsigned int is_empty;
        Vector2 test_point1,test_point2;
};
void mpl_polygon_init(RigidBody *rigid_body,unsigned int index, Vector2 *vertices, unsigned int vertex_count);
void mpl_polygon_transform(RigidBody*rigid_body, unsigned int index);
void mpl_rigid_bodies_init(RigidBody *rigid_body, unsigned int capacity);
void mpl_rigid_body_init(RigidBody *rigid_body,unsigned int index, unsigned int shape,unsigned int width,unsigned int height);
void mpl_update(RigidBody *rigid_bodies, unsigned int body_count,Manifold *manifolds, unsigned int manifold_capacity, int iterations ,float dt, float G);
void mpl_rigid_body_set_position(RigidBody *rigid_body,unsigned int index,float x, float y);
void mpl_rigid_body_set_static(RigidBody *rb, unsigned int index);
void mpl_rigid_body_apply_force(RigidBody *rigid_body, unsigned int index,Vector2 f);