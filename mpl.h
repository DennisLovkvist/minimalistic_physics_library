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
typedef struct Polygon Polygon;
struct Polygon
{
    unsigned int vertex_count;
    Vector2 m_vertices[4];
    Edge m_edges[4];
    Vector2 normals[4];
    RotationMatrix matrix;
    Vector2 position,position_last;
    float orientation,orientation_last,radius;
    Vector2 contact_points[4];
};
typedef struct RigidBody RigidBody;
struct RigidBody
{
    Polygon polygon;
    unsigned int is_polygon,is_colliding,is_static,active,lock_orientation,ignore_gravity;
    float restitution,static_frictionf,dynamic_frictionf;
    Vector2 velocity;
    Vector2 force;
    float angularVelocity,torque,orient;
    float mass, invMass, inertia, invInertia;
    float aprox_point_area,longest_radius,damping,velocity_epsilon;
};
typedef struct Manifold Manifold;
struct Manifold
{
        float e,static_friction,dynamic_friction,penetration;
        RigidBody *A;
        RigidBody *B;
        int contact_count,reference_face_index,incident_face_index;
        Polygon *reference_polygon;
        Polygon *incident_polygon;
        Vector2 normal;
        Vector2 contacts[2];
        unsigned int EMPTY;
        Vector2 test_point1,test_point2;
};
void mpl_polygon_init(Polygon *polygon, Vector2 *vertices, unsigned int vertex_count);
void mpl_polygon_transform(Polygon *polygon);
void mpl_rigid_body_init(RigidBody *rigid_body, unsigned int shape,unsigned int width,unsigned int height);
void mpl_update(RigidBody *bodies, unsigned int body_count,Manifold *manifolds, int iterations ,float dt, float G);
void mpl_rigid_body_set_position(RigidBody *rigid_body,float x, float y);
void mpl_rigid_body_set_static(RigidBody *rb);
void mpl_rigid_body_apply_force(RigidBody *rigid_body,Vector2 f);