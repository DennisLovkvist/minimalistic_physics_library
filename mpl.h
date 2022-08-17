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
    unsigned int is_colliding,active,is_static,lock_orientation,ignore_gravity;
    float restitution,static_frictionf,dynamic_frictionf;
    Vector2 velocity;
    Vector2 force;
    float angularVelocity,torque,orient;
    float mass, invMass, inertia, invInertia;
    float aprox_point_area,longest_radius,damping,velocity_epsilon;
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
    Vector2 edge[2];
};
typedef struct Manifold Manifold;
struct Manifold
{
        float e,static_friction,dynamic_friction,penetration;
        RigidBody A,B;
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