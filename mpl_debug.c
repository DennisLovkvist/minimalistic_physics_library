#include "mpl.h"
#include "delo2d.h"
#include <math.h>

void mpl_debug_draw_polygon(PrimitiveBatch *primitive_batch,unsigned int shader_primitive,Projection projection, RigidBody *rigid_body, unsigned int count)
{
    for (size_t i = 0; i < count; i++)
    {    
        int v_index = i * 4;
        delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_TRIANGLE_LIST);
        if(rigid_body->vertex_count[i] == 3)
        {
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+0],rigid_body->vertices_y[v_index+0],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+1],rigid_body->vertices_y[v_index+1],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+2],rigid_body->vertices_y[v_index+2],1,1,1,0.5);
        }    
        else if(rigid_body->vertex_count[i] == 4)
        {
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+0],rigid_body->vertices_y[v_index+0],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+1],rigid_body->vertices_y[v_index+1],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+2],rigid_body->vertices_y[v_index+2],1,1,1,0.5);
            
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+2],rigid_body->vertices_y[v_index+2],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+3],rigid_body->vertices_y[v_index+3],1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->vertices_x[v_index+0],rigid_body->vertices_y[v_index+0],1,1,1,0.5);
        } 
        else if(rigid_body->vertex_count[i] == 0)
        {        
            Vector2 center = {rigid_body->position_x[i],rigid_body->position_y[i]};    
            
            float r = 0;
            float inc = (3.14*2) /24;
            float radius = rigid_body->radius[i];
            for (int i = 0; i <23; i++)
            {    
                delo2d_primitive_batch_add(primitive_batch,center.x,center.y,1,1,1,0.5);

                delo2d_primitive_batch_add(primitive_batch,center.x + cos(r)*radius,center.y + sin(r)*radius,1,1,1,0.5);
                r += inc;
                delo2d_primitive_batch_add(primitive_batch,center.x + cos(r)*radius,center.y + sin(r)*radius,1,1,1,0.5);         
            } 
                       
        }   
        delo2d_primitive_batch_end(primitive_batch);

/*
        delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_LINE_LIST);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->contact_points_x[v_index+0],rigid_body->contact_points_y[v_index+0],1,0,0,1);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->contact_points_x[v_index+1],rigid_body->contact_points_y[v_index+1],1,0,0,1); 
            delo2d_primitive_batch_add(primitive_batch,rigid_body->contact_points_x[v_index+2],rigid_body->contact_points_y[v_index+2],0,1,0,1);
            delo2d_primitive_batch_add(primitive_batch,rigid_body->contact_points_x[v_index+3],rigid_body->contact_points_y[v_index+3],0,1,0,1);            
        delo2d_primitive_batch_end(primitive_batch);*/
    }
}