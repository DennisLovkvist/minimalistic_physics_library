#include "mpl.h"
#include "delo2d.h"
#include <math.h>

void mpl_debug_draw_polygon(PrimitiveBatch *primitive_batch,unsigned int shader_primitive,Projection projection, Polygon *polygon)
{
    Color color;
    delo2d_color_set_f(&color,0.04,0.6,1,0.5);
    delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_TRIANGLE_LIST);
        if(polygon->vertex_count == 3)
        {
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[0].x,polygon->vertices[0].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[1].x,polygon->vertices[1].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[2].x,polygon->vertices[2].y,color.r,color.g,color.b,color.a);
        }    
        else if(polygon->vertex_count == 4)
        {
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[0].x,polygon->vertices[0].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[1].x,polygon->vertices[1].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[2].x,polygon->vertices[2].y,color.r,color.g,color.b,color.a);
            
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[2].x,polygon->vertices[2].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[3].x,polygon->vertices[3].y,color.r,color.g,color.b,color.a);
            delo2d_primitive_batch_add(primitive_batch,polygon->vertices[0].x,polygon->vertices[0].y,color.r,color.g,color.b,color.a);
        } 
        else if(polygon->vertex_count == 0)
        {

            Vector2 prev;
            float r = -polygon->orientation;
            prev.x = polygon->position.x + cos(r)*polygon->radius;
            prev.y = polygon->position.y + sin(r)*polygon->radius;

            float inc = (3.14*2) /24;

            for (int i = 0; i <23; i++)
            {               
                
                
                delo2d_primitive_batch_add(primitive_batch,polygon->position.x,polygon->position.y,color.r,color.g,color.b,color.a);

                delo2d_primitive_batch_add(primitive_batch,prev.x,prev.y,color.r,color.g,color.b,color.a);
                r += inc;
                float x = polygon->position.x + cos(r)*polygon->radius;
                float y = polygon->position.y + sin(r)*polygon->radius;
                delo2d_primitive_batch_add(primitive_batch,x,y,color.r,color.g,color.b,color.a);


                
                prev.x = x;
                prev.y = y;

                 
            }
            
        }   
    delo2d_primitive_batch_end(primitive_batch);

/*

            delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_LINE_LIST);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[0].x,polygon->contact_points[0].y,1,0,0,1);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[1].x,polygon->contact_points[1].y,1,0,0,1); 
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[2].x,polygon->contact_points[2].y,0,1,0,1);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[3].x,polygon->contact_points[3].y,0,1,0,1);            
            delo2d_primitive_batch_end(primitive_batch);
            */
}