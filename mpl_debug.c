#include "mpl.h"
#include "delo2d.h"

void mpl_debug_draw_polygon(PrimitiveBatch *primitive_batch,unsigned int shader_primitive,Projection projection, Polygon *polygon)
{
    delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_TRIANGLE_LIST);
        if(polygon->vertex_count == 3)
        {
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[0].x,polygon->m_vertices[0].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[1].x,polygon->m_vertices[1].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[2].x,polygon->m_vertices[2].y,1,1,1,0.5);
        }    
        else if(polygon->vertex_count == 4)
        {
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[0].x,polygon->m_vertices[0].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[1].x,polygon->m_vertices[1].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[2].x,polygon->m_vertices[2].y,1,1,1,0.5);
            
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[2].x,polygon->m_vertices[2].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[3].x,polygon->m_vertices[3].y,1,1,1,0.5);
            delo2d_primitive_batch_add(primitive_batch,polygon->m_vertices[0].x,polygon->m_vertices[0].y,1,1,1,0.5);
        }    
    delo2d_primitive_batch_end(primitive_batch);



            delo2d_primitive_batch_begin(primitive_batch,shader_primitive,projection,DELO_LINE_LIST);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[0].x,polygon->contact_points[0].y,1,0,0,1);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[1].x,polygon->contact_points[1].y,1,0,0,1); 
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[2].x,polygon->contact_points[2].y,0,1,0,1);
                delo2d_primitive_batch_add(primitive_batch,polygon->contact_points[3].x,polygon->contact_points[3].y,0,1,0,1);            
            delo2d_primitive_batch_end(primitive_batch);
}