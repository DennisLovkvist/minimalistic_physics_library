/**
 * Author:    Dennis Lövkvist
 * Created:   2022-08-05
 * Version: 1.0
 **/
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "delo2d.h"
#include "stb_image.h"
#include <math.h>

void GLClearError()
{
    while(!glGetError());
}
void GLCheckError()
{
    GLenum error;
    while((error = glGetError()))
    {
        printf("%i",error);
        printf("%c",'\n');
    }
}
void delo2d_rectangle_set(Rectangle *rectengle, int x, int y,int width, int height)
{
    rectengle->x = x;
    rectengle->y = y;
    rectengle->width = width;
    rectengle->height = height;
}
int delo2d_render_setup(GLFWwindow **window, unsigned int width, unsigned int height,const char *title)
{
    if (!glfwInit()){return -1;} 

    *window = glfwCreateWindow(width, height, title, NULL, NULL);

    if (!*window)
    {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(*window);
    
    if(glewInit() != GLEW_OK)
    {
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
    glEnable(GL_TEXTURE_2D);

    return 0;
}
//region render target begin
void delo2d_render_target_create(RenderTarget *rt,float screen_width,float screen_height)
{
    rt->vertices[0] = 1.0f;rt->vertices[1] = -1.0f;rt->vertices[2] = 1.0f;rt->vertices[3] = 0.0f;
    rt->vertices[4] = -1.0f;rt->vertices[5] = -1.0f;rt->vertices[6] = 0.0f;rt->vertices[7] = 0.0f;
    rt->vertices[8] = -1.0f;rt->vertices[9] = 1.0f;rt->vertices[10] = 0.0f;rt->vertices[11] = 1.0f;

    rt->vertices[12] = 1.0f;rt->vertices[13] = 1.0f;rt->vertices[14] = 1.0f;rt->vertices[15] = 1.0f;
    rt->vertices[16] = 1.0f;rt->vertices[17] = -1.0f;rt->vertices[18] = 1.0f;rt->vertices[19] = 0.0f;
    rt->vertices[20] = -1.0f;rt->vertices[21] = 1.0f;rt->vertices[22] = 0.0f;rt->vertices[23] = 1.0f;

	glGenVertexArrays(1, &rt->vao);
	glGenBuffers(1, &rt->vbo);

	glBindVertexArray(rt->vao);
	glBindBuffer(GL_ARRAY_BUFFER, rt->vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(rt->vertices), &rt->vertices, GL_STATIC_DRAW);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));

	glGenFramebuffers(1, &rt->fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, rt->fbo);

	// Create Framebuffer Texture
	glGenTextures(1, &rt->fbt);
	glBindTexture(GL_TEXTURE_2D, rt->fbt);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen_width,  screen_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // Prevents edge bleeding
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Prevents edge bleeding
    
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rt->fbt, 0);

	glGenRenderbuffers(1, &rt->rbo);
	glBindRenderbuffer(GL_RENDERBUFFER, rt->rbo);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, screen_width,  screen_height);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rt->rbo);

    rt->status = glCheckFramebufferStatus(GL_FRAMEBUFFER);

	glBindVertexArray(0);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
    rt->initialized = 1;
    
}
void delo2d_render_target_delete(RenderTarget *rt)
{
    if(rt->initialized == 0)return;
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glDeleteVertexArrays(1, &rt->vao);
    glDeleteBuffers(1,&rt->vbo);
    glDeleteRenderbuffers(1, &rt->rbo);
    glDeleteFramebuffers(2,&rt->fbo);
	glDeleteTextures(1, &rt->fbt);
    rt->initialized = 0;
}
void delo2d_render_target_draw(RenderTarget *render_target, unsigned int shader_id)
{   
    glBindVertexArray(render_target->vao); 
    delo2d_texture_bind(render_target->fbt,0); 
    glUseProgram(shader_id); 
    glDrawArrays(GL_TRIANGLES, 0, 6);
}
void delo2d_render_target_set(unsigned int frame_buffer,float r, float g, float b, float a)
{
    glBindFramebuffer(GL_FRAMEBUFFER, frame_buffer);
    glClearColor(r,g,b,a);
    glClear(GL_COLOR_BUFFER_BIT);  
}
//region rendering end

//texture code begin
void delo2d_texture_load(Texture *texture, char file_path[])
{
    stbi_set_flip_vertically_on_load(0);
    texture->local_buffer = stbi_load(file_path,&texture->width,&texture->height,&texture->bytes_per_pixel,4);
    glGenTextures(1,&texture->renderer_id);
    glBindTexture(GL_TEXTURE_2D,texture->renderer_id);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);

    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,texture->width,texture->height,1,GL_RGBA,GL_UNSIGNED_BYTE,texture->local_buffer);

    glBindTexture(GL_TEXTURE_2D,0);
    if(texture->local_buffer)
    {
        stbi_image_free(texture->local_buffer);
    }
    texture->initialized = 1;
}
void delo2d_texture_delete(Texture *texture)
{
    if(texture->initialized == 0)return;
    glBindTexture(GL_TEXTURE_2D,0);
    glDeleteTextures(1,&texture->renderer_id);
    texture->initialized = 0;    
}
void delo2d_texture_bind(unsigned int texture, unsigned int slot)
{
    glActiveTexture(GL_TEXTURE0 + slot);
    glBindTexture(GL_TEXTURE_2D,texture);
}
void delo2d_texture_unbind()
{    
    glBindTexture(GL_TEXTURE_2D,0);
}
void delo2d_texture_copy(Texture *texture_src,Texture *texture_des)
{
    texture_des->bytes_per_pixel = texture_src->bytes_per_pixel;
    texture_des->width = texture_src->width;
    texture_des->height = texture_src->height;
    texture_des->renderer_id = texture_src->renderer_id;
    texture_des->local_buffer = texture_src->local_buffer;
    texture_des->initialized = texture_src->initialized;
}
//texture code end

//region matrices begin
void delo2d_matrix_mul_vector2fp_matrix33(Vector2fp *vector,float (*R)[3][3])
{
    float x =  ((*vector->x) * (*R)[0][0]) + ((*vector->y) * (*R)[0][1]) + (1 * (*R)[0][2]);
    float y = ((*vector->x) * (*R)[1][0]) + ((*vector->y) * (*R)[1][1]) + (1 * (*R)[1][2]);
    (*vector->x) = x;
    (*vector->y) = y;
}
void delo2d_matrix_orthographic_projection(Projection *projection, float l,float r,float t,float b,float f,float n)
{
    projection->matrix[0][0] = 2.0f/(r-l);  //Scale?  
    projection->matrix[0][1] = 0;              
    projection->matrix[0][2] = 0;                    
    projection->matrix[0][3] = -((r+l)/(r-l));

    projection->matrix[1][0] = 0; 
    projection->matrix[1][1] = 2.0/(t-b); //Scale?             
    projection->matrix[1][2] = 0;                    
    projection->matrix[1][3] = -((t+b)/(t-b));

    projection->matrix[2][0] = 0; 
    projection->matrix[2][1] = 0;              
    projection->matrix[2][2] = 2 / (f-n);                    
    projection->matrix[2][3] = -((f+n)/(f-n));

    projection->matrix[3][0] = 0; 
    projection->matrix[3][1] = 0;              
    projection->matrix[3][2] = 0;                    
    projection->matrix[3][3] = 1;
}
void delo2d_rotation_matrix(float (*R)[3][3],float theta, float tx, float ty)
{
    (*R)[0][0] = cos(theta);
    (*R)[0][1] = sin(theta);
    (*R)[0][2] = tx;

    (*R)[1][0] = -sin(theta);
    (*R)[1][1] = cos(theta);
    (*R)[1][2] = ty;

    (*R)[2][0] = 0;
    (*R)[2][1] = 0;
    (*R)[2][2] = 1;
}
void delo2d_projection_matrix_set(Projection *projection_src,Projection *projection_des)
{
    float *m0 = projection_des->matrix[0];
    float *m1 = projection_src->matrix[0];

    for (size_t i = 0; i < 16; i++)
    {
        m0[i] = m1[i];
    }    
}
//region matrices end

//region quads begin
void delo2d_quad_get(Quad *quad, VertexArray *vertex_array, int element_index)
{    
    int quad_index = element_index*vertex_array->layout_float_count*4;
    int stride = vertex_array->layout_float_count;

    quad->v0.x = &(vertex_array->buffer_position[quad_index]);
    quad->v0.y = &(vertex_array->buffer_position[quad_index+1]);

    quad->v1.x = &(vertex_array->buffer_position[quad_index + (stride)]);
    quad->v1.y = &(vertex_array->buffer_position[quad_index + (stride+1)]); 

    quad->v2.x = &(vertex_array->buffer_position[quad_index + (stride*2)]);
    quad->v2.y = &(vertex_array->buffer_position[quad_index + (stride*2+1)]); 

    quad->v3.x = &(vertex_array->buffer_position[quad_index + (stride*3)]);
    quad->v3.y = &(vertex_array->buffer_position[quad_index + (stride*3+1)]); 
}
void delo2d_quad_define(VertexArray *vertex_array, int quad_index, Rectangle_f *rect_des,Rectangle_f *rect_src, int texture_index,Color color,int flip_horizontally,int flip_vertically)
{
    int n = quad_index * 4;

    delo2d_sprite_vertex_set_element_sprite(vertex_array,n,rect_des->x,rect_des->y,(rect_src->x + rect_src->width*flip_horizontally),(rect_src->y + rect_src->height*flip_vertically),texture_index,color);
    delo2d_sprite_vertex_set_element_sprite(vertex_array,n + 1,rect_des->x + rect_des->width,rect_des->y,rect_src->x + rect_src->width * !flip_horizontally,rect_src->y + rect_src->height * flip_vertically,texture_index,color);
    delo2d_sprite_vertex_set_element_sprite(vertex_array,n + 2,rect_des->x + rect_des->width,rect_des->y + rect_des->height,rect_src->x + rect_src->width * !flip_horizontally,rect_src->y + rect_src->height * !flip_vertically,texture_index,color);    
    delo2d_sprite_vertex_set_element_sprite(vertex_array,n + 3,rect_des->x,rect_des->y + rect_des->height,rect_src->x + rect_src->width*flip_horizontally ,rect_src->y + rect_src->height * !flip_vertically,texture_index,color);

}
void delo2d_quad_translate(Quad *quad,float tx, float ty)
{
    *(quad->v0.x) += tx;
    *(quad->v0.y) += ty;
    *(quad->v1.x) += tx;
    *(quad->v1.y) += ty;
    *(quad->v2.x) += tx;
    *(quad->v2.y) += ty;
    *(quad->v3.x) += tx;
    *(quad->v3.y) += ty;
}
void delo2d_quad_set_position(Quad *quad,int x,int y)
{
    Vector2f center;
    delo2d_quad_get_center(quad,&center);

    Vector2f delta;
    delta.x = x-center.x;
    delta.y = y-center.y;

    delo2d_quad_translate(quad,delta.x,delta.y);
}
void delo2d_quad_get_center(Quad *quad,Vector2f *center)
{
    center->x = (*quad->v0.x + *quad->v1.x + *quad->v2.x + *quad->v3.x)*0.25f;
    center->y = (*quad->v0.y + *quad->v1.y + *quad->v2.y + *quad->v3.y)*0.25f;
}
void delo2d_quad_rotate(Quad *quad, float theta)
{
    Vector2f center;
    center.x = (*(quad->v0.x) + *(quad->v1.x) + *(quad->v2.x) + *(quad->v3.x))/4;
    center.y = (*(quad->v0.y) + *(quad->v1.y) + *(quad->v2.y) + *(quad->v3.y))/4;


    float R[3][3];
    delo2d_rotation_matrix(&R, 0,-center.x,-center.y);

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);

    delo2d_rotation_matrix(&R, theta,0,0);

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);

    delo2d_rotation_matrix(&R, 0,center.x,center.y);

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);
}
void delo2d_quad_rotate_around_point(Quad *quad, float theta,float point_x, float point_y)
{
    float R[3][3];
    delo2d_rotation_matrix(&R, 0,-(point_x),-(point_y));

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);

    delo2d_rotation_matrix(&R, theta,0,0);

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);

    delo2d_rotation_matrix(&R, 0,(point_x),(point_y));

    delo2d_matrix_mul_vector2fp_matrix33(&quad->v0, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v1, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v2, &R);
    delo2d_matrix_mul_vector2fp_matrix33(&quad->v3, &R);
}
void delo2d_quad_scale(Quad *quad,float scale_x,float scale_y)
{
    float dx = *quad->v0.x - *quad->v1.x;
    float diff = (dx-(dx*scale_x))/2;
    *quad->v0.x -= diff;
    *quad->v1.x += diff;
    *quad->v3.x -= diff;
    *quad->v2.x += diff;

    float dy = *quad->v0.y - *quad->v3.y;
    diff = (dy-(dy*scale_y))/2;
    *quad->v0.y -= diff;
    *quad->v3.y += diff;
    *quad->v1.y -= diff;
    *quad->v2.y += diff;

}
void delo2d_quad_skew(Quad *quad,float skew_x,float skew_y,float pivot_point_x,float pivot_point_y)
{
    Vector2f center;
    delo2d_quad_get_center(quad,&center);

    *quad->v0.x -= skew_x * (center.y + pivot_point_y-(*quad->v0.y));
    *quad->v1.x -= skew_x * (center.y + pivot_point_y-(*quad->v1.y));

    *quad->v2.x -= skew_x * (center.y + pivot_point_y-(*quad->v2.y));
    *quad->v3.x -= skew_x * (center.y + pivot_point_y-(*quad->v3.y));

    *quad->v0.y -= skew_y * (center.x + pivot_point_x-(*quad->v0.x));
    *quad->v3.y -= skew_y * (center.x + pivot_point_x-(*quad->v3.x));

    *quad->v1.y += skew_y * (center.x + pivot_point_x-(*quad->v1.x));
    *quad->v2.y += skew_y * (center.x + pivot_point_x-(*quad->v2.x));
}
//region quads end

//vertex array code begin
void delo2d_sprite_vertex_array_create(VertexArray *vertex_array,unsigned int type, unsigned int element_count)
{
    vertex_array->layout_float_count = 9;
    vertex_array->type = type;
    vertex_array->count_elements = element_count;
    vertex_array->count_position = element_count * (type*vertex_array->layout_float_count);
    vertex_array->buffer_position = malloc(vertex_array->count_position * sizeof(float));

    vertex_array->indices_per_element = 6;
    
    int length = vertex_array->count_position;
    for (size_t i = 0; i < length; i++)
    {
        vertex_array->buffer_position[i] = 0.0f;
    } 

    vertex_array->count_index = vertex_array->indices_per_element * element_count;
    vertex_array->buffer_index = malloc(vertex_array->count_index * sizeof(unsigned int));            
            
    for (size_t i = 0; i < element_count; i++)
    {

        int index = i * 6;
        int vertex_index = i*4;
        vertex_array->buffer_index[index + 0] = vertex_index + 0;
        vertex_array->buffer_index[index + 1] = vertex_index + 1;
        vertex_array->buffer_index[index + 2] = vertex_index + 2;
        vertex_array->buffer_index[index + 3] = vertex_index + 2;
        vertex_array->buffer_index[index + 4] = vertex_index + 3;
        vertex_array->buffer_index[index + 5] = vertex_index + 0;
    }   

    glGenVertexArrays(1,&vertex_array->vao); 
    glBindVertexArray(vertex_array->vao); 

    glGenBuffers(1,&vertex_array->vbo);
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);
    glBufferData(GL_ARRAY_BUFFER,vertex_array->count_position * sizeof(float),NULL,GL_DYNAMIC_DRAW);     
    glEnableVertexAttribArray(0);//vertex position float2
    glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE, sizeof(float)*vertex_array->layout_float_count,0);    
    glEnableVertexAttribArray(1);//texture coordinate float2
    glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE, sizeof(float)*vertex_array->layout_float_count,(GLvoid*)(sizeof(float)*2));    
    glEnableVertexAttribArray(2);//texture index float
    glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE, sizeof(GLfloat)*vertex_array->layout_float_count,(GLvoid*)(4 * sizeof(GLfloat)));    
    glEnableVertexAttribArray(3);//color float4
    glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE, sizeof(GLfloat)*vertex_array->layout_float_count,(GLvoid*)(5 * sizeof(GLfloat)));

    glGenBuffers(1,&vertex_array->ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,vertex_array->ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,vertex_array->count_index * sizeof(GLuint),vertex_array->buffer_index,GL_DYNAMIC_DRAW);  

    glBindVertexArray(0); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);

    vertex_array->initialized = 1;
}
void delo2d_sprite_vertex_array_delete(VertexArray *vertex_array)
{    
    if(vertex_array->initialized == 0)return;
    delo2d_sprite_vertex_array_unbind(vertex_array);

    glDeleteVertexArrays(1,&vertex_array->vao); 
    glDeleteBuffers(1,&vertex_array->vbo);
    glDeleteBuffers(1,&vertex_array->ibo);

    free(vertex_array->buffer_position);
    free(vertex_array->buffer_index);
    vertex_array->initialized = 0;
}
void delo2d_sprite_vertex_set_element_sprite(VertexArray *vertex_array, int element,float x, float y, float tex_x,float tex_y,unsigned int texture_slot,Color color)
{
    int index = element * vertex_array->layout_float_count;
    vertex_array->buffer_position[index + 0] = x;
    vertex_array->buffer_position[index + 1] = y;
    vertex_array->buffer_position[index + 2] = tex_x;
    vertex_array->buffer_position[index + 3] = tex_y;
    vertex_array->buffer_position[index + 4] = texture_slot;
    vertex_array->buffer_position[index + 5] = color.r;
    vertex_array->buffer_position[index + 6] = color.g;
    vertex_array->buffer_position[index + 7] = color.b;
    vertex_array->buffer_position[index + 8] = color.a;
}
void delo2d_sprite_vertex_set_tex_data(VertexArray *vertex_array, int element,float tex_x,float tex_y,unsigned int texture_slot)
{
    int index = element * vertex_array->layout_float_count;
    vertex_array->buffer_position[index + 2] = tex_x;
    vertex_array->buffer_position[index + 3] = tex_y;
    vertex_array->buffer_position[index + 4] = texture_slot;
}
void delo2d_sprite_vertex_array_draw(VertexArray *vertex_array,unsigned int count_elements,unsigned int shader_id,Texture *textures,int texture_count,Projection projection)
{
   
    glUseProgram(shader_id); 
    glUniformMatrix4fv(glGetUniformLocation(shader_id,"u_mvp"),1,GL_FALSE,projection.matrix[0]);
    
    int samplers[3] = {0,1,2};
    glUniform1iv(glGetUniformLocation(shader_id,"u_textures"),3,samplers);  

    glBindVertexArray(vertex_array->vao);
     
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);//Bind to update with delo2d_vertex_array_to_graphics_device
    delo2d_sprite_vertex_array_to_graphics_device(vertex_array,0);
    for (size_t i = 0; i < texture_count; i++)
    {      
        delo2d_texture_bind(textures[i].renderer_id,i); 
    }
    
    glUseProgram(shader_id); 

    int draw_index_count = vertex_array->indices_per_element * count_elements;    

    glDrawElements(GL_TRIANGLES,draw_index_count,GL_UNSIGNED_INT,NULL);

    delo2d_texture_unbind();

}
void delo2d_sprite_vertex_array_to_graphics_device(VertexArray *vertex_array, GLintptr offset)
{   
    //delo2d_vertex_array_bind(vertex_array);
    glBufferSubData(GL_ARRAY_BUFFER,offset,vertex_array->count_position * sizeof(float),vertex_array->buffer_position);
}
void delo2d_sprite_vertex_array_bind(VertexArray *vertex_array)
{
    glBindVertexArray(vertex_array->vao);
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,vertex_array->ibo);
}
void delo2d_sprite_vertex_array_unbind()
{
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0); 
}
//vertex array code end

//shader code begin
static unsigned int delo2d_shader_compile(unsigned int type,char *shader_source_code)
{    
    unsigned int id = glCreateShader(type);
    char const* src = shader_source_code;
    glShaderSource(id,1,&src,NULL);
    glCompileShader(id);
    return id;
}
static int delo2d_shader_create(char *vertex_shader_source_code, char *fragment_shader_source_code)
{
    unsigned int program = glCreateProgram();
    unsigned int vs = delo2d_shader_compile(GL_VERTEX_SHADER,vertex_shader_source_code);
    unsigned int fs = delo2d_shader_compile(GL_FRAGMENT_SHADER,fragment_shader_source_code);

    glAttachShader(program,vs);
    glAttachShader(program,fs);

    glLinkProgram(program);
    glValidateProgram(program);

    glDeleteShader(vs);
    glDeleteShader(fs);

    return program;
}
static char* delo2d_shader_load(char *path)
{
    FILE *f = fopen(path, "rb");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *source_code = malloc(fsize + sizeof(char));
    fread(source_code, fsize, 1, f);
    fclose(f);
    source_code[fsize] = '\0';
    return source_code;
}
static int delo2d_find_keyword(char *string, char *sub_string,int tag)
{
    int n = 0;
    size_t length = strlen(string);
    for (size_t i = 0; i < length; i++)
    {
        if(string[i] == sub_string[n])
        {
            n++;
        }
        else
        {
            n = 0;
        }
        if(n == strlen(sub_string))
        {            
            if(tag == 0)
            {
                return i+1;
            }
            else
            {
                return i - (strlen(sub_string)-1);
            }
        }
    }
    return 0;
}
static char* delo2d_shader_parse(char *source_full,char *keyword_begin,char *keyword_end)
{
    int begin = delo2d_find_keyword(source_full,keyword_begin,0);
    int end = delo2d_find_keyword(source_full,keyword_end,1);

    if(end > begin)
    {
        int length = end - begin;
        char *src = malloc(sizeof(char)*length+1);
        memcpy(src,&source_full[begin],sizeof(char)*(length));
        src[length] = '\0'; 
        return src; 
    }
    else
    {
        return NULL;
    }
}
unsigned int delo2d_shader_from_file(char *path_shader)
{
    char *src = delo2d_shader_load(path_shader); 
    char *src_vert = delo2d_shader_parse(src,"#VERT_BEGIN","#VERT_END");
    char *src_frag = delo2d_shader_parse(src,"#FRAG_BEGIN","#FRAG_END");    

    if(src_vert != NULL && src_frag != NULL)
    {
        unsigned int shader = delo2d_shader_create(src_vert,src_frag); 
        
        free(src);
        free(src_vert);
        free(src_frag);
        return shader;
    }
    else
    {
        return 0;
    }    
}
//shader code end

//spritebatch code begin
void delo2d_sprite_batch_create(SpriteBatch *sprite_batch,int capacity)
{

    delo2d_sprite_vertex_array_create(&sprite_batch->vertex_array,DELO_QUAD_LIST,capacity);

    sprite_batch->capacity = capacity;
    sprite_batch->count = 0;
    sprite_batch->rect_des = malloc(sizeof(Rectangle_f)*capacity);
    sprite_batch->rect_src = malloc(sizeof(Rectangle_f)*capacity);
    sprite_batch->texture_index = malloc(sizeof(unsigned int)*capacity);
    sprite_batch->quad_index = malloc(sizeof(unsigned int)*capacity);
    sprite_batch->flip_horizontally = malloc(sizeof(unsigned int)*capacity);
    sprite_batch->flip_vertically = malloc(sizeof(unsigned int)*capacity);
    sprite_batch->updated = malloc(sizeof(unsigned int)*capacity);
    sprite_batch->color = malloc(sizeof(Color)*capacity);
    sprite_batch->scale = malloc(sizeof(Vector2f)*capacity);
    sprite_batch->skew = malloc(sizeof(Vector2f)*capacity);
    sprite_batch->position = malloc(sizeof(Vector2f)*capacity);
    sprite_batch->pivot_point = malloc(sizeof(Vector2f)*capacity);
    sprite_batch->orientation = malloc(sizeof(float)*capacity);
    sprite_batch->textures = malloc(sizeof(Texture)*32);
    sprite_batch->texture_count = 0;
    sprite_batch->called_begin = 0;
    sprite_batch->called_end = 1;
    delo2d_matrix_orthographic_projection(&sprite_batch->projection,0.0f,1920,0.0f,1080,1,-1); 

    for (size_t i = 0; i < capacity; i++)
    {
        sprite_batch->texture_index[i] = 0;
        sprite_batch->quad_index[i] = 0;
        sprite_batch->rect_des[i].x = sprite_batch->rect_des[i].y = 0;
        sprite_batch->rect_des[i].width = sprite_batch->rect_des[i].height = 100;

        sprite_batch->rect_src[i].x = sprite_batch->rect_src[i].y = 0;
        sprite_batch->rect_src[i].width = sprite_batch->rect_src[i].height = 1;
        sprite_batch->color[i].r = sprite_batch->color[i].g = sprite_batch->color[i].b = sprite_batch->color[i].a = 0;
        sprite_batch->flip_horizontally[i] = sprite_batch->flip_vertically[i] = 0;
        sprite_batch->scale[i].x = sprite_batch->scale[i].y = 1;
        sprite_batch->skew[i].x = sprite_batch->skew[i].y = 0;
        sprite_batch->pivot_point[i].x = sprite_batch->pivot_point[i].y = 0;
        sprite_batch->position[i].x = sprite_batch->position[i].y = 0;
        sprite_batch->updated[i] = 1;
    }    
    sprite_batch->initialized = 1;
}
void delo2d_sprite_batch_delete(SpriteBatch *sprite_batch)
{
    if(sprite_batch->initialized == 0)return;

    free(sprite_batch->rect_des);
    free(sprite_batch->rect_src);
    free(sprite_batch->texture_index);
    free(sprite_batch->quad_index);
    free(sprite_batch->flip_horizontally);
    free(sprite_batch->flip_vertically);
    free(sprite_batch->updated);
    free(sprite_batch->color);
    free(sprite_batch->scale);
    free(sprite_batch->skew);
    free(sprite_batch->position);
    free(sprite_batch->pivot_point);
    free(sprite_batch->orientation);
    free(sprite_batch->textures);

    delo2d_sprite_vertex_array_delete(&sprite_batch->vertex_array);

    sprite_batch->initialized = 0;
}
void delo2d_sprite_batch_begin(SpriteBatch *sprite_batch,unsigned int shader, Projection projection)
{
    if(sprite_batch->called_end == 0)return;

    sprite_batch->shader_id = shader;
    delo2d_projection_matrix_set(&projection,&sprite_batch->projection);

    sprite_batch->called_begin = 1;
}
void delo2d_sprite_batch_end(SpriteBatch *sprite_batch)
{
    if(sprite_batch->called_begin == 0)return;
    
    delo2d_sprite_batch_to_vertex_array(sprite_batch,&sprite_batch->vertex_array);     

    delo2d_sprite_vertex_array_draw(&sprite_batch->vertex_array,sprite_batch->count,sprite_batch->shader_id,sprite_batch->textures,sprite_batch->texture_count,sprite_batch->projection);

    sprite_batch->called_end = 1;
    sprite_batch->texture_count = 0;
    sprite_batch->count = 0;
    
}
void delo2d_sprite_define(Sprite *sprite, int dx, int dy,int dw, int dh,int sx, int sy,int sw, int sh,unsigned int texture_index, unsigned int texture_width, unsigned int texture_height, unsigned int stride,unsigned int frames, float duration, Color color,float scale_x,float scale_y,float skew_x,float skew_y,unsigned int flip_horizontally,unsigned int flip_vertically)
{
    sprite->frame = 0;
    sprite->time = 0;
    sprite->updated_tex_coords = 0;
    sprite->duration = duration;
    sprite->rect_des.x = dx;
    sprite->rect_des.y = dy;
    sprite->rect_des.width = dw;
    sprite->rect_des.height = dh;

    sprite->position.x = dx + dw*0.5f;
    sprite->position.y = dy + dh*0.5f;

    sprite->rect_src.x = sx;
    sprite->rect_src.y = sy;
    sprite->rect_src.width = sw;
    sprite->rect_src.height = sh;

    sprite->quad_index = 0;
    sprite->stride = stride;
    sprite->frames = frames;
    sprite->texture_index = texture_index;
    sprite->texture_width = texture_width;
    sprite->texture_height = texture_height;
    sprite->color.r = color.r;
    sprite->color.g = color.g;
    sprite->color.b = color.b;
    sprite->color.a = color.a;
    sprite->flip_horizontally = flip_horizontally;
    sprite->flip_vertically = flip_vertically;
    sprite->updated_tex_coords = 0;
    sprite->orientation = 0;
    sprite->scale.x = scale_x;
    sprite->scale.y = scale_y;
    sprite->skew.x = skew_x;
    sprite->skew.y = skew_y;
    sprite->offset.x = 0;
    sprite->offset.y = 0;
    sprite->pivot_point.x = sprite->pivot_point.y = 0;
    sprite->loop = 1;
}
void delo2d_sprite_batch_add(SpriteBatch *sprite_batch, Sprite *sprite,Texture *texture)
{
    int texture_index = -1;
    delo2d_sprite_batch_add_texture(sprite_batch,texture, &texture_index);

    int index = sprite_batch->count;
    sprite_batch->count ++;
    
    sprite_batch->rect_src[index].x = ((float)sprite->rect_src.x);
    sprite_batch->rect_src[index].y = ((float)sprite->rect_src.y);
    sprite_batch->rect_src[index].width = ((float)sprite->rect_src.width);
    sprite_batch->rect_src[index].height = ((float)sprite->rect_src.height);
   

    sprite_batch->rect_des[index].x = sprite->rect_des.x;
    sprite_batch->rect_des[index].y = sprite->rect_des.y;
    sprite_batch->rect_des[index].width = sprite->rect_des.width;
    sprite_batch->rect_des[index].height = sprite->rect_des.height;

    sprite->batch_index = index;
    sprite->quad_index = index;
    sprite_batch->texture_index[index] = texture_index; 
    sprite_batch->color[index].r = sprite->color.r;
    sprite_batch->color[index].g = sprite->color.g;
    sprite_batch->color[index].b = sprite->color.b;
    sprite_batch->color[index].a = sprite->color.a;
    sprite_batch->flip_horizontally[index] = sprite->flip_horizontally;
    sprite_batch->flip_vertically[index] = sprite->flip_vertically;

    sprite_batch->scale[index].x = sprite->scale.x;
    sprite_batch->scale[index].y = sprite->scale.y;

    sprite_batch->skew[index].x = sprite->skew.x;
    sprite_batch->skew[index].y = sprite->skew.y;

    sprite_batch->pivot_point[index].x = sprite->pivot_point.x;
    sprite_batch->pivot_point[index].y = sprite->pivot_point.y;

    sprite_batch->position[index].x = sprite->position.x;
    sprite_batch->position[index].y = sprite->position.y;

    sprite_batch->orientation[index] = sprite->orientation;
    sprite_batch->updated[index] = 1;
}
void delo2d_sprite_batch_add_texture(SpriteBatch *sprite_batch,Texture *texture, int *texture_index)
{
    
    *texture_index = -1;
    for (size_t i = 0; i < sprite_batch->texture_count; i++)
    {
        if(sprite_batch->textures[i].renderer_id == texture->renderer_id)
        {
            *texture_index = i;
            break;
        }
    }
    if(*texture_index == -1)
    { 
        delo2d_texture_copy(texture,&sprite_batch->textures[sprite_batch->texture_count]);
        *texture_index = sprite_batch->texture_count;
        sprite_batch->texture_count ++;
    }
}
void delo2d_sprite_batch_to_vertex_array(SpriteBatch *sprite_batch,VertexArray *vertex_array)
{    
    int length = sprite_batch->count;
    for(int i = 0; i < length; i++)
    {
        delo2d_quad_define(vertex_array,i,&(sprite_batch->rect_des)[i],&(sprite_batch->rect_src)[i],sprite_batch->texture_index[i],sprite_batch->color[i],sprite_batch->flip_horizontally[i],sprite_batch->flip_vertically[i]);
        Quad quad;
        delo2d_quad_get(&quad,vertex_array,i);

        delo2d_quad_set_position(&quad,sprite_batch->position[i].x,sprite_batch->position[i].y);
        
        sprite_batch->updated[i] = 0;
        //delo2d_quad_skew_top(&quad,sprite_batch->skew[i].x);
        
        delo2d_quad_skew(&quad,sprite_batch->skew[i].x,sprite_batch->skew[i].y,sprite_batch->pivot_point[i].x,sprite_batch->pivot_point[i].y);

        delo2d_quad_scale(&quad,sprite_batch->scale[i].x,sprite_batch->scale[i].y);

        Vector2f center;
        delo2d_quad_get_center(&quad,&center);
        delo2d_quad_rotate_around_point(&quad,sprite_batch->orientation[i], center.x + sprite_batch->pivot_point[i].x, center.y + sprite_batch->pivot_point[i].y);          
    }    
}
void delo2d_sprite_scale_dest_rect(Sprite *sprite, float scale_x, float scale_y)
{
    sprite->rect_des.width *= scale_x;
    sprite->rect_des.height *= scale_y;
    sprite->position.x = sprite->rect_des.x+sprite->rect_des.width*scale_x;
    sprite->position.y = sprite->rect_des.y+sprite->rect_des.height*scale_y;
} 
void delo2d_sprite_rotate(Sprite *sprite,float rotation,VertexArray *vertex_array)
{
    Quad quad;
    delo2d_quad_get(&quad,vertex_array,sprite->quad_index);
    delo2d_quad_rotate(&quad,rotation);
    sprite->orientation += rotation;
}
void delo2d_sprite_rotate_around_point(Sprite *sprite,float rotation,float point_x, float point_y,VertexArray *vertex_array)
{
    Quad quad;
    delo2d_quad_get(&quad,vertex_array,sprite->quad_index);
    delo2d_quad_rotate_around_point(&quad,rotation, point_x, point_y);
    sprite->orientation += rotation;
}
void delo2d_sprite_set_orientation_around_point(Sprite *sprite,float orientation,float point_x, float point_y,VertexArray *vertex_array)
{
     float delta = orientation - sprite->orientation;
    delo2d_sprite_rotate_around_point(sprite,delta,point_x,point_y,vertex_array);
}
void delo2d_sprite_set_orientation(Sprite *sprite,float orientation,VertexArray *vertex_array)
{
    float delta = orientation - sprite->orientation;
    delo2d_sprite_rotate(sprite,delta,vertex_array);
}
void delo2d_sprite_translate(Sprite *sprite,float tx,float ty,VertexArray *vertex_array)
{
    Quad quad;
    delo2d_quad_get(&quad,vertex_array,sprite->quad_index);
    delo2d_quad_translate(&quad,tx,ty);

    sprite->position.x += tx;
    sprite->position.y += ty;

    //delo2d_quad_get_center(&quad,&sprite->position);    
}
void delo2d_sprite_animate(Sprite *sprite,float dt)
{
    sprite->time +=dt;

    if(sprite->time > sprite->duration)
    {
        if(sprite->loop == 1)
        {
            sprite->time = 0;
        }
        else
        {
            return;
        }
    }

    sprite->frame = (sprite->time / sprite->duration)*(float)sprite->frames;

    sprite->rect_src.x = sprite->offset.x + (sprite->frame % sprite->stride) * sprite->rect_src.width;
    
    sprite->rect_src.y = sprite->offset.y + (sprite->frame/sprite->stride) * sprite->rect_src.height;

    sprite->updated_tex_coords = 1;      
}
//spritebatch code end

//camera code begin
void delo2d_camera_move(Projection *projection, float tx, float ty, float screen_width, float screen_height)
{
    *projection->matrix[3] += tx/screen_width;
    *projection->matrix[7] -= ty/screen_height;
}
void delo2d_camera_set_position(Projection *projection, float x, float y, float screen_width, float screen_height)
{
    *projection->matrix[3] = x/screen_width;
    *projection->matrix[7] = y/screen_height;
}
void delo2d_camera_set_zoom(Projection *projection, float z, float screen_width, float screen_height)
{    
    *projection->matrix[0] = z/screen_width;
    *projection->matrix[5] = z/(-screen_height);   
}
//camera code end
 
//color code begin
void delo2d_color_set_f(Color *color,float r, float g, float b, float a)
{
    color->r = r;
    color->g = g;
    color->b = b;
    color->a = a;
}
void delo2d_color_set_i(Color *color,int r, int g, int b, int a)
{
    color->r = (float)r/255.0f;
    color->g = (float)g/255.0f;
    color->b = (float)b/255.0f;
    color->a = (float)a/255.0f;
}
void delo2d_color_lerp(Color *result, Color *color_a,Color *color_b, float facor)
{
    result->r = (color_a->r + color_b->r ) / 2; 
    result->g = (color_a->g + color_b->g ) / 2; 
    result->b = (color_a->b + color_b->b ) / 2; 
    result->a = (color_a->a + color_b->a ) / 2; 
}
//color code end

void delo2d_input_update(GLFWwindow *window, KeyboardInput *ki,KeyboardInput *ki_prev)
{
    ki->move_up = glfwGetKey(window, GLFW_KEY_W);
    ki->move_l = glfwGetKey(window, GLFW_KEY_A);
    ki->move_dn = glfwGetKey(window, GLFW_KEY_S);
    ki->move_r = glfwGetKey(window, GLFW_KEY_D);

    ki_prev->move_up = ki->move_up;
    ki_prev->move_l = ki->move_l;
    ki_prev->move_dn = ki->move_dn;
    ki_prev->move_r = ki->move_r;
}   
void delo2d_input_init(KeyboardInput *ki,KeyboardInput *ki_prev)
{
    ki_prev->move_up = GLFW_RELEASE;
    ki_prev->move_l = GLFW_RELEASE;
    ki_prev->move_dn = GLFW_RELEASE;
    ki_prev->move_r = GLFW_RELEASE;
} 

//vertex array primitives code begin
void delo2d_primitive_vertex_set_element(VertexArrayPrimitives *vertex_array, int element,float x, float y,float r,float g, float b, float a)
{
    int index = element * vertex_array->layout_float_count;
    vertex_array->buffer_position[index + 0] = x;
    vertex_array->buffer_position[index + 1] = y;
    vertex_array->buffer_position[index + 2] = r;
    vertex_array->buffer_position[index + 3] = g;
    vertex_array->buffer_position[index + 4] = b;
    vertex_array->buffer_position[index + 5] = a;

}
void delo2d_primitive_vertex_array_draw(VertexArrayPrimitives *vertex_array,unsigned int shader_id,Projection projection)
{   
    glUseProgram(shader_id); 
    glUniformMatrix4fv(glGetUniformLocation(shader_id,"u_mvp"),1,GL_FALSE,projection.matrix[0]);

    glBindVertexArray(vertex_array->vao);
     
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);//Bind to update with delo2d_vertex_array_to_graphics_device
    delo2d_primitive_vertex_array_to_graphics_device(vertex_array,0);
 
    if(vertex_array->type == DELO_TRIANGLE_LIST)
    {
        glDrawArrays(GL_TRIANGLES,0,vertex_array->count);
    }
    else if (vertex_array->type == DELO_LINE_LIST)
    {
        glDrawArrays(GL_LINES,0,vertex_array->count);
    }  
}
void delo2d_primitive_vertex_array_to_graphics_device(VertexArrayPrimitives *vertex_array, GLintptr offset)
{   
    delo2d_primitive_vertex_array_bind(vertex_array);
    glBufferSubData(GL_ARRAY_BUFFER,offset,vertex_array->count * vertex_array->layout_float_count * sizeof(float),vertex_array->buffer_position);
}
void delo2d_primitive_vertex_array_bind(VertexArrayPrimitives *vertex_array)
{
    glBindVertexArray(vertex_array->vao);
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);
}
void delo2d_primitive_vertex_array_unbind()
{
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
}
void delo2d_primitive_vertex_array_create(VertexArrayPrimitives *vertex_array,unsigned int vertex_capacity, unsigned int shader)
{
    vertex_array->count = 0;
    vertex_array->capacity = vertex_capacity;
    vertex_array->layout_float_count = 6;
    vertex_array->buffer_position = malloc(vertex_capacity * sizeof(float));
    
    int length = vertex_capacity;
    for (size_t i = 0; i < length; i++)
    {
        vertex_array->buffer_position[i] = 0.0f;
    } 
   
    glGenVertexArrays(1,&vertex_array->vao); 
    glBindVertexArray(vertex_array->vao); 

    glGenBuffers(1,&vertex_array->vbo);
    glBindBuffer(GL_ARRAY_BUFFER,vertex_array->vbo);
    glBufferData(GL_ARRAY_BUFFER,vertex_capacity * sizeof(float),NULL,GL_DYNAMIC_DRAW);   

    glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE, sizeof(GLfloat)*vertex_array->layout_float_count,0);   
    glEnableVertexAttribArray(0);//vertex position float2    
    glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE, sizeof(GLfloat)*vertex_array->layout_float_count,(GLvoid*)(2 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);//color float4

    glUseProgram(shader);

    glBindBuffer(GL_ARRAY_BUFFER,0);
    glBindVertexArray(0); 
    vertex_array->initialized = 1;
}
void delo2d_primitive_vertex_array_delete(VertexArrayPrimitives *vertex_array)
{
    if(vertex_array->initialized == 0)return;
    delo2d_sprite_vertex_array_unbind(vertex_array);

    glDeleteVertexArrays(1,&vertex_array->vao); 
    glDeleteBuffers(1,&vertex_array->vbo);

    free(vertex_array->buffer_position);
    vertex_array->initialized = 0;
}
//vertex array primitives code end
//primitive batch code begin
void delo2d_primitive_batch_create(PrimitiveBatch *primitive_batch,int capacity)
{
    primitive_batch->vertex_array.capacity = capacity;
    primitive_batch->vertex_array.count = 0;
    delo2d_primitive_vertex_array_create(&primitive_batch->vertex_array,capacity,primitive_batch->shader_id);
    primitive_batch->initialized = 1;
}
void delo2d_primitive_batch_delete(PrimitiveBatch *primitive_batch)
{
    if(primitive_batch->initialized == 0)return;    
    primitive_batch->vertex_array.capacity = 0;
    primitive_batch->vertex_array.count = 0;
    delo2d_primitive_vertex_array_delete(&primitive_batch->vertex_array);
    primitive_batch->initialized = 0;    
}
void delo2d_primitive_batch_begin(PrimitiveBatch *primitive_batch,unsigned int shader,Projection projection,unsigned int primitive_type)
{
    if(primitive_batch->called_end == 0)return;
    primitive_batch->vertex_array.type = primitive_type;
    primitive_batch->shader_id = shader;
    delo2d_projection_matrix_set(&projection,&primitive_batch->projection);

    primitive_batch->called_begin = 1;
}
void delo2d_primitive_batch_add(PrimitiveBatch *primitive_batch,int x, int y,float r,float g, float b, float a)
{
    int index = primitive_batch->vertex_array.count;
    primitive_batch->vertex_array.count ++;
    delo2d_primitive_vertex_set_element(&primitive_batch->vertex_array,index,x,y,r,g,b,a);    
}
void delo2d_primitive_batch_end(PrimitiveBatch *primitive_batch)
{
    if(primitive_batch->called_begin == 0)return;  

    delo2d_primitive_vertex_array_to_graphics_device(&primitive_batch->vertex_array,0);
    delo2d_primitive_vertex_array_draw(&primitive_batch->vertex_array,primitive_batch->shader_id,primitive_batch->projection);

    primitive_batch->called_end = 1;
    primitive_batch->vertex_array.count = 0;
}
//primitive batch code end
