/**
 * Author:    Dennis Lövkvist
 * Created:   2022-08-05
 * Version: 1.0
 **/

#include <GL/glew.h>//glew fetches function pointers from graphics card
#include <GLFW/glfw3.h>//glfw cross-platform windows
#include <stdio.h>
#include <stdlib.h>
#include "delo2d.h"
#include <sys/time.h>
#include "mpl.h"
#include "delo_math.h"
#include "mpl_debug.h"
#include <math.h>
#include <time.h>

#define GLEW_STATIC 1
#define COUNT_SHADERS 7
#define COUNT_TEXTURES 3
#define WINDOW_TITLE "delo2d boilerplate"

const unsigned int screen_width = 800;
const unsigned int screen_height = 600;

int main(void)
{ 
    GLFWwindow *window;
    unsigned int shader_primitive;
    Projection projection;
    PrimitiveBatch primitive_batch;

    if(delo2d_render_setup(&window, screen_width, screen_height,WINDOW_TITLE) == -1){return -1;}//setup and initialization for opengl
    
    delo2d_matrix_orthographic_projection(&projection,0.0f,(float)screen_width,0.0f,(float)screen_height,1,-1);//creates an orthographic_projection to be used as our camera
    
    glEnable(GL_BLEND);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE);
    glBlendEquation(GL_FUNC_ADD);

    glfwSwapInterval(0);

    shader_primitive = delo2d_shader_from_file("shaders/delo2d_primitive_default.glsl");
    
    Color color_white;
    delo2d_color_set_f(&color_white,1,1,1,1);

    float dt = 0;
    struct timeval t1, t2;
    double elapsedTime;

    delo2d_primitive_batch_create(&primitive_batch,4000);

    int rbs = 100;
    int mfs = 1000;
    Manifold *manifolds = malloc(sizeof(Manifold)*mfs);


    RigidBody rb;
    mpl_rigid_bodies_init(&rb,100);

    mpl_rigid_body_init(&rb,0,2,20,screen_height);
    mpl_rigid_body_init(&rb,1,2,20,screen_height);
    mpl_rigid_body_init(&rb,2,2,screen_width-20,20);

    for (size_t i = 3; i < rbs; i++)
    {
        float x = rand()%(((screen_width-40)+1)-40) + 40;
        float y = (rand()%((200+1)-40) + 40)-500;
        float s = rand()%((3)-0) + 0;
        mpl_rigid_body_init(&rb,i,2,20,20);
        mpl_rigid_body_set_position(&rb,i,x,y);
    }
    

    mpl_rigid_body_set_position(&rb,0,0,screen_height/2);
    mpl_rigid_body_set_position(&rb,1,screen_width,screen_height/2);
    mpl_rigid_body_set_position(&rb,2,screen_width/2,400);

    srand(time(NULL));
    for (size_t i = 3; i < rbs; i++)
    {        
        float x = rand()%(((screen_width-40)+1)-40) + 40;
        float y = (rand()%((200+1)-40) + 40);
        float s = rand()%((3)-0) + 0;
        mpl_rigid_body_init(&rb,i,s,20,20);
        mpl_rigid_body_set_position(&rb,i,x,y);
    }

    mpl_rigid_body_set_static(&rb,0); 
    mpl_rigid_body_set_static(&rb,1); 
    mpl_rigid_body_set_static(&rb,2); 
            
      
    while (!glfwWindowShouldClose(window))
    { 
        gettimeofday(&t1, NULL);

        mpl_update(&rb, rbs,manifolds,mfs, 1 ,dt, 600.0f);  
            
        
        delo2d_render_target_set(0,0,0,0,1);//mpl


        mpl_debug_draw_polygon(&primitive_batch,shader_primitive,projection,&rb,rbs);


        glfwSwapBuffers(window);

        gettimeofday(&t2, NULL);
        
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;// sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;// us to ms        

        dt =  elapsedTime *= 0.001; 


        glfwPollEvents();
    }

    delo2d_primitive_batch_delete(&primitive_batch);
    glDeleteProgram(shader_primitive);
    glfwTerminate();

    return 0;
}