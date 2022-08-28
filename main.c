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
#include <time.h>
#include "game.h"
#include <math.h>

#define GLEW_STATIC 1
#define WINDOW_TITLE "delo2d boilerplate"

const unsigned int screen_width = 800;
const unsigned int screen_height = 600;

int main(void)
{ 
    GLFWwindow *window;
    unsigned int shader_primitive;
    unsigned int shader_circle;
    Projection projection;
    PrimitiveBatch primitive_batch;
    SpriteBatch sprite_batch;

    if(delo2d_render_setup(&window, screen_width, screen_height,WINDOW_TITLE) == -1){return -1;}//setup and initialization for opengl
    
    delo2d_matrix_orthographic_projection(&projection,0.0f,(float)screen_width,0.0f,(float)screen_height,1,-1);//creates an orthographic_projection to be used as our camera
    
    glEnable(GL_BLEND);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE);
    glBlendEquation(GL_FUNC_ADD);

    glfwSwapInterval(0);

    shader_primitive = delo2d_shader_from_file("shaders/delo2d_primitive_default.glsl");
    shader_circle = delo2d_shader_from_file("shaders/delo2d_sprite_circle.glsl");
        
    float dt = 0;
    struct timeval t1, t2;
    double elapsedTime;

    delo2d_primitive_batch_create(&primitive_batch,4*7*4);
    
    delo2d_sprite_batch_create(&sprite_batch,CAPACITY_RIGID_BODIES);

    game_init();

    float total = 0;
    int steps = 0;

    while (!glfwWindowShouldClose(window))
    {  
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        gettimeofday(&t1, NULL);  
        
        game_update(xpos,ypos,dt);

        gettimeofday(&t2, NULL);
        
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;// sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;// us to ms        

        dt =  elapsedTime *= 0.001; 

        total += dt;
        steps += 1;
                     
        game_draw(window,&sprite_batch,&primitive_batch,projection,shader_primitive,shader_circle);
        
        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    delo2d_primitive_batch_delete(&primitive_batch);
    glDeleteProgram(shader_primitive);
    glfwTerminate();

    return 0;
}