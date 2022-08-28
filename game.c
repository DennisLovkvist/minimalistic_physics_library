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
#include "mpl.h"
#include "mpl_debug.h"
#include "game.h"
#include <math.h>

Manifold *manifolds;    
RigidBody *rigid_bodies;
RigidBodyMeta *rigid_bodies_meta;

Sprite sprite_circle;

int game_init()
{ 
  
    manifolds = malloc(sizeof(Manifold)*CAPACITY_MANIFOLDS);    
    rigid_bodies = malloc(sizeof(RigidBody)*CAPACITY_RIGID_BODIES);
    rigid_bodies_meta = malloc(sizeof(RigidBodyMeta)*CAPACITY_RIGID_BODIES);

    int screen_width = 800;
    int screen_height = 600;    

    mpl_rigid_body_init(&rigid_bodies[0],&rigid_bodies_meta[0],2,20,screen_height,100);
    mpl_rigid_body_init(&rigid_bodies[1],&rigid_bodies_meta[1],2,20,screen_height,100);
    mpl_rigid_body_init(&rigid_bodies[2],&rigid_bodies_meta[2],2,screen_width-20,20,100);
    mpl_rigid_body_init(&rigid_bodies[3],&rigid_bodies_meta[3],2,screen_width-20,20,100);

    mpl_rigid_body_init(&rigid_bodies[4],&rigid_bodies_meta[4],2,50,50,100);

    mpl_rigid_body_set_position(&rigid_bodies[0],0,screen_height/2);
    mpl_rigid_body_set_position(&rigid_bodies[1],screen_width,screen_height/2);
    mpl_rigid_body_set_position(&rigid_bodies[2],screen_width/2,600);
    mpl_rigid_body_set_position(&rigid_bodies[3],screen_width/2,0);

    mpl_rigid_body_set_position(&rigid_bodies[4],screen_width/2,screen_height/2);
    rigid_bodies[4].lock_orientation = 1;
    rigid_bodies_meta[4].ignore_gravity = 1;

    srand(time(NULL));
    for (size_t i = 5; i < CAPACITY_RIGID_BODIES; i++)
    {        
        float x = rand()%(((screen_width-40)+1)-40) + 40;
        float y = (rand()%((400+1)-300) + 300);
        mpl_rigid_body_init(&rigid_bodies[i],&rigid_bodies_meta[i],0,20,20,6);
        mpl_rigid_body_set_position(&rigid_bodies[i],x,y);
    }  

    mpl_rigid_body_set_static(&rigid_bodies[0],&rigid_bodies_meta[0]);
    mpl_rigid_body_set_static(&rigid_bodies[1],&rigid_bodies_meta[1]); 
    mpl_rigid_body_set_static(&rigid_bodies[2],&rigid_bodies_meta[2]); 
    mpl_rigid_body_set_static(&rigid_bodies[3],&rigid_bodies_meta[3]);          

    int radius = 20;
    Color color;
    delo2d_color_set_f(&color,0.04,0.6,1,0.5);
    delo2d_sprite_define(&sprite_circle,0,0,radius,radius,0,0,radius,radius,0,radius,radius,0,0,0,color,1,1,0,0,0,0);
    
}
void game_update(double xpos, double ypos, float dt)
{
    Vector2 delta = {(xpos-rigid_bodies[4].polygon.position.x),(ypos-rigid_bodies[4].polygon.position.y)};
    float r = atan2(delta.y,delta.x);
    float d = sqrt(delta.x * delta.x + delta.y * delta.y);
    Vector2 v = {cos(r)*d*20,sin(r)*d*20};
    rigid_bodies[4].velocity.x = v.x;
    rigid_bodies[4].velocity.y = v.y;
  
    mpl_update(rigid_bodies,rigid_bodies_meta, CAPACITY_RIGID_BODIES,manifolds,CAPACITY_MANIFOLDS, 2 ,0.0016f, 12000.0f);            
}

void game_draw(GLFWwindow *window,SpriteBatch *sprite_batch,PrimitiveBatch *primitive_batch,Projection projection,unsigned int shader_primitive,unsigned int shader_circle)
{   
    delo2d_render_target_set(0,0,0,0,1);

    delo2d_sprite_batch_begin(sprite_batch,shader_circle,projection);
    glUseProgram(shader_circle); 
    glUniform1f(glGetUniformLocation(shader_circle,"radius"),10);

    for (size_t i = 5; i < CAPACITY_RIGID_BODIES; i++)
    {
        sprite_circle.position.x = rigid_bodies[i].polygon.position.x;
        sprite_circle.position.y = rigid_bodies[i].polygon.position.y;
        sprite_circle.orientation = rigid_bodies[i].polygon.orientation;
        delo2d_sprite_batch_add_no_texture(sprite_batch,&sprite_circle);
    }      
    delo2d_sprite_batch_end(sprite_batch);

    for (size_t i = 0; i < 5; i++)
    {
        mpl_debug_draw_polygon(primitive_batch,shader_primitive,projection,&rigid_bodies[i].polygon);
    }

}
