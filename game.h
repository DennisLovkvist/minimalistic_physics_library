#pragma once

#define CAPACITY_RIGID_BODIES 200
#define CAPACITY_MANIFOLDS 200*24

int game_init();
void game_update(double xpos, double ypos, float dt);
void game_draw(GLFWwindow *window,SpriteBatch *sprite_batch,PrimitiveBatch *primitive_batch,Projection projection,unsigned int shader_primitive,unsigned int shader_circle);
