#VERT_BEGIN
#version 330 core

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 tex_coord;

out vec2 v_tex_coord;

void main()
{
    gl_Position = vec4(position.x, position.y, 0.0, 1.0); 
    v_tex_coord = tex_coord;
}  
#VERT_END
#FRAG_BEGIN
#version 330 core
out vec4 color;
in vec2 v_tex_coord;

uniform sampler2D u_texture;

void main()
{
    color = texture(u_texture,v_tex_coord);   
}
#FRAG_END