#VERT_BEGIN
#version 330 core

layout(location = 0) in vec4 position;
layout(location = 1) in vec2 tex_coord;
layout(location = 2) in float tex_index;
layout(location = 3) in vec4 color;

out vec2 v_tex_coord;
out vec4 v_color;

uniform mat4 u_mvp;

void main()
{     
    gl_Position = position * u_mvp;
    v_tex_coord = tex_coord;
    v_color = color;
}

#VERT_END
#FRAG_BEGIN
#version 330 core

layout(location = 0) out vec4 color;
in vec2 v_tex_coord;
in vec4 v_color;
uniform float radius;

void main()
{ 
    
    float dx = radius - v_tex_coord.x;
    float dy = radius - v_tex_coord.y;
    float dist = sqrt(dx*dx+dy*dy);

    if(dist > radius)
    {
        color = vec4(0,0,0,0);
    }
    else
    {
        color = vec4(1,1,1,1)*v_color;
    }   
    if(abs(v_tex_coord.y - radius) < 1 && v_tex_coord.x > radius)
    {
        color.r = 0;
        color.g = 0;
        color.b = 0;
        color.a = 1;
    }
}
#FRAG_END