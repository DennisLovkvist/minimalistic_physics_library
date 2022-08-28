#VERT_BEGIN
#version 330 core

layout(location = 0) in vec4 position;
layout(location = 1) in vec4 color;

out vec4 v_color;

uniform mat4 u_mvp;

void main()
{     
    gl_Position = position * u_mvp;
    v_color = color;
}

#VERT_END
#FRAG_BEGIN
#version 330 core

layout(location = 0) out vec4 color;
in vec4 v_color;

void main()
{ 
    color = v_color;
}
#FRAG_END