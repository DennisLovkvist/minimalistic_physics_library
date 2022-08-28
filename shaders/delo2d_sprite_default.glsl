#VERT_BEGIN
#version 330 core

layout(location = 0) in vec4 position;
layout(location = 1) in vec2 tex_coord;
layout(location = 2) in float tex_index;
layout(location = 3) in vec4 color;

out vec2 v_tex_coord;
out float v_tex_index;
out vec4 v_color;

uniform mat4 u_mvp;

void main()
{     
    gl_Position = position * u_mvp;
    v_tex_coord = tex_coord;
    v_tex_index = tex_index;
    v_color = color;
}

#VERT_END
#FRAG_BEGIN
#version 330 core

uniform sampler2D u_textures[3];
layout(location = 0) out vec4 color;
in vec2 v_tex_coord;
in float v_tex_index;
in vec4 v_color;

void main()
{ 
    int index = int(v_tex_index);
    //vec4 tex_color = texture2D(u_textures[index],v_tex_coord); 
    vec4 tex_color = texelFetch(u_textures[index], ivec2(v_tex_coord), 0);

    


    color = tex_color*v_color;
}
#FRAG_END