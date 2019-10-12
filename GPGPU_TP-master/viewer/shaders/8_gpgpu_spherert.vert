#version 410

// Standard parameters of a VAO
in vec4 vertex;
in vec4 normal;
in vec4 color;
in vec2 texcoords;

out vec4 vertColor;
out vec4 vertNormal;
out vec2 textCoords;
out vec4 position;

void main( void )
{
    // You need to use color/normal/textcoords. Bad things can happen otherwise
    vertColor = color;
    vertNormal = normal;
    textCoords = texcoords;

    position = vertex;
    gl_Position = vertex;
}
