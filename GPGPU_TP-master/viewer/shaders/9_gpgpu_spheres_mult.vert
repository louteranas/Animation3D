#version 410

// Standard parameters of a VAO
in vec4 vertex;
in vec4 normal;
in vec4 color;
in vec2 texcoords;

out vec2 textCoords;
out vec4 position;

uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;
uniform bool noColor;
uniform vec3 lightPosition;

// Camera-space coordinates
out vec4 eyeVector;
out vec4 lightVector;
out vec4 lightSpace; // placeholder for shadow mapping
out vec4 vertColor;
out vec4 vertNormal;


void main( void )
{
    // You need to use color/normal/textcoords. Bad things can happen otherwise
    vertColor = color;
    vertNormal = normal;
    textCoords = texcoords;

    position = vertex;
    gl_Position = vertex;
}
