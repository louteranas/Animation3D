#version 410

uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;
uniform bool noColor;
uniform vec3 lightPosition;

in vec4 vertex;
in vec4 normal;
in vec4 color;

out vec4 eyeVector;
out vec4 lightVector;
out vec4 vertColor;
out vec4 vertNormal;
out vec4 vertPos;

void main( void )
{
    if (noColor) vertColor = vec4(0.2, 0.6, 0.7, 1.0 );
    else vertColor = color;
    vertNormal.xyz = normalize(normalMatrix * normal.xyz);
    vertNormal.w = 0.0;

    gl_Position = perspective * matrix * vertex;
}
