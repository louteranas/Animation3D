#version 410

in vec4 vertNormal;
in vec4 vertColor;

out vec4 fragColor;

void main( void )
{
    fragColor = vertColor;
}
