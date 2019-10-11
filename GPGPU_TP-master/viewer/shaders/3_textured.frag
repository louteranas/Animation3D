#version 410

uniform float lightIntensity;
uniform sampler2D colorTexture;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec2 textCoords;
in vec4 lightSpace;

out vec4 fragColor;

void main( void )
{
   fragColor = vertColor;
}
