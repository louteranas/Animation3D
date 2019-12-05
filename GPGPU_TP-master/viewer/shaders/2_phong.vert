#version 410

uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;
uniform bool noColor;
uniform vec3 lightPosition;

// Shadow Mapping
uniform mat4 matrixLight;
uniform mat4 perspectiveLight;

// World coordinates
in vec4 vertex;
in vec4 normal;
in vec4 color;

// Camera-space coordinates
out vec4 eyeVector;
out vec4 lightVector;
out vec4 lightSpace; // placeholder for shadow mapping
out vec4 positionScreen;
out vec4 vertColor;
out vec4 vertNormal;

void main( void )
{
    if (noColor) vertColor = vec4(0.2, 0.6, 0.7, 1.0 );
    else vertColor = color;
    vertNormal.xyz = normalize(normalMatrix * normal.xyz);
    vertNormal.w = 0.0;
    
    // lightPosCS is light position in camera space coords
    vec4 lightPosCS = matrix * vec4(lightPosition,1);
    vec4 vertexPositionCS = matrix * vertex;

    // light vector and eye vector in camera space coordinates
    lightVector = normalize(lightPosCS - vertexPositionCS);
    vec4 eyePosition = vec4(0.0, 0.0, 0.0, 1.0);
    eyeVector = normalize(eyePosition - vertexPositionCS);

    // shadow mapping: position in screen space and in light space
    lightSpace = perspectiveLight * matrixLight * vertex;
    positionScreen = perspective * matrix * vertex;


    // VERY IMPORTANT : only GL position has to be converted to
    // screen space, so all other calculus do not GET the perspective 
    // matrix

    // gl_position is the vertex position in screen space
    gl_Position = perspective * matrix * vertex;
}
