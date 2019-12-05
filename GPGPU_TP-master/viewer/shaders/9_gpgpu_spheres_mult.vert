#version 410

/****************************************************************************************************************************************************/
/***************************************************************** PARAMETERS ***********************************************************************/
/****************************************************************************************************************************************************/

/* 
* Parameters UNIFORM of the VAO
*/
uniform bool noColor;

/* 
* Parameters IN of the VAO
*/
in vec4 vertex;
in vec4 normal;
in vec4 color;
in vec2 texcoords;

/* 
* Parameters OUT of the VAO
*/
out vec2 textCoords;
out vec4 position;

// Camera-space coordinates
out vec4 eyeVector;
out vec4 lightVector;
out vec4 lightSpace; // placeholder for shadow mapping
out vec4 vertColor;
out vec4 vertNormal;


/****************************************************************************************************************************************************/
/********************************************************************** MAIN ************************************************************************/
/****************************************************************************************************************************************************/

void main( void )
{
    // You need to use color/normal/textcoords. Bad things can happen otherwise
    vec4 temp = /*vec4(1,1,1,1); /*/vec4(0.5, 0.75, 0.81, 1.0 );
    if (noColor) vertColor = temp;
    else vertColor = temp;

    vertNormal = normal;
    textCoords = texcoords;

    // position
    position = vertex;
    gl_Position = vertex;
}
