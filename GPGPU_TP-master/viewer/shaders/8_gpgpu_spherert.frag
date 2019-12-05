#version 410
#define M_PI 3.14159265358979323846
#define EPS 0.000001

/****************************************************************************************************************************************************/
/***************************************************************** PARAMETERS ***********************************************************************/
/****************************************************************************************************************************************************/

/* 
* Parameters UNIFORM
*/
uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;
uniform int numBounds;
uniform bool transparent;
uniform float shininess;
uniform float eta;

/* 
* Parameters IN: position
*/
in vec4 position;

/* 
* Parameters OUT: result color
*/
out vec4 fragColor;


/****************************************************************************************************************************************************/
/***************************************************************** ENV MAP **************************************************************************/
/****************************************************************************************************************************************************/

/* 
* compute color from the env map
*/
vec4 getColorFromEnvironnement(in vec3 direction){
    // the env map is like an infinite sphere arround the world
    // so the first thing to do is to transform the direction coords into 
    //sheprique coords, then we know that thi is in -pi, pi 
    //and Theta between 0 ans 2*Pi, so we transform this to get a value between 0 and 1 
    // in order to be able to use the function Texture 2D that return the frag color coresponding
    //point of the ray's impact with the env Map
    vec2 coord2D; // our interpolated coords between 0 and 1  
    float thi; // the angle between the X axis and the direction vector's pro111jection on the plane XY
    float theta; // the angle between the Z axis and the direction vector 
    // 
    float scalar = dot(normalize(direction), vec3(0, -1, 0));
    theta = acos(scalar); // we change the direction of the vector to not get a flipped reflection
    thi = atan(direction.x, direction.z); 
    coord2D.x = 0.5 + thi/(2*M_PI); // projection on the interval 0-1
    coord2D.y = theta/M_PI; // projection on the interval 0-1

    return texture2D(envMap,coord2D);
}

/****************************************************************************************************************************************************/
/***************************************************************** RAY TRACING **********************************************************************/
/****************************************************************************************************************************************************/

/* 
* Compute delta = 4(CP dot u)² - 4(CP - r²)
*/
float computeDelta(in vec3 cp, in float cpU){
    // CP
    float cpModule = length(cp);
    // CP²
    float cpModuleSquare = cpModule * cpModule;
    // r²
    float radiusSquare = radius*radius;
    // (CP dot u)²
    float cpUSquare = cpU * cpU;
    return 4.*cpUSquare-4.*(cpModuleSquare-radiusSquare);
}

/* 
* find ray sphere intersection with start (eye), direction (u) and intersection (to be the result
*    return true if intersect, false if not
*/
bool raySphereIntersect(in vec3 start, in vec3 direction,in bool inside, out vec3 intersection) {
    // CP
    vec3 cp = start - center;
    // CP dot u
    float cpU = dot(cp,direction);
    // delta
    float delta = computeDelta(cp, cpU);
    // if no intersection
    if(delta < 0){
        return false;
    }
    // compute lambda: (-b-sqrt(delta))/2a
    
    float lambda = (-2.*cpU-sqrt(delta))/2.;
    if(inside){
        lambda = (-2.*cpU+sqrt(delta))/2.;
    }
    // compute intersection = P+lambda*u
    intersection = start + lambda*direction;
    return true;
}

/* 
* compute cos theta d with u and n
*/
float getCosThetha(vec3 intersection, vec3 u, bool inside){
    // normal
    vec3 normal = normalize(intersection - center);
    // if inside the sphere
    if(inside){
        normal = normalize(center - intersection);
    }
    // cos theta = abs((n.u))
    return abs(dot(normal, normalize(u)));
}

/* 
* compute fresnel coefficients
*/
float fresnelCoeff(float cosThethaD, float etaN){
     // Ci coeff 
     float Ci = pow(max(EPS, (pow(etaN, 2) - (1 - pow(cosThethaD, 2)))), 0.5);
     // Fs coef
     float fracFs = (cosThethaD - Ci) / (cosThethaD + Ci);
     float Fs = pow(abs(fracFs), 2);
     // Fp coeff
     float fracFp = (pow(etaN,2)*cosThethaD - Ci) / (pow(etaN,2)*cosThethaD + Ci);
     float Fp = pow(abs(fracFp), 2);
     ///// Fresnel coeff
     float F = (Fs + Fp)/2;
     return F;
}

/* 
* compute reflected and refracted rays, using the glsl functions reflect and refract
*/
void computeReflectedRefractedRays(in vec3 intersection, in vec3 u, in float etaN, in bool inside, out vec3 reflectedRay, out vec3 refractedRay){
    // normal = intersection-center */
    vec3 normal = normalize(intersection - center);
    if(inside){
        normal = normalize(center - intersection);
    }
    // reflected
    reflectedRay = reflect(u, normal);
    // refracted using eta
    refractedRay = refract(u, normal, etaN);
    
}

/* 
* compute result color with n bounds
*/
vec4 computeResultColor(vec3 u, vec3 eye, int n){
    // parameters
    vec3 intersection;
    bool hasIntersect = raySphereIntersect(eye, u, false, intersection);
    vec4 resultColor;
    vec3 refractedRay;
    vec3 reflectedRay;
    
    // if intersection for the first ray
    if(hasIntersect) {
        vec4 result;
        // compute reflected and refracted rays
        computeReflectedRefractedRays(intersection, u, 1./eta, false, reflectedRay,refractedRay);
        // compute costheta
        float cosThetha = getCosThetha(intersection, u, false);
        // compute fresnel coeffs
        float fresnelReflexion = fresnelCoeff(cosThetha, eta);
        float fresnelTrans = 1 - fresnelReflexion;
        // keep last transmission coeff
        float lastCoeff = fresnelTrans;
        // if bounds and transparency:
        if(n != 0 && transparent && !(fresnelReflexion>1.)){
            // result = coeff fresnel * env map, and continue
            result = fresnelReflexion * getColorFromEnvironnement(normalize(reflectedRay));
        }
        else{
            // result = env map and stop
            return getColorFromEnvironnement(normalize(reflectedRay));
        }
        // new u
        u = normalize(refractedRay);
        // if bounds and transparency
        if(n >0 && transparent){
            // for each bounds
            for(int i = 0; i<n; i++){
                // intersection
                hasIntersect = raySphereIntersect(intersection, u, true, intersection);
                // reflected/refracted rays
                computeReflectedRefractedRays(intersection, u, eta, true, reflectedRay,refractedRay);
                // costheta
                cosThetha = getCosThetha(intersection, u, true);
                // fresnel coeff
                fresnelReflexion = lastCoeff * fresnelCoeff(cosThetha, 1./eta);
                fresnelTrans = lastCoeff - fresnelReflexion;
                lastCoeff = fresnelReflexion;
                // new u
                u = normalize(reflectedRay);
                // if last coeff is not equal to 1
                if(!(lastCoeff>1.)){
                    // result += coeff * envMap
                    result = result + fresnelTrans * getColorFromEnvironnement(refractedRay);
                }
            }
        }
        resultColor = result;
    } else {
        // if no intersection: env map
        resultColor = getColorFromEnvironnement(u);
    }
    return resultColor;
}


/****************************************************************************************************************************************************/
/********************************************************************** MAIN ************************************************************************/
/****************************************************************************************************************************************************/

void main(void)
{
    // Step 1: I need pixel coordinates. Division by w?
    vec4 worldPos = position;
    worldPos.z = 1; // near clipping plane
    worldPos = persp_inverse * worldPos;
    worldPos /= worldPos.w;
    worldPos.w = 0;
    worldPos = normalize(worldPos);

    // Step 2: ray direction:
    vec3 u = normalize((mat_inverse * worldPos).xyz);
    vec3 eye = (mat_inverse * vec4(0, 0, 0, 1)).xyz;

    // Step 3: compute frag color
    // n = number of rebounds
    int n = numBounds;
    fragColor = computeResultColor(u, eye, n);
}