#version 410
#define M_PI 3.14159265358979323846

uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;
uniform int numBounds;
uniform bool transparent;
uniform float shininess;
uniform float eta;

in vec4 position;

out vec4 fragColor;


#define EPS                 0.000001

/* compute the color of the pixel of impact between the ray and the env Map*/
vec4 getColorFromEnvironment(in vec3 direction)
{

    vec2 coord2D; 
    float thi; 
    float theta;  
    float scalar = dot(normalize(direction), vec3(0, -1, 0));
    theta = acos(scalar); 
    thi = atan(direction.x, direction.z); 
    coord2D.x = 0.5 + thi/(2*M_PI); 
    coord2D.y = theta/M_PI; 

    return texture2D(envMap,coord2D);
}


/* Compute delta = 4(CP dot u)² - 4(CP - r²) */
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

/* find ray sphere intersection with start (eye), direction (u) and intersection (to be the result
    return true if intersect, false if not */
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

/* compute cos theta d with u and n */
float getCosThetha(vec3 intersection, vec3 u, bool inside){
    vec3 normal = normalize(intersection - center);
    if(inside){
        normal = normalize(center - intersection);
    }
    return abs(dot(normal, normalize(u)));
}

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


/* compute result color with n spheres */
vec4 computeResultColor(vec3 u, vec3 eye, int n){

    // TODO
    vec3 intersection;
    bool hasIntersect = raySphereIntersect(eye, u, false, intersection);
    vec4 resultColor;
    vec3 refractedRay;
    vec3 reflectedRay;
    
    if(hasIntersect) {
        vec4 result;
        // Step 4: compute reflected and refracted rays
        computeReflectedRefractedRays(intersection, u, 1./eta, false, reflectedRay,refractedRay);
        float cosThetha = getCosThetha(intersection, u, false);
        float fresnelReflexion = fresnelCoeff(cosThetha, eta);
        float fresnelTrans = 1 - fresnelReflexion;
        float lastCoeff = fresnelTrans;
        if(n != 0 && transparent && !(fresnelReflexion>1.)){
            result = fresnelReflexion * getColorFromEnvironment(normalize(reflectedRay));
        }
        else{
            return getColorFromEnvironment(normalize(reflectedRay));
        }
        u = normalize(refractedRay);
        if(n >0 && transparent){
            for(int i = 0; i<n; i++){
                hasIntersect = raySphereIntersect(intersection, u, true, intersection);
                computeReflectedRefractedRays(intersection, u, eta, true, reflectedRay,refractedRay);
                cosThetha = getCosThetha(intersection, u, true);
                fresnelReflexion = lastCoeff * fresnelCoeff(cosThetha, 1./eta);
                fresnelTrans = lastCoeff - fresnelReflexion;
                lastCoeff = fresnelReflexion;
                u = normalize(reflectedRay);
                if(!(lastCoeff>1.)){
                result = result + fresnelTrans * getColorFromEnvironment(refractedRay);
                }
            }
        }
        resultColor = result;
    } else {
        resultColor = getColorFromEnvironment(u);
    }
    return resultColor;
}

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