#version 410
#define M_PI 3.14159265358979323846

uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;

uniform bool transparent;
uniform float shininess;
uniform float eta;

in vec4 position;

out vec4 fragColor;


vec4 getColorFromEnvironment(in vec3 direction)
{
    vec2 coord2D;
    float thi;
    float theta;
    // 
    theta = M_PI + acos(dot(normalize(direction), vec3(0, 1, 0)));
    thi = atan(direction.x, direction.z);
    coord2D.x = 0.5 + thi/(2*M_PI);
    coord2D.y = theta/M_PI;
    return texture2D(envMap,coord2D);
}

/* compute reflected and refracted rays of u */
void computeReflectedRefractedRays(in vec3 intersection, in vec3 u, out vec3 reflectedRay, out vec3 refractedRay){
    // normal = intersection-center */
    vec3 normal = normalize(intersection - center);
    // reflected
    reflectedRay = reflect(u, normal);
    // refracted using eta
    refractedRay = refract(u, normal, eta);
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
    return 4*cpUSquare-4*(cpModuleSquare-radiusSquare);
}

/* find ray sphere intersection with start (eye), direction (u) and intersection (to be the result
    return true if intersect, false if not */
bool raySphereIntersect(in vec3 start, in vec3 direction, out vec3 intersection) {
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
    float lambda = -cpU-sqrt(delta);
    // compute intersection = P+lambda*u
    intersection = start + lambda*direction;
    return true;
}

/**
* this fuction calculate the Fresnet Coefficient
* For more details about the variable names, check the TP page
**/
float fresnelCoeff(float cosThethaD){
     // Ci coeff 
     float Ci = pow((pow(eta, 2) - (1 - pow(cosThethaD, 2))), 0.5);
     // Fs coef
     float fracFs = (cosThethaD - Ci) / (cosThethaD + Ci);
     float Fs = pow(abs(fracFs), 2);
     // Fp coeff
     float fracFp = (pow(eta,2)*cosThethaD - Ci) / (pow(eta,2)*cosThethaD + Ci);
     float Fp = pow(abs(fracFp), 2);
     ///// Fresnel coeff
     float F = (Fs + Fp)/2;
     if(F > 1.){
         return 1.0;
     }
     return F;
}

float getCosThetha(vec3 intersect, vec3 u){
    vec3 normal = normalize(intersect - center);
    return dot(normal, normalize(u));

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

    // Step 3: ray intersection
    vec3 intersection;
    bool intersect = raySphereIntersect(eye, u, intersection);

    vec4 resultColor;
    float fresnelRefract[4];
    float fresnelReflex[4];
    if(intersect){
        // Step 4: compute reflected and refracted rays
        vec3 reflectedRay;
        vec3 refractedRay;
        vec4 result;
        computeReflectedRefractedRays(intersection, u,reflectedRay,refractedRay);
        
        float cosThetha = getCosThetha(intersection, u);
        float fresnelReflexion = fresnelCoeff(cosThetha);
        float fresnelTrans = 1 - fresnelReflexion;
        result = fresnelReflexion * getColorFromEnvironment(reflectedRay);
        // we neef to change the reflected and refracted rays for inside the sphere;
        vec3 temp = reflectedRay;
        reflectedRay = refractedRay; // rayon qui reste dans la shpere
        refractedRay = temp; // rayon qio sort
        for(int i = 0; i<5; i++){
            intersect = raySphereIntersect(intersection, reflectedRay, intersection);
            if(intersect){
                computeReflectedRefractedRays(intersection, reflectedRay, reflectedRay,refractedRay);
                cosThetha = getCosThetha(intersection, reflectedRay);
                fresnelReflexion = fresnelCoeff(cosThetha);
                fresnelTrans = 1 - fresnelReflexion;
                result += fresnelTrans * (fresnelReflexion * getColorFromEnvironment(reflectedRay);
            }    
        }
    } else {
        resultColor = getColorFromEnvironment(u);
    }
    

    fragColor = resultColor;
}
