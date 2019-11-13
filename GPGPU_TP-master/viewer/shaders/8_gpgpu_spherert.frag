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

/* define others spheres */
vec3 centers[1] = vec3[1](vec3(center.x, center.y, center.z));
float radiuss[1] = float[1](radius);


/* compute the color of the pixel of impact between the ray and the env Map*/
vec4 getColorFromEnvironment(in vec3 direction)
{
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
    theta = acos(dot(normalize(direction), vec3(0, -1, 0))); // we change the direction of the vector to not get a flipped reflection
    thi = atan(direction.x, direction.z); 
    coord2D.x = 0.5 + thi/(2*M_PI); // projection on the interval 0-1
    coord2D.y = theta/M_PI; // projection on the interval 0-1
    return texture2D(envMap,coord2D);
}

/* compute reflected and refracted rays of u */
/* n2: outgoing milieu */
void computeReflectedRefractedRays(in int index, in vec3 intersection, in vec3 u, in float etaN, out vec3 reflectedRay, out vec3 refractedRay){
    // normal = intersection-center */
    vec3 normal = normalize(intersection - centers[index]);
    // reflected
    reflectedRay = reflect(u, normal);
    // refracted using eta
    refractedRay = refract(u, normal, etaN);
}

/* Compute delta = 4(CP dot u)² - 4(CP - r²) */
float computeDelta(in int index, in vec3 cp, in float cpU){
    // CP
    float cpModule = length(cp);
    // CP²
    float cpModuleSquare = cpModule * cpModule;
    // r²
    float radiusSquare = radiuss[index]*radiuss[index];
    // (CP dot u)²
    float cpUSquare = cpU * cpU;
    return 4.*cpUSquare-4.*(cpModuleSquare-radiusSquare);
}

/* find ray sphere intersection with start (eye), direction (u) and intersection (to be the result
    return true if intersect, false if not */
bool raySphereIntersect(in int index, in vec3 start, in vec3 direction, out vec3 intersection) {
    // CP
    vec3 cp = start - centers[index];
    // CP dot u
    float cpU = dot(cp,direction);
    // delta
    float delta = computeDelta(index, cp, cpU);
    // if no intersection
    if(delta < 0){
        return false;
    }
    // compute lambda: (-b-sqrt(delta))/2a
    float lambda = (-2.*cpU-sqrt(delta))/2.;
    // compute intersection = P+lambda*u
    intersection = start + lambda*direction;
    return true;
}

/**
* this fuction calculate the Fresnet Coefficient
* see https://en.wikipedia.org/wiki/Fresnel_equations
* n1: incoming milieu
* n2: outgoing milieu
* For more details about the variable names, check the TP page
**/
float fresnelCoeff(float cosThethaD, float etaN){
     // Ci coeff 
     float Ci = sqrt(1.- etaN*etaN*(1. - cosThethaD*cosThethaD));
     // Fs coef
     float fracFs = abs(etaN*cosThethaD - Ci) / abs(etaN*cosThethaD + Ci);
     float Fs = fracFs*fracFs;
     // Fp coeff
     float fracFp = abs(etaN*etaN*Ci - cosThethaD) / abs(etaN*etaN*Ci + cosThethaD);
     float Fp = fracFp*fracFp;
     ///// Fresnel coeff
     float F = (Fs + Fp)/2.;
     if(F > 1.){
         return 1.;
     }
     return F;

}

/* compute cos theta d with u and n */
float getCosThetha(in int index, vec3 intersection, vec3 u, bool inside){
    vec3 normal = normalize(intersection - centers[index]);
    if(inside){
        normal = normalize(centers[index] - intersection);
    }
    return dot(normal, normalize(-1.*u));
}

/* compute result color with n spheres */
vec4 computeResultColor(vec3 u, vec3 eye, int n){
    vec4 resultColor;
    bool rayIntersected = false;

    for(int j = 0; j<n; j++){
        // Step 3: ray intersection
        vec3 intersection;
        bool intersect = raySphereIntersect(j, eye, u, intersection);

        if(intersect){
            rayIntersected = true;
            // Step 4: compute reflected and refracted rays
            vec3 reflectedRay;
            vec3 refractedRay;
            vec4 result;
            int numberOfRebounds = 1;
            computeReflectedRefractedRays(j, intersection, u, 1./eta, reflectedRay,refractedRay);
            float cosThetha = getCosThetha(j, intersection, u, false);
            float fresnelReflexion = fresnelCoeff(cosThetha, 1./eta);
            float fresnelTrans = 1 - fresnelReflexion;
            float lastCoeff = fresnelTrans;
            // We have just calculated the first reflected ray, after this the ray that we need is the refracted one since we are inside the sphere,so we need the ray that get out of the sphere.
            if(numberOfRebounds  != 0 || transparent){
                result.rgb = fresnelReflexion * getColorFromEnvironment(normalize(reflectedRay)).rgb;
            }
            else{
                result.rgb = getColorFromEnvironment(normalize(reflectedRay)).rgb;
            }
            // now the incoming ray if the last refracted from the first raying coming from the eye
            u = normalize(refractedRay);
            if(numberOfRebounds >0 && transparent){
                for(int i = 0; i<numberOfRebounds; i++){
                    // the first intersection param is the last intersection and it represents our start point of the ray the second param is u, the direction of the ray with is the normalised last reflected/refracted ray 
                    //depending on the iteration the last param is again interection which give us the new intersection point
                    intersect = raySphereIntersect(j, intersection, u, intersection);
                    // we check if there is an intersection but in out case it's useless since we are in a sphere so we will always have an intersection computing reflected and refracted ray
                    computeReflectedRefractedRays(j, intersection, u, 1./eta, reflectedRay,refractedRay);
                    // computing the angle between the direction and normal in intersection point 
                    cosThetha = getCosThetha(j, intersection, u, true);
                    // we multiply with the last coeff from the last calculated ray
                    fresnelReflexion = lastCoeff * fresnelCoeff(cosThetha, eta);
                    fresnelTrans = lastCoeff - fresnelReflexion;
                    // we update the last coeff which is the the reflexion one because we follow the ray that stays inside the shpere
                    lastCoeff = fresnelTrans;
                    // the new ray to trace is the one staying inside the sphere aka the reflected one
                    u = normalize(reflectedRay);
                    // and we add to the result the color of the pixel in which the ray that got out of the sphere aka the refracted one intersects the envMap
                    result.rgb = result.rgb + fresnelReflexion * getColorFromEnvironment(refractedRay).rgb;   
                }
            }
            resultColor.rgb = result.rgb;
            
        }
    }
    if(!rayIntersected){
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
    // n = number of spheres
    int n = 1;
    int NumRebounds = 0;
    fragColor = computeResultColor(u, eye, n);
}
