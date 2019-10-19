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
    // TODO
    return vec4(1);
}

/* compute reflected and refracted rays of u */
void computeReflectedRefractedRays(in vec3 intersection, in vec3 eye,in vec3 u, out vec3 reflectedRay, out vec3 refractedRay){
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
    float cpModule = abs(cp);
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

    if intersect{
        // Step 4: compute reflected and refracted rays
        vec3 reflectedRay;
        vec3 refractedRay;
        computeReflectedRefractedRays(intersection, eye,u,reflectedRay,refractedRay);
    }

    vec4 resultColor = vec4(0,0,0,1);
    fragColor = resultColor;
}
