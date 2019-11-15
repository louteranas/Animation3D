#version 410
#define M_PI 3.14159265358979323846


uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;

// light position
uniform vec3 lightPosition;

uniform bool transparent;
uniform float shininess;
uniform float eta;

in vec4 position;

out vec4 fragColor;

/* define others spheres */
const int numberOfSpheres = 2;
vec3 centers[numberOfSpheres] = vec3[](
    vec3(center.x, center.y, center.z),
    vec3(center.x+2*radius, center.y+2*radius, center.z+2*radius)
);
float radiuss[numberOfSpheres] = float[](
    radius,
    radius
);

/* compute ambient lighting when the intersection is in the shadow of a light */
vec4 computeAmbientLighting(){
    // TODO
    return vec4(0,1,0,1);
}

/* compute color source lighting */
vec4 computeColorFromLightSource(in vec3 start, in vec3 normal){
    // TODO
    return vec4(1,0,0,1);
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

bool raySphereIntersectOne(in vec3 start, in vec3 direction, in int indexSphere, out float lambda){
    vec3 cp = start - centers[indexSphere];
    // CP dot u
    float cpU = dot(cp,direction);
    // delta
    float delta = computeDelta(indexSphere, cp, cpU);
    // if intersection
    if(delta >= 0){
        // compute lambda: (-b-sqrt(delta))/2a
        lambda = (-2.*cpU-sqrt(delta))/2.;
        return true;
    } else {
        return false;
    }
}

/* compute the color of the pixel of impact between the ray and the color source*/
vec4 getColorFromLightSource(in vec3 start, in vec3 normal, in int indexSphere)
{
    // 1) VERIFY IF INTERSECT A SPHERE
    bool intersect = false;
    vec3 direction = normalize(start - lightPosition);
    float lambda = 0;
    for(int i = 0; i< numberOfSpheres; i++){
        if(indexSphere != i){
            intersect = intersect || raySphereIntersectOne(start,direction,i,lambda);
        }
    }
    
    // 2) COLOR
    if(intersect){
        // ambient lighting
        return computeAmbientLighting();
    } else {
        // color source
        return computeColorFromLightSource(start,normal);
    }
}

/* find ray sphere intersection with start (eye), direction (u) and intersection (to be the result
    return true if intersect, false if not */
bool raySphereIntersect(in vec3 start, in vec3 direction, out int indexSphere, out vec3 intersection, out vec3 normal) {
    bool intersect = false;
    int counterSphere = 0;
    float lambda; 
    float lambdaMin;

    // compute intersection, computing lambda minimum
    for(int i = 0; i< numberOfSpheres; i++){
        // verify intersection with sphere i
        intersect = raySphereIntersectOne(start, direction, i, lambda);
        if(intersect){
            // save lambda or not
            if(counterSphere == 0){
                lambdaMin = lambda;
                indexSphere = i;
                counterSphere += 1;
            } else {
                if(lambda < lambdaMin){
                    lambdaMin = lambda;
                    indexSphere = i;
                    counterSphere += 1;
                }
            }
        }
    }
    // if intersect nobody
    if(counterSphere == 0){
        return false;
    } else {
        // compute intersection = P+lambdaMin*u
        intersection = start + lambdaMin*direction;
        // compute normal
        normal = intersection - centers[indexSphere];
        return true;
    }
}

/* compute result color with n spheres */
vec4 computeResultColor(vec3 u, vec3 eye){
    vec4 resultColor;

    vec3 intersection;
    vec3 reflectedRay;
    vec3 normal;
    int indexSphere;

    // limit for number of bounds
    int limit = 4;

    // 1) INIT STACKS
    vec3 stackOfIntersections[4] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    vec3 stackOfNormals[4] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    vec4 stackOfColors[4] = vec4[](
        vec4(0.0,0.0,0.0,0.0),
        vec4(0.0,0.0,0.0,0.0),
        vec4(0.0,0.0,0.0,0.0),
        vec4(0.0,0.0,0.0,0.0)
    );
    vec3 stackOfIncoming[4] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    int counter = 0;

    // 2) FOLLOW THE RAY 
    // until it reaches the limit or it leaves the scene
    bool intersect = raySphereIntersect(eye, u, indexSphere,intersection,normal);
    while(counter < limit && intersect){
        // update information about the intersection, normal, color and incoming direction
        stackOfIntersections[counter].xyz = intersection.xyz;
        stackOfNormals[counter].xyz = normal.xyz;
        stackOfIncoming[counter].xyz = u.xyz;
        stackOfColors[counter].rgb = getColorFromLightSource(intersection, normal, indexSphere).rgb;
        counter += 1;

        // new reflected ray
        reflectedRay = reflect(u, normal);
        u = normalize(reflectedRay);
        intersect = raySphereIntersect(intersection,u,indexSphere,intersection,normal);
    }

    // NO INTERSECTION
    if(counter == 0){
        resultColor.rgb = vec4(0,0,0,1).rgb;
        return resultColor;
    } else {
        // NO MORE INTERSECTION: FOLLOW THE RAY STARTING FROM THE LAST INTERSECTION
        counter -= 1;
        while(counter >= 0){
            resultColor.rgb +=  stackOfColors[counter].rgb;
            counter -= 1;
        }
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
    fragColor = computeResultColor(u, eye);
}
