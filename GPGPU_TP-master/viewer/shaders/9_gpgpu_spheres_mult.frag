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
int numberOfSpheres = 1;
vec3 centers[1] = vec3[1](vec3(center.x, center.y, center.z));
float radiuss[1] = float[1](radius);

/* compute the color of the pixel of impact between the ray and the color source*/
vec4 getColorFromLightSource(in vec3 start, in int indexSphere)
{
    // VERIFY IF INTERSECT A SPHERE
    bool intersect = false;
    vec3 direction = normalize(start - sourceCenter);
    for(int i = 0; i< numberOfSpheres; i++){
        // CP
        vec3 cp = start - centers[i];
        // CP dot u
        float cpU = dot(cp,direction);
        // delta
        float delta = computeDelta(i, cp, cpU);
        if(delta >= 0 && indexSphere != i){
            intersect = true;
        }
    
    if(intersect){
        // AMBIENT LIGHTING
        // TODO
        return vec4(0,0,0,1);
    }
    // COLOR FROM THE SOURCE
    return vec4(0,0,0,1);
}

/* compute reflected and refracted rays of u */
/* n2: outgoing milieu */
void computeReflected(in int index, in vec3 intersection, in vec3 u,out vec3 reflectedRay){
    // normal = intersection-center */
    vec3 normal = normalize(intersection - centers[index]);
    // reflected
    reflectedRay = reflect(u, normal);
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
bool raySphereIntersect(in vec3 start, in vec3 direction, out int indexSphere, out vec3 intersection, out vec3 normal) {
    bool intersect = false;
    int counterSphere = 0;
    float lambda; 
    float lambdaMin;

    // compute intersection, computing lambda minimum
    for(int i = 0; i< numberOfSpheres; i++){
        // CP
        vec3 cp = start - centers[i];
        // CP dot u
        float cpU = dot(cp,direction);
        // delta
        float delta = computeDelta(i, cp, cpU);
        // if intersection
        if(delta >= 0){
            intersect = true;
            // compute lambda: (-b-sqrt(delta))/2a
            lambda = (-2.*cpU-sqrt(delta))/2.;
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
    // compute intersection = P+lambdaMin*u
    intersection = start + lambdaMin*direction;
    // compute normal
    normal = intersection - centers[indexSphere];
    return true;
}

/* compute result color with n spheres */
vec4 computeResultColor(vec3 u, vec3 eye){
    vec4 resultColor;
    bool rayIntersected = false;
    vec3 intersection;
    vec3 reflectedRay;
    vec3 normal;
    int indexSphere;

    // limit for number of bounds
    int limit = 2;

    // STACKS
    vec3 stackOfIntersections[limit];
    vec3 stackOfNormals[limit];
    vec4 stackOfColors[limit];
    vec3 stackOfUncoming[limit];
    int counter = 0;

    // FOLLOW THE RAY UNTIL IT REACHES THE LIMIT
    bool intersect = raySphereIntersect(eye, u, indexSphere,intersection,normal);
    while(counter < limit && intersect == true){
        // update information about the intersection, normal, color and incoming direction
        stackOfIntersections[counter] = intersection;
        stackOfNormals[counter] = normal;
        stackOfIncoming[counter] = u;
        // stackOfColors[counter] = getColorFromLightSource();
        counter += 1;

        // new reflected ray
        computeReflected(indexSphere,intersection, u,reflectedRay);
        u = normalize(reflectedRay);
        intersect = raySphereIntersect(intersection, u,indexSphere,intersection,normal);
    }

    // NO INTERSECTION
    if(counter == 0){
        resultColor.rgb = vec4(0,0,0,1).rgb;
        return resultColor;
    } else {
        // NO MORE INTERSECTION: FOLLOW THE RAY STARTING FROM THE LAST INTERSECTION
        counter -= 1;
        while(counter >= 0){
            resultColor.rgb +=  vec4(0,0,0,1).rgb; //stackOfColors[counter].rgb;
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
