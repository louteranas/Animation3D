#version 410
#define M_PI 3.14159265358979323846

// uniform color models
uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform bool noColor;

// matrix perspective
uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;

// env map
uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;

// light position
uniform vec3 lightPosition;

uniform bool transparent;
uniform float eta;

in vec4 position;
in vec4 vertColor;

out vec4 fragColor;

/* define others spheres */
const int numberOfSpheres = 5;
vec3 centers[numberOfSpheres] = vec3[](
    vec3(center.x, center.y, center.z),
    vec3(center.x+2*radius, center.y+2*radius, center.z+2*radius),
    vec3(center.x-2*radius, center.y-2*radius, center.z-2*radius),
    vec3(center.x-2*radius, center.y+2*radius, center.z+2*radius),
    vec3(center.x-2*radius, center.y+2*radius, center.z-2*radius)
);
float radiuss[numberOfSpheres] = float[](
    radius,
    radius,
    radius,
    radius,
    radius
);

/* compute color from plane */
vec4 getColorFromPlane(in vec3 direction){
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

/* compute ambient lighting when the intersection is in the shadow of a light */
vec4 computeAmbientLighting(){
    // ambient reflection param 
     float Ka = 0.7;
     // setting the ambiantLighting - Ca
     vec4 ambientLight = Ka * lightIntensity * vertColor;
    return ambientLight;
}

/* compute diffuse lighting */
vec4 computeDiffuseLighting(in vec4 vertNormal, in vec4 lightVector){
    vec4 diffuseLighting;
    // Diffuse reflection param 
    float Kd = 0.7;
    // setting the Diffuse lighting - Cd
    diffuseLighting = Kd * vertColor * lightIntensity * max(0, dot(vertNormal, lightVector));
    return diffuseLighting;
}


/**
* this fuction calculate the microfacet normal distribution D(thetaH)
* For more details about the variable names, check the TP page
**/
bool indicatrice(float cosThetaH){
     if(cosThetaH < 0){
          return false;
     }
     return true;
     // maybe it's okay because theta is positive 
     // if not, so the graph and lights are reversed and we have theta positive
}

float NormalDistrib(float cosThetaH, float alpha){
     if(!indicatrice(cosThetaH)){
          return 0;
     }
     float frac1 = 1 / (pow(cosThetaH, 4) * M_PI);
     float tanThetaSquare = (1 - pow(cosThetaH, 2))/pow(cosThetaH, 2);
     float frac2 = pow(alpha/100, 2)/ pow((pow(alpha/100, 2) + pow(tanThetaSquare, 2)), 2);
     return frac1 * frac2;
}

/**
* this fuction calculate the GGX distribution G1(thetaI/thetaO)
* For more details about the variable names, check the TP page
**/
float GGXDistrib(float cosTheta, float alpha){
     float tanThetaSquare = (1 - pow(cosTheta, 2))/pow(cosTheta, 2);
     float base = 1 + sqrt(1 + tanThetaSquare * pow(alpha/100, 2));
     return 2 / base;
}

/**
* this fuction calculate the specular lighting using the blinn-phing model
* For more details about the variable names, check the TP page
**/
vec4 specularLightingBP(float cosThethaD, vec4 halfVector, vec4 normal){
    return fresnelCoeff(cosThethaD,1./eta) * vertColor * pow(max(0, dot(normal, halfVector)), shininess) * lightIntensity;
}

/**
* this fuction calculate the specular lighting using the cook-torrance model
* For more details about the variable names, check the TP page
**/
vec4 specularLightingCT(float cosThethaD, vec4 halfVector, float alpha, vec4 normal, vec4 lightVector, vec4 eyeVector){
     float cosThetaH = dot(normal,halfVector);
     float cosThetaI = dot(normal,lightVector);
     float cosThetaO = dot(normal,eyeVector);
     float top = fresnelCoeff(cosThethaD,1./eta) * NormalDistrib(cosThetaH, alpha) * GGXDistrib(cosThetaI, alpha) * GGXDistrib(cosThetaO, alpha);
     float bottom = 4 * cosThetaI * cosThetaO;
     return (top/bottom)*vertColor*lightIntensity;
}

/* compute specular lighting */
vec4 computeSpecularLighting(in vec4 normal, in vec4 lightVector, in vec4 eyeVector){
    vec4 specularLighting;

    // shininess and alpha reversed
    float alpha = (2.-shininess/100);
    
    // ThetaD, angle between normal & lightVector we only need its cos value
    vec4 halfVector = normalize(lightVector + eyeVector);
    float cosThethaD = dot(halfVector, lightVector); // /(length(halfVector)*length(lightVector)) ? what is the correct formula? OK because normalize :);
    if(blinnPhong)
        // using the blinn phong model
        specularLighting = specularLightingBP(cosThethaD, halfVector, normal);
    else {
        // using the cook torrance model
        specularLighting = specularLightingCT(cosThethaD, halfVector, alpha, normal, lightVector, eyeVector);
    }
    return specularLighting; 
}

/* compute color source lighting */
vec4 computeColorFromLightSource(in bool intersect, in vec3 start, in vec3 u, in vec3 normal){
    /* parameters
        - vertNormal = normal
        - lightVector = position light - start
        - eyeVector = -1*u
    */
    vec4 vertNormal = normalize(vec4(normal.xyz,1.));
    vec4 lightVector = normalize(vec4(lightPosition,1)-vec4(start,1));
    vec4 eyeVector = normalize(vec4(-1*u,1));

    // ambient lighting
    vec4 ambientLighting = computeAmbientLighting();

    if(intersect){
        return ambientLighting;
    } else {
        // diffuse lighting
        vec4 diffuseLighting = computeDiffuseLighting(vertNormal, lightVector);
        // specular lighting
        vec4 specularLighting = computeSpecularLighting(vertNormal,lightVector, eyeVector);
        return ambientLighting + diffuseLighting + specularLighting;
    }
}

/* compute cos theta d with u and n */
float getCosThetha(vec3 normal, vec3 u){
    return dot(normal, normalize(1.*u));
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
vec4 getColorFromLightSource(in vec3 u, in vec3 start, in vec3 normal, in int indexSphere)
{
    // 1) VERIFY IF INTERSECT A SPHERE
    bool intersect = false;
    vec3 direction = normalize(lightPosition.xyz-start);
    float lambda = 0;
    for(int i = 0; i< numberOfSpheres; i++){
        if(indexSphere != i){
            intersect = intersect || raySphereIntersectOne(start,direction,i,lambda);
        }
    }
    
    // 2) COLOR
    return computeColorFromLightSource(intersect, start, u, normal);
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

/* compute if u intersect with the plane defined by point Q and normal n */
bool rayPlaneIntersect(in vec3 start, in vec3 u, out vec3 intersection){
    // define the plane
    vec3 Q = vec3(0,-100,0);
    vec3 n = vec3(0,1,0);

    // compute intersection
    vec3 PQ = Q-start;
    float lambda = -1*(dot(PQ,n)/dot(u,n));
    if(lambda > 0){
        intersection = start+lambda*u;
        return true;
    } else {
        return false;
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
    const int limit = 10;

    // 1) INIT STACKS
    vec3 stackOfIntersections[limit] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    vec3 stackOfNormals[limit] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    vec4 stackOfColors[limit] = vec4[](
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0),
        vec4(0.0,0.0,0.0,1.0)
    );
    vec3 stackOfIncoming[limit] = vec3[](
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0),
        vec3(0.0,0.0,0.0)
    );
    float stackOfFresnel[limit] = float[](0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

    int counter = 0;

    // 2) FOLLOW THE RAY 
    // until it reaches the limit or it leaves the scene
    bool intersectSphere = raySphereIntersect(eye, u, indexSphere,intersection,normal);
    while(counter < limit && intersectSphere){
        // update information about the intersection, normal, color and incoming direction
        stackOfIntersections[counter].xyz = intersection.xyz;
        stackOfNormals[counter].xyz = normal.xyz;
        stackOfIncoming[counter].xyz = u.xyz;
        stackOfColors[counter] = getColorFromLightSource(u,intersection, normal, indexSphere);
        float cosThethaD = getCosThetha(normal, u);
        float F = fresnelCoeff(cosThethaD, 1./eta);
        stackOfFresnel[counter] = F;
        counter += 1;

        // new reflected ray
        reflectedRay = reflect(u, normal);
        u = normalize(reflectedRay);
        intersectSphere = raySphereIntersect(intersection,u,indexSphere,intersection,normal);
    }

    resultColor = vec4(0.,0.,0.,1.);

    if(counter != 0){
        counter-=1;
        if(rayPlaneIntersect(stackOfIntersections[counter],stackOfIncoming[counter],intersection)){
            resultColor = getColorFromPlane(u);
        }
        for(int j = counter; j>=0; j--){
            resultColor = stackOfFresnel[j]*resultColor + stackOfColors[j];
        }
    } else {
        if(rayPlaneIntersect(eye,u,intersection)){
            resultColor = getColorFromPlane(u);
        }
    }
    

    // // NO INTERSECTION
    // if(counter == 0){
    //     if(rayPlaneIntersect(eye,u,intersection)){
    //         resultColor = getColorFromPlane(u);
    //     } else {
    //         resultColor = vec4(0.0,0.0,0.0,1.0);
    //     }
    //     return resultColor;
    // } else {
    //     // NO MORE INTERSECTION: FOLLOW THE RAY STARTING FROM THE LAST INTERSECTION
    //     counter -= 1;
    //     while(counter >= 0){
    //         resultColor +=  stackOfColors[counter];
    //         counter -= 1;
    //     }
    // }
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
