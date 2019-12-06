#version 410
#define EPS 0.000001
#define PI 3.1415926538

/****************************************************************************************************************************************************/
/***************************************************************** PARAMETERS ***********************************************************************/
/****************************************************************************************************************************************************/

/* 
* Parameters UNIFORM
*/
uniform float lightIntensity;
uniform sampler2D colorTexture;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;
uniform bool shadowMapping;

/* 
* Parameters IN
*/
in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec2 textCoords;
in vec4 lightSpace;

/* 
* Parameters OUT: result color
*/
out vec4 fragColor;

/* 
* Normalize
*/
vec4 vertNormalN = normalize(vertNormal);
vec4 eyeVectorN = normalize(eyeVector);
vec4 lightVectorN = normalize(lightVector);
vec4 vertColorN = texture2D(colorTexture,textCoords);

/**
* this fuction calculate the Fresnet Coefficient
* For more details about the variable names, check the TP page
**/
float fresnetCoeffRl(float cosThethaD){
     // Ci coeff 
     float Ci = pow(max(EPS,(pow(eta, 2) - (1 - pow(cosThethaD, 2)))), 0.5);
     // Fs coef
     float fracFs = (cosThethaD - Ci) / (cosThethaD + Ci);
     float Fs = pow(abs(fracFs), 2);
     // Fp coeff
     float fracFp = (pow(eta,2)*cosThethaD - Ci) / (pow(eta,2)*cosThethaD + Ci);
     float Fp = pow(abs(fracFp), 2);
     ///// Fresnel coeff
     float F = (Fs + Fp)/2;
     return F;
}

/**
* this function computes the indicatrice of the cosThetaH (false if <0, true else)
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

/**
* this fuction calculates the microfacet normal distribution D(thetaH)
* For more details about the variable names, check the TP page
**/
float NormalDistrib(float cosThetaH, float alpha){
     if(!indicatrice(cosThetaH)){
          return 0;
     }
     float frac1 = 1 / (pow(cosThetaH, 4) * PI);
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
vec4 specularLightingBP(float cosThethaD, vec4 halfVector){
   return fresnetCoeffRl(cosThethaD) * vertColorN * pow(max(EPS, dot(vertNormalN, halfVector)), shininess) * lightIntensity;
}

/**
* this fuction calculate the specular lighting using the cook-torrance model
* For more details about the variable names, check the TP page
**/
vec4 specularLightingCT(float cosThethaD, vec4 halfVector, float alpha){
     float cosThetaH = dot(vertNormalN,halfVector);
     float cosThetaI = dot(vertNormalN,lightVectorN);
     float cosThetaO = dot(vertNormalN,eyeVectorN);
     float top = fresnetCoeffRl(cosThethaD) * NormalDistrib(cosThetaH, alpha) * GGXDistrib(cosThetaI, alpha) * GGXDistrib(cosThetaO, alpha);
     float bottom = 4 * cosThetaI * cosThetaO;
     return (top/bottom)*vertColorN*lightIntensity;
}


/****************************************************************************************************************************************************/
/********************************************************************** MAIN ************************************************************************/
/****************************************************************************************************************************************************/

void main( void )
{

     /********** alpha ***********************/
     // alpha = 1-(shininess/200) in order to be between 0 and 1
     // shininess and alpha reversed
     float alpha = (2.-shininess/100);

     /********** Ambient light setup *********/
     // ambient reflection param 
     float Ka = 0.1;
     // setting the ambiantLighting - Ca
     vec4 ambientLight = Ka * vertColorN * lightIntensity;

     /********** Diffuse light setup *********/
     // Diffuse reflection param 
     float Kd = 0.5;
     // setting the Diffuse lighting - Cd
     vec4 diffuseLighting = Kd * vertColorN * lightIntensity * max(EPS, dot(vertNormalN, lightVectorN));
     
     /********** Specular light setup *********/
     // half Vector
     vec4 halfVector = normalize(lightVectorN + eyeVectorN);
     // ThetaD, angle between halfVector & lightVectorN we only need its cos value
     float cosThethaD = dot(halfVector, lightVectorN); // /(length(halfVector)*length(lightVectorN)) ? what is the correct formula? OK because normalize :);
     // setting the specular lighting - Cs
     vec4 specularLighting;
     // float temp = specularLightingCT(cosThethaD, halfVector);
     if(blinnPhong){
          // using the blinn phong model
          specularLighting = specularLightingBP(cosThethaD, halfVector);
     } else {
          // using the cook torrance model
          specularLighting = specularLightingCT(cosThethaD, halfVector, alpha);
     }

     /********** Shadow mapping *********/
     if(shadowMapping)
     {
          ////////////////////////////// TO TEST
          vec3 coord = lightSpace.xyz/lightSpace.w;
          // transform to [-1,1] then [0,1]
          coord = coord*0.5+0.5;
          // depth value in shadow map
          float depthValue = (texture(shadowMap, coord.xy).r);
          fragColor = vec4(vec3(depthValue),1);
          // // distance between pixel position and light position
          // float distanceLightSource = coord.z;
          // // check if depth < distance
          // if(depthValue - distanceLightSource < EPS){
          //      // object between the pixel and the light => color = shadow
          //      fragColor = ambientLight;
          // } else {
          //      // no object between the pixel and the light => color = phong model
          //      fragColor = ambientLight + diffuseLighting + specularLighting;
          // }
     } else {
        // object color is the sum of ambient, specular and diffuse lights with the object base color 
        fragColor = ambientLight + diffuseLighting + specularLighting;
     }

}
