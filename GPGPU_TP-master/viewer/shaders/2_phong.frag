#version 410
#define PI 3.1415926538

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform float etaComplex;
uniform sampler2D shadowMap;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

/**
* this fuction calculate the Fresnet Coefficient when eta is complexe
* For more details about the variable names, check the TP page
**/
float fresnetCoeffCmp(float cosThethaD){
     // Ci coeff 
     float etaCarreR = pow(eta, 2) - pow(etaComplex, 2);
     float etaCarreI = 2 * eta * etaComplex;
     float CiRcarre = etaCarreR - (1 - pow(cosThethaD, 2));
     float CiImcarre = etaCarreI;
     float CiRaison = pow(CiRcarre + CiImcarre, 0.25);
     float CiArgum = 0.5 * atan(CiImcarre, CiRcarre);
     // calcul Fs
     float fracTopFsRaison = pow(pow(cosThethaD - (CiRaison*cos(CiArgum)),2) + pow(CiRaison * sin(CiArgum), 2), 0.5);
     float fracBottomFsRaison = pow(pow(cosThethaD + (CiRaison*cos(CiArgum)),2) + pow(CiRaison * sin(CiArgum), 2), 0.5);
     float Fs = pow(fracTopFsRaison/fracBottomFsRaison, 2);
     // calcul Fp 
     float fracTopFpRaison = pow(pow(etaCarreR*cosThethaD - CiRaison*cos(CiArgum),2) + pow(etaCarreI*cosThethaD - CiRaison * sin(CiArgum),2), 0.5);
     float fracBottomFpRaison = pow(pow(etaCarreR*cosThethaD + CiRaison*cos(CiArgum),2) + pow(etaCarreI*cosThethaD + CiRaison * sin(CiArgum),2), 0.5);
     float Fp = pow(fracTopFpRaison/fracBottomFpRaison, 2);
     ///// Fresnel coeff
     float F = (Fs + Fp)/2;
     return F;
}




/**
* this fuction calculate the Fresnet Coefficient
* For more details about the variable names, check the TP page
**/
float fresnetCoeffRl(float cosThethaD){
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
     return F;
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
     if(etaComplex > 0.1)
          return fresnetCoeffCmp(cosThethaD) * vertColor * pow(max(0, dot(vertNormal, halfVector)), shininess) * lightIntensity;
     else
          return fresnetCoeffRl(cosThethaD) * vertColor * pow(max(0, dot(vertNormal, halfVector)), shininess) * lightIntensity;
}

/**
* this fuction calculate the specular lighting using the cook-torrance model
* For more details about the variable names, check the TP page
**/
vec4 specularLightingCT(float cosThethaD, vec4 halfVector, float alpha){
     float cosThetaH = dot(vertNormal,halfVector);
     float cosThetaI = dot(vertNormal,lightVector);
     float cosThetaO = dot(vertNormal,eyeVector);
     float top = fresnetCoeffRl(cosThethaD) * NormalDistrib(cosThetaH, alpha) * GGXDistrib(cosThetaI, alpha) * GGXDistrib(cosThetaO, alpha);
     if(etaComplex > 0.1)
          top = fresnetCoeffCmp(cosThethaD) * NormalDistrib(cosThetaH, alpha) * GGXDistrib(cosThetaI, alpha) * GGXDistrib(cosThetaO, alpha);
     float bottom = 4 * cosThetaI * cosThetaO;
     return (top/bottom)*vertColor*lightIntensity;
}


void main( void )
{
     // This is the place where there's work to be done
     /// we must normalise vectors before using them in case they change

     /********** Normalize *******************/
     normalize(vertNormal);
     normalize(eyeVector);
     normalize(lightVector);

     /********** alpha ***********************/
     // alpha = 1-(shininess/200) in order to be between 0 and 1
     // shininess and alpha reversed
     float alpha = (2.-shininess/100);

     /********** Ambient light setup *********/

     // ambient reflection param 
     float Ka = 0.7;
     // setting the ambiantLighting - Ca
     vec4 ambientLight = Ka * vertColor * lightIntensity;


     /********** Diffuse light setup *********/

     // Diffuse reflection param 
     float Kd = 0.7;
     // setting the Diffuse lighting - Cd
     vec4 diffuseLighting = Kd * vertColor * lightIntensity * max(0, dot(vertNormal, lightVector));
     
     
     /********** Specular light setup *********/

     // half Vector
     vec4 halfVector = normalize(lightVector + eyeVector);
     // ThetaD, angle between halfVector & lightVector we only need its cos value
     float cosThethaD = dot(halfVector, lightVector); // /(length(halfVector)*length(lightVector)) ? what is the correct formula? OK because normalize :);
     // setting the specular lighting - Cs
     vec4 specularLighting;
     // float temp = specularLightingCT(cosThethaD, halfVector);
     if(blinnPhong)
          // using the blinn phong model
          specularLighting = specularLightingBP(cosThethaD, halfVector);
     else
          // using the cook torrance model
          specularLighting = specularLightingCT(cosThethaD, halfVector, alpha);


     /************** 
     object color is the sum of ambient, specular and diffuse lights with the object base color 
     ***************/
     fragColor = ambientLight + diffuseLighting + specularLighting;
}
