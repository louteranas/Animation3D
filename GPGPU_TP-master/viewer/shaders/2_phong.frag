#version 410
#define PI 3.1415926538

uniform float lightIntensity;
uniform bool blinnPhong;
uniform float shininess;
uniform float eta;
uniform sampler2D shadowMap;

in vec4 eyeVector;
in vec4 lightVector;
in vec4 vertColor;
in vec4 vertNormal;
in vec4 lightSpace;

out vec4 fragColor;

/**
* this fuction calculate the Fresnet Coefficient
* For more details about the variable names, check the TP page
**/
float fresnetCoeff(float cosThethaD,  float eta){
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
     return true;
}

float NormalDistrib(float cosThetaH, float alpha){
     if(!indicatrice(cosThetaH)){
          return 0;
     }
     float frac1 = 1 / (pow(cosThetaH, 4) * PI);
     float tanThetaSquare = (1 - pow(cosThetaH, 2))/pow(cosThetaH, 2);
     float frac2 = pow(alpha, 2)/ pow((pow(alpha, 2) + pow(tanThetaSquare, 2)), 2);
     return frac1 * frac2;
}

/**
* this fuction calculate the GGX distribution G1(thetaI/thetaO)
* For more details about the variable names, check the TP page
**/
float GGXDistrib(float cosTheta, float alpha){
     float tanThetaSquare = (1 - pow(cosTheta, 2))/pow(cosTheta, 2);
     float base = 1 + sqrt(1 + tanThetaSquare * pow(alpha, 2));
     return 2 / base;
}



void main( void )
{
     // This is the place where there's work to be done
     /// we must normalise vectors before using them in case they change
     

     /***** AMBIENT LIGHT SETUP *******/

     // ambient reflection param 
     float Ka = 0.7;
     // setting the ambiantLighting - Ca
     vec4 ambientLight = Ka * vertColor * lightIntensity;

     /********** Diffuse light setup *********/

     // Diffuse reflection param 
     float Kd = 0.7;
     // setting the Diffuse lighting - Cd
     vec4 diffuseLighting = Kd * vertColor * lightIntensity * max(0, dot(vertNormal, normalize(lightVector)));
     
     
     /********** Specular light setup *********/

     // half Vector
     vec4 halfVector = normalize(lightVector + eyeVector);
     // ThetaD, angle between halfVector & lightVector we only need its cos value
     float cosThethaD = dot(halfVector, lightVector); // /(length(halfVector)*length(lightVector)) ? what is the correct formula?;
     // setting the specular lighting - Cs
     vec4 specularLighting = fresnetCoeff(cosThethaD, eta) * vertColor * pow(max(0, dot(vertNormal, halfVector)), shininess) * lightIntensity;

     // ====> complexe eta ?


     /********** Cook-Torrance model *********/
     // Angles: thetaH, thetaI, thetaO
     float thetaHcos = dot(vertNormal,halfVector);
     float thetaIcos = dot(vertNormal,normalize(lightVector));
     float thetaOcos = dot(vertNormal,normalize(eyeVector));

     /************** 
     object color is the sum of ambient, specular and diffuse lights with the object base color 
     ***************/
     fragColor = ambientLight + diffuseLighting + specularLighting;


}
