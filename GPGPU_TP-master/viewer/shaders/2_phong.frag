#version 410

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

void main( void )
{
     // This is the place where there's work to be done

     /// we must normalise vectors before using them 
     //eyeVector = normalize(eyeVector);
     //lightVector = normalize(lightVector);
     

     /***** AMBIENT LIGHT SETUP *******/

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
     float cosThethaD = dot(halfVector, lightVector);
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
     // setting the specular lighting - Cs
     vec4 specularLighting = F * vertColor * pow(max(0, dot(vertNormal, halfVector)), shininess) * lightIntensity;
     fragColor = ambientLight + diffuseLighting + specularLighting;


}
