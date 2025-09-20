//
// Sky quad shader
//

$input v_param0

#include "lib/sky.sh"

uniform vec4 u_sunColour;
uniform vec4 u_sunDirection;
uniform vec4 u_groundColour;

void main()
{
    vec3 dir = normalize(v_param0.xyz);  // world space direction of ray through pixel

    vec3 skyColour = skyColourTable(dir, u_sunDirection.xyz);
    vec3 sunColour = u_sunColour.xyz;

    float overcast = u_sunColour.w;
    float ss = overcast * overcast;

    // Apply sun
    float cosSA = 0.999957;  // cos(sunAngle) = cos(0.5332)
    float cosAA = 1e-5;  // tiny anti-aliasing factor
    float cosSA1 = cosSA + cosAA;  // outer rim
    float cosSA0 = cosSA - 0.01 * saturate(ss * ss - 0.125);  // inner rim -- smaller as cloud cover increases

    float cosSun = dot(dir, u_sunDirection.xyz);

    // smooth out sun term as it gets overcast
    float sunWeight = smoothstep(cosSA0 - 2e-5, cosSA1 + 2e-5f, cosSun);

    vec3 colour = skyColour + sunColour * sunWeight;

    gl_FragColor.xyz = toneMap(colour);
    gl_FragColor.a = 1.0;
}
