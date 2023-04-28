#include "Common.glsl"

//#define USE_DEFAULT

#ifndef USE_DEFAULT
//max level is 4
uniform vec4 uSHCoeffs[9];

#endif //USE_DEFAULT

in vec2 vTexCoord0;
in vec3 vWorldNormal;
in vec3 vWorldPosition;

out vec4 Target;

void get_lighting_SH(out vec3 L00,
    out vec3 L1_1, out vec3 L10, out vec3 L11,
    out vec3 L2_2, out vec3 L2_1, out vec3 L20, out vec3 L21, out vec3 L22)
{
#ifdef USE_DEFAULT
    vec3 grace00  = vec3( 0.078908,  0.043710,  0.054161);
    vec3 grace1_1 = vec3( 0.039499,  0.034989,  0.060488);
    vec3 grace10  = vec3(-0.033974, -0.018236, -0.026940);
    vec3 grace11  = vec3(-0.029213, -0.005562,  0.000944);
    vec3 grace2_2 = vec3(-0.011141, -0.005090, -0.012231);
    vec3 grace2_1 = vec3(-0.026240, -0.022401, -0.047479);
    vec3 grace20  = vec3(-0.015570, -0.009471, -0.014733);
    vec3 grace21  = vec3( 0.056014,  0.021444,  0.013915);
    vec3 grace22  = vec3( 0.021205, -0.005432, -0.030374);

    L00 = grace00;
    L1_1 = grace1_1, L10 = grace10, L11 = grace11;
    L2_2 = grace2_2, L2_1 = grace2_1, L20 = grace20, L21 = grace21, L22 = grace22;

#else //!USE_DEFAULT
    L00 = uSHCoeffs[0].xyz;
    L1_1 = uSHCoeffs[1].xyz, L10  = uSHCoeffs[2].xyz, L11 = uSHCoeffs[3].xyz;
    L2_2 = uSHCoeffs[4].xyz, L2_1 = uSHCoeffs[5].xyz, L20 = uSHCoeffs[6].xyz, L21 = uSHCoeffs[7].xyz, L22 = uSHCoeffs[8].xyz;
#endif //USE_DEFAULT
}

// vec3 irradiance_coeffs(vec3 n)
// {
// 	float c1 = 0.429043;
// 	float c2 = 0.511664;
// 	float c3 = 0.743125;
// 	float c4 = 0.886227;
// 	float c5 = 0.247708;

//     vec3 L00;
//     vec3 L1_1, L10, L11;
//     vec3 L2_2, L2_1, L20, L21, L22;
//     get_lighting_SH(L00, L1_1, L10, L11, L2_2, L2_1, L20, L21, L22);

//     float x = n.x, y = n.y, z = n.z;
// 	float x2 = x*x, y2 = y*y, z2 = z*z;
// 	float xy = x*y, yz = y*z, xz = x*z;

// 	return c1*L22*(x2-y2) + c3*L20*z2 + c4*L00 - c5*L20 
//             + 2*c1*(L2_2*xy + L21*xz + L2_1*yz) 
//             + 2*c2*(L11*x+L1_1*y+L10*z) ;
// }

vec3 irradiance_coeffs(vec3 n)
{
    vec3 L00;
    vec3 L1_1, L10, L11;
    vec3 L2_2, L2_1, L20, L21, L22;
    get_lighting_SH(L00, L1_1, L10, L11, L2_2, L2_1, L20, L21, L22);

    float x = n.x, y = n.y, z = n.z;
    float x2 = x*x, y2 = y*y, z2 = z*z;
    vec3 result = L00;

    result += L1_1 * -sqrt(3.0f/(4.0f*PI))*y + 
              L10 *   sqrt(3.0f/(4.0f*PI))*z + 
              L11 *  -sqrt(3.0f/(4.0f*PI))*x;

    result += L2_2 *  sqrt(15.0f/( 4.0f*PI))*y*x +
              L2_1 * -sqrt(15.0f/( 4.0f*PI))*y*z +
              L20  *  sqrt( 5.0f/(16.0f*PI))*(3.0f*z2-1.0f) +
              L21  * -sqrt(15.0f/( 4.0f*PI))*x*z +
              L22  *  sqrt(15.0f/(16.0f*PI))*(x2-y2);

    return result;

    // if (L >= 1)
    // {
    //     result[i++] = -sqrt(3.0f/(4.0f*pi))*y;
    //     result[i++] =  sqrt(3.0f/(4.0f*pi))*z;
    //     result[i++] = -sqrt(3.0f/(4.0f*pi))*x;
    // }

    // if (L >= 2)
    // {
    //     result[i++] =  sqrt(15.0f/(4.0f*pi))*y*x;
    //     result[i++] = -sqrt(15.0f/(4.0f*pi))*y*z;
    //     result[i++] =  sqrt(5.0f/(16.0f*pi))*(3.0f*z2-1.0f);
    //     result[i++] = -sqrt(15.0f/(4.0f*pi))*x*z;
    //     result[i++] =  sqrt(15.0f/(16.0f*pi))*(x2-y2);
    // }
}

void main()
{
	vec3 albedo = vec3(1.0);
	vec3 normal = normalize(vWorldNormal);
	vec3 irradiance = irradiance_coeffs(normal) / PI;
	vec3 color = albedo * irradiance;

	color = tonemapLinear(color, uExposure);
	color = applyDithering(color, gl_FragCoord.xy / uResolution, uElapsedTime);
	Target = vec4(color, 1.0);
}
