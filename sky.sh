//
// Functions for evaluating sky model for sky box quad and lighting
//
// See https://github.com/andrewwillmott/sun-sky
//

#ifndef SKY_LIB_SH
#define SKY_LIB_SH

uniform vec4 u_skyInfo[3];
// C++ bind code:
//    const SkyHosek& hk = ...;
//
//    Vec3f cH(hk.mCoeffsXYZ[0][7], hk.mCoeffsXYZ[1][7], hk.mCoeffsXYZ[2][7]);
//    Vec3f cI(hk.mCoeffsXYZ[0][2] - 1.0f, hk.mCoeffsXYZ[1][2] - 1.0f, hk.mCoeffsXYZ[2][2] - 1.0f);
//    Vec3f rd(hk.mRadXYZ);
//
//    Vec4f skyInfo[3] =
//    {
//        { cH, lumScale },
//        { cI, mRoughness },
//        { rd, mParams.skyboxScale }
//    };
//
//    SetUniform("skyInfo", skyInfo);


uniform sampler2D s_skyBRDFTable;
// C++ bind code:
//    SkyHosek& skyHosek = ...;
//    SkyTable& skyTable = ...;
//    SkyBRDF&  skyBRDF  = ...;
//
//    skyHosek.Update(sunDir, turbidity, groundAlbedo, overcast);
//    skyTable.FindThetaGammaTables(skyHosek);
//    skyBRDF.FindBRDFTables(skyTable, skyHosek);
//
//    int width  =     SkyBRDF::kTableSize;
//    int height = 4 * SkyBRDF::kBRDFSamples;
//
//    imageData.resize(sizeof(Vec4f) * width * height);
//    float* f = imageData.data();
//    skyBRDF.FillBRDFTexture(width, height, *(float(*) [][4]) f);
//    UpdateTexture("skyBRDFTable", imageData.size(), imageData.data());
//

vec3 XYZToRGB(vec3 c)
{
    const vec3 XYZ_TO_R = vec3( 3.2404542, -1.5371385, -0.4985314);
    const vec3 XYZ_TO_G = vec3(-0.9692660,  1.8760108,  0.0415560);
    const vec3 XYZ_TO_B = vec3( 0.0556434, -0.2040259,  1.0572252);

    return vec3
    (
        dot(XYZ_TO_R, c),
        dot(XYZ_TO_G, c),
        dot(XYZ_TO_B, c)
    );
}

vec3 skyColourTable(float cosTheta, float cosGamma, float r)
{
    // remap roughness to table rows
    float row_t  = (0.5 + 7.0 * r) / 8.0 / 4.0;
    float row_g  = row_t + 0.25;
    float row_fh = row_t + 0.50;

    cosTheta = cosTheta >= 0.0f ? sqrt(cosTheta) : -sqrt(-cosTheta);        // remap to concentrate around horizon (0): REMAP_THETA
    float t = 0.5 * (cosTheta + 1.0);

    float g = 0.5 * (1.0 - cosGamma);  // remap to concentrate around sun
    g = sqrt(g);

    vec3 F = texture2D(s_skyBRDFTable, vec2(t, row_t)).xyz;
    vec3 G = texture2D(s_skyBRDFTable, vec2(g, row_g)).xyz;

    vec3 cH = u_skyInfo[0].xyz;
    vec3 cI = u_skyInfo[1].xyz;

    vec4 FHandH = texture2D(s_skyBRDFTable, vec2(t, row_fh));
    vec3 H  = cH * vec3_splat(FHandH.w);
    vec3 FH = cH * FHandH.xyz;

    H  +=     cI;
    FH += F * cI;

    // (1 - F(theta)) * (1 + G(phi) + H(theta))
    vec3 XYZ = (vec3_splat(1.0) - F) * (vec3_splat(1.0) + G) + H - FH;

    XYZ = max(XYZ, vec3_splat(0.0));
    XYZ *= u_skyInfo[2].xyz;

    float lumScale = u_skyInfo[0].w;

    return XYZToRGB(XYZ * lumScale);
}

vec3 skyColourTable(vec3 viewDir, vec3 sunDir)
{
    float cosTheta = viewDir.z;
    float cosGamma = dot(viewDir, sunDir);

    return skyColourTable(cosTheta, cosGamma, 0.0);
}

vec3 skyColourTable(vec3 viewDir, vec3 sunDir, float roughness)
{
    float cosTheta = viewDir.z;
    float cosGamma = dot(viewDir, sunDir);

    return skyColourTable(cosTheta, cosGamma, roughness);
}

#endif
