//
//  File:       SunSky.cpp
//
//  Function:   Implements SunSky.h
//
//  Author:     Andrew Willmott, Preetham sun/sky based on code by Brian Smits,
//              Hosek code used L. Hosek & A. Wilkie code Version 1.4a as reference
//
//  Copyright:  Andrew Willmott, see HosekDataXYZ.h for Hosek data copyright
//

#include "VL234f.h"

#include "SunSky.h"

#include "HosekDataXYZ.h"

#include <stdint.h>
#include <float.h>
#include <assert.h>

using namespace SSLib;

// #define SIM_CLAMP                // emulate normalised integer texture, for CPU-side checking
// #define SUPPORT_OVERCAST_CLAMP   // support overcast functionality via tables. (Without overcast, the theta table is within 0-1 anyway.)
// #define HOSEK_G_FIX              // experimental fix for hue ringing in rare sunset occasions with BRDF version of Hosek
// #define HOSEK_BRDF_ANALYTIC_H    // for A/B'ing H/FH tables vs. direct ZH evaluation
// #define LOCAL_DEBUG              // dump debug/tuning info

#ifdef LOCAL_DEBUG
    #include <stdio.h>
#endif

namespace
{
    // XYZ/RGB for sRGB primaries
    const Vec3f kXYZToR( 3.2404542f, -1.5371385f, -0.4985314f);
    const Vec3f kXYZToG(-0.9692660f,  1.8760108f,  0.0415560f);
    const Vec3f kXYZToB( 0.0556434f, -0.2040259f,  1.0572252f);

    const Vec3f kRGBToX(0.4124564f,  0.3575761f,  0.1804375f);
    const Vec3f kRGBToY(0.2126729f,  0.7151522f,  0.0721750f);
    const Vec3f kRGBToZ(0.0193339f,  0.1191920f,  0.9503041f);

    const Vec2f kClearChroma       (2.0f / 3.0f, 1.0f / 3.0f);
    const Vec2f kOvercastChroma    (1.0f / 3.0f, 1.0f / 3.0f);
    const Vec2f kPartlyCloudyChroma(1.0f / 3.0f, 1.0f / 3.0f);

    inline Vec3f xyYToXYZ(const Vec3f& c)
    {
        return Vec3f(c.x, c.y, 1.0f - c.x - c.y) * c.z / c.y;
    }

    inline Vec3f xyYToRGB(const Vec3f& xyY)
    {
        Vec3f XYZ = xyYToXYZ(xyY);

        return Vec3f
        (
            dot(kXYZToR, XYZ),
            dot(kXYZToG, XYZ),
            dot(kXYZToB, XYZ)
        );
    }

    inline Vec3f XYZToRGB(const Vec3f& XYZ)
    {
        return Vec3f
        (
            dot(kXYZToR, XYZ),
            dot(kXYZToG, XYZ),
            dot(kXYZToB, XYZ)
        );
    }

    inline Vec3f RGBToXYZ(const Vec3f& rgb)
    {
        return Vec3f
        (
            dot(kRGBToX, rgb),
            dot(kRGBToY, rgb),
            dot(kRGBToZ, rgb)
        );
    }

    inline float ClampUnit(float s)
    {
        if (s <= 0.0f)
            return 0.0f;
        if (s >= 1.0f)
            return 1.0f;
        return s;
    }

    inline float ClampPositive(float s)
    {
        if (s <= 0.0f)
            return 0.0f;
        return s;
    }
    inline Vec3f ClampPositive(Vec3f v)
    {
        return { ClampPositive(v.x), ClampPositive(v.y), ClampPositive(v.z) };
    }

#ifdef SIM_CLAMP
    inline Vec3f ClampUnit(Vec3f v)
    {
        return { ClampUnit(v.x), ClampUnit(v.y), ClampUnit(v.z) };
    }
#endif

    const float    kFToIBiasF32 = (float)(uint32_t(3) << 22);
    const int32_t  kFToIBiasS32 = 0x4B400000;

    inline int32_t FloorToInt32(float s) // force int 'fast' version of this
    {
        union
        {
            float f;
            int32_t i;
        } fi;

        fi.f = floorf(s) + kFToIBiasF32;
        return fi.i - kFToIBiasS32;
    }

    inline float LerpClamp(float s)
    {
        if (s <= 0.0f)         return 0.0f;
        if (s >= 1.0f - 1e-6f) return 1.0f - 1e-6f;
        return s;
    }

    template<class T_V> T_V Lerp(float s, int n, const T_V c[])
    {
        s = LerpClamp(s);
        assert(s >= 0.0f && s < 1.0f);

        s *= n - 1;

        int si0 = FloorToInt32(s);
        int si1 = (si0 + 1);

        float sf = s - si0;

        return c[si0] * (1 - sf)
             + c[si1] *    sf   ;
    }

    template<class T_V> T_V BiLerp(float s, float t, int w, int h, const T_V c[])
    {
        s = LerpClamp(s);
        t = LerpClamp(t);
        assert(s >= 0.0f && t >= 0.0f);
        assert(s  < 1.0f && t  < 1.0f);

        s *= w - 1;
        t *= h - 1;

        int si0 = FloorToInt32(s);
        int ti0 = FloorToInt32(t);

        int si1 = (si0 + 1);
        int ti1 = (ti0 + 1);

        float sf = s - si0;
        float tf = t - ti0;

        ti0 *= w;
        ti1 *= w;

        return c[si0 + ti0] * (1 - sf) * (1 - tf)
             + c[si1 + ti0] *    sf    * (1 - tf)
             + c[si0 + ti1] * (1 - sf) *    tf
             + c[si1 + ti1] *    sf    *    tf;
    }

    inline float DegreesToRadians(float d)
    {
        return d * (vl_twoPi / 360.0f);
    }
}

Vec3f SSLib::SunDirection(float timeOfDay, float timeZone, int julianDay, float latitude, float longitude)
{
    float solarTime = timeOfDay
        +  (0.170f * sinf(4 * vl_pi * (julianDay - 80) / 373)
          - 0.129f * sinf(2 * vl_pi * (julianDay -  8) / 355))
        + (longitude / 15 - timeZone);

    float solarDeclination = (0.4093f * sinf(2 * vl_pi * (julianDay - 81) / 368));

    float latRads = DegreesToRadians(latitude);

    float sx = cosf(solarDeclination)                 * sinf(vl_pi * solarTime / 12);
    float sy = cosf(latRads) * sinf(solarDeclination)
             + sinf(latRads) * cosf(solarDeclination) * cosf(vl_pi * solarTime / 12);
    float sz = sinf(latRads) * sinf(solarDeclination)
             - cosf(latRads) * cosf(solarDeclination) * cosf(vl_pi * solarTime / 12);

    return Vec3f(sx, sy, sz);
}

const float SSLib::kSunDiameter   = 1.392f;
const float SSLib::kSunDistance   = 149.6f;
const float SSLib::kSunCosAngle   = sqrtf(1.0f - sqr(0.5f * kSunDiameter / kSunDistance));   // = 0.999989
const float SSLib::kSunSolidAngle = vl_twoPi * (1.0f - kSunCosAngle);  // = 6.79998e-05 steradians


//------------------------------------------------------------------------------
// Preetham sun model
//------------------------------------------------------------------------------

namespace
{
    // From Preetham. Was used by Hosek as sun source, so can be used in conjunction
    // with either.
    const float kSunRadiance[16][16][3] =
    {
        { {39.4028, 1.98004, 5.96046e-08}, {68821.4, 29221.3, 3969.28}, {189745, 116333, 43283.4}, {284101, 199843, 103207}, {351488, 265139, 161944}, {400584, 315075, 213163}, {437555, 353806, 256435}, {466261, 384480, 292823}, {489140, 409270, 323569}, {507776, 429675, 349757}, {523235, 446739, 372260}, {536260, 461207, 391767}, {547379, 473621, 408815}, {556978, 484385, 423827}, {565348, 493805, 437137}, {572701, 502106, 449002} },
        { {34.9717, 0.0775114, 0}, {33531, 11971.9, 875.627}, {127295, 71095, 22201.3}, {216301, 142827, 66113.9}, {285954, 205687, 115900}, {339388, 256990, 163080}, {380973, 298478, 205124}, {414008, 332299, 241816}, {440780, 360220, 273675}, {462869, 383578, 301382}, {481379, 403364, 325586}, {497102, 420314, 346848}, {510615, 434983, 365635}, {522348, 447795, 382333}, {532628, 459074, 397255}, {541698, 469067, 410647} },
        { {10.0422, 0, 0.318865}, {16312.8, 4886.47, 84.98}, {85310.4, 43421.5, 11226.2}, {164586, 102046, 42200.5}, {232559, 159531, 82822.4}, {287476, 209581, 124663}, {331656, 251771, 163999}, {367569, 287173, 199628}, {397168, 317025, 231420}, {421906, 342405, 259652}, {442848, 364181, 284724}, {460784, 383030, 307045}, {476303, 399483, 326987}, {489856, 413955, 344876}, {501789, 426774, 360988}, {512360, 438191, 375548} },
        { {2.3477, 5.96046e-08, 0.129991}, {117.185, 30.0648, 0}, {57123.3, 26502.1, 5565.4}, {125170, 72886.2, 26819.8}, {189071, 123708, 59081.9}, {243452, 170892, 95209.2}, {288680, 212350, 131047}, {326303, 248153, 164740}, {357842, 278989, 195638}, {384544, 305634, 223657}, {407381, 328788, 248954}, {427101, 349038, 271779}, {444282, 366866, 292397}, {459372, 382660, 311064}, {472723, 396734, 328012}, {484602, 409337, 343430} },
        { {0.383395, 0, 0.027703}, {58.0534, 12.8383, 0}, {38221.6, 16163.6, 2681.55}, {95147.4, 52043, 16954.8}, {153669, 95910.9, 42062}, {206127, 139327, 72640.8}, {251236, 179082, 104653}, {289639, 214417, 135896}, {322383, 245500, 165343}, {350467, 272796, 192613}, {374734, 296820, 217644}, {395864, 318050, 240533}, {414400, 336900, 261440}, {430773, 353719, 280544}, {445330, 368800, 298027}, {458337, 382374, 314041} },
        { {0.0560895, 0, 0.00474608}, {44.0061, 8.32402, 0}, {25559, 9849.99, 1237.01}, {72294.8, 37148.7, 10649}, {124859, 74345.6, 29875.8}, {174489, 113576, 55359.1}, {218617, 151011, 83520.3}, {257067, 185252, 112054}, {290413, 216016, 139698}, {319390, 243473, 165842}, {344686, 267948, 190241}, {366896, 289801, 212852}, {386513, 309371, 233736}, {403942, 326957, 252998}, {419513, 342823, 270764}, {433487, 357178, 287149} },
        { {0.00811136, 0, 0.000761211}, {38.0318, 6.09287, 0}, {17083.4, 5996.83, 530.476}, {54909.7, 26508.7, 6634.5}, {101423, 57618.7, 21163.3}, {147679, 92573, 42135.2}, {190207, 127327, 66606.4}, {228134, 160042, 92352.6}, {261593, 190061, 117993}, {291049, 217290, 142758}, {317031, 241874, 166258}, {340033, 264051, 188331}, {360490, 284081, 208945}, {378771, 302212, 228135}, {395184, 318667, 245976}, {409974, 333634, 262543} },
        { {0.00118321, 0, 0.000119328}, {34.5228, 4.62524, 0}, {11414.1, 3646.94, 196.889}, {41690.9, 18909.8, 4091.39}, {82364.6, 44646.9, 14944.8}, {124966, 75444.4, 32024.3}, {165467, 107347, 53075.4}, {202437, 138252, 76076.7}, {235615, 167214, 99627}, {265208, 193912, 122858}, {291580, 218327, 145272}, {315124, 240580, 166611}, {336208, 260851, 186761}, {355158, 279331, 205696}, {372256, 296206, 223440}, {387729, 311636, 240030} },
        { {0.000174701, 0, 1.84774e-05}, {31.4054, 3.4608, 0}, {7624.24, 2215.02, 48.0059}, {31644.8, 13484.4, 2490.1}, {66872.4, 34589.1, 10515}, {105728, 61477.4, 24300.5}, {143926, 90494.6, 42256.1}, {179617, 119420, 62635.3}, {212200, 147105, 84088.4}, {241645, 173041, 105704}, {268159, 197064, 126911}, {292028, 219187, 147374}, {313550, 239512, 166913}, {333008, 258175, 185447}, {350650, 275321, 202953}, {366683, 291081, 219433} },
        { {2.61664e-05, 0, 2.86102e-06}, {27.3995, 2.42835, 5.96046e-08}, {391.889, 104.066, 0}, {24013.1, 9611.97, 1489.37}, {54282.4, 26792.1, 7366.53}, {89437, 50090, 18406.3}, {125174, 76280.7, 33609.8}, {159354, 103145, 51538.2}, {191098, 129407, 70945.4}, {220163, 154409, 90919.4}, {246607, 177864, 110847}, {270613, 199690, 130337}, {292410, 219912, 149156}, {312229, 238614, 167173}, {330289, 255902, 184328}, {346771, 271876, 200589} },
        { {3.93391e-06, 0, 4.76837e-07}, {21.8815, 1.51091, 0}, {106.645, 26.2423, 0}, {18217.8, 6848.77, 869.811}, {44054, 20748.7, 5134.5}, {75644.5, 40807, 13913.2}, {108852, 64293.6, 26704.2}, {141364, 89082.8, 42380.1}, {172081, 113831, 59831.4}, {200579, 137777, 78179.7}, {226776, 160529, 96794.7}, {250759, 181920, 115250}, {272686, 201910, 133270}, {292739, 220530, 150685}, {311103, 237847, 167398}, {327934, 253933, 183349} },
        { {6.55651e-07, 0, 1.19209e-07}, {15.4347, 0.791314, 0}, {67.98, 15.4685, 0}, {13818.5, 4877.71, 490.832}, {35746.5, 16065.3, 3556.94}, {63969.8, 33240.3, 10492.5}, {94648, 54185.5, 21192.5}, {125394, 76932.4, 34825.1}, {154946, 100125, 50435.6}, {182726, 122930, 67203.7}, {208530, 144877, 84504.4}, {232352, 165726, 101891}, {254283, 185376, 119059}, {274458, 203811, 135807}, {293024, 221062, 152009}, {310113, 237169, 167579} },
        { {5.96046e-08, 0, 0}, {9.57723, 0.336247, 0}, {52.9113, 11.1074, 0}, {10479.8, 3472.19, 262.637}, {29000.9, 12436.5, 2445.87}, {54089.5, 27073.4, 7891.84}, {82288.3, 45662.7, 16796.5}, {111218, 66434.7, 28595.3}, {139508, 88064, 42494.5}, {166453, 109678, 57749.2}, {191743, 130747, 73756.6}, {215288, 150968, 90064.3}, {237114, 170191, 106348}, {257311, 188355, 122384}, {275989, 205455, 138022}, {293255, 221507, 153152} },
        { {0, 0, 0}, {5.37425, 0.109694, 0}, {44.9811, 8.68891, 5.96046e-08}, {7946.76, 2470.32, 128.128}, {23524.7, 9625.27, 1666.58}, {45729.5, 22047.9, 5917.85}, {71535.2, 38477.1, 13293.2}, {98636.4, 57365.7, 23460.6}, {125598, 77452, 35785}, {151620, 97851, 49607}, {176299, 117990, 64359}, {199469, 137520, 79594.4}, {221098, 156245, 94979.6}, {241228, 174066, 110274}, {259937, 190947, 125309}, {277307, 206875, 139956} },
        { {0, 0, 0}, {2.83079, 0.0199037, 0}, {40.0718, 7.10214, 0}, {6025.35, 1756.45, 51.1916}, {19080.1, 7447.79, 1122.67}, {38657, 17952.9, 4422.16}, {62181.1, 32419.5, 10503.8}, {87471.2, 49531.4, 19230.6}, {113069, 68115.1, 30117.9}, {138102, 87295.1, 42596.4}, {162092, 106474, 56143.2}, {184805, 125266, 70327.1}, {206156, 143438, 84812.9}, {226144, 160857, 99349.8}, {244814, 177459, 113755}, {262220, 193206, 127887} },
        { {0, 0, 0}, {1.43779, 0, 0.00738072}, {36.6245, 5.93644, 0}, {4568.17, 1248.02, 9.13028}, {15473.4, 5761.51, 745.266}, {32674.7, 14616.6, 3291.16}, {54045.1, 27313.1, 8284.85}, {77563.8, 42764.4, 15747.9}, {101783, 59900.8, 25332.8}, {125782, 77874.7, 36561.6}, {149022, 96078.4, 48962}, {171213, 114101, 62125.3}, {192218, 131678, 75721.7}, {211998, 148648, 89495.8}, {230564, 164920, 103255}, {247950, 180437, 116847}}
    };
}

Vec3f SSLib::SunRGB(float cosTheta, float turbidity)
{
    if (cosTheta < 0.0f)
        return vl_0;

    float s = cosTheta;
    float t = (turbidity - 2.0f) / 10.0f; // useful range is 2-12
    Vec3f sun = BiLerp(s, t, 16, 16, (const Vec3f*) kSunRadiance);

    // 683 converts from watts to candela at 540e12 hz. Really, we should be weighting by luminous efficiency curve rather than CIE_Y
    sun *=  683.0f;
    sun *= kSunSolidAngle;

    return sun;
}


//------------------------------------------------------------------------------
// CIE models
//------------------------------------------------------------------------------

float SSLib::CIEOvercastSkyLuminance(const Vec3f& v, float Lz)
{
    float cosTheta = v.z;

    return Lz * (1.0f + 2.0f * cosTheta) / 3.0f;
}

float SSLib::CIEClearSkyLuminance(const Vec3f& v, const Vec3f& toSun, float Lz)
{
    float cosThetaV = v.z;

    if (cosThetaV < 0.0f)
        cosThetaV = 0.0f;

    float cosThetaS = toSun.z;
    float thetaS    = acosf(cosThetaS);

    float cosGamma = dot(toSun, v);
    float gamma    = acosf(cosGamma);

    float top1 = 0.91f + 10 * expf(-3 * gamma ) + 0.45f * sqr(cosGamma);
    float bot1 = 0.91f + 10 * expf(-3 * thetaS) + 0.45f * sqr(cosThetaS);

    float top2 = 1 - expf(-0.32f / (cosThetaV + 1e-6f));
    float bot2 = 1 - expf(-0.32f);

    return Lz * (top1 * top2) / (bot1 * bot2);
}

float SSLib::CIEPartlyCloudySkyLuminance(const Vec3f& v, const Vec3f& toSun, float Lz)
{
    float cosThetaV = v.z;

    float cosThetaS = toSun.z;
    float thetaS    = acosf(cosThetaS);

    float cosGamma = dot(toSun, v);
    float gamma    = acosf(cosGamma);

    float top1 = 0.526f + 5 * expf(-1.5f * gamma);
    float bot1 = 0.526f + 5 * expf(-1.5f * thetaS);

    float top2 = 1 - expf(-0.8f / (cosThetaV + 1e-6f));
    float bot2 = 1 - expf(-0.8f);

    return Lz * (top1 * top2) / (bot1 * bot2);
}

// More modern setup, see Darula & Kittler 2002
namespace
{
    const float kCIEStandardSkyCoeffs[15][5] =
    {
        4.0, -0.70,  0, -1.0, 0,        // Overcast. When normalised this is a fit of the older (1 + 2cos(theta)) / 3 formula found in CIEOvercastSkyLuminance.
        4.0, -0.70,  2, -1.5, 0.15,     // Overcast, with steep luminance gradation and slight brightening towards the sun

        1.1, -0.80,  0, -1.0, 0,        // Overcast, moderately graded with azimuthal uniformity
        1.1, -0.80,  2, -1.5, 0.15,     // Overcast, moderately graded and slight brightening towards the sun

        0.0, -1.00,  0, -1.0, 0,        // Sky of uniform luminance
        0.0, -1.00,  2, -1.5, 0.15,     // Partly cloudy sky, no gradation towards zenith, slight brightening towards the sun
        0.0, -1.00,  5, -2.5, 0.30,     // Partly cloudy sky, no gradation towards zenith, brighter circumsolar region
        0.0, -1.00, 10, -3.0, 0.45,     // Partly cloudy sky, no gradation towards zenith, distinct solar corona

       -1.0, -0.55,  2, -1.5, 0.15,     // Partly cloudy, with the obscured sun
       -1.0, -0.55,  5, -2.5, 0.30,     // Partly cloudy, with brighter circumsolar region
       -1.0, -0.55, 10, -3.0, 0.45,     // White-blue sky with distinct solar corona

       -1.0, -0.32, 10, -3.0, 0.45,     // CIE Standard Clear Sky, low illuminance turbidity. T <= 2.45
       -1.0, -0.32, 16, -3.0, 0.30,     // CIE Standard Clear Sky, polluted atmosphere

       -1.0, -0.15, 16, -3.0, 0.30,     // Cloudless turbid sky with broad solar corona
       -1.0, -0.15, 24, -2.8, 0.15,     // White-blue turbid sky with broad solar corona
    };

    inline float CIELumRatio(const Vec3f& v, const Vec3f& toSun, float a, float b, float c, float d, float e)
    {
        float cosThetaV = v.z;

        float cosThetaS = toSun.z;
        float thetaS    = acosf(cosThetaS);

        float cosGamma = dot(toSun, v);
        float gamma    = acosf(cosGamma);

        float top1 = 1.0f + a * expf(b / (cosThetaV + 1e-6f));
        float bot1 = 1.0f + a * expf(b);

        float top2 = 1.0f + c * expf(d * gamma ) + e * sqr(cosGamma);
        float bot2 = 1.0f + c * expf(d * thetaS) + e * sqr(cosThetaS);

        return (top1 * top2) / (bot1 * bot2);
    }
}

float SSLib::CIEStandardSky(int type, const Vec3f& v, const Vec3f& toSun, float Lz)
{
    const float* c = kCIEStandardSkyCoeffs[type];
    return CIELumRatio(v, toSun, c[0], c[1], c[2], c[3], c[4]) * Lz;
}

float SSLib::ZenithLuminance(float thetaS, float T)
{
    float chi = (4.0f / 9.0f - T / 120.0f) * (vl_pi - 2.0f * thetaS);
    float Lz = (4.0453f * T - 4.9710f) * tanf(chi) - 0.2155f * T + 2.4192f;
    Lz *= 1000.0;   // conversion from kcd/m^2 to cd/m^2
    return Lz;
}





//------------------------------------------------------------------------------
// - Preetham
//------------------------------------------------------------------------------

namespace
{
    //
    // The Perez function is as follows:
    //
    //    P(t, g) =  (1 + A e ^ (B / cos(t))) (1 + C e ^ (D g)  + E cos(g)  ^ 2)
    //               --------------------------------------------------------
    //               (1 + A e ^ B)            (1 + C e ^ (D ts) + E cos(ts) ^ 2)
    //

    // A: sky
    // B: sky tightness
    // C: sun
    // D: sun tightness, higher = tighter
    // E: rosy hue around sun

    inline float PerezUpper(const float* lambdas, float cosTheta, float gamma, float cosGamma)
    {
        return  (1.0f + lambdas[0] * expf(lambdas[1] / (cosTheta + 1e-6f)))
              * (1.0f + lambdas[2] * expf(lambdas[3] * gamma) + lambdas[4] * sqr(cosGamma));
    }

    inline float PerezLower(const float* lambdas, float cosThetaS, float thetaS)
    {
        return  (1.0f + lambdas[0] * expf(lambdas[1]))
              * (1.0f + lambdas[2] * expf(lambdas[3] * thetaS) + lambdas[4] * sqr(cosThetaS));
    }
}

cSunSkyPreetham::cSunSkyPreetham() :
    mToSun(vl_z),

    // mPerez_x[5] = {}
    // mPerez_y[5] = {}
    // mPerez_Y[5] = {}
    mZenith(vl_0),
    mPerezInvDen(vl_1)
{
}

void cSunSkyPreetham::Update(const Vec3f& sun, float turbidity, float overcast, float horizCrush)
{
    mToSun = sun;

    float   T = turbidity;

    mPerez_Y[0] =  0.17872f * T - 1.46303f;
    mPerez_Y[1] = -0.35540f * T + 0.42749f;
    mPerez_Y[2] = -0.02266f * T + 5.32505f;
    mPerez_Y[3] =  0.12064f * T - 2.57705f;
    mPerez_Y[4] = -0.06696f * T + 0.37027f;

    mPerez_x[0] = -0.01925f * T - 0.25922f;
    mPerez_x[1] = -0.06651f * T + 0.00081f;
    mPerez_x[2] = -0.00041f * T + 0.21247f;
    mPerez_x[3] = -0.06409f * T - 0.89887f;
    mPerez_x[4] = -0.00325f * T + 0.04517f;

    mPerez_y[0] = -0.01669f * T - 0.26078f;
    mPerez_y[1] = -0.09495f * T + 0.00921f;
    mPerez_y[2] = -0.00792f * T + 0.21023f;
    mPerez_y[3] = -0.04405f * T - 1.65369f;
    mPerez_y[4] = -0.01092f * T + 0.05291f;

    float cosTheta = mToSun.z;
    float theta  = acosf(cosTheta);    // angle from zenith rather than horizon
    float theta2 = sqr(theta);
    float theta3 = theta2 * theta;
    float T2 = sqr(T);

    // mZenith stored as xyY
    mZenith.z = ZenithLuminance(theta, T);

    mZenith.x =
        ( 0.00165f * theta3 - 0.00374f * theta2 + 0.00208f * theta + 0.00000f) * T2 +
        (-0.02902f * theta3 + 0.06377f * theta2 - 0.03202f * theta + 0.00394f) * T  +
        ( 0.11693f * theta3 - 0.21196f * theta2 + 0.06052f * theta + 0.25885f);

    mZenith.y =
        ( 0.00275f * theta3 - 0.00610f * theta2 + 0.00316f * theta + 0.00000f) * T2 +
        (-0.04214f * theta3 + 0.08970f * theta2 - 0.04153f * theta + 0.00515f) * T  +
        ( 0.15346f * theta3 - 0.26756f * theta2 + 0.06669f * theta + 0.26688f);

    // Adjustments (extensions)

    if (cosTheta < 0.0f)    // Handle sun going below the horizon
    {
        float s = ClampUnit(1.0f + cosTheta * 50.0f);   // goes from 1 to 0 as the sun sets

        // Take C/E which control sun term to zero
        mPerez_x[2] *= s;
        mPerez_y[2] *= s;
        mPerez_Y[2] *= s;
        mPerez_x[4] *= s;
        mPerez_y[4] *= s;
        mPerez_Y[4] *= s;
    }

    if (overcast != 0.0f)      // Handle overcast term
    {
        float invOvercast = 1.0f - overcast;

        // lerp back towards unity
        mPerez_x[0] *= invOvercast;  // main sky chroma -> base
        mPerez_y[0] *= invOvercast;

        // sun flare -> 0 strength/base chroma
        mPerez_x[2] *= invOvercast;
        mPerez_y[2] *= invOvercast;
        mPerez_Y[2] *= invOvercast;
        mPerez_x[4] *= invOvercast;
        mPerez_y[4] *= invOvercast;
        mPerez_Y[4] *= invOvercast;

        // lerp towards a fit of the CIE cloudy sky model: 4, -0.7
        mPerez_Y[0] = lerp(mPerez_Y[0],  4.0f, overcast);
        mPerez_Y[1] = lerp(mPerez_Y[1], -0.7f, overcast);

        // lerp base colour towards white point
        mZenith.x = mZenith.x * invOvercast + 0.333f * overcast;
        mZenith.y = mZenith.y * invOvercast + 0.333f * overcast;
    }

    if (horizCrush != 0.0f)
    {
        // The Preetham sky model has a "muddy" horizon, which can be objectionable in
        // typical game views. We allow artistic control over it.
        mPerez_Y[1] *= horizCrush;
        mPerez_x[1] *= horizCrush;
        mPerez_y[1] *= horizCrush;
    }

    // initialize sun-constant parts of the Perez functions

    Vec3f perezLower    // denominator terms are dependent on sun angle only
    (
        PerezLower(mPerez_x, cosTheta, theta),
        PerezLower(mPerez_y, cosTheta, theta),
        PerezLower(mPerez_Y, cosTheta, theta)
    );

    mPerezInvDen = mZenith / perezLower;
}

Vec3f cSunSkyPreetham::SkyRGB(const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(mToSun, v);
    float gamma    = acosf(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    Vec3f xyY
    (
        PerezUpper(mPerez_x, cosTheta, gamma, cosGamma),
        PerezUpper(mPerez_y, cosTheta, gamma, cosGamma),
        PerezUpper(mPerez_Y, cosTheta, gamma, cosGamma)
    );

    xyY *= mPerezInvDen;

    return xyYToRGB(xyY);
}

float cSunSkyPreetham::SkyLuminance(const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(mToSun, v);
    float gamma    = acosf(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    return PerezUpper(mPerez_Y, cosTheta, gamma, cosGamma) * mPerezInvDen.z;
}

Vec2f cSunSkyPreetham::SkyChroma(const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(mToSun, v);
    float gamma    = acosf(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    return Vec2f
    (
        PerezUpper(mPerez_x, cosTheta, gamma, cosGamma) * mPerezInvDen.x,
        PerezUpper(mPerez_y, cosTheta, gamma, cosGamma) * mPerezInvDen.y
    );
}


//------------------------------------------------------------------------------
// cSunSkyHosek
//------------------------------------------------------------------------------

namespace
{
    inline float EvalQuintic(const float w[6], const float data[6])
    {
        return    w[0] * data[0]
                + w[1] * data[1]
                + w[2] * data[2]
                + w[3] * data[3]
                + w[4] * data[4]
                + w[5] * data[5];
    }

    inline void EvalQuintic(const float w[6], const float data[6][9], float coeffs[9])
    {
        for (int i = 0; i < 9; i++)
        {
            coeffs[i] =   w[0] * data[0][i]
                        + w[1] * data[1][i]
                        + w[2] * data[2][i]
                        + w[3] * data[3][i]
                        + w[4] * data[4][i]
                        + w[5] * data[5][i];
        }
    }

    inline void FindQuinticWeights(float s, float w[6])
    {
        float s1 = s;
        float s2 = s1 * s1;
        float s3 = s1 * s2;
        float s4 = s2 * s2;
        float s5 = s2 * s3;

        float is1 = 1.0f - s1;
        float is2 = is1 * is1;
        float is3 = is1 * is2;
        float is4 = is2 * is2;
        float is5 = is2 * is3;

        w[0] = is5;
        w[1] = is4 * s1 *  5.0f;
        w[2] = is3 * s2 * 10.0f;
        w[3] = is2 * s3 * 10.0f;
        w[4] = is1 * s4 *  5.0f;
        w[5] =       s5;
    }


    float FindHosekCoeffs
    (
        const float   dataset9[2][10][6][9],    // albedo x 2, turbidity x 10, quintics x 6, weights x 9
        const float   datasetR[2][10][6],       // albedo x 2, turbidity x 10, quintics x 6
        float         turbidity,
        float         albedo,
        float         solarElevation,
        float         coeffs[9]
    )
    {
        int tbi = int(floorf(turbidity));

        if (tbi < 1)
            tbi = 1;
        else if (tbi > 9)
            tbi = 9;

        float tbf = turbidity - tbi;

        const float s = powf(solarElevation / vl_halfPi, (1.0f / 3.0f));

        float quinticWeights[6];
        FindQuinticWeights(s, quinticWeights);

        float ic[4][9];

        EvalQuintic(quinticWeights, dataset9[0][tbi - 1], ic[0]);
        EvalQuintic(quinticWeights, dataset9[1][tbi - 1], ic[1]);
        EvalQuintic(quinticWeights, dataset9[0][tbi    ], ic[2]);
        EvalQuintic(quinticWeights, dataset9[1][tbi    ], ic[3]);

        float ir[4] =
        {
            EvalQuintic(quinticWeights, datasetR[0][tbi - 1]),
            EvalQuintic(quinticWeights, datasetR[1][tbi - 1]),
            EvalQuintic(quinticWeights, datasetR[0][tbi    ]),
            EvalQuintic(quinticWeights, datasetR[1][tbi    ]),
        };

        float cw[4] =
        {
            (1.0f - albedo) * (1.0f - tbf),
            albedo          * (1.0f - tbf),
            (1.0f - albedo) * tbf,
            albedo          * tbf,
        };

        for (int i = 0; i < 9; i++)
            coeffs[i] = cw[0] * ic[0][i]
                      + cw[1] * ic[1][i]
                      + cw[2] * ic[2][i]
                      + cw[3] * ic[3][i];

        return cw[0] * ir[0] + cw[1] * ir[1] + cw[2] * ir[2] + cw[3] * ir[3];
    }

    // Hosek:
    // (1 + A e ^ (B / cos(t))) (1 + C e ^ (D g) + E cos(g) ^ 2   + F mieM(g, G)  + H cos(t)^1/2 + (I - 1))
    //
    // These bits are the same as Preetham, but do different jobs in some cases
    // A: sky gradient, carries white -> blue gradient
    // B: sky tightness
    // C: sun, carries most of sun-centred blue term
    // D: sun tightness, higher = tighter
    // E: rosy hue around sun
    //
    // Hosek-specific:
    // F: mie term, does most of the heavy lifting for sunset glow
    // G: mie tuning
    // H: zenith gradient
    // I: constant term balanced with H

    // Notes:
    // A/B still carries some of the "blue" base of sky, but much comes from C/D
    // C/E minimal effect in sunset situations, carry bulk of sun halo in sun-overhead
    // F/G sunset glow, but also takes sun halo from yellowish to white overhead

    float EvalHosekCoeffs(const float coeffs[9], float cosTheta, float gamma, float cosGamma)
    {
        // Current coeffs ordering is AB I CDEF HG
        //                            01 2 3456 78
        const float expM   = expf(coeffs[4] * gamma);   // D g
        const float rayM   = cosGamma * cosGamma;       // Rayleigh scattering
        const float mieM   = (1.0f + rayM) / powf((1.0f + coeffs[8] * coeffs[8] - 2.0f * coeffs[8] * cosGamma), 1.5f);  // G
        const float zenith = sqrtf(cosTheta);           // vertical zenith gradient

        return (1.0f
                     + coeffs[0] * expf(coeffs[1] / (cosTheta + 0.01f))     // A, B
               )
             * (1.0f
                     + coeffs[3] * expM     // C
                     + coeffs[5] * rayM     // E
                     + coeffs[6] * mieM     // F
                     + coeffs[7] * zenith   // H
                     + (coeffs[2] - 1.0f)   // I
               );
    }
}


cSunSkyHosek::cSunSkyHosek() :
    mToSun(vl_z)
{
}

void cSunSkyHosek::Update(const Vec3f& sun, float turbidity, Vec3f rgbAlbedo, float overcast)
{
    mToSun = sun;

    float solarElevation = mToSun.z > 0.0f ? asinf(mToSun.z) : 0.0f;    // altitude rather than zenith, so sin rather than cos

    Vec3f albedo = RGBToXYZ(rgbAlbedo);

    // Note that the hosek coefficients change with time of day, vs. Preetham where the 'upper' coefficients stay the same,
    // and only the scaler mPerezInvDen, consisting of time-dependent normalisation and zenith luminnce factors, changes.
    mRadXYZ[0] = FindHosekCoeffs(kHosekCoeffsX, kHosekRadX, turbidity, albedo.x, solarElevation, mCoeffsXYZ[0]);
    mRadXYZ[1] = FindHosekCoeffs(kHosekCoeffsY, kHosekRadY, turbidity, albedo.y, solarElevation, mCoeffsXYZ[1]);
    mRadXYZ[2] = FindHosekCoeffs(kHosekCoeffsZ, kHosekRadZ, turbidity, albedo.z, solarElevation, mCoeffsXYZ[2]);

    mRadXYZ *= 683; // convert to luminance in lumens

    if (mToSun.z < 0.0f)   // sun below horizon?
    {
        float s = ClampUnit(1.0f + mToSun.z * 50.0f);   // goes from 1 to 0 as the sun sets
        float is = 1.0f - s;

        // Emulate Preetham's zenith darkening
        float darken = ZenithLuminance(acosf(mToSun.z), turbidity) / ZenithLuminance(vl_halfPi, turbidity);

        // Take C/E/F which control sun term to zero
        for (int j = 0; j < 3; j++)
        {
            mCoeffsXYZ[j][3] *= s;
            mCoeffsXYZ[j][5] *= s;
            mCoeffsXYZ[j][6] *= s;

            // Take horizon term H to zero, as it's an orange glow at this point
            mCoeffsXYZ[j][7] *= s;

            // Take I term back to 1
            mCoeffsXYZ[j][2] *= s;
            mCoeffsXYZ[j][2] += is;
        }

        mRadXYZ *= darken;
    }

    if (overcast != 0.0f)      // Handle overcast term
    {
        float is = overcast;
        float s = 1.0f - overcast;     // goes to 0 as we go to overcast

        // Hosek isn't self-normalising, unlike Preetham/CIE, which divides by PreethamLower().
        // Thus when we lerp to the CIE overcast model, we get some non-linearities.
        // We deal with this by using ratios of normalisation terms to balance.
        // Another difference is that Hosek is relative to the average radiance,
        // whereas CIE is the zenith radiance, so rather than taking the zenith
        // as normalising as in CIE, we average over the zenith and two horizon
        // points.
        float cosGammaZ = mToSun.z;
        float gammaZ    = acosf(cosGammaZ);
        float cosGammaH = mToSun.y;
        float gammaHP   = acosf(+mToSun.y);
        float gammaHN   = vl_pi - gammaHP;

        float sc0 = EvalHosekCoeffs(mCoeffsXYZ[1], 1.0f, gammaZ, cosGammaZ) * 2.0f
                  + EvalHosekCoeffs(mCoeffsXYZ[1], 0.0f, gammaHP, +cosGammaH)
                  + EvalHosekCoeffs(mCoeffsXYZ[1], 0.0f, gammaHN, -cosGammaH);

        for (int j = 0; j < 3; j++)
        {
            // sun flare -> 0 strength/base chroma
            // Take C/E/F which control sun term to zero
            mCoeffsXYZ[j][3] *= s;
            mCoeffsXYZ[j][5] *= s;
            mCoeffsXYZ[j][6] *= s;

            // Take H back to 0
            mCoeffsXYZ[j][7] *= s;

            // Take I term back to 1
            mCoeffsXYZ[j][2] *= s;
            mCoeffsXYZ[j][2] += is;

            // Take A/B to  CIE cloudy sky model: 4, -0.7
            mCoeffsXYZ[j][0] = lerp(mCoeffsXYZ[j][0],  4.0f, is);
            mCoeffsXYZ[j][1] = lerp(mCoeffsXYZ[j][1], -0.7f, is);
        }

        float sc1 = EvalHosekCoeffs(mCoeffsXYZ[1], 1.0f, gammaZ, cosGammaZ) * 2.0f
                  + EvalHosekCoeffs(mCoeffsXYZ[1], 0.0f, gammaHP, +cosGammaH)
                  + EvalHosekCoeffs(mCoeffsXYZ[1], 0.0f, gammaHN, -cosGammaH);

        float rescale = sc0 / sc1;
        mRadXYZ *= rescale;

        // move back to white point
        mRadXYZ.x = lerp(mRadXYZ.x, mRadXYZ.y, is);
        mRadXYZ.z = lerp(mRadXYZ.z, mRadXYZ.y, is);
    }

#ifdef LOCAL_DEBUG
    for (int j = 0; j < 3; j++)
    {
        printf("%c: Rad = %8.2f, C[9] = ", 'X' + j, mRadXYZ[j]);
        for (int i = 0; i < 9; i++)
            printf("%5.2f ", mCoeffsXYZ[j][i]);
        printf("\n");
    }
#endif
}

float cSunSkyHosek::SkyLuminance(const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(mToSun, v);
    float gamma    = acosf(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    return EvalHosekCoeffs(mCoeffsXYZ[1], cosTheta, gamma, cosGamma) * mRadXYZ.y;
}

Vec3f cSunSkyHosek::SkyXYZ(const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(mToSun, v);
    float gamma    = acosf(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    Vec3f XYZ =
    {
        EvalHosekCoeffs(mCoeffsXYZ[0], cosTheta, gamma, cosGamma),
        EvalHosekCoeffs(mCoeffsXYZ[1], cosTheta, gamma, cosGamma),
        EvalHosekCoeffs(mCoeffsXYZ[2], cosTheta, gamma, cosGamma)
    };

    XYZ *= mRadXYZ;

    return XYZ;
}

Vec3f cSunSkyHosek::SkyRGB(const Vec3f& v) const
{
    return XYZToRGB(SkyXYZ(v));
}


//------------------------------------------------------------------------------
// cSunSkyTable
//------------------------------------------------------------------------------

// Pre-calculated table version, particularly useful for shader use

namespace
{
    // remap cosGamma to concentrate table around the sun location
    inline float MapGamma(float g)
    {
        return sqrtf(0.5f * (1.0f - g));
    }
    inline float UnmapGamma(float g)
    {
        return 1 - 2 * sqr(g);
    }
}

void cSunSkyTable::FindThetaGammaTables(const cSunSkyPreetham& pt)
{
    float dt = 1.0f / (kTableSize - 1);
    float t = dt * 1e-6f;    // epsilon to avoid divide by 0, which can lead to NaN when m_perez_[1] = 0

    for (int i = 0; i < kTableSize; i++)
    {
        float cosTheta = t;

        mThetaTable[i][0] =  -pt.mPerez_x[0] * expf(pt.mPerez_x[1] / cosTheta);
        mThetaTable[i][1] =  -pt.mPerez_y[0] * expf(pt.mPerez_y[1] / cosTheta);
        mThetaTable[i][2] =  -pt.mPerez_Y[0] * expf(pt.mPerez_Y[1] / cosTheta);

        mMaxTheta = mMaxTheta >= mThetaTable[i][2] ? mMaxTheta : mThetaTable[i][2];

        float cosGamma = UnmapGamma(t);
        float gamma = acosf(cosGamma);

        mGammaTable[i][0] =  pt.mPerez_x[2] * expf(pt.mPerez_x[3] * gamma) + pt.mPerez_x[4] * sqr(cosGamma);
        mGammaTable[i][1] =  pt.mPerez_y[2] * expf(pt.mPerez_y[3] * gamma) + pt.mPerez_y[4] * sqr(cosGamma);
        mGammaTable[i][2] =  pt.mPerez_Y[2] * expf(pt.mPerez_Y[3] * gamma) + pt.mPerez_Y[4] * sqr(cosGamma);

        mMaxGamma = mMaxGamma >= mGammaTable[i][2] ? mMaxGamma : mGammaTable[i][2];

        t += dt;
    }

    mXYZ = false;

#ifdef SUPPORT_OVERCAST_CLAMP
    mMaxTheta = -pt.mPerez_Y[0] * expf(pt.mPerez_Y[1]);
#endif

#ifdef SIM_CLAMP
    for (int i = 0; i < kTableSize; i++)
    {
        mThetaTable[i].z /= mMaxTheta;
        mGammaTable[i].z /= mMaxGamma;
        mThetaTable[i] = ClampUnit(mThetaTable[i]);
        mGammaTable[i] = ClampUnit(mGammaTable[i]);
    }
#endif

#ifdef LOCAL_DEBUG
    printf("PTx: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].x);
    printf("\n");
    printf("PTy: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].y);
    printf("\n");
    printf("PTY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].z);
    printf("\n");

    printf("PGx: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].x);
    printf("\n");
    printf("PGy: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].y);
    printf("\n");
    printf("PGY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].z);
    printf("\n");
#endif
}

void cSunSkyTable::FindThetaGammaTables(const cSunSkyHosek& hk)
{
    float dt = 1.0f / (kTableSize - 1);
    float t = 0.0f;

    const float (&coeffsXYZ)[3][9] = hk.mCoeffsXYZ;

    mMaxGamma = 1.0f;

    for (int i = 0; i < kTableSize; i++)
    {
        float cosTheta = t;
        float cosGamma = UnmapGamma(t);
        float gamma = acosf(cosGamma);

        float rayM = cosGamma * cosGamma;

        for (int j = 0; j < 3; j++)
        {
            const float (&coeffs)[9] = coeffsXYZ[j];

            mThetaTable[i][j] = -coeffs[0] * expf(coeffs[1] / (cosTheta + 0.01f));

            float expM   = expf(coeffs[4] * gamma);
            float mieM   = (1.0f + rayM) / powf((1.0f + coeffs[8] * coeffs[8] - 2.0f * coeffs[8] * cosGamma), 1.5f);

            mGammaTable[i][j] = coeffs[3] * expM + coeffs[5] * rayM + coeffs[6] * mieM;

            mMaxGamma = mMaxGamma >= mGammaTable[i][j] ? mMaxGamma : mGammaTable[i][j];
        }

        t += dt;
    }

    mXYZ = true;

#ifdef SUPPORT_OVERCAST_CLAMP
    mMaxTheta = -coeffsXYZ[1][0] * expf(coeffsXYZ[1][1]);
#endif

#ifdef SIM_CLAMP
    for (int i = 0; i < kTableSize; i++)
    {
        mThetaTable[i] = ClampUnit(mThetaTable[i] / mMaxTheta);
        mGammaTable[i] = ClampUnit(mGammaTable[i] / mMaxGamma);
    }
#endif

#ifdef LOCAL_DEBUG
    printf("HTY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].x);
    printf("\n");
    printf("HTY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].y);
    printf("\n");
    printf("HTY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mThetaTable[i].z);
    printf("\n");

    printf("HGX: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].x);
    printf("\n");
    printf("HGY: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].y);
    printf("\n");
    printf("HGZ: ");
    for (int i = 0; i < kTableSize; i+= 4)
        printf("%5.2f ", mGammaTable[i].z);
    printf("\n");
#endif
}

Vec3f cSunSkyTable::SkyRGB(const cSunSkyPreetham& pt, const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(pt.mToSun, v);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    float t = cosTheta;
    float g = MapGamma(cosGamma);

    Vec3f F = Lerp(t, kTableSize, mThetaTable);
    Vec3f G = Lerp(g, kTableSize, mGammaTable);

#ifdef SIM_CLAMP
    F.z *= mMaxTheta;
    G.z *= mMaxGamma;
#endif

    Vec3f xyY = (Vec3f(vl_1) - F) * (Vec3f(vl_1) + G);

    xyY *= pt.mPerezInvDen;

    return xyYToRGB(xyY);
}

Vec3f cSunSkyTable::SkyRGB(const cSunSkyHosek& hk, const Vec3f& v) const
{
    float cosTheta = v.z;
    float cosGamma = dot(hk.mToSun, v);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    float t = cosTheta;
    float g = MapGamma(cosGamma);

    Vec3f F = Lerp(t, kTableSize, mThetaTable);
    Vec3f G = Lerp(g, kTableSize, mGammaTable);

    const float zenith = sqrtf(cosTheta);
    Vec3f H
    (
        hk.mCoeffsXYZ[0][7] * zenith + (hk.mCoeffsXYZ[0][2] - 1.0f),
        hk.mCoeffsXYZ[1][7] * zenith + (hk.mCoeffsXYZ[1][2] - 1.0f),
        hk.mCoeffsXYZ[2][7] * zenith + (hk.mCoeffsXYZ[2][2] - 1.0f)
    );

#ifdef SIM_CLAMP
    F *= mMaxTheta;
    G *= mMaxGamma;
#endif

    // (1 - F(theta)) * (1 + G(phi) + H(theta))
    Vec3f XYZ = (Vec3f(vl_1) - F) * (Vec3f(vl_1) + G + H);

    XYZ *= hk.mRadXYZ;

    return XYZToRGB(XYZ);
}

namespace
{
    inline uint8_t ToU8(float f)
    {
        if (f <= 0.0f)
            return 0;
        if (f >= 1.0f)
            return 255;

        return uint8_t(f * 255.0f + 0.5f);
    }
}

void cSunSkyTable::FillTexture(int width, int height, uint8_t image[][4]) const
{
    CL_ASSERT(width == kTableSize);
    CL_ASSERT(height == 2);

    (void) width;
    (void) height;

    for (int i = 0; i < width; i++)
    {
        Vec3f c = mThetaTable[i];

        if (mXYZ)
            c /= mMaxTheta;
        else
            c.z /= mMaxTheta;

        image[i][2] = ToU8(c.x);
        image[i][1] = ToU8(c.y);
        image[i][0] = ToU8(c.z);
        image[i][3] = 255;
    }

    image += width;

    for (int i = 0; i < width; i++)
    {
        Vec3f c = mGammaTable[i];

        if (mXYZ)
            c /= mMaxGamma;
        else
            c.z /= mMaxGamma;

        image[i][2] = ToU8(c.x);
        image[i][1] = ToU8(c.y);
        image[i][0] = ToU8(c.z);
        image[i][3] = 255;
    }
}

void cSunSkyTable::FillTexture(int width, int height, float image[][4]) const
{
    CL_ASSERT(width == kTableSize);
    CL_ASSERT(height == 2);

    (void) width;
    (void) height;

    for (int i = 0; i < width; i++)
        *(Vec4f*) image[i] = Vec4f(mThetaTable[i], 1.0f);

    image += width;

    for (int i = 0; i < width; i++)
        *(Vec4f*) image[i] = Vec4f(mGammaTable[i], 1.0f);
}

//------------------------------------------------------------------------------
// cSunSkyBRDF
//------------------------------------------------------------------------------

// Extended version that enables lerping between perfect mirror and fully
// diffuse.

namespace
{
    // ZH routines from SHLib
    const float kZH_Y_0 = sqrtf( 1 / (   4 * vl_pi));  //         1
    const float kZH_Y_1 = sqrtf( 3 / (   4 * vl_pi));  //         z
    const float kZH_Y_2 = sqrtf( 5 / (  16 * vl_pi));  // 1/2     (3 z^2 - 1)
    const float kZH_Y_3 = sqrtf( 7 / (  16 * vl_pi));  // 1/2     (5 z^3 - 3 z)
    const float kZH_Y_4 = sqrtf( 9 / ( 256 * vl_pi));  // 1/8     (35 z^4 - 30 z^2 + 3)
    const float kZH_Y_5 = sqrtf(11 / ( 256 * vl_pi));  // 1/8     (63 z^5 - 70 z^3 + 15 z)
    const float kZH_Y_6 = sqrtf(13 / (1024 * vl_pi));  // 1/16    (231 z^6 - 315 z^4 + 105 z^2 - 5)

    void CalcCosPowerSatZH7(float n, float zcoeffs[7])
    {
        zcoeffs[0] =   1.0f / (n + 1);
        zcoeffs[1] =   1.0f / (n + 2);
        zcoeffs[2] =   3.0f / (n + 3) -   1.0f / (n + 1);
        zcoeffs[3] =   5.0f / (n + 4) -   3.0f / (n + 2);
        zcoeffs[4] =  35.0f / (n + 5) -  30.0f / (n + 3) +   3.0f / (n + 1);
        zcoeffs[5] =  63.0f / (n + 6) -  70.0f / (n + 4) +  15.0f / (n + 2);
        zcoeffs[6] = 231.0f / (n + 7) - 315.0f / (n + 5) + 105.0f / (n + 3) - 5.0f / (n + 1);

        // apply norm constants
        zcoeffs[0] *= vl_twoPi * kZH_Y_0;
        zcoeffs[1] *= vl_twoPi * kZH_Y_1;
        zcoeffs[2] *= vl_twoPi * kZH_Y_2;
        zcoeffs[3] *= vl_twoPi * kZH_Y_3;
        zcoeffs[4] *= vl_twoPi * kZH_Y_4;
        zcoeffs[5] *= vl_twoPi * kZH_Y_5;
        zcoeffs[6] *= vl_twoPi * kZH_Y_6;

        // [0]: 2pi sqrtf( 1 / (   4 * vl_pi)) / 1
        // we'll multiply by alpha = sqrtf(4.0f * vl_pi / (2 * i + 1)) in convolution, leaving 2pi.
    }

    template<class T> void ConvolveZH7WithZH7Norm(const float brdfCoeffs[7], const T zhCoeffsIn[7], T zhCoeffsOut[7])
    {
        zhCoeffsOut[0] = zhCoeffsIn[0];

        for (int i = 1; i < 7; i++)
        {
            float invAlpha = sqrtf(2.0f * i + 1);

            zhCoeffsOut[i] = zhCoeffsIn[i] * (brdfCoeffs[i] / (invAlpha * brdfCoeffs[0]));
        }
    }

    template<class T> void AddZH7Sample(float z, T c, T zhCoeffs[7])
    {
        float z2 = z * z;
        float z3 = z2 * z;
        float z4 = z2 * z2;
        float z5 = z2 * z3;
        float z6 = z3 * z3;

        zhCoeffs[0] += c * kZH_Y_0;
        zhCoeffs[1] += c * kZH_Y_1 * z;
        zhCoeffs[2] += c * kZH_Y_2 * (3 * z2 - 1);
        zhCoeffs[3] += c * kZH_Y_3 * (5 * z3 - 3 * z);
        zhCoeffs[4] += c * kZH_Y_4 * (35 * z4 - 30 * z2 + 3);
        zhCoeffs[5] += c * kZH_Y_5 * (63 * z5 - 70 * z3 + 15 * z);
        zhCoeffs[6] += c * kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);
    }

    template<class T> T EvalZH7(float z, const T zhCoeffs[7])
    {
        float z2 = z * z;
        float z3 = z2 * z;
        float z4 = z2 * z2;
        float z5 = z2 * z3;
        float z6 = z3 * z3;

        T c;

        c  = zhCoeffs[0] * kZH_Y_0;
        c += zhCoeffs[1] * kZH_Y_1 * z;
        c += zhCoeffs[2] * kZH_Y_2 * (3 * z * z - 1);
        c += zhCoeffs[3] * kZH_Y_3 * (5 * z3 - 3 * z);
        c += zhCoeffs[4] * kZH_Y_4 * (35 * z4 - 30 * z2 + 3);
        c += zhCoeffs[5] * kZH_Y_5 * (63 * z5 - 70 * z3 + 15 * z);
        c += zhCoeffs[6] * kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);

        return c;
    }

    // Windowing a la Peter-Pike Sloan
    inline float WindowScale(int n, float gamma)
    {
        float nt = float(n * (n + 1));
        return 1.0f / (1.0f + gamma * nt * nt);
    }

    template<class T> void ApplyZH7Windowing(float gamma, T coeffs[7])
    {
        for (int i = 0; i < 7; i++)
            coeffs[i] *= WindowScale(i, gamma);
    }
    // End

    inline Vec3f Bias_xyY(Vec3f c)  // effectively make delta lum proportional to real lum, and scale xy by lum
    {
        c.z += 1.0f;
        c.x *= c.z;
        c.y *= c.z;
        return c;
    }

    inline Vec3f Unbias_xyY(Vec3f c)  // Return to delta xyY form
    {
        float y = c.z;
        if (y < 1e-2f)
            y = 1e-2f;

        c.x /= y;
        c.y /= y;
        c.z -= 1.0f;
        return c;
    }

    template<class T> void FindZH7FromThetaTable(int tableSize, const T table[], T zhCoeffs[7])
    {
        float dt = 1.0f / (tableSize - 1);
        float t = 0.0f;

        float w = vl_twoPi * 2 * dt; // 2pi dz = 2pi 2 dt

        for (int i = 0; i < 7; i++)
            zhCoeffs[i] = vl_0;

        for (int i = 0; i < tableSize; i++)
        {
            AddZH7Sample(2 * t - 1, table[i] * w, zhCoeffs);
            t += dt;
        }
    }

    template<class T> void GenerateThetaTableFromZH7(const T zhCoeffs[7], int tableSize, T table[])
    {
        float dt = 1.0f / (tableSize - 1);
        float t = 0.0f;

        for (int i = 0; i < tableSize; i++)
        {
            table[i] = EvalZH7(2 * t - 1, zhCoeffs);
            t += dt;
        }
    }

    template<class T> void FindZH7FromGammaTable(int tableSize, const T table[], T zhCoeffs[7])
    {
        float dg = 1.0f / (tableSize - 1);
        float g = 0.0f;

        float w = vl_twoPi * 4 * dg; // 2pi dz = 2pi -4 g dg

        for (int i = 0; i < 7; i++)
            zhCoeffs[i] = vl_0;

        for (int i = 0; i < tableSize; i++)
        {
            AddZH7Sample(UnmapGamma(g), table[i] * w * g, zhCoeffs);
            g += dg;
        }
    }

    template<class T> void GenerateGammaTableFromZH7(T zhCoeffs[7], int tableSize, T table[])
    {
        float dg = 1.0f / (tableSize - 1);
        float g = 0.0f;

        for (int i = 0; i < tableSize; i++)
        {
            table[i] = EvalZH7(UnmapGamma(g), zhCoeffs);
            g += dg;
        }
    }

    const float kThetaW = 0.01f;    // windowing gamma to use for theta table
    const float kGammaW = 0.002f;   // windowing gamma to use for gamma table

    const float kThetaWHosek = 0.05f;   // windowing gamma to use for Hosek theta table
    const float kGammaWHosek = 0.005f;  // windowing gamma to use for gamma table

#ifdef COMPACT_BRDF_TABLE
    const float kRowPowers[cSunSkyBRDF::kBRDFSamples - 1] = { 1.0f, 8.0f, 256.0f };
#else
    const float kRowPowers[cSunSkyBRDF::kBRDFSamples - 1] = { 1.0f, 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 96.0f };
#endif
}

void cSunSkyBRDF::FindBRDFTables(const cSunSkyTable& table, const cSunSkyPreetham&)
{
    // The BRDF tables cover the entire sphere, so we must resample theta from the Perez/Hosek tables which cover a hemisphere.
    Vec3f thetaTable[cSunSkyTable::kTableSize];
    const Vec3f* gammaTable = table.mGammaTable;

    // Fill top hemisphere of table
    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[kTableSize / 2 + i] = table.mThetaTable[2 * i];

    // Fill lower hemisphere with term that evaluates close to 0, to avoid below-ground luminance leaking in.
    // TODO: modify this using ground albedo
    Vec3f lowerHemi = Vec3f(0.0f, 0.0f, 0.999f);  // 0.999 as theta table is used as (1 - F)

    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[i] = lowerHemi;

    // Project tables into ZH coefficients
    Vec3f zhCoeffsTheta[7];
    Vec3f zhCoeffsGamma[7];

    // theta table works better if we operate on something proportional to real luminance
    Vec3f biasedThetaTable[kTableSize];
    for (int i = 0; i < kTableSize; i++)
        biasedThetaTable[i] = Bias_xyY(thetaTable[i]);

    FindZH7FromThetaTable(kTableSize, biasedThetaTable, zhCoeffsTheta);
    FindZH7FromGammaTable(kTableSize,       gammaTable, zhCoeffsGamma);

    // row 0 is the original unconvolved signal

    // Firstly, fill -ve z with reflected +ve z, ramped to 0 at z = -1. This avoids discontinuities at the horizon.
    for (int i = 0; i < kTableSize / 2; i++)
    {
        thetaTable[i] = table.mThetaTable[kTableSize - 1 - 2 * i];

        // Ramp luminance down
        float s = (i + 0.5f) / (kTableSize / 2);
        s = sqrtf(s);

        thetaTable[i].z = thetaTable[i].z * s + (1 - s);
    }

    for (int i = 0; i < kTableSize; i++)
    {
        mBRDFThetaTable[0][i] = thetaTable[i];
        mBRDFGammaTable[0][i] = gammaTable[i];
    }

    // rows 1..n-1 are successive convolutions
    Vec3f zhCoeffsThetaConv[7];
    Vec3f zhCoeffsGammaConv[7];

    for (int r = 1; r < kBRDFSamples; r++)
    {
        int rs = kBRDFSamples - r - 1;
        float s = kRowPowers[rs];

        float csCoeffs[7];
        CalcCosPowerSatZH7(s, csCoeffs);

        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsTheta, zhCoeffsThetaConv);
        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsGamma, zhCoeffsGammaConv);

        // Scale up to full windowing at full spec power...
        float rw = sqrtf(rs / float(kBRDFSamples - 2));

        ApplyZH7Windowing(kThetaW * rw, zhCoeffsThetaConv);
        ApplyZH7Windowing(kGammaW * rw, zhCoeffsGammaConv);

        GenerateThetaTableFromZH7(zhCoeffsThetaConv, kTableSize, mBRDFThetaTable[r]);
        GenerateGammaTableFromZH7(zhCoeffsGammaConv, kTableSize, mBRDFGammaTable[r]);

        for (int i = 0; i < kTableSize; i++)
            mBRDFThetaTable[r][i] = Unbias_xyY(mBRDFThetaTable[r][i]);
    }

    mMaxTheta = table.mMaxTheta;
    mMaxGamma = table.mMaxGamma;
    mXYZ = false;
    mHasHTerm = false;

#ifdef SIM_CLAMP
    for (int r = 1; r < kBRDFSamples; r++)
        for (int i = 0; i < kTableSize; i++)
        {
            mBRDFThetaTable  [r][i] = ClampUnit(mBRDFThetaTable  [r][i]);
            mBRDFGammaTable  [r][i] = ClampUnit(mBRDFGammaTable  [r][i]);
            mBRDFThetaTableH [r][i] = ClampUnit(mBRDFThetaTableH [r][i]);
            mBRDFThetaTableFH[r][i] = ClampUnit(mBRDFThetaTableFH[r][i]);
        }
#endif
}

void cSunSkyBRDF::FindBRDFTables(const cSunSkyTable& table, const cSunSkyHosek& )
{
    // The BRDF tables cover the entire sphere, so we must resample theta from the Perez/Hosek tables which cover a hemisphere.
    Vec3f thetaTable[cSunSkyTable::kTableSize];
    const Vec3f* gammaTable = table.mGammaTable;

    // Fill top hemisphere of table
    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[kTableSize / 2 + i] = table.mThetaTable[2 * i];

    // Fill lower hemisphere with term that evaluates close to 0, to avoid below-ground luminance leaking in.
    // TODO: modify this using ground albedo
    Vec3f lowerHemi = Vec3f(0.999f);  // 0.999 as theta table is used as (1 - F)

    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[i] = lowerHemi;

    // Project tables into ZH coefficients
    Vec3f zhCoeffsTheta[7];
    Vec3f zhCoeffsGamma[7];

    // theta table works better if we operate on something proportional to real luminance
    Vec3f biasedThetaTable[kTableSize];

    for (int i = 0; i < kTableSize; i++)
        biasedThetaTable[i] = thetaTable[i] + vl_one;

    FindZH7FromThetaTable(kTableSize, biasedThetaTable, zhCoeffsTheta);
    FindZH7FromGammaTable(kTableSize,       gammaTable, zhCoeffsGamma);

    // row 0 is the original unconvolved signal

    // Firstly, fill -ve z with reflected +ve z, ramped to 0 at z = -1. This avoids discontinuities at the horizon.
    for (int i = 0; i < kTableSize / 2; i++)
    {
        thetaTable[i] = table.mThetaTable[kTableSize - 1 - 2 * i];

        // Ramp luminance down
        float s = (i + 0.5f) / (kTableSize / 2);
        s = sqrtf(s);

        thetaTable[i] = thetaTable[i] * s + Vec3f(1 - s);
    }

    for (int i = 0; i < kTableSize; i++)
    {
        mBRDFThetaTable[0][i] = thetaTable[i];
        mBRDFGammaTable[0][i] = gammaTable[i];
    }

    // Construct H term table -- just the zenith part as we can potentially store as 4th component
    for (int i = 0; i < kTableSize / 2; i++)
        mBRDFThetaTableH[0][i] = vl_0;

    for (int i = 0; i < kTableSize / 2; i++)
    {
        float cosTheta = i / float(kTableSize / 2 - 1); // XXX
        float zenith = sqrtf(cosTheta);

        mBRDFThetaTableH[0][kTableSize / 2 + i] = zenith;
    }

    // Calculate FH as we get slightly better results convolving the full term
    // rather than approximating, at the cost of an extra table
    for (int i = 0; i < kTableSize; i++)
        mBRDFThetaTableFH[0][i] = mBRDFThetaTableH[0][i] * thetaTable[i];

    float zhCoeffsH[7];
    FindZH7FromThetaTable(kTableSize, mBRDFThetaTableH [0], zhCoeffsH);

    Vec3f zhCoeffsFH[7];
    FindZH7FromThetaTable(kTableSize, mBRDFThetaTableFH[0], zhCoeffsFH);

    // rows 1..n-1 are successive convolutions

    for (int r = 1; r < kBRDFSamples; r++)
    {
        int rs = kBRDFSamples - r - 1;
        float s = kRowPowers[rs];

        float csCoeffs[7];
        CalcCosPowerSatZH7(s, csCoeffs);

        Vec3f zhCoeffsThetaConv[7];
        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsTheta, zhCoeffsThetaConv);
        Vec3f zhCoeffsGammaConv[7];
        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsGamma, zhCoeffsGammaConv);

        float zhCoeffsThetaConvH [7];
        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsH , zhCoeffsThetaConvH);
        Vec3f zhCoeffsThetaConvFH[7];
        ConvolveZH7WithZH7Norm(csCoeffs, zhCoeffsFH, zhCoeffsThetaConvFH);

        // Scale up to full windowing at full spec power...
        float rw = sqrtf(rs / float(kBRDFSamples - 2));

        ApplyZH7Windowing(kThetaWHosek * rw, zhCoeffsThetaConv);
        ApplyZH7Windowing(kGammaWHosek * rw, zhCoeffsGammaConv);
        ApplyZH7Windowing(kThetaWHosek * rw, zhCoeffsThetaConvH);
        ApplyZH7Windowing(kThetaWHosek * rw, zhCoeffsThetaConvFH);

        // Generate convolved tables from ZH
        GenerateThetaTableFromZH7(zhCoeffsThetaConv,   kTableSize, mBRDFThetaTable  [r]);
        GenerateGammaTableFromZH7(zhCoeffsGammaConv,   kTableSize, mBRDFGammaTable  [r]);
        GenerateThetaTableFromZH7(zhCoeffsThetaConvH,  kTableSize, mBRDFThetaTableH [r]);
        GenerateThetaTableFromZH7(zhCoeffsThetaConvFH, kTableSize, mBRDFThetaTableFH[r]);

        for (int i = 0; i < kTableSize; i++)
        {
            mBRDFThetaTable[r][i] -= vl_one;

        #ifdef HOSEK_G_FIX
            // Ringing on the Hosek G term leads to blue spots opposite the sun
            // in sunset situatons. Trying to solve this completely via windowing
            // leads to excessive blurring, so we also compensate by scaling down
            // the far pole in this situation.
            // TODO: investigate further, not sure this is worth it.
            float g = i / float(kTableSize - 1);
            float cosGamma = UnmapGamma(g);
            if (cosGamma < -0.6f)
                mBRDFGammaTable[r][i] *= ClampUnit(1.0f - sqr(-0.6f - cosGamma) * 1.5f * rw);
        #endif
        }
    }

    mMaxTheta = table.mMaxTheta;
    mMaxGamma = table.mMaxGamma;
    mXYZ = true;
    mHasHTerm = true;

#ifdef SIM_CLAMP
    for (int r = 1; r < kBRDFSamples; r++)
        for (int i = 0; i < kTableSize; i++)
        {
            mBRDFThetaTable  [r][i] = ClampUnit(mBRDFThetaTable  [r][i]);
            mBRDFGammaTable  [r][i] = ClampUnit(mBRDFGammaTable  [r][i]);
            mBRDFThetaTableH [r][i] = ClampUnit(mBRDFThetaTableH [r][i]);
            mBRDFThetaTableFH[r][i] = ClampUnit(mBRDFThetaTableFH[r][i]);
        }
#endif
}

Vec3f cSunSkyBRDF::ConvolvedSkyRGB(const cSunSkyPreetham& pt, const Vec3f& v, float r) const
{
    CL_ASSERT(!mXYZ);

    float cosTheta = v.z;
    float cosGamma = dot(pt.mToSun, v);

    float t = 0.5f * (cosTheta + 1);
    float g = MapGamma(cosGamma);

    Vec3f F = BiLerp(t, r, kTableSize, kBRDFSamples, mBRDFThetaTable[0]);
    Vec3f G = BiLerp(g, r, kTableSize, kBRDFSamples, mBRDFGammaTable[0]);

#ifdef SIM_CLAMP
    F.z *= mMaxTheta;
    G.z *= mMaxGamma;
#endif

    // (1 - F(theta)) * (1 + G(phi))
    Vec3f xyY = (Vec3f(vl_1) - F) * (Vec3f(vl_1) + G);

    xyY *= pt.mPerezInvDen;

    return xyYToRGB(xyY);
}

Vec3f cSunSkyBRDF::ConvolvedSkyRGB(const cSunSkyHosek& hk, const Vec3f& v, float r) const
{
    CL_ASSERT(mXYZ);

    float cosTheta = v.z;
    float cosGamma = dot(hk.mToSun, v);

    float t = 0.5f * (cosTheta + 1);
    float g = MapGamma(cosGamma);

    Vec3f F = BiLerp(t, r, kTableSize, kBRDFSamples, mBRDFThetaTable[0]);
    Vec3f G = BiLerp(g, r, kTableSize, kBRDFSamples, mBRDFGammaTable[0]);

#ifdef SIM_CLAMP
    F *= mMaxTheta;
    G *= mMaxGamma;
#endif

#ifdef HOSEK_BRDF_ANALYTIC_H
    // Evaluate H term on the fly, for cross-checking.
    // This could be done in shader to avoid the additional table lookup,
    // but it's a fair bit of math, and using the table variant ensures better
    // consistency with the F/G parts.
    float zhZ[7];
    CalcCosPowerSatZH7(0.5f, zhZ);

    float zhR[7];
    CalcCosPowerSatZH7(1.0f / (r + 0.001f), zhR);

    float zhZR[7];
    ConvolveZH7WithZH7Norm(zhR, zhZ, zhZR);

    float zenith = EvalZH7(cosTheta, zhZR);

    Vec3f H
    (
        hk.mCoeffsXYZ[0][7] * zenith + (hk.mCoeffsXYZ[0][2] - 1.0f),
        hk.mCoeffsXYZ[1][7] * zenith + (hk.mCoeffsXYZ[1][2] - 1.0f),
        hk.mCoeffsXYZ[2][7] * zenith + (hk.mCoeffsXYZ[2][2] - 1.0f)
    );
    Vec3f FH = F * H;
#else
    Vec3f cH(hk.mCoeffsXYZ[0][7], hk.mCoeffsXYZ[1][7], hk.mCoeffsXYZ[2][7]);
    Vec3f cI(hk.mCoeffsXYZ[0][2] - 1.0f, hk.mCoeffsXYZ[1][2] - 1.0f, hk.mCoeffsXYZ[2][2] - 1.0f);

    Vec3f H  = BiLerp(t, r, kTableSize, kBRDFSamples, mBRDFThetaTableH [0]) * cH;
    Vec3f FH = BiLerp(t, r, kTableSize, kBRDFSamples, mBRDFThetaTableFH[0]) * cH;

#ifdef SIM_CLAMP
    FH *= mMaxTheta;
#endif

    H  +=     cI;
    FH += F * cI;
#endif

    // (1 - F(theta)) * (1 + G(phi) + H(theta))
    Vec3f XYZ = (Vec3f(vl_1) - F) * (Vec3f(vl_1) + G) + H - FH;

    XYZ = ClampPositive(XYZ);
    XYZ *= hk.mRadXYZ;

    return XYZToRGB(XYZ);
}

void cSunSkyBRDF::FillBRDFTexture(int width, int height, uint8_t image[][4]) const
{
    CL_ASSERT(width == kTableSize);
    CL_ASSERT((height == 2 * kBRDFSamples) || (mHasHTerm && height == 4 * kBRDFSamples));

    (void) width;
    (void) height;

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
        {
            Vec3f c = mBRDFThetaTable[j][i];

            if (mXYZ)
                c /= mMaxTheta;
            else
                c.z /= mMaxTheta;

            image[i][2] = ToU8(c.x);
            image[i][1] = ToU8(c.y);
            image[i][0] = ToU8(c.z);

            if (mHasHTerm)
                image[i][3] = ToU8(mBRDFThetaTableH[j][i]);
            else
                image[i][3] = 255;
        }

        image += width;
    }

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
        {
            Vec3f c = mBRDFGammaTable[j][i];

            if (mXYZ)
                c /= mMaxGamma;
            else
                c.z /= mMaxGamma;

            image[i][2] = ToU8(c.x);
            image[i][1] = ToU8(c.y);
            image[i][0] = ToU8(c.z);
            image[i][3] = 255;
        }

        image += width;
    }

    if (height <= 2)
        return;

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
        {
            Vec3f c = mBRDFThetaTableFH[j][i];

            c /= mMaxTheta;

            image[i][2] = ToU8(c.x);
            image[i][1] = ToU8(c.y);
            image[i][0] = ToU8(c.z);
            image[i][3] = ToU8(mBRDFThetaTableH[j][i]);
        }

        image += width;
    }
}

void cSunSkyBRDF::FillBRDFTexture(int width, int height, float image[][4]) const
{
    CL_ASSERT(width == kTableSize);
    CL_ASSERT(height == 2 * kBRDFSamples);

    (void) width;
    (void) height;

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
            *(Vec4f*) image[i] = Vec4f(mBRDFThetaTable[j][i], mHasHTerm ? mBRDFThetaTableH[j][i] : 1.0f);

        image += width;
    }

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
            *(Vec4f*) image[i] = Vec4f(mBRDFGammaTable[j][i], 1.0f);

        image += width;
    }

    if (height <= 2)
        return;

    for (int j = 0; j < kBRDFSamples; j++)
    {
        for (int i = 0; i < width; i++)
            *(Vec4f*) image[i] = Vec4f(mBRDFThetaTableFH[j][i], mBRDFThetaTableH[j][i]);

        image += width;
    }
}


//------------------------------------------------------------------------------
// cSunSky -- composite class for easier comparison
//------------------------------------------------------------------------------

cSunSky::cSunSky() :
    mSkyType(kPreetham),
    mToSun(vl_0),
    mTurbidity(2.5f),
    mAlbedo(vl_0),
    mOvercast(0.0f),
    mRoughness(0.0f),
    mZenithY(0.0f)
{
}

void cSunSky::SetSkyType(tSkyType skyType)
{
    mSkyType = skyType;
}

tSkyType cSunSky::SkyType() const
{
    return mSkyType;
}

void cSunSky::SetSunDir(const Vec3f& v)
{
    mToSun = v;
}

void cSunSky::SetTurbidity(float turbidity)
{
    mTurbidity = turbidity;
}

void cSunSky::SetAlbedo(Vec3f rgb)
{
    mAlbedo = rgb;
}

void cSunSky::SetOvercast(float overcast)
{
    mOvercast = overcast;
}

void cSunSky::SetRoughness(float roughness)
{
    mRoughness = roughness;
}

void cSunSky::Update()
{
    mZenithY = ZenithLuminance(acosf(mToSun.z), mTurbidity);
    
    mPreetham.Update(mToSun, mTurbidity, mOvercast);
    mHosek   .Update(mToSun, mTurbidity, mAlbedo, mOvercast);

    if (mSkyType == kPreethamTable || mSkyType == kPreethamBRDF)
        mTable.FindThetaGammaTables(mPreetham);
    else if (mSkyType == kHosekTable || mSkyType == kHosekBRDF)
        mTable.FindThetaGammaTables(mHosek);

    if (mSkyType == kPreethamBRDF)
        mBRDF.FindBRDFTables(mTable, mPreetham);
    if (mSkyType == kHosekBRDF)
        mBRDF.FindBRDFTables(mTable, mHosek);
}

Vec3f cSunSky::SkyRGB(const Vec3f& v) const
{
    switch (mSkyType)
    {
    case kPreetham:
        return mPreetham.SkyRGB(v);

    case kPreethamTable:
        return mTable.SkyRGB(mPreetham, v);
    case kPreethamBRDF:
        return mBRDF.ConvolvedSkyRGB(mPreetham, v, mRoughness);

    case kHosek:
        return mHosek.SkyRGB(v);
    case kHosekTable:
        return mTable.SkyRGB(mHosek, v);
    case kHosekBRDF:
        return mBRDF.ConvolvedSkyRGB(mHosek, v, mRoughness);

    case kCIEClear:
        return Vec3f(CIEClearSkyLuminance       (v, mToSun, mZenithY));
    case kCIEOvercast:
        return Vec3f(CIEOvercastSkyLuminance    (v,         mZenithY));
    case kCIEPartlyCloudy:
        return Vec3f(CIEPartlyCloudySkyLuminance(v, mToSun, mZenithY));

    default:
        return vl_0;
    }
}

float cSunSky::SkyLuminance(const Vec3f& v) const
{
    if (v.z < 0.0)
        return 0.0;

    switch (mSkyType)
    {
    case kPreetham:
        return mPreetham.SkyLuminance(v);
    case kCIEClear:
        return CIEClearSkyLuminance       (v, mToSun, mZenithY);
    case kCIEOvercast:
        return CIEOvercastSkyLuminance    (v,         mZenithY);
    case kCIEPartlyCloudy:
        return CIEPartlyCloudySkyLuminance(v, mToSun, mZenithY);
    case kHosek:
        return mHosek.SkyLuminance(v);
    default:
        return 0;
    }
}

Vec2f cSunSky::SkyChroma(const Vec3f& v) const
{
    if (v.z < 0.0)
        return Vec2f(0, 0);

    switch (mSkyType)
    {
    case kPreetham:
        return mPreetham.SkyChroma(v);
    case kCIEClear:
        return kClearChroma;
    case kCIEOvercast:
        return kOvercastChroma;
    case kCIEPartlyCloudy:
        return kPartlyCloudyChroma;
    case kHosek:
        {
            Vec3f xyz = mHosek.SkyXYZ(v);
            return Vec2f(xyz.x / xyz.z, xyz.y / xyz.z);
        }
    default:
        return kOvercastChroma;
    }
}

float cSunSky::AverageLuminance() const
{
    switch (mSkyType)
    {
    case kPreetham:
    case kPreethamTable:
    case kPreethamBRDF:
        return mPreetham.mPerezInvDen.z;

    case kCIEClear:
        return CIEClearSkyLuminance(vl_0, mToSun, mZenithY);
    case kCIEOvercast:
        return mZenithY;
    case kCIEPartlyCloudy:
        return CIEPartlyCloudySkyLuminance(Vec3f(0.0f, 0.0f, 0.0f), mToSun, mZenithY);

    case kHosek:
    case kHosekTable:
    case kHosekBRDF:
        return mHosek.mRadXYZ.y;

    default:
        return 1.0f;
    }
}

