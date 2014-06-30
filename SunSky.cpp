/*
    File:       SunSky.cpp
 
    Function:   Implements a sky radiance function for a clear sky,
                and for various overcast/clear skies, using the Utah sun sky
                model and the CIE standard sky models respectively.

    Author:     Andrew Willmott, originally based on code by Brian Smits.
 
    Notes:      
*/

#include "SunSky.h"

namespace
{
    // XYZ -> Gamut conversion for a modern-day Dell monitor
    const Vec3f kXYZToR( 2.80298,   -1.18735,  -0.437286);
    const Vec3f kXYZToG(-1.07909,    1.97927,   0.0423073);
    const Vec3f kXYZToB( 0.0746015, -0.248513,  1.08196);

    const Vec2f kClearChroma       (2.0 / 3.0, 1.0 / 3.0);
    const Vec2f kOvercastChroma    (1.0 / 3.0, 1.0 / 3.0);
    const Vec2f kPartlyCloudyChroma(1.0 / 3.0, 1.0 / 3.0);

    float kSunRadiance[16][16][3] =
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

    inline float DegreesToRadians(float d)
    {
        return d * (vl_twoPi / 360.0f);
    }

    inline float PerezFunctionFromDensity
    (
        const float*     lambdas,
        float            invDenLum,
        float            cosTheta,
        float            gamma,
        float            cosGamma
    )
    {
        float num =  (1.0 + lambdas[0] * exp(lambdas[1] / cosTheta))
                  * (1.0 + lambdas[2] * exp(lambdas[3] * gamma) + lambdas[4] * sqr(cosGamma));

        return num * invDenLum;
    }

    inline Vec3f xyYToXYZ(const Vec3f& c)
    {
        return Vec3f(c[0], c[1], 1.0 - c[0] - c[1]) * c[2] / c[1];
    }

    inline Vec3f xyYToRGB(const Vec3f& xyY)
    {
        Vec3f cieXYZ = xyYToXYZ(xyY);
        Vec3f rgb;

        rgb[0] = dot(kXYZToR, cieXYZ);
        rgb[1] = dot(kXYZToG, cieXYZ);
        rgb[2] = dot(kXYZToB, cieXYZ);

        return rgb;
    }
}

//
// NOTES
//  All times in decimal form (6.25 = 6:15 AM)
//  All angles in radians apart from longitude/latitude, theta angles measured
//  from vertical (zenith) rather than from horizon.
//

cSunSky::cSunSky() :
    mInited(false),
    mSkyType(kUtahClearSky),

    mLatitude(0.0f),
    mLongitude(0.0f),
    mJulianDay(180.0f),
    mTimeOfDay(12.0f),
    mStandardMeridian(0.0f),
    mTurbidity(4.0f),

    mToSun(vl_z),
    mThetaS(0),
    mPhiS(0),
    mSunSolidAngle(0),

    mZenith(vl_0),

    // mPerez_x[5] = 0
    // mPerez_y[5] = 0
    // mPerez_Y[5] = 0

    mPerezInvDen(vl_1),

    // mThetaTable[64];
    // mGammaTable[64];

    mRoughness(0.0f)
//    mBRDFThetaTable[64][4] = {}
//    mBRDFGammaTable[64][4] = {}

{
}

void cSunSky::SetAllParams
(
    float       timeOfDay,          
    float       timeZone,
    int         julianDay,          
    float       latitude,
    float       longitude,
    float       turbidity           
)
{
    mLatitude = latitude;
    mLongitude = longitude;
    mJulianDay = julianDay;
    mTimeOfDay = timeOfDay;
    mStandardMeridian = timeZone * 15.0;
    mTurbidity = turbidity;

    InitTurbidityConstants();
}

void cSunSky::SetTime(float   timeOfDay)
{
    mTimeOfDay = timeOfDay;
    InitTimeConstants();
}

void cSunSky::InitTimeConstants()
// init time-dependent constants
{
    FindSunDirection();
    
//    mToSun = Vec3f(cos(mPhiS) * sin(mThetaS), sin(mPhiS) * sin(mThetaS), cos(mThetaS));

    // Units: W cm^-2 Sr^-1.
    // Sun diameter: 1.392 10^9 m, avg. distance: 149.6 10^9
    mSunSolidAngle = 0.25 * vl_pi * sqr(1.392) / sqr(149.6);   // = 6.7443e-05

    float   theta2 = sqr(mThetaS);
    float   theta3 = theta2 * mThetaS;
    float   T = mTurbidity;
    float   T2 = sqr(T);

    float   chi = (4.0 / 9.0 - T / 120.0) * (vl_pi - 2.0 * mThetaS);

    // mZenith stored as xyY
    mZenith[2] = (4.0453 * T - 4.9710) * tan(chi) - .2155 * T + 2.4192;
    mZenith[2] *= 1000.0;   // conversion from kcd/m^2 to cd/m^2

    mZenith[0] =
        ( 0.00165 * theta3 - 0.00374 * theta2 + 0.00208 * mThetaS + 0.0)     * T2 +
        (-0.02902 * theta3 + 0.06377 * theta2 - 0.03202 * mThetaS + 0.00394) * T  +
        ( 0.11693 * theta3 - 0.21196 * theta2 + 0.06052 * mThetaS + 0.25885);

    mZenith[1] =
        ( 0.00275 * theta3 - 0.00610 * theta2 + 0.00316 * mThetaS + 0.0)     * T2 +
        (-0.04214 * theta3 + 0.08970 * theta2 - 0.04153 * mThetaS + 0.00515) * T  +
        ( 0.15346 * theta3 - 0.26756 * theta2 + 0.06669 * mThetaS + 0.26688);

    // initialize sun-constant parts of the Perez functions
    mPerezInvDen[0] = mZenith[0] / ((1.0 + mPerez_x[0] * exp(mPerez_x[1])) *
       (1.0 + mPerez_x[2] * exp(mPerez_x[3] * mThetaS) + mPerez_x[4] * sqr(cos(mThetaS))));
    mPerezInvDen[1] = mZenith[1] / ((1.0 + mPerez_y[0] * exp(mPerez_y[1])) *
       (1.0 + mPerez_y[2] * exp(mPerez_y[3] * mThetaS) + mPerez_y[4] * sqr(cos(mThetaS))));
    mPerezInvDen[2] = mZenith[2] / ((1.0 + mPerez_Y[0] * exp(mPerez_Y[1])) *
       (1.0 + mPerez_Y[2] * exp(mPerez_Y[3] * mThetaS) + mPerez_Y[4] * sqr(cos(mThetaS))));

    FindThetaGammaTables();
    FindBRDFTables();
}

void cSunSky::InitTurbidityConstants()
// initialize turbidity-dependent constants
{
    mInited = true;

    float   T = mTurbidity;

    mPerez_Y[0] =  0.17872 * T - 1.46303;
    mPerez_Y[1] = -0.35540 * T + 0.42749;
    mPerez_Y[2] = -0.02266 * T + 5.32505;
    mPerez_Y[3] =  0.12064 * T - 2.57705;
    mPerez_Y[4] = -0.06696 * T + 0.37027;

    mPerez_x[0] = -0.01925 * T - 0.25922;
    mPerez_x[1] = -0.06651 * T + 0.00081;
    mPerez_x[2] = -0.00041 * T + 0.21247;
    mPerez_x[3] = -0.06409 * T - 0.89887;
    mPerez_x[4] = -0.00325 * T + 0.04517;

    mPerez_y[0] = -0.01669 * T - 0.26078;
    mPerez_y[1] = -0.09495 * T + 0.00921;
    mPerez_y[2] = -0.00792 * T + 0.21023;
    mPerez_y[3] = -0.04405 * T - 1.65369;
    mPerez_y[4] = -0.01092 * T + 0.05291;

    InitTimeConstants();    
}

void cSunSky::SunThetaPhi(float* stheta, float* sphi) const
{
    *sphi = mPhiS;
    *stheta = mThetaS;
}

void cSunSky::FindSunDirection()
{
    float solarTime = mTimeOfDay +
         (0.170 * sin(4 * vl_pi * (mJulianDay - 80) / 373) 
        - 0.129 * sin(2 * vl_pi * (mJulianDay - 8) / 355)) +
        (mLongitude - mStandardMeridian) / 15.0;

    float solarDeclination = (0.4093 * sin(2 * vl_pi * (mJulianDay - 81) / 368));

    float latRads = DegreesToRadians(mLatitude);

    float solarAltitude = 
        asin(sin(latRads) * sin(solarDeclination) -
             cos(latRads) * cos(solarDeclination) * cos(vl_pi * solarTime / 12));


    float sx = cos(solarDeclination) * sin(vl_pi * solarTime / 12);
    float sy = cos(latRads) * sin(solarDeclination) +
               sin(latRads) * cos(solarDeclination) * cos(vl_pi * solarTime / 12);

    float solarAzimuth = atan2(sy, sx);

    mPhiS = solarAzimuth;
    mThetaS = vl_halfPi - solarAltitude;

    mToSun = Vec3f(sx, sy, cos(mThetaS));
}

Vec3f cSunSky::SunDirection() const
{
    return mToSun;
}

float cSunSky::SunSolidAngle() const
{
    return mSunSolidAngle;
}

//
// The Perez function is as follows:
// 
//    P(t, g) =  (1 + l1 e ^ (l2 / cos(t))) (1 + l3 e ^ (l4 g)  + l5 cos(g)  ^ 2)
//               --------------------------------------------------------
//               (1 + l1 e ^ l2)            (1 + l3 e ^ (l4 ts) + l5 cos(ts) ^ 2)
//    

float cSunSky::PerezFunction(const float* lambdas, float theta, float gamma, float lvz) const
{
    float num =  (1.0 + lambdas[0] * exp(lambdas[1] / cos(theta))) 
               * (1.0 + lambdas[2] * exp(lambdas[3] * gamma) + lambdas[4] * sqr(cos(gamma)));

    float den =  (1.0 + lambdas[0] * exp(lambdas[1]))
               * (1.0 + lambdas[2] * exp(lambdas[3] * mThetaS) + lambdas[4] * sqr(cos(mThetaS)));

    return lvz * num / den;
}

Vec3f cSunSky::SkyRGB(const Vec3f& v) const
{
    switch (mSkyType)
    {
    case kUtahClearSky:
        return ClearSkyRGB(v);

    case kUtahClearSkyTable:
        return SkyRGBFromTable(v);
    case kUtahClearSkyBRDF:
        return ConvolvedSkyRGB(v, mRoughness);

    case kCIEClearSky:
        return CIEClearSkyLuminance(v)        * Vec3f(0.0f, 0.0f, 0.3f);
    case kCIECloudySky:
        return CIEOvercastSkyLuminance(v)     * Vec3f(0.3f);
    case kCIEPartlyCloudySky:
        return CIEPartlyCloudySkyLuminance(v) * Vec3f(0.3f);
    default:
        return Vec3f(1.0f, 0.0f, 0.0f);
    }
}

float cSunSky::SkyLuminance(const Vec3f& v) const
{
    if (v[2] < 0.0)
        return 0.0;

    switch (mSkyType)
    {
    case kUtahClearSky:
        return ClearSkyLuminance(v);
    case kCIEClearSky:
        return CIEClearSkyLuminance(v);
    case kCIECloudySky:
        return CIEOvercastSkyLuminance(v);
    case kCIEPartlyCloudySky:
        return CIEPartlyCloudySkyLuminance(v);
    default:
        return 0;
    }
}

Vec2f cSunSky::SkyChroma(const Vec3f& v) const
{
    if (v[2] < 0.0)
        return Vec2f(0, 0);

    switch (mSkyType)
    {
    case kUtahClearSky:
        return ClearSkyChroma(v);
    case kCIEClearSky:
        return kClearChroma;
    case kCIECloudySky:
        return kOvercastChroma;
    case kCIEPartlyCloudySky:
        return kPartlyCloudyChroma;
    default:
        return Vec2f(0, 0);
    }
}

Vec3f cSunSky::SunRGB() const
{
    // Look up table data from original sun.pic
    float cosTheta = mToSun[2];

    if (cosTheta < 0.0f)
        return vl_0;

    const float kTableSizeMinusOne = sizeof(kSunRadiance) / sizeof(kSunRadiance[0]);

    float s = kTableSizeMinusOne * cosTheta;
    float t = kTableSizeMinusOne * (mTurbidity - 2.0f) / 10.0f; // useful range is 2-12

    int is = floor(s);
    int it = floor(t);

    s -= is;
    t -= it;

    // find bilinear weights
    float b00 = (1 - s) * (1 - t);
    float b01 = s       * (1 - t);
    float b10 = (1 - s) * t;
    float b11 = s       * t;

    const float* t00 = kSunRadiance[it + 0][is + 0];
    const float* t01 = kSunRadiance[it + 0][is + 1];
    const float* t10 = kSunRadiance[it + 1][is + 0];
    const float* t11 = kSunRadiance[it + 1][is + 1];

    Vec3f c00(t00[0], t00[1], t00[2]);
    Vec3f c01(t01[0], t01[1], t01[2]);
    Vec3f c10(t10[0], t10[1], t10[2]);
    Vec3f c11(t11[0], t11[1], t11[2]);

    Vec3f sun;
    sun  = c00 * b00;
    sun += c01 * b01;
    sun += c10 * b10;
    sun += c11 * b11;

    // 683 converts from watts to candela at 540e12 hz. Really, we should be weighting by luminous efficiency curve rather than CIE_Y
    sun *=  683.0f;

    return sun;
}


Vec3f cSunSky::ClearSkyRGB(const Vec3f& v) const
{
    float cosTheta = v[2];
    float cosGamma = dot(mToSun, v);
    float gamma = acos(cosGamma);

    if (cosTheta < 0.0f)
        cosTheta = 0.0f;

    Vec3f xyY;
    xyY[0] = PerezFunctionFromDensity(mPerez_x, mPerezInvDen[0], cosTheta, gamma, cosGamma);
    xyY[1] = PerezFunctionFromDensity(mPerez_y, mPerezInvDen[1], cosTheta, gamma, cosGamma);
    xyY[2] = PerezFunctionFromDensity(mPerez_Y, mPerezInvDen[2], cosTheta, gamma, cosGamma);

    Vec3f cieXYZ = xyYToXYZ(xyY);
    Vec3f rgb;

    rgb[0] = dot(kXYZToR, cieXYZ);
    rgb[1] = dot(kXYZToG, cieXYZ);
    rgb[2] = dot(kXYZToB, cieXYZ);

    return rgb;
}

float cSunSky::ClearSkyLuminance(const Vec3f& v) const
{
    float cosTheta = v[2];
    float cosGamma = dot(mToSun, v);
    float gamma = acos(cosGamma);

    float luminance;

    luminance = PerezFunctionFromDensity(mPerez_Y, mPerezInvDen[2], cosTheta, gamma, cosGamma);

    return luminance;
}

Vec2f cSunSky::ClearSkyChroma(const Vec3f& v) const
{
    float cosTheta = v[2];
    float cosGamma = dot(mToSun, v);
    float gamma = acos(cosGamma);

    Vec2f chroma;

    chroma[0] = PerezFunctionFromDensity(mPerez_x, mPerezInvDen[0], cosTheta, gamma, cosGamma);
    chroma[1] = PerezFunctionFromDensity(mPerez_y, mPerezInvDen[1], cosTheta, gamma, cosGamma);

    return chroma;
}

float cSunSky::CIEOvercastSkyLuminance(const Vec3f& v) const
{
    float sinTheta = v[2];
    float zenithLuminance = mZenith[2];

    return zenithLuminance * (1.0 + 2.0 * sinTheta) / 3.0;
}

float cSunSky::CIEClearSkyLuminance(const Vec3f& v) const
{
    float zenithLuminance = mZenith[2];
    
    float sinTheta = v[2];
    float cosTheta = sqrt(1 - sqr(sinTheta));
    
    float sinS  = mToSun[2];
    float cosS2 = 1 - sqr(sinS);
    float S     = asin(sinS);
    
    float cosGamma = dot(mToSun, v);    
    float gamma    = acos(cosGamma);
    
    float top1 = 0.91 + 10 * exp(-3 * gamma) + 0.45 * sqr(cosGamma);
    float bot1 = 0.91 + 10 * exp(-3 * S)     + 0.45 * cosS2;
    
    float top2 = 1 - exp(-0.32 / cosTheta);
    float bot2 = 1 - exp(-0.32);

    return zenithLuminance * (top1 * top2) / (bot1 * bot2);
}

float cSunSky::CIEPartlyCloudySkyLuminance(const Vec3f& v) const
{
    float zenithLuminance = mZenith[2];
    
    float sinS = mToSun[2];
    float S = asin(sinS);

    float cosGamma = dot(mToSun, v);    
    float gamma = acos(cosGamma);
    
    float top1 = 0.526 + 5 * exp(-1.5 * gamma);
    float bot1 = 0.526 + 5 * exp(-1.5 * S);
    
    float top2 = 1 - exp(-0.8 / cosGamma);
    float bot2 = 1 - exp(-0.8);

    return zenithLuminance * (top1 * top2) / (bot1 * bot2);
}

void cSunSky::SetSkyType(tSkyType skyType)
{
    mSkyType = skyType;
}

tSkyType cSunSky::SkyType() const
{
    return mSkyType;
}

void cSunSky::SetRoughness(float roughness)
{
    mRoughness = roughness;
}

float cSunSky::Roughness() const
{
    return mRoughness;
}

// Pre-calculated table version, particularly useful for shader use

void cSunSky::FindThetaGammaTables()
{
    const int tableSize = sizeof(mThetaTable) / sizeof(mThetaTable[0]);

    float dt = 1.0 / (tableSize - 1);
    float cosTheta = 0.0;

    for (int i = 0; i < tableSize; i++)
    {
        mThetaTable[i][0] =  mPerez_x[0] * expf(mPerez_x[1] / cosTheta);
        mThetaTable[i][1] =  mPerez_y[0] * expf(mPerez_y[1] / cosTheta);
        mThetaTable[i][2] =  mPerez_Y[0] * expf(mPerez_Y[1] / cosTheta);

        float cosGamma = 1.0 - 2.0 * sqr(cosTheta);   // remap cosGamma to concentrate table around the sun
        float gamma = acosf(cosGamma);
        
        mGammaTable[i][0] =  mPerez_x[2] * expf(mPerez_x[3] * gamma) + mPerez_x[4] * sqr(cosGamma);
        mGammaTable[i][1] =  mPerez_y[2] * expf(mPerez_y[3] * gamma) + mPerez_y[4] * sqr(cosGamma);
        mGammaTable[i][2] =  mPerez_Y[2] * expf(mPerez_Y[3] * gamma) + mPerez_Y[4] * sqr(cosGamma);

        cosTheta += dt;
    }
}

namespace
{
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

    inline Vec3f PerezFunctionTable(float cosTheta, float cosGamma, int tableSize, const Vec3f thetaTable[], const Vec3f gammaTable[])
    {
        int ti = FloorToInt32(cosTheta * (tableSize - 1));

        if (ti < 0)
            ti = 0;

        float g = sqrtf(0.5f * (1.0f - cosGamma));  // remap to table parameterisation

        int gi = FloorToInt32(g * (tableSize - 1));

        // (1 + F(theta)) * (1 + G(phi))
        return (Vec3f(vl_1) + thetaTable[ti]) * (Vec3f(vl_1) + gammaTable[gi]);
    }
}

Vec3f cSunSky::SkyRGBFromTable(const Vec3f& v) const
{
    float cosTheta = v[2];
    float cosGamma = dot(mToSun, v);

    Vec3f perez = ::PerezFunctionTable(cosTheta, cosGamma, sizeof(mThetaTable) / sizeof(mThetaTable[0]), mThetaTable, mGammaTable);
    Vec3f xyY = perez * mPerezInvDen;

    Vec3f cieXYZ = xyYToXYZ(xyY);
    Vec3f rgb;

    rgb[0] = dot(kXYZToR, cieXYZ);
    rgb[1] = dot(kXYZToG, cieXYZ);
    rgb[2] = dot(kXYZToB, cieXYZ);

    return rgb;
}


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
    }

    void ConvolveZHWithCosPow(float n, int numZCoeffs, const Vec3f zhCoeffsIn[], Vec3f zhCoeffsOut[])
    {
        float brdfCoeffs[7];
        CalcCosPowerSatZH7(n, brdfCoeffs);
        
        for (int i = 0; i < numZCoeffs; i++)
            brdfCoeffs[i] *= sqrt(4 * vl_pi / (2 * i + 1));
        
        for (int i = 1; i < numZCoeffs; i++)
            brdfCoeffs[i] /= brdfCoeffs[0]; // normalise power

        zhCoeffsOut[0] = zhCoeffsIn[0];

        for (int i = 1; i < numZCoeffs; i++)
            zhCoeffsOut[i] = zhCoeffsIn[i] * brdfCoeffs[i];
    }

    inline void AddZHSample(float z, Vec3f c, Vec3f zhCoeffs[7])
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

    Vec3f SampleZH(float z, const Vec3f zhCoeffs[7])
    {
        float z2 = z * z;
        float z3 = z2 * z;
        float z4 = z2 * z2;
        float z5 = z2 * z3;
        float z6 = z3 * z3;

        Vec3f c;
        
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
        float nt = n * (n + 1);
        return 1.0f / (1.0f + gamma * nt * nt);
    }

    void ApplyZHWindowing(float gamma, int n, Vec3f* coeffs)
    {
        for (int i = 0; i < n; i++)
            coeffs[i] *= WindowScale(i, gamma);
    }
    // End

    void FindZHFromThetaTable(int tableSize, const Vec3f table[], int numZCoeffs, Vec3f zhCoeffs[])
    {
        float dt = 1.0 / (tableSize - 1);
        float t = 0.0;
        
        for (int i = 0; i < numZCoeffs; i++)
            zhCoeffs[i] = vl_0;
        
        for (int i = 0; i < tableSize; i++)
        {
            Vec3f wt = table[i];

            wt[2] += 1.0f;
            wt[0] *= wt[2];
            wt[1] *= wt[2];

            AddZHSample(2 * t - 1, wt, zhCoeffs);

            t += dt;
        }

        for (int i = 0; i < numZCoeffs; i++)
            zhCoeffs[i] *= (4.0f * vl_pi) / tableSize; // dw = 2pi dz
    }

    void GenerateThetaTableFromZH(int numZCoeffs, const Vec3f zhCoeffs[], int tableSize, Vec3f table[])
    {
        float dt = 1.0 / (tableSize - 1);
        float t = 0.0;

        for (int i = 0; i < tableSize; i++)
        {
            Vec3f c = SampleZH(2 * t - 1, zhCoeffs);

            float y = c[2];
            if (y < 1e-2f)
                y = 1e-2f;

            c[0] /= y;
            c[1] /= y;
            c[2] -= 1.0f;

            table[i] = c;

            t += dt;
        }
    }
    
    void FindZHFromGammaTable(int tableSize, const Vec3f table[], int numZCoeffs, Vec3f zhCoeffs[])
    {        
        float dt = 1.0 / (tableSize - 1);
        float t = 0.0;
        
        for (int i = 0; i < numZCoeffs; i++)
            zhCoeffs[i] = vl_0;
        
        for (int i = 0; i < tableSize; i++)
        {
            AddZHSample(1 - 2 * sqr(t), table[i] * (4 * t), zhCoeffs);  // 4t is dw weight

            t += dt;
        }
        
        // integral of -4t dt over [0, 1] is 2 as you'd expect, but, discrete case above.
        for (int i = 0; i < numZCoeffs; i++)
            zhCoeffs[i] *= vl_twoPi / tableSize; // dw = 2pi dz
    }

    void GenerateGammaTableFromZH(int numZCoeffs, Vec3f zhCoeffs[], int tableSize, Vec3f table[])
    {
        float dt = 1.0 / (tableSize - 1);
        float t = 0.0;

        for (int i = 0; i < tableSize; i++)
        {
            table[i] = SampleZH(1 - 2 * sqr(t), zhCoeffs);
            t += dt;
        }
    }

    const float kThetaW = 0.01f;
    const float kGammaW = 0.002f;
}

void cSunSky::FindBRDFTables()
{
    static const float kRowPowers[] = { 1.0f, 16.0f, 128.0f };
    const int numRows = 4;

    const int kTableSize = sizeof(mThetaTable) / sizeof(mThetaTable[0]);

    // The BRDF tables cover the entire sphere, so we must resample theta from the Perez tables which cover a hemisphere.
    Vec3f thetaTable[kTableSize];
    const Vec3f* gammaTable = mGammaTable;

    // Fill top hemisphere of table
    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[kTableSize / 2 + i] = mThetaTable[2 * i];

    // Fill lower hemisphere with term that evaluates close to 0, to avoid below-ground luminance leaking in.
    for (int i = 0; i < kTableSize / 2; i++)
        thetaTable[i] = Vec3f(0.0f, 0.0f, -0.999f);

    // Project tables into ZH coefficients
    Vec3f zhCoeffsTheta[7];
    Vec3f zhCoeffsGamma[7];

    FindZHFromThetaTable(kTableSize, thetaTable, 7, zhCoeffsTheta);
    FindZHFromGammaTable(kTableSize, gammaTable, 7, zhCoeffsGamma);

    // Fill -ve z with reflected +ve z, ramped to 0 at z = -1. This avoids discontinuities at the horizon.
    for (int i = 0; i < kTableSize / 2; i++)
    {
        thetaTable[i] = mThetaTable[kTableSize - 1 - 2 * i];

        // Ramp luminance down
        float s = (i + 0.5f) / (kTableSize / 2);
        s = sqrtf(s);
        thetaTable[i][2] = thetaTable[i][2] * s + (s - 1);
    }

    for (int i = 0; i < kTableSize; i++)
    {
        mBRDFThetaTable[0][i] = thetaTable[i];
        mBRDFGammaTable[0][i] = gammaTable[i];
    }

    Vec3f zhCoeffsThetaConv[7];
    Vec3f zhCoeffsGammaConv[7];

    for (int r = 1; r < numRows; r++)
    {
        float s = kRowPowers[numRows - r - 1];

        ConvolveZHWithCosPow(s, 7, zhCoeffsTheta, zhCoeffsThetaConv);
        ConvolveZHWithCosPow(s, 7, zhCoeffsGamma, zhCoeffsGammaConv);

        ApplyZHWindowing(kThetaW, 7, zhCoeffsThetaConv);
        ApplyZHWindowing(kGammaW, 7, zhCoeffsGammaConv);

        GenerateThetaTableFromZH(7, zhCoeffsThetaConv, kTableSize, mBRDFThetaTable[r]);
        GenerateGammaTableFromZH(7, zhCoeffsGammaConv, kTableSize, mBRDFGammaTable[r]);
    }
}


Vec3f cSunSky::ConvolvedSkyRGB(const Vec3f& v, float roughness) const
{
    const int brdfSlices = sizeof(mBRDFThetaTable) / sizeof(mBRDFThetaTable[0]);
    const int tableSize = sizeof(mBRDFThetaTable[0]) / sizeof(mBRDFThetaTable[0][0]);

    int level = lrintf(roughness * (brdfSlices - 1));

    const Vec3f* thetaTable = mBRDFThetaTable[level];
    const Vec3f* gammaTable = mBRDFGammaTable[level];

    float cosTheta = v[2];
    float cosGamma = dot(mToSun, v);

    float c = 0.5f * (cosTheta + 1);
    float g = 0.5f * (1.0f - cosGamma);
    g = sqrtf(g);

    int ti = FloorToInt32(c * (tableSize - 1));
    int gi = FloorToInt32(g * (tableSize - 1));

    // (1 + F(theta)) * (1 + G(phi))
    Vec3f perez = (Vec3f(vl_1) + thetaTable[ti]) * (Vec3f(vl_1) + gammaTable[gi]);

    Vec3f xyY = perez * mPerezInvDen;

    return xyYToRGB(xyY);
}



