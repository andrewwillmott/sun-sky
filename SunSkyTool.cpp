//
//  File:       SunSkyTool.cpp
//
//  Function:   Test tool for SunSky.*
//
//  Copyright:  Andrew Willmott
//

#define  _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "SunSky.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef _MSC_VER
    #include <unistd.h>
    #include <strings.h>
#else
    #include <string.h>
    #define strcasecmp _stricmp
#endif

#include "stb_image_mini.h"

// #define EDGE_FILL

using namespace SSLib;

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

namespace
{
    // latitude, longitude
    const Vec2f kLondon       (+51.5f,  0.0f);
    const Vec2f kAuckland     (-37.0f,  174.8f);
    const Vec2f kPittsburgh   ( 40.5f, -80.22f);
    const Vec2f kOakland      ( 37.8f, -122.2f);
    const Vec2f kSanFrancisco ( 37.8f, -122.4f);
    const Vec2f kJakarta      (-6.21f,  106.85f);

    inline float saturate(float s)
    {
        if (s < 0.0f)
            return 0.0f;
        if (s >= 1.0f)
            return 1.0f;
        return s;
    }

    inline float Max(float a, float b)
    {
        return b < a ? a : b;
    }

    inline Vec3f MaxElts(const Vec3f& a, const Vec3f& b)
    {
        return Vec3f
        (
            Max(a.x, b.x),
            Max(a.y, b.y),
            Max(a.z, b.z)
        );
    }

    inline uint32_t RGBFToU32(Vec3f rgb)
    {
        return
              0xFF000000
            | (int) lrintf(saturate(rgb.x) * 255.0f) <<  0
            | (int) lrintf(saturate(rgb.y) * 255.0f) <<  8
            | (int) lrintf(saturate(rgb.z) * 255.0f) << 16;
    }

    bool ArgCountError(const char* opt, int expected, int argc)
    {
        if (argc < expected)
        {
            fprintf(stderr, "Not enough arguments for %s: expected %d, have %d\n", opt, expected, argc);
            return true;
        }

        return false;
    }

    inline int HemiInset(float y2, int width)
    {
        float maxX2 = 1.0f - y2;
        float maxX = sqrtf(maxX2);

        return (int) lrintf(ceil((1.0f - maxX) * width / 2.0f));
    }

    Vec3f pow(Vec3f v, float n)
    {
        return Vec3f(powf(v.x, n), powf(v.y, n), powf(v.z, n));
    }

    Vec3f toneMapLinear(Vec3f c, float weight)
    {
        return c * weight;
    }

    Vec3f toneMapExp(Vec3f c, float weight)
    {
        return Vec3f(vl_one) - Vec3f(exp(-weight * c.x), exp(-weight * c.y), exp(-weight * c.z));
    }

    Vec3f toneMapReinhard(Vec3f c, float weight)
    {
        c *= weight;
        return c / (Vec3f(vl_one) + c);
    }

    typedef Vec3f ToneMapFunc(Vec3f c, float weight);

    enum kToneMapType
    {
        kToneMapLinear,
        kToneMapExponential,
        kToneMapReinhard,
        kNumToneMapTypes
    };

    ToneMapFunc* kToneMapFuncs[kNumToneMapTypes + 1] =
    {
        toneMapLinear,
        toneMapExp,
        toneMapReinhard,
        nullptr
    };

    struct cMapInfo
    {
        float weight    = 5e-5f;
        float gamma     = 2.2f;
        float hemiSign  = 1.0f;
        bool  fisheye   = false;

        ToneMapFunc* toneMap = toneMapLinear;
    };

    bool PFMWrite(const char* filename, int width, int height, Vec3f* image)
    {
        FILE* f = fopen(filename, "w");

        if (!f)
            return false;

        fprintf(f, "PF\n");
        fprintf(f, "%d %d\n", width, height);
        fprintf(f, "-1.0\n");   // -ve = little endian

        image += width * height;

        for (int i = 0; i < height; i++)
        {
            image -= width;
            fwrite(image, sizeof(float) * 3, width, f);
        }

        fclose(f);

        return true;
    }
}


//------------------------------------------------------------------------------
// Projected (or fisheye) hemisphere in LDR (png) and HDR (pfm)
//------------------------------------------------------------------------------

namespace
{
    /// Fill top-down projection of upper or lower hemisphere
    void SkyToHemisphere(const cSunSky& sunSky, int width, int height, uint8_t* data, int stride, const cMapInfo& mi)
    {
        float invGamma = 1.0;

        if (mi.gamma > 0.0)
            invGamma = 1.0f / mi.gamma;

        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            uint32_t* row = (uint32_t*) data;

            float y = 2.0f * (i + 0.5f) / height - 1.0f;
            float y2 = y * y;

            int sw = HemiInset(y2, width);

            for (int j = sw; j < width - sw; j++)
            {
                float x = 2.0f * (j + 0.5f) / width - 1.0f;
                float x2 = x * x;
                float h2 = x2 + y2;

                Vec3f v;

                if (mi.fisheye)
                {
                    float theta = vl_halfPi - vl_halfPi * sqrtf(h2);
                    float phi = atan2f(y, x);
                    v = Vec3f(cos(phi) * cos(theta), sin(phi) * cos(theta), sin(theta));
                }
                else
                    v = Vec3f(x, y, mi.hemiSign * sqrtf(1.0f - h2));

                Vec3f c = sunSky.SkyRGB(v);

                c = mi.toneMap(c, mi.weight);
                c = pow(c, invGamma);

                row[j] = RGBFToU32(c);
            }

            // fill in surrounds
        #ifdef EDGE_FILL
            for (int j = 0; j < sw; j++)
                row[j] = row[sw];
            for (int j = width - sw; j < width; j++)
                row[j] = row[width - sw - 1];
        #else
            for (int j = 0; j < sw; j++)
                row[j] = 0xFF000000;
            for (int j = width - sw; j < width; j++)
                row[j] = 0xFF000000;
        #endif

            data -= stride;
        }
    }

    struct cStats
    {
        Vec3f avg;
        Vec3f max;
        Vec3f dev;
    };

    void SkyToHemisphere(const cSunSky& sunSky, int width, int height, Vec3f* data, const cMapInfo& mi, cStats* stats)
    {
        Vec3f maxElts(vl_0);
        Vec3f sumElts(vl_0);
        Vec3f varElts(vl_0);
        int samples = 0;
        int stride = width;

        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            Vec3f* row = data;

            float y = 2.0f * (i + 0.5f) / height - 1.0f;
            float y2 = y * y;

            int sw = HemiInset(y2, width);

            for (int j = sw; j < width - sw; j++)
            {
                float x = 2.0f * (j + 0.5f) / width - 1.0f;
                float x2 = x * x;
                float h2 = x2 + y2;

                Vec3f v(x, y, mi.hemiSign * sqrtf(1.0f - h2));

                Vec3f c = sunSky.SkyRGB(v);

                if (stats)
                {
                    maxElts = MaxElts(maxElts, c);
                    sumElts += c;
                    varElts += c * c;
                    samples++;
                }

                c *= mi.weight;

                row[j] = c;
            }

        #ifdef EDGE_FILL
            // fill in surrounds by replicating edge texels
            for (int j = 0; j < sw; j++)
                row[j] = row[sw];
            for (int j = width - sw; j < width; j++)
                row[j] = row[width - sw - 1];
        #else
            // fill in surrounds by replicating edge texels
            for (int j = 0; j < sw; j++)
                row[j] = vl_0;
            for (int j = width - sw; j < width; j++)
                row[j] = vl_0;
        #endif

            data -= stride;
        }

        if (stats)
        {
            stats->avg = sumElts / float(samples);
            stats->max = maxElts;
            varElts = varElts / float(samples) - sqr(stats->avg);
            stats->dev = Vec3f(sqrtf(varElts.x), sqrtf(varElts.y), sqrtf(varElts.z));
        }
    }
}


//------------------------------------------------------------------------------
// Cubemap generation in LDR (png) and HDR (pfm)
//------------------------------------------------------------------------------

namespace
{
    const int kFaceIndices[6][3] =
    {
        0, 2, 1,
        2, 0, 1,
        0, 2, 1,
        2, 0, 1,
        0, 1, 2,
        0, 1, 2,
    };
    const float kFaceSigns[6][3] =
    {
        +1.0f, +1.0f, +1.0f,
        +1.0f, -1.0f, +1.0f,
        -1.0f, -1.0f, +1.0f,
        -1.0f, +1.0f, +1.0f,
        +1.0f, -1.0f, +1.0f,
        +1.0f, +1.0f, -1.0f,
    };

    void SkyToCubeFace(const cSunSky& sunSky, int face, int width, int height, uint8_t* data, int stride, const cMapInfo& mi)
    {
        float invGamma = 1.0;

        if (mi.gamma > 0.0)
            invGamma = 1.0f / mi.gamma;

        const float* signs   = kFaceSigns  [face];
        const int*   indices = kFaceIndices[face];

        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            uint32_t* row = (uint32_t*) data;

            for (int j = 0; j < width; j++)
            {
                Vec3f facePos(2 * (j + 0.5f) / width - 1, 2 * (i + 0.5f) / height - 1, 1.0f);

                Vec3f faceDir
                (
                    signs[0] * facePos[indices[0]],
                    signs[1] * facePos[indices[1]],
                    signs[2] * facePos[indices[2]]
                );

                faceDir = norm(faceDir);
                Vec3f faceColour = sunSky.SkyRGB(faceDir);

                faceColour = mi.toneMap(faceColour, mi.weight);
                faceColour = pow(faceColour, invGamma);

                row[j] = RGBFToU32(faceColour);
            }

            data -= stride;
        }
    }

    void SkyToCubeFace(const cSunSky& sunSky, int face, int width, int height, Vec3f* data, const cMapInfo& mi)
    {
        const float* signs   = kFaceSigns  [face];
        const int*   indices = kFaceIndices[face];

        int stride = width;
        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            Vec3f* row = data;

            for (int j = 0; j < width; j++)
            {
                Vec3f facePos(2 * (j + 0.5f) / width - 1, 2 * (i + 0.5f) / height - 1, 1.0f);

                Vec3f faceDir
                (
                    signs[0] * facePos[indices[0]],
                    signs[1] * facePos[indices[1]],
                    signs[2] * facePos[indices[2]]
                );

                faceDir = norm(faceDir);

                row[j] = sunSky.SkyRGB(faceDir) * mi.weight;
            }

            data -= stride;
        }
    }
}


//------------------------------------------------------------------------------
// Panorama generation in LDR (png) and HDR (pfm)
//------------------------------------------------------------------------------

namespace
{
    void SkyToPanoramic(const cSunSky& sunSky, int height, uint8_t* data, int stride, const cMapInfo& mi)
    {
        float invGamma = 1.0;

        if (mi.gamma > 0.0)
            invGamma = 1.0f / mi.gamma;

        int width = 2 * height;
        if (stride == 0)
            stride = 4 * width;

        float da = vl_pi / height;
        float phi = vl_pi - 0.5f * da;

        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            uint32_t* row = (uint32_t*) data;

            float theta = 0.5f * da;
            float sp = sinf(phi);
            float cp = cosf(phi);

            for (int j = 0; j < width; j++)
            {
                float st = sinf(theta);
                float ct = cosf(theta);

                // middle of image is north, east to right, west to left, edges are south
                Vec3f dir(-st * sp, -ct * sp, cp);

                Vec3f c = sunSky.SkyRGB(dir);

                c = mi.toneMap(c, mi.weight);
                c = pow(c, invGamma);

                row[j] = RGBFToU32(c);
                theta += da;
            }

            data -= stride;
            phi -= da;
        }
    }

    void SkyToPanoramic(const cSunSky& sunSky, int height, Vec3f* data, const cMapInfo& mi)
    {
        int width = 2 * height;
        int stride = width;

        float da = vl_pi / height;
        float phi = vl_pi - 0.5f * da;

        data += (height - 1) * stride;

        for (int i = 0; i < height; i++)
        {
            Vec3f* row = data;

            float theta = 0.5f * da;
            float sp = sinf(phi);
            float cp = cosf(phi);

            for (int j = 0; j < width; j++)
            {
                float st = sinf(theta);
                float ct = cosf(theta);

                // middle of image is north, east to right, west to left, edges are south
                Vec3f dir(-st * sp, -ct * sp, cp);

                row[j] = sunSky.SkyRGB(dir) * mi.weight;
                theta += da;
            }

            data -= stride;
            phi -= da;
        }
    }

    const float kMinAutoLum    = 2000.0f;
    const float kAutoLumTarget = 0.4f;
}



//------------------------------------------------------------------------------
// Main program
//------------------------------------------------------------------------------

namespace
{
    struct cEnumInfo
    {
        const char* mName;
        const char* mShort;
        int         mValue;
    };

    cEnumInfo kSkyTypeEnum[] =
    {
        "Preetham",         "pt",   kPreetham,
        "PreethamTable",    "ptt",  kPreethamTable,
        "PreethamBRDF",     "ptb",  kPreethamBRDF,
        "Hosek",            "hk",   kHosek,
        "HosekTable",       "hkt",  kHosekTable,
        "HosekBRDF",        "hkb",  kHosekBRDF,
        "cieClear",         "cc",   kCIEClear,
        "cieOvercast",      "co",   kCIEOvercast,
        "ciePartlyCloudy",  "cp",   kCIEPartlyCloudy,
        nullptr, nullptr, 0
    };

    cEnumInfo kToneMapTypeEnum[] =
    {
        "linear",       "l",   kToneMapLinear,
        "exponential",  "ex",  kToneMapExponential,
        "reinhard",     "rh",  kToneMapReinhard,
        nullptr, nullptr, 0
    };

    int ArgEnum(const cEnumInfo info[], const char* name, int defaultValue = -1)
    {
        for ( ; info->mName; info++)
            if (strcasecmp(info->mName, name) == 0 || strcasecmp(info->mShort, name) == 0)
            {
                return info->mValue;
                break;
            }

        return defaultValue;
    }

    int Help(const char* command)
    {
        printf
        (
            "%s <options>\n"
            "\n"
            "Options:\n"
            "  -h : this help\n"
            "  -s <skyType> : use given sky type (default: Preetham.)\n"
            "  -t <time>          : 0 - 24\n"
            "  -d <day of year>   : 0 - 365\n"
            "  -b <turbidity>     : 2 - 12\n"
            "  -x <l>|<r g b>     : 0 - 1, specify ground albedo for Hosek\n"
            "  -r <roughness>     : 0 - 1, specify roughness for BRDF types\n"
            "  -l <latitude> <longitude>\n"
            "  -w <luminance scale>\n"
            "  -g <gamma>\n"
            "  -e <tonemapType> : use given tonemap operator (default: linear)\n"
            "  -a : autoscale intensity\n"
            "  -i : invert hemisphere\n"
            "  -f : fisheye rather than cos projection\n"
            "  -c : output cubemap instead\n"
            "  -p : output panorama instead\n"
            "  -m : output movie, record day as sky.mp4, requires ffmpeg\n"
            "  -v : verbose\n"
            , command
        );

        printf("\n" "skyType:\n");
        for (const cEnumInfo* info = kSkyTypeEnum; info->mName; info++)
            printf("  %-16s (%s)\n", info->mName, info->mShort);

        printf("\n" "toneMapType:\n");
        for (const cEnumInfo* info = kToneMapTypeEnum; info->mName; info++)
            printf("  %-16s (%s)\n", info->mName, info->mShort);

        return 0;
    }
}

int main(int argc, const char* argv[])
{
    const char* command = argv[0];
    argv++; argc--; // chomp

    if (argc == 0)
        return Help(command);

    // Set up defaults
    time_t unixTime;
    ::time(&unixTime);

    struct tm* theTime = localtime(&unixTime);

    int julianDay = theTime->tm_yday;

    float localTime = theTime->tm_hour + theTime->tm_min / 60.0f + theTime->tm_sec / 3600.0f;
    bool dst = (theTime->tm_isdst != 0);
    Vec2f latLong   = kLondon;
    float turbidity = 2.5;
    Vec3f albedo    = vl_0;
    float overcast  = 0.0f;

    cMapInfo mi;
    mi.weight = -1.0f;

    float roughness = -1.0f;
    bool autoscale  = false;
    bool cubeMap    = false;
    bool panoramic  = false;
    bool movie      = false;
    bool verbose    = false;
    tSkyType skyType = kPreetham;

    // Options
    while (argc > 0 && argv[0][0] == '-')
    {
        const char* option = argv[0] + 1;
        argv++; argc--;

        switch (option[0])
        {
        case 'h':
        case '?':
            return Help(command);

        case 't':
            if (ArgCountError(option, 1, argc))
                return -1;
            localTime = (float) atof(argv[0]);
            argv++; argc--;
            break;

        case 'd':
            if (ArgCountError(option, 1, argc))
                return -1;
            julianDay = atoi(argv[0]);
            argv++; argc--;
            break;

        case 'b':
            if (ArgCountError(option, 1, argc))
                return -1;
            turbidity = (float) atof(argv[0]);
            argv++; argc--;
            break;

        case 'o':
            overcast = (float) atof(argv[0]);
            argv++; argc--;
            break;

        case 'x':
            if (ArgCountError(option, 1, argc))
                return -1;
            albedo.x = (float) atof(argv[0]);
            argv++; argc--;

            if (argc >= 1 && argv[0][0] != '-')
            {
                albedo.y = (float) atof(argv[0]);
                argv++; argc--;
            }
            else
                albedo.y = albedo.x;

            if (argc >= 1 && argv[0][0] != '-')
            {
                albedo.z = (float) atof(argv[0]);
                argv++; argc--;
            }
            else
                albedo.z = albedo.y;

            break;

        case 'l':
            if (ArgCountError(option, 2, argc))
                return -1;
            latLong[0] = (float) atof(argv[0]);
            argv++; argc--;
            latLong[1] = (float) atof(argv[0]);
            argv++; argc--;

            dst = false;    // don't take dst from local time info
            break;

        case 'w':
            if (ArgCountError(option, 1, argc))
                return -1;
            mi.weight = (float) atof(argv[0]);
            argv++; argc--;

            break;

        case 'g':
            if (ArgCountError(option, 1, argc))
                return -1;
            mi.gamma = (float) atof(argv[0]);
            argv++; argc--;
            break;

        case 'a':
            autoscale = !autoscale;
            break;
        case 'c':
            cubeMap = !cubeMap;
            break;
        case 'p':
            panoramic = !panoramic;
            break;
        case 'm':
            movie = !movie;
            break;
        case 'i':
            mi.hemiSign = -mi.hemiSign;
            break;
        case 'f':
            mi.fisheye = !mi.fisheye;
            break;
        case 'v':
            verbose = !verbose;
            break;

        case 's':
            {
                if (ArgCountError(option, 1, argc))
                    return -1;

                const char* typeName = argv[0];
                argv++; argc--;

                skyType = (tSkyType) ArgEnum(kSkyTypeEnum, typeName, kNumSkyTypes);

                if (skyType == kNumSkyTypes)
                {
                    fprintf(stderr, "Unknown sky type: %s\n", typeName);
                    return -1;
                }
            }
            break;

        case 'e':
            {
                if (ArgCountError(option, 1, argc))
                    return -1;

                const char* typeName = argv[0];
                argv++; argc--;

                mi.toneMap = kToneMapFuncs[ArgEnum(kToneMapTypeEnum, typeName, kNumToneMapTypes)];

                if (!mi.toneMap)
                {
                    fprintf(stderr, "Unknown tone map type: %s\n", typeName);
                    return -1;
                }
            }
            break;

        case 'r':
            if (ArgCountError(option, 1, argc))
                return -1;
            roughness = saturate((float) atof(argv[0]));
            argv++; argc--;
            break;


        default:
            fprintf(stderr, "Unrecognised option: %s\n", option);
            return -1;
        }
    }

    if (argc > 0)
    {
        fprintf(stderr, "Unrecognised arguments starting with %s\n", argv[0]);
        return -1;
    }

    float timeZone = rintf(latLong[1] / 15.0f);    // estimate for now

    if (dst)
        timeZone += 1.0;

    Vec3f sunDir = SunDirection(localTime, timeZone, julianDay, latLong[0], latLong[1]);

    cSunSky sunSky;
    sunSky.SetSkyType(skyType);

    sunSky.SetSunDir(sunDir);
    sunSky.SetTurbidity(turbidity);
    sunSky.SetAlbedo(albedo);
    sunSky.SetOvercast(overcast);

    if (roughness >= 0.0f)
        sunSky.SetRoughness(roughness);

    sunSky.Update();

    if (verbose)
    {
        printf("Time: %g, time zone: %g, day: %d, latitude: %g, longitude: %g, turbidity: %g, albedo: %g\n", localTime, timeZone, julianDay, latLong[0], latLong[1], turbidity, albedo.y);

        float theta = asinf (sunDir.z);
        float phi   = atan2f(sunDir.y, sunDir.x);

        theta *= 180.0f / vl_pi;
        phi   *= 180.0f / vl_pi;

        phi   = 90.0f - phi;    // make relative to North rather than East, and clockwise.

        if (phi < 0.0f)
            phi += 360.0f;

        printf("Sun elevation      : %g\n", theta);
        printf("Sun compass heading: %g\n", phi  );
    }

    if (mi.weight < 0.0f)
        switch (skyType)
        {
        case kHosek:
        case kHosekTable:
        case kHosekBRDF:
            mi.weight = 8e-5f;
            break;
        default:
            mi.weight = 5e-5f;
        }

    if (!movie && autoscale)
    {
        float avgLum = sunSky.AverageLuminance();

        if (verbose)
            printf("Average luminance: %g\n", avgLum);

        // Once we get dark enough (sun below horizon), stop auto-scaling, so we don't snap to black
        if (avgLum < kMinAutoLum)
            avgLum = kMinAutoLum;

        float lumScale = kAutoLumTarget / avgLum;

        if (verbose)
            printf("Autoscaling luminance by: %g\n", lumScale);

        mi.weight = lumScale;
    }

    if (verbose)
        printf("Ouput: weight = %g, gamma = %g\n", mi.weight, mi.gamma);

    char fileName[32];

    if (panoramic)
    {
        uint32_t image   [256][512];
        Vec3f    imageHDR[256][512];

        SkyToPanoramic(sunSky, 256, (uint8_t*) image, 0, mi);

        snprintf(fileName, 32, "sky-panoramic.png");

        if (stbi_write_png(fileName, 512, 256, 4, image[0], 0) != 0)
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);

        SkyToPanoramic(sunSky, 256, imageHDR[0], mi);

        snprintf(fileName, 32, "sky-panoramic.pfm");

        if (PFMWrite(fileName, 512, 256, imageHDR[0]))
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);
    }
    else if (cubeMap)
    {
        uint32_t image   [256][256];
        Vec3f    imageHDR[256][256];

        for (int i = 0; i < 6; i++)
        {
            SkyToCubeFace(sunSky, i, 256, 256, (uint8_t*) image, 4 * 256, mi);

            snprintf(fileName, 32, "sky-cube-%d.png", i);

            if (stbi_write_png(fileName, 256, 256, 4, image[0], 0) != 0)
                printf("wrote %s\n", fileName);
            else
                printf("failed to write %s\n", fileName);

            SkyToCubeFace(sunSky, i, 256, 256, imageHDR[0], mi);

            snprintf(fileName, 32, "sky-cube-%d.pfm", i);

            if (PFMWrite(fileName, 256, 256, imageHDR[0]))
                printf("wrote %s\n", fileName);
            else
                printf("failed to write %s\n", fileName);
        }
    }
#ifndef _MSC_VER
    else if (movie)
    {
        uint32_t image[256][256];

        // crf = constant rate factor, 0 - 51, 0 is lossless, 51 worst
        // -preset = veryfast/faster/fast/medium/slow/slower/veryslow
        const char* cmd = "/usr/local/bin/ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s 256x256 -i - -threads 0 -preset medium -y -pix_fmt yuv420p -crf 10 sky.mp4";

        // open pipe to ffmpeg's stdin in binary write mode
        FILE* ffmpeg = popen(cmd, "w");

        if (!ffmpeg)
        {
            perror("cmd");
            return -1;
        }

        for (float time = 6.0f; time <= 21.0f; time += 0.1f)
        {
            sunSky.SetSunDir(SunDirection(time, timeZone, julianDay, latLong[0], latLong[1]));
            sunSky.Update();

            if (autoscale)
            {
                float avgLum = sunSky.AverageLuminance();

                if (verbose)
                    printf("Average luminance: %g\n", avgLum);

                // Once we get dark enough (sun below horizon), stop auto-scaling, so we don't snap to black
                if (avgLum < kMinAutoLum)
                    avgLum = kMinAutoLum;

                float lumScale = kAutoLumTarget / avgLum;

                if (verbose)
                    printf("Autoscaling luminance by: %g\n", lumScale);

                mi.weight = lumScale;
            }

            SkyToHemisphere(sunSky, 256, 256, (uint8_t*) image, 1024, mi);

            fwrite(image, sizeof(image), 1, ffmpeg);
        }

        if (pclose(ffmpeg) == 0)
            printf("wrote sky.mp4\n");
        else
            printf("failed to write sky.mp4\n");
    }
#endif
    else
    {
        uint32_t image   [256][256];
        Vec3f    imageHDR[256][256];

        SkyToHemisphere(sunSky, 256, 256, (uint8_t*) image, 1024, mi);

        snprintf(fileName, 32, "sky-hemi.png");

        if (stbi_write_png(fileName, 256, 256, 4, image[0], 0) != 0)
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);

        cStats stats;
        SkyToHemisphere(sunSky, 256, 256, imageHDR[0], mi, verbose ? &stats : nullptr);

        snprintf(fileName, 32, "sky-hemi.pfm");

        if (PFMWrite(fileName, 256, 256, imageHDR[0]))
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);

        if (verbose)
        {
            printf("avg: %8.2f, %8.2f, %8.2f\n", stats.avg.x, stats.avg.y, stats.avg.z);
            printf("max: %8.2f, %8.2f, %8.2f\n", stats.max.x, stats.max.y, stats.max.z);
            printf("dev: %8.2f, %8.2f, %8.2f\n", stats.dev.x, stats.dev.y, stats.dev.z);
        }
    }

    return 0;
}
