/*
    File:       SunSkyTest.cpp

    Function:   Simple test harness for various sky models implemented by cSunSky.

    Author:     Andrew Willmott

    Copyright:  (c) 1999-2014, Andrew Willmott

    Notes:      
*/

#include "SunSky.h"

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

extern "C"
{
    #include "targa.c"
}

namespace
{
    // latitude, longitude
    const Vec2f kAuckland     (-37.0,  174.8);
    const Vec2f kPittsburgh   ( 40.5, -80.22);
    const Vec2f kOakland      ( 37.8, -122.2);
    const Vec2f kSanFrancisco ( 37.8, -122.4);
    const Vec2f kLondon       (+51.5,  0.0);

    // TGA byte order BGRA
    inline float saturate(float s)
    {
        if (s < 0.0f)
            return 0.0f;
        if (s >= 1.0f)
            return 1.0f;
        return s;
    }

    inline uint32_t RGBFToU32(Vec3f rgb)
    {
        return
              0xFF000000
            | lrintf(saturate(rgb[0]) * 255.0f) << 16
            | lrintf(saturate(rgb[1]) * 255.0f) <<  8
            | lrintf(saturate(rgb[2]) * 255.0f);
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

        return lrintf(ceil((1.0f - maxX) * width / 2.0f));
    }

    struct cEnumInfo
    {
        const char* mName;
        int         mValue;
    };

    cEnumInfo kSkyTypeEnum[] =
    {
        "clearUtah",        kUtahClearSky,
        "clearUtahTable",   kUtahClearSkyTable,
        "clearUtahBRDF",    kUtahClearSkyBRDF,
        "clearCIE",         kCIEClearSky,
        "cloudyCIE",        kCIECloudySky,
        "partlyCloudyCIE",  kCIEPartlyCloudySky,
        0, 0
    };
}

/// Fill top-down projection of upper or lower hemisphere
void SkyToHemisphere(const cSunSky& sunSky, int width, int height, uint8_t* data, int stride, float weight = 1.0f, float gamma = 2.2f, float hemiSign = 1.0f)
{
    float invGamma = 0.0;

    if (gamma > 0.0)
        invGamma = 1.0f / gamma;

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

            Vec3f v(x, y, hemiSign * sqrtf(1.0f - h2));

            Vec3f c = sunSky.SkyRGB(v);

            c *= weight;

            if (invGamma > 0.0)
            {
                c[0] = powf(c[0], invGamma);
                c[1] = powf(c[1], invGamma);
                c[2] = powf(c[2], invGamma);
            }

            row[j] = RGBFToU32(c);
        }

        // fill in surrounds by replicating edge texels
        for (int j = 0; j < sw; j++)
            row[j] = row[sw];
        for (int j = width - sw; j < width; j++)
            row[j] = row[width - sw - 1];

        data += stride;
    }
}

void SkyToCubeFace(const cSunSky& sunSky, int face, int width, int height, uint8_t* data, int stride, float weight = 1.0f, float gamma = 0.0f)
{
    float     invGamma = 0.0;
    Vec3f     groundColour(0.0f, 1.0f, 0.0f);

    if (gamma != 0.0 && gamma != 1.0f)
        invGamma = 1.0 / gamma;

    for (int i = 0; i < height; i++)
    {
        uint32_t* row = (uint32_t*) data;

        for (int j = 0; j < width; j++)
        {
            Vec3f faceDir;
            Vec3f faceColour;
            Vec2f facePos((j + 0.5) / width, (i + 0.5) / height);

            facePos = (2.0 * facePos - vl_1);

            switch (face)
            {
            case 0:
                // +ve Y
                faceDir = Vec3f(facePos[0], 1.0, facePos[1]);
                break;
            case 1:
                // +ve X
                faceDir = Vec3f(1.0, -facePos[0], facePos[1]);
                break;
            case 2:
                // -ve Y
                faceDir = Vec3f(-facePos[0], -1.0, facePos[1]);
                break;
            case 3:
                // -ve X
                faceDir = Vec3f(-1.0, facePos[0], facePos[1]);
                break;
            case 4:
                // +ve Z
                faceDir = Vec3f(facePos[0], -facePos[1], 1.0);
                break;
            case 5:
                // -ve Z
                faceDir = Vec3f(facePos[0], facePos[1], -1.0);
                break;
            }                            

            faceDir = norm(faceDir);
            faceColour = sunSky.SkyRGB(faceDir) * weight;

            if (invGamma != 0.0)
            {
                faceColour[0] = powf(faceColour[0], invGamma);
                faceColour[1] = powf(faceColour[1], invGamma);
                faceColour[2] = powf(faceColour[2], invGamma);
            }
            
            row[j] = RGBFToU32(faceColour);
        }

        data += stride;
    }
}

namespace
{
    int Help(const char* command)
    {
        printf
        (
            "%s <options>\n"
            "\n"
            "Options:\n"
            "  -h : this help\n"
            "  -t <time>\n"
            "  -d <day of year> : 0-365\n"
            "  -b <tubidity> : 2 - 30\n"
            "  -l <longitude> <latitude>\n"
            "  -w <normalisation weight>\n"
            "  -g <gamma>\n"
            "  -c <cubeMap>\n"
            , command
        );

        return 0;
    }
}

int main(int argc, const char* argv[])
{
    const char* command = argv[0];
    argv++; argc--; // chomp

    // Set up defaults
    time_t unixTime;
    ::time(&unixTime);

    struct tm* theTime = localtime(&unixTime);

    float julianDay = theTime->tm_yday;

    float localTime = theTime->tm_hour + theTime->tm_min / 60.0 + theTime->tm_sec / 3600.0;

    float timeZone = 0.0f;
    if (theTime->tm_isdst != 0)
        timeZone += 1.0;

    Vec2f latLong = kLondon;
    float turbidity = 2.5;

    float gamma = 1.0f;
    float hemiSign = 1.0f; // upper
    float weight = 5e-5f;

    float roughness = -1.0f;
    bool cubeMap = false;
    tSkyType skyType = kUtahClearSky;

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
            localTime = atof(argv[0]);
            argv++; argc--;
            break;

        case 'd':
            if (ArgCountError(option, 1, argc))
                return -1;
            julianDay = atof(argv[0]);
            argv++; argc--;
            break;

        case 'b':
            if (ArgCountError(option, 1, argc))
                return -1;
            turbidity = atof(argv[0]);
            argv++; argc--;
            break;

        case 'l':
            if (ArgCountError(option, 2, argc))
                return -1;
            latLong[0] = atof(argv[0]);
            argv++; argc--;
            latLong[1] = atof(argv[0]);
            argv++; argc--;
            break;

        case 'w':
            if (ArgCountError(option, 1, argc))
                return -1;
            weight = atof(argv[0]);
            argv++; argc--;

            break;

        case 'g':
            if (ArgCountError(option, 1, argc))
                return -1;
            gamma = atof(argv[0]);
            argv++; argc--;
            break;

        case 'c':
            cubeMap = !cubeMap;
            break;
        case 'f':
            hemiSign = -hemiSign;
            break;

        case 's':
            {
            if (ArgCountError(option, 1, argc))
                return -1;

                const char* typeName = argv[0];
                argv++; argc--;

                tSkyType newSkyType = kNumSkyTypes;

                for (const cEnumInfo* info = kSkyTypeEnum; info->mName; info++)
                    if (strcasecmp(info->mName, typeName) == 0)
                    {
                        newSkyType = tSkyType(info->mValue);
                        break;
                    }

                if (newSkyType == kNumSkyTypes)
                {
                    fprintf(stderr, "Unknown sky type: %s\n", typeName);
                    return -1;
                }

                skyType = newSkyType;
            }
            break;

        case 'r':
            if (ArgCountError(option, 1, argc))
                return -1;
            roughness = saturate(atof(argv[0]));
            skyType = kUtahClearSkyBRDF;
            argv++; argc--;
            break;


        default:
            fprintf(stderr, "Unrecognised option: %s", option);
            return -1;
        }
    }

    if (argc > 0)
    {
        fprintf(stderr, "Unrecognised arguments starting with %s\n", argv[0]);
        return -1;
    }

    cSunSky sunSky;

    sunSky.SetSkyType(skyType);
    if (roughness >= 0.0f)
        sunSky.SetRoughness(roughness);
    sunSky.SetAllParams(localTime, timeZone, julianDay, latLong[0], latLong[1], turbidity);

    uint32_t image[256][256];

    tga_image tga;
    init_tga_image(&tga, (uint8_t*) image, 256, 256, 32);
    tga.image_type = TGA_IMAGE_TYPE_BGR;
    tga.image_descriptor &= ~TGA_T_TO_B_BIT;

    if (!cubeMap)
    {
        SkyToHemisphere(sunSky, 256, 256, (uint8_t*) image, 1024, weight, gamma, hemiSign);

        char fileName[32];
        snprintf(fileName, 32, "sky-hemi.tga");
        tga_write(fileName, &tga);
    }
    else
    {
        for (int i = 0; i < 6; i++)
        {
            SkyToCubeFace(sunSky, i, 256, 256, (uint8_t*) image, 1024, weight, gamma);

            char fileName[32];
            snprintf(fileName, 32, "sky-cube-%d.tga", i);
            tga_write(fileName, &tga);
        }
    }

    return 0;
}
