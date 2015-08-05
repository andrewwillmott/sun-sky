/*
    File:       SunSkyTest.cpp

    Function:   Simple test harness for various sky models implemented by cSunSky.

    Author:     Andrew Willmott

    Copyright:  (c) 2014, Andrew Willmott

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

void SkyToHemisphere(const cSunSky& sunSky, int width, int height, Vec3f* data, float weight = 1.0f, float hemiSign = 1.0f)
{
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

            Vec3f v(x, y, hemiSign * sqrtf(1.0f - h2));

            Vec3f c = sunSky.SkyRGB(v);

            c *= weight;

            row[j] = c;
        }

        // fill in surrounds by replicating edge texels
        for (int j = 0; j < sw; j++)
            row[j] = row[sw];
        for (int j = width - sw; j < width; j++)
            row[j] = row[width - sw - 1];

        data += width;
    }
}

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
}

void SkyToCubeFace(const cSunSky& sunSky, int face, int width, int height, uint8_t* data, int stride, float weight = 1.0f, float gamma = 0.0f)
{
    float invGamma = 0.0;

    if (gamma != 0.0 && gamma != 1.0f)
        invGamma = 1.0 / gamma;

    const float* signs   = kFaceSigns  [face];
    const int*   indices = kFaceIndices[face];

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
            Vec3f faceColour = sunSky.SkyRGB(faceDir) * weight;

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

void SkyToCubeFace(const cSunSky& sunSky, int face, int width, int height, Vec3f* data, float weight = 1.0f)
{
    const float* signs   = kFaceSigns  [face];
    const int*   indices = kFaceIndices[face];

    for (int i = 0; i < height; i++)
    {
        Vec3f* row = data;

        for (int j = 0; j < width; j++)
        {
            Vec3f facePos(2 * (j + 0.5f) / width - 1, 2 * (i + 0.5f) / height - 1, 1.0f);

            Vec3f faceDir;
            faceDir[0] = signs[0] * facePos[indices[0]];
            faceDir[1] = signs[1] * facePos[indices[1]];
            faceDir[2] = signs[2] * facePos[indices[2]];

            faceDir = norm(faceDir);

            row[j] = sunSky.SkyRGB(faceDir) * weight;
        }

        data += width;
    }
}

void SkyToPanoramic(const cSunSky& sunSky, int height, Vec3f* data, float weight = 1.0f)
{
    int width = 2 * height;
    float da = vl_pi / height;
    float phi = vl_pi - 0.5f * da;

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

            row[j] = sunSky.SkyRGB(dir) * weight;

            theta += da;
        }

        data += width;
        phi -= da;
    }
}

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
            "  -t <time> : 0-24\n"
            "  -d <day of year> : 0-365\n"
            "  -b <tubidity> : 2 - 30\n"
            "  -l <longitude> <latitude>\n"
            "  -w <normalisation weight>\n"
            "  -g <gamma>\n"
            "  -a : autoscale intensity\n"
            "  -f : flip hemisphere\n"
            "  -c : output cubemap instead\n"
            "  -p : output panorama instead\n"
            "  -v : verbose\n"
            "  -s <skyType> : use given sky type\n"
            "  -r <roughness> : specify roughness for clearUtahBRDF\n"
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

    Vec2f latLong = kLondon;
    float turbidity = 2.5;

    float gamma = 1.0f;
    float hemiSign = 1.0f; // upper
    float weight = 5e-5f;

    float roughness = -1.0f;
    bool autoscale  = false;
    bool cubeMap    = false;
    bool panoramic  = false;
    bool verbose    = false;
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

        case 'a':
            autoscale = !autoscale;
            break;
        case 'c':
            cubeMap = !cubeMap;
            break;
        case 'p':
            panoramic = !panoramic;
            break;
        case 'f':
            hemiSign = -hemiSign;
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
            fprintf(stderr, "Unrecognised option: %s\n", option);
            return -1;
        }
    }

    if (argc > 0)
    {
        fprintf(stderr, "Unrecognised arguments starting with %s\n", argv[0]);
        return -1;
    }

    float timeZone = lrintf(latLong[1] / 15.0f);    // estimate for now

    if (theTime->tm_isdst != 0)
        timeZone += 1.0;

    cSunSky sunSky;

    sunSky.SetSkyType(skyType);

    if (roughness >= 0.0f)
        sunSky.SetRoughness(roughness);

    sunSky.SetAllParams(localTime, timeZone, julianDay, latLong[0], latLong[1], turbidity);

    if (verbose)
    {
        printf("Time: %g, timezone: %g, day: %g, latitude: %g, longitude: %g, turbidity: %g\n", localTime, timeZone, julianDay, latLong[0], latLong[1], turbidity);

        float theta, phi;
        sunSky.SunThetaPhi(&theta, &phi);

        theta *= 180.0 / M_PI;
        phi   *= 180.0 / M_PI;

        theta = 90.0f - theta;
        phi   = 90.0f - phi;    // make relative to North rather than East, and clockwise.

        if (phi < 0.0f)
            phi += 360.0f;

        printf("Sun elevation, direction: %g, %g\n", theta, phi);
    }

    if (autoscale)
    {
        float horizLum = sunSky.SkyHorizonLuminance();

        printf("Horizon luminance: %g\n", horizLum);

        if (horizLum > 0.0f)
        {
            float lumScale = 8300.0f / horizLum;

            if (verbose)
                printf("Autoscaling luminance by: %g\n", lumScale);

            weight *= lumScale;
        }
    }

    if (verbose)
    {
        printf("Ouput: weight = %g, gamma = %g\n", weight, gamma);
    }

    uint32_t image   [256][256];
    Vec3f    imageHDR[256][256];

    char fileName[32];

    tga_image tga;
    init_tga_image(&tga, (uint8_t*) image, 256, 256, 32);
    tga.image_type = TGA_IMAGE_TYPE_BGR;
    tga.image_descriptor &= ~TGA_T_TO_B_BIT;

    if (!cubeMap)
    {
        SkyToHemisphere(sunSky, 256, 256, (uint8_t*) image, 1024, weight, gamma, hemiSign);

        snprintf(fileName, 32, "sky-hemi.tga");

        if (tga_write(fileName, &tga) == TGA_NOERR)
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);

        SkyToHemisphere(sunSky, 256, 256, imageHDR[0], weight, hemiSign);

        snprintf(fileName, 32, "sky-hemi.pfm");

        if (PFMWrite(fileName, 256, 256, imageHDR[0]))
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);
    }
    else
    {
        Vec3f imageHDR[6][256][256];

        for (int i = 0; i < 6; i++)
        {
            SkyToCubeFace(sunSky, i, 256, 256, (uint8_t*) image, 1024, weight, gamma);

            snprintf(fileName, 32, "sky-cube-%d.tga", i);

            if (tga_write(fileName, &tga) == TGA_NOERR)
                printf("wrote %s\n", fileName);
            else
                printf("failed to write %s\n", fileName);

            SkyToCubeFace(sunSky, i, 256, 256, imageHDR[i][0], weight);
        }

        snprintf(fileName, 32, "sky-cube.pfm");

        if (PFMWrite(fileName, 256, 256 * 6, imageHDR[0][0]))
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);
    }

    if (panoramic)
    {
        Vec3f imageHDR[256][512];

        SkyToPanoramic(sunSky, 256, imageHDR[0], weight);

        snprintf(fileName, 32, "sky-panoramic.pfm");

        if (PFMWrite(fileName, 512, 256, imageHDR[0]))
            printf("wrote %s\n", fileName);
        else
            printf("failed to write %s\n", fileName);
    }

    return 0;
}
