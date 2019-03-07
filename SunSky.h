//
//  File:       SunSky.h
//
//  Function:   Implements various sky models -- Preetham, Hosek, CIE* -- as
//              well as a separable table-based model for fast BRDF convolutions
//
//  Author:     Andrew Willmott, Preetham sun/sky based on code by Brian Smits,
//              Hosek code used L. Hosek & A. Wilkie code Version 1.4a as reference
//
//  Copyright:  Andrew Willmott, see HosekDataXYZ.h for Hosek data copyright
//

#ifndef SUNSKY_H
#define SUNSKY_H

#include "VL234f.h"

namespace SSLib
{
    extern const float kSunDiameter;
    extern const float kSunDistance;
    extern const float kSunCosAngle;
    extern const float kSunSolidAngle;

    Vec3f SunDirection
    (
        float       timeOfDay,         // 24-hour, decimal: 0.0-23.99
        float       timeZone,          // relative to UTC: west -ve, east +ve
        int         julianDay,         // day of year: 1-365
        float       latitude,          // Degrees: north is positive
        float       longitude          // Degrees: east is positive
    );
    // Returns the local sun direction at the given time/location. +Y = north, +X = east, +Z = up.

    // Utilities
    Vec3f SunRGB(float cosTheta, float turbidity = 3.0f);   // Returns RGB for given sun elevation
    float ZenithLuminance(float thetaS, float T);           // Returns luminance estimate for given solar altitude and turbidity

    float CIEOvercastSkyLuminance    (const Vec3f& v, float Lz);                        // CIE standard overcast sky
    float CIEClearSkyLuminance       (const Vec3f& v, const Vec3f& toSun, float Lz);    // CIE standard clear sky
    float CIEPartlyCloudySkyLuminance(const Vec3f& v, const Vec3f& toSun, float Lz);    // CIE standard partly cloudy sky
    float CIEStandardSky   (int type, const Vec3f& v, const Vec3f& toSun, float Lz);    // Returns one of 15 standard skies: type = 0-14. See kCIEStandardSkyCoeffs


    //--------------------------------------------------------------------------
    // cSunSkyPreetham
    //--------------------------------------------------------------------------

    class cSunSkyPreetham
    {
    public:
        cSunSkyPreetham();

        void        Update(const Vec3f& sun, float turbidity, float overcast = 0.0f, float hCrush = 0.0f); // update model with given settings

        Vec3f       SkyRGB      (const Vec3f &v) const;     // Returns luminance/chroma converted to RGB
        float       SkyLuminance(const Vec3f &v) const;     // Returns the luminance of the sky in direction v. v must be normalized. Luminance is in Nits = cd/m^2 = lumens/sr/m^2 */
        Vec2f       SkyChroma   (const Vec3f &v) const;     // Returns the chroma of the sky in direction v. v must be normalized.

        // Data
        Vec3f       mToSun;

        float       mPerez_x[5];
        float       mPerez_y[5];
        float       mPerez_Y[5];

        Vec3f       mZenith;
        Vec3f       mPerezInvDen;
    };


    //--------------------------------------------------------------------------
    // cSunSkyHosek
    //--------------------------------------------------------------------------

    class cSunSkyHosek
    {
    public:
        cSunSkyHosek();

        void        Update(const Vec3f& sun, float turbidity, Vec3f albedo = vl_0, float overcast = 0.0f); // update model with given settings

        Vec3f       SkyXYZ      (const Vec3f &v) const;     // Returns CIE XYZ
        Vec3f       SkyRGB      (const Vec3f &v) const;     // Returns luminance/chroma converted to RGB
        float       SkyLuminance(const Vec3f &v) const;     // Returns CIE XYZ

        // Data
        Vec3f       mToSun;
        float       mCoeffsXYZ[3][9];   // Hosek 9-term distribution coefficients
        Vec3f       mRadXYZ;            // Overall average radiance
    };


    //--------------------------------------------------------------------------
    // cSunSkyTable
    //--------------------------------------------------------------------------

    class cSunSkyTable
    {
    public:
        // Table-based version - faster than per-sample Perez/Hosek function evaluation, suitable for shader use via 64 x 2 texture
        // For a fixed time, Preetham can be expressed in the form
        //   K (1 + F(theta))(1 + G(gamma))
        // where theta is the zenith angle of v, and gamma the angle between v and the sun direction.
        // Hosek can be expressed in the form
        //   K (1 + F(theta))(1 + G(gamma) + H(theta))
        // where H is trivial to evaluate in a shader, involving a constant term and sqrt(v.z).
        // Note: the F term is generally negative, so we use F' = -F in the tables.

        void        FindThetaGammaTables(const cSunSkyPreetham& pt);
        void        FindThetaGammaTables(const cSunSkyHosek& pt);

        Vec3f       SkyRGB(const cSunSkyPreetham& pt, const Vec3f& v) const;  // Use precalculated table to return fast sky colour on CPU
        Vec3f       SkyRGB(const cSunSkyHosek& hk,    const Vec3f& v) const;  // Use precalculated table to return fast sky colour on CPU

        void        FillTexture(int width, int height, uint8_t image[][4]) const;  // Fill kTableSize x 2 BGRA8 texture with tables
        void        FillTexture(int width, int height, float   image[][4]) const;  // Fill kTableSize x 2 RGBAF32 texture with tables

        // Table acceleration
        enum { kTableSize = 64 };
        Vec3f       mThetaTable[kTableSize];
        Vec3f       mGammaTable[kTableSize];
        float       mMaxTheta = 1.0f;       // To avoid clipping when using non-float textures. Currently only necessary if overcast is being used.
        float       mMaxGamma = 1.0f;       // To avoid clipping when using non-float textures.
        bool        mXYZ      = false;      // Whether tables are storing xyY (Preetham) or XYZ (Hosek)
    };


    //--------------------------------------------------------------------------
    // cSunSkyBRDF
    //--------------------------------------------------------------------------

    // #define COMPACT_BRDF_TABLE   // enable for smaller half-size table

    class cSunSkyBRDF
    {
    public:
        // Extended version of cSunSkyTable that uses zonal harmonics to produce table rows
        // convolved with increasing cosine powers. This is an approximation, because
        // conv(AB) != conv(A) conv(B), but, because the Perez-form evaluation is
        // (1 + F(theta)) (1 + G(gamma)), the approximation is only for the small
        // order-2 FG term.
        void        FindBRDFTables(const cSunSkyTable& table, const cSunSkyPreetham& pt);
        void        FindBRDFTables(const cSunSkyTable& table, const cSunSkyHosek& hk);

        Vec3f       ConvolvedSkyRGB(const cSunSkyPreetham& pt, const Vec3f& v, float roughness) const; // return sky term convolved with roughness, 1 = fully diffuse
        Vec3f       ConvolvedSkyRGB(const cSunSkyHosek& pt,    const Vec3f& v, float roughness) const; // return sky term convolved with roughness, 1 = fully diffuse

        void        FillBRDFTexture(int width, int height, uint8_t image[][4]) const; // Fill kTableSize x (kBRDFSamples x 2|4) BGRA8 texture with tables
        void        FillBRDFTexture(int width, int height, float   image[][4]) const; // Fill kTableSize x (kBRDFSamples x 2|4) RGBAF32 texture with tables
                    // Note: for Hosek, the H term will be in the 'w' component of the theta section, and if a kBRDFSamples x 4 size texture is supplied,
                    // the two additional sections will contain the FH term table. (Using this improves accuracy but can be skipped.)

        enum { kTableSize = cSunSkyTable::kTableSize };
    #ifdef COMPACT_BRDF_TABLE
        enum { kBRDFSamples = 4 };
    #else
        enum { kBRDFSamples = 8 };
    #endif
        Vec3f       mBRDFThetaTable[kBRDFSamples][kTableSize];
        Vec3f       mBRDFGammaTable[kBRDFSamples][kTableSize];

        // Additional tables for 'H' term in Hosek.
        float       mBRDFThetaTableH [kBRDFSamples][kTableSize];
        Vec3f       mBRDFThetaTableFH[kBRDFSamples][kTableSize];
        bool        mHasHTerm = false;

        float       mMaxTheta = 1.0f;       // To avoid clipping when using non-float textures. Currently only necessary if overcast is being used.
        float       mMaxGamma = 1.0f;       // To avoid clipping when using non-float textures.
        bool        mXYZ      = false;      // Whether tables are storing xyY (Preetham) or XYZ (Hosek)
    };


    //--------------------------------------------------------------------------
    // cSunSky
    // Composite sun/sky model for easy comparisons
    //--------------------------------------------------------------------------

    enum tSkyType
    {
        kPreetham,
        kPreethamTable, // Table accelerated version, for sanity checking -- really designed for texture/shader use
        kPreethamBRDF,  // Sky model convolved with roughness (sat(cos^n) brdf)
        kHosek,
        kHosekTable,
        kHosekBRDF,
        kCIEClear,
        kCIEOvercast,
        kCIEPartlyCloudy,
        kNumSkyTypes
    };

    class cSunSky
    {
    public:
        cSunSky();

        void        SetSkyType(tSkyType skyType);
        tSkyType    SkyType() const;

        void        SetSunDir   (const Vec3f& sun);
        void        SetTurbidity(float turbidity);
        void        SetAlbedo   (Vec3f rgb);  // Set ground-bounce factor
        void        SetOvercast (float overcast);   // 0 = clear, 1 = completely overcast
        void        SetRoughness(float roughness);  // Set roughness for BRDF tables

        void        Update();                       // update model given above settings

        float       SkyLuminance(const Vec3f &v) const;     // Returns the luminance of the sky in direction v. v must be normalized. Luminance is in Nits = cd/m^2 = lumens/sr/m^2 */
        Vec2f       SkyChroma   (const Vec3f &v) const;     // Returns the chroma of the sky in direction v. v must be normalized.
        Vec3f       SkyRGB      (const Vec3f &v) const;     // Returns luminance/chroma converted to RGB

        float       AverageLuminance() const;

    protected:
        // Data
        tSkyType    mSkyType;

        Vec3f       mToSun;
        float       mTurbidity;
        Vec3f       mAlbedo;
        float       mOvercast;
        float       mRoughness;

        // Various models
        float           mZenithY;   // for CIE functions

        cSunSkyPreetham mPreetham;
        cSunSkyHosek    mHosek;
        cSunSkyTable    mTable;
        cSunSkyBRDF     mBRDF;
    };
}

#endif 
