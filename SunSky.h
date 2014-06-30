/*
    File:       SunSky.h

    Function:   Implements a sky radiance function for a clear sky

    Author:     Andrew Willmott
                Utah sun/sky based on code by Brian Smits (bes@phoenix.cs.utah.edu)

    Copyright:    
*/

#ifndef SUNSKY_H
#define SUNSKY_H

#include "VL234f.h"

enum tSkyType
{
    kUtahClearSky,
    kUtahClearSkyTable,     ///< Table accelerated version, for sanity checking -- really designed for texture/shader use
    kUtahClearSkyBRDF,      ///< Sky model convolved with roughness (sat(cos^n) brdf)
    kCIEClearSky,
    kCIECloudySky,
    kCIEPartlyCloudySky,
    kNumSkyTypes
};

class cSunSky
{
public:
    cSunSky();

    void SetAllParams
    (
        float       timeOfDay,         // 24-hour, decimal: 0.0-23.99
        float       timeZone,          // relative to UTC: west -ve, east +ve
        int         julianDay,         // day of year: 1-365
        float       latitude,          // Degrees: north is positive
        float       longitude,         // Degrees: east is positive
        float       turbidity          // 2-30
    );

    void        SetTime(float timeOfDay);
    
    Vec3f       SunDirection() const;
    void        SunThetaPhi(float* theta, float* phi) const;    ///< theta = angle from zenith, phi = heading

    float       SunSolidAngle() const;

    float       SkyLuminance(const Vec3f &v) const;     ///< Returns the luminance of the sky in direction v. v must be normalized. Luminance is in Nits = cd/m^2 = lumens/sr/m^2 */
    Vec2f       SkyChroma   (const Vec3f &v) const;     ///< Returns the chroma of the sky in direction v. v must be normalized.
    Vec3f       SkyRGB      (const Vec3f &v) const;     ///< Returns luminance/chroma converted to RGB

    Vec3f       SunRGB() const; ///< Returns the RGB "luminance" of the sun. You must still map from Nits to an appropriate [0, 1] range. */

    void        SetSkyType(tSkyType skyType);
    tSkyType    SkyType() const;

    void        SetRoughness(float roughness);  ///< Set roughness for BRDF tables
    float       Roughness() const;

    // Direct access, though really, you want to be pulling routines out wholesale for non-research code.
    Vec3f       SkyRGBFromTable(const Vec3f& v) const;  ///< Use precalculated table to return fast sky colour
    Vec3f       ConvolvedSkyRGB(const Vec3f& v, float roughness) const; ///< return sky term convolved with roughness, 1 = fully diffuse

protected:
    void        FindSunDirection(); ///< Compute the local sun position from time values and longitude/latitude

    void        InitTimeConstants();
    void        InitTurbidityConstants();
    float       PerezFunction(const float lambdas[5], float theta, float phi, float lvz) const;
    
    Vec3f       ClearSkyRGB      (const Vec3f& v) const;
    float       ClearSkyLuminance(const Vec3f &v) const;
    Vec2f       ClearSkyChroma   (const Vec3f &v) const;

    // Table-based version - faster than per-sample Perez function evaluation
    void        FindThetaGammaTables();

    void        FindBRDFTables();

    float       CIEOvercastSkyLuminance(const Vec3f& v) const;
    float       CIEClearSkyLuminance(const Vec3f& v) const;
    float       CIEPartlyCloudySkyLuminance(const Vec3f& v) const;

    // Data
    bool        mInited;
    tSkyType    mSkyType;
    float       mLatitude;
    float       mLongitude;
    int         mJulianDay;
    float       mTimeOfDay;
    float       mStandardMeridian;
    float       mTurbidity;

    Vec3f       mToSun;
    float       mThetaS;
    float       mPhiS;
    float       mSunSolidAngle;

    Vec3f       mZenith;

    float       mPerez_x[5];
    float       mPerez_y[5];
    float       mPerez_Y[5];

    Vec3f       mPerezInvDen;

    Vec3f       mThetaTable[64];
    Vec3f       mGammaTable[64];

    float       mRoughness;
    Vec3f       mBRDFThetaTable[4][64];
    Vec3f       mBRDFGammaTable[4][64];
};

#endif 
