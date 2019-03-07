// VL234 monolithic header


typedef float Real;

#ifdef VL_DEBUG
    #include <assert.h>

    #define CL_ASSERT(M_X) assert(M_X);
    #define CL_ASSERT_MSG(M_X, M_MSG) assert(M_X);
    #define CL_ERROR(M_MSG) assert(false);
    #define CL_WARNING(M_X)
    #define CL_WARNING_MSG(M_X, M_MSG)
    #define CL_EXPECT(M_X) assert(M_X)
    #define CL_EXPECT_MSG(M_X, M_MSG) assert(M_X)
    #define CL_RANGE(M_X, M_0, M_1) assert((M_X) >= M_0 && (M_X) < M_1)
    #define CL_RANGE_MSG(M_X, M_0, M_1, M_MSG) assert((M_X) >= M_0 && (M_X) < M_1)
    #define CL_INDEX(M_X, M_N) assert((M_X) >= 0 && (M_X) < M_N)
    #define CL_INDEX_MSG(M_X, M_N, M_MSG) assert((M_X) >= 0 && (M_X) < M_N)
#else
    #define CL_ASSERT(M_X)

    #define CL_ASSERT_MSG(M_X, M_MSG)

    #define CL_ERROR(M_X)

    #define CL_WARNING(M_X)
    #define CL_WARNING_MSG(M_X, M_MSG)
    #define CL_EXPECT(M_X) (void) (M_X)
    #define CL_EXPECT_MSG(M_X, M_MSG) (void) (M_X)

    #define CL_RANGE(M_X, M_0, M_1)

    #define CL_RANGE_MSG(M_X, M_0, M_1, M_MSG)
    #define CL_INDEX(M_X, M_N)
    #define CL_INDEX_MSG(M_X, M_N, M_MSG)
#endif


#ifndef VL_NEW
    #define VL_NEW new

    #define VL_DELETE delete

#endif
#ifndef VL_RESTRICT
    #define VL_RESTRICT
#endif

/*
    File:       VLConstants.h

    Function:   Contains various constants for VL.

    Copyright:  Andrew Willmott
*/

#ifndef VLConstants_H
#define VLConstants_H

#include <cmath>
#include <cstdint>


// --- Mathematical constants -------------------------------------------------

const Real          vl_pi        = Real(3.14159265358979323846);
const Real          vl_halfPi    = Real(vl_pi / 2.0);
const Real          vl_quarterPi = Real(vl_pi / 4.0);
const Real          vl_twoPi     = Real(2.0 * vl_pi);

#ifdef HUGE_VAL
    const Real vl_inf = Real(HUGE_VAL);
#endif

enum    ZeroOrOne   { vl_zero = 0, vl_0 = 0, vl_one = 1, vl_I = 1, vl_1 = 1 };
enum    Block       { vl_Z = 0, vl_B = 1, vl_block = 1 };
enum    Axis        { vl_x, vl_y, vl_z, vl_w };
typedef Axis        vl_axis;

const uint32_t      VL_REF_FLAG = uint32_t(1) << (sizeof(uint32_t) * 8 - 1);
const uint32_t      VL_REF_MASK = (~VL_REF_FLAG);
enum                { VL_SV_END = -1 };
const uint32_t      VL_SV_MAX_INDEX = (1 << 30);
const int           VL_MAX_STEPS = 10000;

#endif
/*
    File:           VLMath.h

    Function:       Various math definitions for VL
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_MATH_H
#define VL_MATH_H

#include <stdlib.h>

// --- Inlines ----------------------------------------------------------------

// additions to arithmetic functions

#ifdef VL_HAS_IEEEFP
    #include <ieeefp.h>
    #define vl_is_finite(X) finite(X)
#elif defined (__GNUC__) && defined(VL_USE_MISC)
    #define vl_is_finite(X) finite(X)
#else
    #define vl_is_finite(X) (1)
#endif

#ifdef VL_HAS_DRAND
    inline double vl_rand()
    { return drand48()); }
    inline double vl_srand(Uint seed)
    { srand48(seed); }
#else
    inline double vl_rand()
    { return rand() / (RAND_MAX + 1.0); }
    inline void vl_srand(uint32_t seed)
    { srand(seed); }
#endif

inline float len(float x)
{ return (fabsf(x)); }
inline double len(double x)
{ return (fabs(x)); }

template<class T_Value> inline T_Value sqr(T_Value x)
{ return x * x; }

template<class T_Value> inline T_Value cube(T_Value x)
{ return x * x * x; }

inline float sqrlen(float r)
{ return sqr(r); }
inline double sqrlen(double r)
{ return sqr(r); }

inline float lerp(float a, float b, float s)
{ return (1.0f - s) * a + s * b; }
inline double lerp(double a, double b, double s)
{ return (1.0 - s) * a + s * b; }

template<class T_Value> inline T_Value lerp(T_Value x, T_Value y, Real s)
{ return x + (y - x) * s; }
    
template<> inline float lerp<float>(float x, float y, Real s)
{ return float(x + (y - x) * s); }

template<class T> inline T vl_min(T a, T b)
{
    return a < b ? a : b;
}

template<class T> inline T vl_max(T a, T b)
{
    return a > b ? a : b;
}

inline double sign(double d)
{
    if (d < 0)
        return -1.0;
    else
        return +1.0;
}

inline double sinc(double x) 
{
    if (fabs(x) < 1.0e-6) 
        return 1.0;
    else 
        return sin(x) / x;
}

inline float sincf(float x) 
{
    if (fabsf(x) < 1.0e-6f) 
        return 1.0f;
    else 
        return sinf(x) / x;
}

#ifndef VL_HAVE_SINCOS
    inline void sincosf(double phi, double* sinv, double* cosv)
    {
        *sinv = sin(phi);
        *cosv = cos(phi);
    }
    inline void sincosf(float phi, float* sinv, float* cosv)
    {
        *sinv = sinf(phi);
        *cosv = cosf(phi);
    }
#endif

template<class T_Value> inline T_Value clip(T_Value x, T_Value min, T_Value max)
{
    if (x <= min)
        return min;
    else if (x >= max)
        return max;
    else
        return x;
}

// useful routines

inline void SetReal(float& a, double b)
{ a = float(b); }
inline void SetReal(double& a, double b)
{ a = b; }

template <class S, class T> inline void ConvertVec(const S& u, T& v)
{
    for (int i = 0; i < u.Elts(); i++)
        v[i] = u[i];
}

template <class T> inline void ConvertVec(const T& u, T& v)
{ v = u; }

template <class S, class T> inline void ConvertMat(const S& m, T& n)
{
    for (int i = 0; i < m.Rows(); i++)
        ConvertVec(m[i], n[i]);
}

template <class T> inline void ConvertMat(const T& m, T& n)
{ n = m; }


#endif

/*
    File:           VLVec2.h

    Function:       Defines a length-2 vector.
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_Vec2_H
#define VL_Vec2_H




// --- Vec2 Class -------------------------------------------------------------


class Vec2f
{
public:

    // Constructors
    Vec2f();
    Vec2f(float x, float y);      // (x, y)
    Vec2f(const Vec2f& v);          // Copy constructor
    Vec2f(ZeroOrOne k);             // v[i] = vl_zero
    Vec2f(Axis k);                  // v[k] = 1
    explicit Vec2f(float s);
    explicit Vec2f(const float v[]);

    // Accessor functions
    int           Elts() const { return 2; };   // Element count

    float&       operator [] (int i);          // Indexing by row
    const float& operator [] (int i) const;    // Indexing by row

    float*       Ref();                        // Return pointer to data
    const float* Ref() const;                  // Return pointer to data

    // Assignment operators
    Vec2f&    operator =  (const Vec2f& a);                      
    Vec2f&    operator =  (ZeroOrOne k);                     
    Vec2f&    operator =  (Axis k);                      
    
    Vec2f&    operator += (const Vec2f& a);
    Vec2f&    operator -= (const Vec2f& a);
    Vec2f&    operator *= (const Vec2f& a);
    Vec2f&    operator *= (float s);
    Vec2f&    operator /= (const Vec2f& a);
    Vec2f&    operator /= (float s);

    // Comparison operators
    bool     operator == (const Vec2f& a) const; // v == a?
    bool     operator != (const Vec2f& a) const; // v != a?

    // Arithmetic operators
    Vec2f    operator + (const Vec2f& a) const;  // v + a
    Vec2f    operator - (const Vec2f& a) const;  // v - a
    Vec2f    operator - () const;                // -v
    Vec2f    operator * (const Vec2f& a) const;  // v * a (vx * ax, ...)
    Vec2f    operator * (float s) const;        // v * s
    Vec2f    operator / (const Vec2f& a) const;  // v / a (vx / ax, ...)
    Vec2f    operator / (float s) const;        // v / s

    // Initialisers
    Vec2f&    MakeZero();                        // Zero vector
    Vec2f&    MakeUnit(int i, float k = vl_one);// I[i]
    Vec2f&    MakeBlock(float k = vl_one);      // All-k vector
    
    Vec2f&    Normalize();                       // normalize vector

    // Data
    float  x;
    float  y;
};


// --- Vec operators ----------------------------------------------------------

Vec2f    operator * (float s, const Vec2f& v);// s * v

float   dot      (const Vec2f& a, const Vec2f& b); // v . a
float   len      (const Vec2f& v);                 // || v ||
float   sqrlen   (const Vec2f& v);                 // v . v
Vec2f    norm     (const Vec2f& v);                 // v / || v ||
Vec2f    norm_safe(const Vec2f& v);                 // v / || v ||, handles || v || = 0
Vec2f    cross    (const Vec2f& v);                 // cross prod.
Vec2f    inv      (const Vec2f& v);                 // inverse: 1 / v
Vec2f    abs      (const Vec2f& v);                 // abs(v_i)

    
// --- Inlines ----------------------------------------------------------------


inline float& Vec2f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 2, "(Vec2::[i]) index out of range");
    return (&x)[i];
}

inline const float& Vec2f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 2, "(Vec2::[i]) index out of range");
    return (&x)[i];
}

inline Vec2f::Vec2f()
{
}

inline Vec2f::Vec2f(float x_in, float y_in) :
    x(x_in),
    y(y_in)
{}

inline Vec2f::Vec2f(const Vec2f& v) :
    x(v.x),
    y(v.y)
{}

inline Vec2f::Vec2f(float s) :
    x(s),
    y(s)
{
}

inline Vec2f::Vec2f(const float v[]) :
    x(v[0]),
    y(v[1])
{
}

inline float* Vec2f::Ref()
{
    return &x;
}

inline const float* Vec2f::Ref() const
{
    return &x;
}

inline Vec2f& Vec2f::operator = (const Vec2f& v)
{
    x = v.x;
    y = v.y;
    
    return *this;
}

inline Vec2f& Vec2f::operator += (const Vec2f& v)
{
    x += v.x;
    y += v.y;
    
    return *this;
}

inline Vec2f& Vec2f::operator -= (const Vec2f& v)
{
    x -= v.x;
    y -= v.y;
    
    return *this;
}

inline Vec2f& Vec2f::operator *= (const Vec2f& v)
{
    x *= v.x;
    y *= v.y;
    
    return *this;
}

inline Vec2f& Vec2f::operator *= (float s)
{
    x *= s;
    y *= s;
    
    return *this;
}

inline Vec2f& Vec2f::operator /= (const Vec2f& v)
{
    x /= v.x;
    y /= v.y;
    
    return *this;
}

inline Vec2f& Vec2f::operator /= (float s)
{
    x /= s;
    y /= s;
    
    return *this;
}

inline Vec2f Vec2f::operator + (const Vec2f& a) const
{
    Vec2f result;
    
    result.x = x + a.x;
    result.y = y + a.y;
    
    return result;
}

inline Vec2f Vec2f::operator - (const Vec2f& a) const
{
    Vec2f result;
    
    result.x = x - a.x;
    result.y = y - a.y;
    
    return result;
}

inline Vec2f Vec2f::operator - () const
{
    Vec2f result;
    
    result.x = -x;
    result.y = -y;
    
    return result;
}

inline Vec2f Vec2f::operator * (const Vec2f& a) const           
{
    Vec2f result;
    
    result.x = x * a.x;
    result.y = y * a.y;
    
    return result;
}

inline Vec2f Vec2f::operator * (float s) const
{
    Vec2f result;
    
    result.x = x * s;
    result.y = y * s;
    
    return result;
}

inline Vec2f operator * (float s, const Vec2f& v)
{
    return v * s;
}

inline Vec2f Vec2f::operator / (const Vec2f& a) const
{
    Vec2f result;
    
    result.x = x / a.x;
    result.y = y / a.y;
    
    return result;
}

inline Vec2f Vec2f::operator / (float s) const
{
    Vec2f result;
    
    result.x = x / s;
    result.y = y / s;
    
    return result;
}

inline float dot(const Vec2f& a, const Vec2f& b)
{
    return a.x * b.x + a.y * b.y;
}

inline Vec2f cross(const Vec2f& a)
{
    Vec2f result;
    
    result.x = -a.y;
    result.y =  a.x;
    
    return result;
}

inline float len(const Vec2f& v)
{
    return sqrt(dot(v, v));
}

inline float sqrlen(const Vec2f& v)
{
    return dot(v, v);
}

inline Vec2f norm(const Vec2f& v)   
{
    CL_ASSERT_MSG(sqrlen(v) > float(vl_zero), "normalising length-zero vector");
    return v / len(v);
}

inline Vec2f norm_safe(const Vec2f& v)   
{
    return v / (len(v) + float(1e-8));
}

inline Vec2f inv(const Vec2f& v)
{
    return Vec2f(float(1) / v.x, float(1) / v.y);
}

inline Vec2f abs(const Vec2f& v) 
{
    return Vec2f(abs(v.x), abs(v.y));
}

inline Vec2f& Vec2f::MakeUnit(int i, float k)
{
    if (i == 0)
    { x = k; y = vl_zero; }
    else if (i == 1)
    { x = vl_zero; y = k; }
    else 
        CL_ERROR("(Vec2::Unit) illegal unit vector");

    return *this;
}

inline Vec2f& Vec2f::MakeZero()
{
    x = vl_zero; y = vl_zero;
    return *this;
}

inline Vec2f& Vec2f::MakeBlock(float k)
{
    x = k; y = k;
    return *this;
}

inline Vec2f& Vec2f::Normalize()
{
    CL_ASSERT_MSG(sqrlen(*this) > float(vl_zero), "normalising length-zero vector");
    *this /= len(*this);
    return *this;
}


inline Vec2f::Vec2f(ZeroOrOne k) 
{
    x = float(k);
    y = float(k);
}

inline Vec2f::Vec2f(Axis k) 
{
    MakeUnit(k, vl_one);
}

inline Vec2f& Vec2f::operator = (ZeroOrOne k)
{
    x = float(k); y = float(k);
    
    return *this;
}

inline Vec2f& Vec2f::operator = (Axis k)
{
    MakeUnit(k, vl_1);
    
    return *this;
}

inline bool Vec2f::operator == (const Vec2f& a) const
{
    return x == a.x && y == a.y;
}

inline bool Vec2f::operator != (const Vec2f& a) const
{
    return x != a.x || y != a.y;
}

#endif

/*
    File:           VLVec3.h

    Function:       Defines a length-3 vector.
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_Vec3_H
#define VL_Vec3_H


 
 
// --- Vec3 Class -------------------------------------------------------------


class Vec2f;

class Vec3f
{
public:
    // Constructors
    Vec3f();
    Vec3f(float x, float y, float z);// [x, y, z]
    Vec3f(const Vec3f& v);              // Copy constructor
    Vec3f(const Vec2f& v, float w);    // Hom. 2D vector
    Vec3f(ZeroOrOne k);
    Vec3f(Axis a);
    explicit Vec3f(float s);
    explicit Vec3f(const float v[]);
    
    // Accessor functions
    int           Elts() const { return 3; };   // Element count

    float&       operator [] (int i);          // Indexing by row
    const float& operator [] (int i) const;    // Indexing by row

    float*       Ref();                        // Return pointer to data
    const float* Ref() const;                  // Return pointer to data

    // Assignment operators
    Vec3f&    operator =  (const Vec3f& a);
    Vec3f&    operator =  (ZeroOrOne k);
    Vec3f&    operator += (const Vec3f& a);
    Vec3f&    operator -= (const Vec3f& a);
    Vec3f&    operator *= (const Vec3f& a);
    Vec3f&    operator *= (float s);
    Vec3f&    operator /= (const Vec3f& a);
    Vec3f&    operator /= (float s);
    
    // Comparison operators
    bool      operator == (const Vec3f& a) const; // v == a?
    bool      operator != (const Vec3f& a) const; // v != a?
    bool      operator <  (const Vec3f& a) const; // v <  a?
    bool      operator >= (const Vec3f& a) const; // v >= a?

    // Arithmetic operators
    Vec3f     operator + (const Vec3f& a) const;  // v + a
    Vec3f     operator - (const Vec3f& a) const;  // v - a
    Vec3f     operator - () const;                // -v
    Vec3f     operator * (const Vec3f& a) const;  // v * a (vx * ax, ...)
    Vec3f     operator * (float s) const;        // v * s
    Vec3f     operator / (const Vec3f& a) const;  // v / a (vx / ax, ...)
    Vec3f     operator / (float s) const;        // v / s    

    // Initialisers
    Vec3f&    MakeZero();                        // Zero vector
    Vec3f&    MakeUnit(int i, float k = vl_one);// I[i]
    Vec3f&    MakeBlock(float k = vl_one);      // All-k vector

    Vec3f&    Normalize();                       // normalize vector

    // Conversion
    Vec2f&          AsVec2();
    const Vec2f&    AsVec2() const;

    // Data
    float x;
    float y;
    float z;
};


// --- Vec operators ----------------------------------------------------------

Vec3f    operator * (float s, const Vec3f& v); // s * v
float   dot(const Vec3f& a, const Vec3f& b);   // v . a

float   len      (const Vec3f& v);             // || v ||
float   sqrlen   (const Vec3f& v);             // v . v
Vec3f    norm     (const Vec3f& v);             // v / || v ||
Vec3f    norm_safe(const Vec3f& v);             // v / || v ||, handles || v || = 0
Vec3f    cross    (const Vec3f& a, const Vec3f& b); // a x b
Vec3f    cross_x  (const Vec3f& v);             // v x e_x
Vec3f    cross_y  (const Vec3f& v);             // v x e_y
Vec3f    cross_z  (const Vec3f& v);             // v x e_z
Vec3f    inv      (const Vec3f& v);             // inverse: 1 / v
Vec2f    proj     (const Vec3f& v);             // homogeneous projection
Vec3f    abs      (const Vec3f& v);             // abs(v_i)


// --- Inlines ----------------------------------------------------------------



inline float& Vec3f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 3, "(Vec3::[i]) index out of range");
    return (&x)[i];
}

inline const float& Vec3f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 3, "(Vec3::[i]) index out of range");
    return (&x)[i];
}

inline Vec3f::Vec3f()
{
}

inline Vec3f::Vec3f(float x_in, float y_in, float z_in) :
    x(x_in),
    y(y_in),
    z(z_in)
{}

inline Vec3f::Vec3f(const Vec3f& v)  :
    x(v.x),
    y(v.y),
    z(v.z)
{}

inline Vec3f::Vec3f(const Vec2f& v, float z_in) :
    x(v.x),
    y(v.y),
    z(z_in)
{}

inline Vec3f::Vec3f(float s) :
    x(s),
    y(s),
    z(s)
{}

inline Vec3f::Vec3f(const float v[]) :
    x(v[0]),
    y(v[1]),
    z(v[2])
{}

inline float* Vec3f::Ref()
{
    return &x;
}

inline const float* Vec3f::Ref() const
{
    return &x;
}

inline Vec3f& Vec3f::operator = (const Vec3f& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    
    return *this;
}

inline Vec3f& Vec3f::operator += (const Vec3f& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    
    return *this;
}

inline Vec3f& Vec3f::operator -= (const Vec3f& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    
    return *this;
}

inline Vec3f& Vec3f::operator *= (const Vec3f& a)
{
    x *= a.x;
    y *= a.y;
    z *= a.z;
    
    return *this;
}

inline Vec3f& Vec3f::operator *= (float s)
{
    x *= s;
    y *= s;
    z *= s;
    
    return *this;
}

inline Vec3f& Vec3f::operator /= (const Vec3f& a)
{
    x /= a.x;
    y /= a.y;
    z /= a.z;
    
    return *this;
}

inline Vec3f& Vec3f::operator /= (float s)
{
    x /= s;
    y /= s;
    z /= s;
    
    return *this;
}

inline Vec3f Vec3f::operator + (const Vec3f& a) const
{
    Vec3f result;
    
    result.x = x + a.x;
    result.y = y + a.y;
    result.z = z + a.z;
    
    return result;
}

inline Vec3f Vec3f::operator - (const Vec3f& a) const
{
    Vec3f result;
    
    result.x = x - a.x;
    result.y = y - a.y;
    result.z = z - a.z;
    
    return result;
}

inline Vec3f Vec3f::operator - () const
{
    Vec3f result;
    
    result.x = -x;
    result.y = -y;
    result.z = -z;
    
    return result;
}

inline Vec3f Vec3f::operator * (const Vec3f& a) const
{
    Vec3f result;
    
    result.x = x * a.x;
    result.y = y * a.y;
    result.z = z * a.z;
    
    return result;
}

inline Vec3f Vec3f::operator * (float s) const
{
    Vec3f result;
    
    result.x = x * s;
    result.y = y * s;
    result.z = z * s;
    
    return result;
}

inline Vec3f Vec3f::operator / (const Vec3f& a) const
{
    Vec3f result;
    
    result.x = x / a.x;
    result.y = y / a.y;
    result.z = z / a.z;
    
    return result;
}

inline Vec3f Vec3f::operator / (float s) const
{
    Vec3f result;
    
    result.x = x / s;
    result.y = y / s;
    result.z = z / s;
    
    return result;
}

inline Vec3f operator * (float s, const Vec3f& v)
{
    return v * s;
}

inline Vec3f& Vec3f::MakeUnit(int n, float k)
{
    switch (n)
    {
    case 0:
        { x = k; y = vl_zero; z = vl_zero; } break;
    case 1:
        { x = vl_zero; y = k; z = vl_zero; } break;
    case 2:
        { x = vl_zero; y = vl_zero; z = k; } break;
    default:
        CL_ERROR("(Vec3::Unit) illegal unit vector");
    }

    return *this;
}

inline Vec3f& Vec3f::MakeZero()
{
    x = vl_zero; y = vl_zero; z = vl_zero;
    return *this;
}

inline Vec3f& Vec3f::MakeBlock(float k)
{
    x = k; y = k; z = k;
    return *this;
}

inline Vec3f& Vec3f::Normalize()
{
    CL_ASSERT_MSG(sqrlen(*this) > vl_zero, "normalising length-zero vector");
    *this /= len(*this);
    return *this;
}


inline Vec3f::Vec3f(ZeroOrOne k) 
{
    x = float(k); y = float(k); z = float(k);
}

inline Vec3f& Vec3f::operator = (ZeroOrOne k)
{
    x = float(k); y = float(k); z = float(k);
    
    return *this;
}

inline Vec3f::Vec3f(Axis a)
{
    MakeUnit(a, vl_one);
}


inline bool Vec3f::operator == (const Vec3f& a) const
{
    return x == a.x && y == a.y && z == a.z;
}

inline bool Vec3f::operator != (const Vec3f& a) const
{
    return x != a.x || y != a.y || z != a.z;
}

inline bool Vec3f::operator < (const Vec3f& a) const
{
    return x < a.x && y < a.y && z < a.z;
}

inline bool Vec3f::operator >= (const Vec3f& a) const
{
    return x >= a.x && y >= a.y && z >= a.z;
}

inline Vec2f& Vec3f::AsVec2()
{
    return (Vec2f&) *this;
}

inline const Vec2f& Vec3f::AsVec2() const
{
    return (const Vec2f&) *this;
}


inline float dot(const Vec3f& a, const Vec3f& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline float len(const Vec3f& v)
{
    return sqrt(dot(v, v));
}

inline float sqrlen(const Vec3f& v)
{
    return dot(v, v);
}

inline Vec3f norm(const Vec3f& v)   
{
    CL_ASSERT_MSG(sqrlen(v) > float(vl_zero), "normalising length-zero vector");
    return v / len(v);
}

inline Vec3f norm_safe(const Vec3f& v)
{
    return v / (len(v) + float(1e-8));
}

inline Vec3f cross(const Vec3f& a, const Vec3f& b)
{
    Vec3f result;

    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;

    return result;
}

inline Vec3f cross_x(const Vec3f& v)
{ return Vec3f(float(vl_zero), v.z, -v.y); }

inline Vec3f cross_y(const Vec3f& v)
{ return Vec3f(-v.z, float(vl_zero), v.x); }

inline Vec3f cross_z(const Vec3f& v)
{ return Vec3f(v.y, -v.x, float(vl_zero)); }

inline Vec3f inv(const Vec3f& v)
{
    return Vec3f(float(1) / v.x, float(1) / v.y, float(1) / v.z);
}

inline Vec2f proj(const Vec3f& v) 
{
    Vec2f result;
    
    CL_ASSERT_MSG(v.z != float(vl_zero), "(Vec3/proj) last elt. is zero");
    
    result.x = v.x / v.z;
    result.y = v.y / v.z;
    
    return result;
}

inline Vec3f abs(const Vec3f& v) 
{
    return Vec3f(abs(v.x), abs(v.y), abs(v.z));
}

#endif

/*
    File:           VLVec4.h

    Function:       Defines a length-4 vector.
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_Vec4_H
#define VL_Vec4_H




// --- Vec4 Class -------------------------------------------------------------


class Vec2f;
class Vec3f;

class Vec4f
{
public:
    // Constructors
    Vec4f();
    Vec4f(float x, float y, float z, float w);  // [x, y, z, w]
    Vec4f(const Vec4f& v);                      // Copy constructor
    Vec4f(const Vec3f& v, float w);            // Hom. 3D vector
    Vec4f(ZeroOrOne k);
    Vec4f(Axis k);
    explicit Vec4f(float s);
    explicit Vec4f(const float v[]);
    
    // Accessor functions
    int           Elts() const { return 4; };   // Element count

    float&       operator [] (int i);          // Indexing by row
    const float& operator [] (int i) const;    // Indexing by row

    float*       Ref();                        // Return pointer to data
    const float* Ref() const;                  // Return pointer to data

    // Assignment operators
    Vec4f&        operator =  (const Vec4f& a);  
    Vec4f&        operator =  (ZeroOrOne k);
    Vec4f&        operator =  (Axis k);
    Vec4f&        operator += (const Vec4f& a);
    Vec4f&        operator -= (const Vec4f& a);
    Vec4f&        operator *= (const Vec4f& a);
    Vec4f&        operator *= (float s);
    Vec4f&        operator /= (const Vec4f& a);
    Vec4f&        operator /= (float s);
    
    // Comparison operators
    bool          operator == (const Vec4f& a) const; // v == a ?
    bool          operator != (const Vec4f& a) const; // v != a ?

    // Arithmetic operators
    Vec4f         operator + (const Vec4f& a) const;  // v + a
    Vec4f         operator - (const Vec4f& a) const;  // v - a
    Vec4f         operator - () const;                // -v
    Vec4f         operator * (const Vec4f& a) const;  // v * a (vx * ax, ...)
    Vec4f         operator * (float s) const;        // v * s
    Vec4f         operator / (const Vec4f& a) const;  // v / a (vx / ax, ...)
    Vec4f         operator / (float s) const;        // v / s

    // Initialisers
    Vec4f&        MakeZero();                        // Zero vector
    Vec4f&        MakeUnit(int i, float k = vl_one);// kI[i]
    Vec4f&        MakeBlock(float k = vl_one);      // All-k vector

    Vec4f&        Normalize();                       // normalize vector

    // Conversion
    Vec2f&          AsVec2();
    const Vec2f&    AsVec2() const;
    Vec3f&          AsVec3();
    const Vec3f&    AsVec3() const;

    // Data
    float x;
    float y;
    float z;
    float w;
};


// --- Vec operators ----------------------------------------------------------

Vec4f   operator * (float s, const Vec4f& v);  // Left mult. by s
float  dot  (const Vec4f& a, const Vec4f& b);  // v . a
Vec4f   cross(const Vec4f& a, const Vec4f& b, const Vec4f& c); // 4D cross prod.

float  len      (const Vec4f& v);              // || v ||
float  sqrlen   (const Vec4f& v);              // v . v
Vec4f   norm     (const Vec4f& v);              // v / || v ||
Vec4f   norm_safe(const Vec4f& v);              // v / || v ||, handles || v || = 0
Vec4f   inv      (const Vec4f& v);              // 1 / v
Vec3f   proj     (const Vec4f& v);              // hom. projection
Vec4f   abs      (const Vec4f& v);              // abs(v_i)

void    normalize(Vec4f& v);                    // v = norm(v)


// --- Inlines ----------------------------------------------------------------



inline float& Vec4f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 4, "(Vec4::[i]) index out of range");
    return (&x)[i];
}

inline const float& Vec4f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 4, "(Vec4::[i]) index out of range");
    return (&x)[i];
}


inline Vec4f::Vec4f()
{
}

inline Vec4f::Vec4f(float x_in, float y_in, float z_in, float w_in) :
    x(x_in),
    y(y_in),
    z(z_in),
    w(w_in)
{}

inline Vec4f::Vec4f(const Vec4f& v) :
    x(v.x),
    y(v.y),
    z(v.z),
    w(v.w)
{}

inline Vec4f::Vec4f(const Vec3f& v, float w_in) :
    x(v.x),
    y(v.y),
    z(v.z),
    w(w_in)
{
}

inline Vec4f::Vec4f(float s) :
    x(s),
    y(s),
    z(s),
    w(s)
{}

inline Vec4f::Vec4f(const float v[]) :
    x(v[0]),
    y(v[1]),
    z(v[2]),
    w(v[3])
{
}

inline float* Vec4f::Ref()
{
    return &x;
}

inline const float* Vec4f::Ref() const
{
    return &x;
}

inline Vec4f& Vec4f::operator = (const Vec4f& v)
{
    x = v.x;
    y = v.y;
    z = v.z;
    w = v.w;
    
    return *this;
}

inline Vec4f& Vec4f::operator += (const Vec4f& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;
    
    return *this;
}

inline Vec4f& Vec4f::operator -= (const Vec4f& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;
    
    return *this;
}

inline Vec4f& Vec4f::operator *= (const Vec4f& v)
{
    x *= v.x;
    y *= v.y;
    z *= v.z;
    w *= v.w;
    
    return *this;
}

inline Vec4f& Vec4f::operator *= (float s)
{
    x *= s;
    y *= s;
    z *= s;
    w *= s;
    
    return *this;
}

inline Vec4f& Vec4f::operator /= (const Vec4f& v)
{
    x /= v.x;
    y /= v.y;
    z /= v.z;
    w /= v.w;
    
    return *this;
}

inline Vec4f& Vec4f::operator /= (float s)
{
    x /= s;
    y /= s;
    z /= s;
    w /= s;
    
    return *this;
}


inline Vec4f Vec4f::operator + (const Vec4f& a) const
{
    Vec4f result;
    
    result.x = x + a.x;
    result.y = y + a.y;
    result.z = z + a.z;
    result.w = w + a.w;
    
    return result;
}

inline Vec4f Vec4f::operator - (const Vec4f& a) const
{
    Vec4f result;
    
    result.x = x - a.x;
    result.y = y - a.y;
    result.z = z - a.z;
    result.w = w - a.w;
    
    return result;
}

inline Vec4f Vec4f::operator - () const
{
    Vec4f result;
    
    result.x = -x;
    result.y = -y;
    result.z = -z;
    result.w = -w;
    
    return result;
}

inline Vec4f Vec4f::operator * (const Vec4f& a) const           
{
    Vec4f result;

    result.x = x * a.x;
    result.y = y * a.y;
    result.z = z * a.z;
    result.w = w * a.w;
    
    return result;
}

inline Vec4f Vec4f::operator * (float s) const
{
    Vec4f result;
    
    result.x = x * s;
    result.y = y * s;
    result.z = z * s;
    result.w = w * s;
    
    return result;
}

inline Vec4f Vec4f::operator / (const Vec4f& a) const
{
    Vec4f result;
    
    result.x = x / a.x;
    result.y = y / a.y;
    result.z = z / a.z;
    result.w = w / a.w;
    
    return result;
}

inline Vec4f Vec4f::operator / (float s) const
{
    Vec4f result;
    
    result.x = x / s;
    result.y = y / s;
    result.z = z / s;
    result.w = w / s;
    
    return result;
}

inline Vec4f operator * (float s, const Vec4f& v)
{
    return v * s;
}

inline Vec4f& Vec4f::MakeZero()
{
    x = vl_zero; y = vl_zero; z = vl_zero; w = vl_zero;
    return *this;
}

inline Vec4f& Vec4f::MakeBlock(float k)
{
    x = k; y = k; z = k; w = k;
    return *this;
}

inline Vec4f& Vec4f::Normalize()
{
    CL_ASSERT_MSG(sqrlen(*this) > vl_zero, "normalising length-zero vector");
    *this /= len(*this);
    return *this;
}

inline Vec4f::Vec4f(ZeroOrOne k)    
{
    MakeBlock(float(k));
}

inline Vec4f::Vec4f(Axis k)     
{
    MakeUnit(k, vl_1);
}

inline Vec4f& Vec4f::operator = (ZeroOrOne k)
{
    MakeBlock(float(k));
    
    return *this;
}

inline Vec4f& Vec4f::operator = (Axis k)
{
    MakeUnit(k, vl_1);
    
    return *this;
}

inline Vec2f& Vec4f::AsVec2()
{
    return (Vec2f&) *this;
}

inline const Vec2f& Vec4f::AsVec2() const
{
    return (const Vec2f&) *this;
}

inline Vec3f& Vec4f::AsVec3()
{
    return (Vec3f&) *this;
}

inline const Vec3f& Vec4f::AsVec3() const
{
    return (const Vec3f&) *this;
}



inline float dot(const Vec4f& a, const Vec4f& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline float len(const Vec4f& v)
{
    return sqrt(dot(v, v));
}

inline float sqrlen(const Vec4f& v)
{
    return dot(v, v);
}

inline Vec4f norm(const Vec4f& v)   
{
    CL_ASSERT_MSG(sqrlen(v) > vl_zero, "normalising length-zero vector");
    return v / len(v);
}

inline Vec4f norm_safe(const Vec4f& v)   
{
    return v / (len(v) + float(1e-8));
}

inline Vec4f inv(const Vec4f& v)
{
    return Vec4f(float(1) / v.x, float(1) / v.y, float(1) / v.z, float(1) / v.w);
}

inline Vec4f abs(const Vec4f& v)
{
    return Vec4f(abs(v.x), abs(v.y), abs(v.z), abs(v.w));
}


#endif

/*
    File:           VLMat2.h

    Function:       Defines a 2 x 2 matrix.
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_Mat2_H
#define VL_Mat2_H





// --- Mat2 Class -------------------------------------------------------------

class Mat2f 
{
public:
    // Constructors
    Mat2f();
    Mat2f(float a, float b, float c, float d);
    Mat2f(const Mat2f& m);
    Mat2f(const Vec2f& v0, const Vec2f& v1);
    Mat2f(ZeroOrOne k);
    Mat2f(Block k);
    
    // Accessor functions
    int           Elts() const { return 4; };
    int           Rows() const { return 2; };
    int           Cols() const { return 2; };
        
    Vec2f&       operator [] (int i);
    const Vec2f& operator [] (int i) const;

    float&       operator () (int i, int j);
    float        operator () (int i, int j) const;
    
    float*       Ref();
    const float* Ref() const;

    // Assignment operators
    Mat2f&        operator =  (const Mat2f& m); 
    Mat2f&        operator =  (ZeroOrOne k); 
    Mat2f&        operator =  (Block k); 
    Mat2f&        operator += (const Mat2f& m);     
    Mat2f&        operator -= (const Mat2f& m);      
    Mat2f&        operator *= (const Mat2f& m);      
    Mat2f&        operator *= (float s);                
    Mat2f&        operator /= (float s);                
    
    // Comparison operators    
    bool          operator == (const Mat2f& m) const; // M == N?
    bool          operator != (const Mat2f& m) const; // M != N?
    
    // Arithmetic operators    
    Mat2f         operator + (const Mat2f& m) const;  // M + N
    Mat2f         operator - (const Mat2f& m) const;  // M - N
    Mat2f         operator - () const;                // -M
    Mat2f         operator * (const Mat2f& m) const;  // M * N
    Mat2f         operator * (float s) const;        // M * s
    Mat2f         operator / (float s) const;        // M / s
    
    // Initialisers
    void          MakeZero();                         // Zero matrix
    void          MakeIdentity();                     // Identity matrix
    void          MakeDiag(float k = vl_one);        // I
    void          MakeBlock(float k = vl_one);       // all elts=k
    
    // Vector Transformations
    Mat2f&        MakeRot(Real theta);
    Mat2f&        MakeScale(const Vec2f& s);

    // Data
    Vec2f x;
    Vec2f y;
};


// --- Matrix operators -------------------------------------------------------


Vec2f&  operator *= (Vec2f& v, const Mat2f& m);       // v *= m
Vec2f   operator *  (const Mat2f& m, const Vec2f& v); // m * v
Vec2f   operator *  (const Vec2f& v, const Mat2f& m); // v * m
Mat2f    operator *  (float s, const Mat2f& m);        // s * m

Mat2f    trans(const Mat2f& m);              // Transpose
float   trace(const Mat2f& m);              // Trace
Mat2f    adj  (const Mat2f& m);              // Adjoint
float   det  (const Mat2f& m);              // Determinant
Mat2f    inv  (const Mat2f& m);              // Inverse
Mat2f    abs  (const Mat2f& m);              // abs(m_ij)
Mat2f    oprod(const Vec2f& a, const Vec2f& b); // Outer product

// The xform functions help avoid dependence on whether row or column
// vectors are used to represent points and vectors.
Vec2f    xform(const Mat2f& m, const Vec2f& v); // Transform of v by m
Mat2f    xform(const Mat2f& m, const Mat2f& n); // Xform v -> m(n(v))


// --- Inlines ----------------------------------------------------------------


inline Mat2f::Mat2f()
{
}

inline Mat2f::Mat2f(float a, float b, float c, float d)
{
    x.x = a;  x.y = b;
    y.x = c;  y.y = d;
}

inline Mat2f::Mat2f(const Mat2f& m)
{
    x = m.x;
    y = m.y;
}

inline Mat2f::Mat2f(const Vec2f& v0, const Vec2f& v1)
{
    x = v0;
    y = v1;
}

inline Vec2f& Mat2f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 2, "(Mat2::[i]) index out of range");
    return (&x)[i];
}

inline const Vec2f& Mat2f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 2, "(Mat2::[i]) index out of range");
    return (&x)[i];
}

inline float& Mat2f::operator () (int i, int j)
{
    CL_RANGE_MSG(i, 0, 2, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 2, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float Mat2f::operator () (int i, int j) const
{
    CL_RANGE_MSG(i, 0, 2, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 2, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float* Mat2f::Ref()
{
    return &x.x;
}

inline const float* Mat2f::Ref() const
{
    return &x.x;
}

inline void Mat2f::MakeZero()
{
    x.x = vl_zero; x.y = vl_zero;
    y.x = vl_zero; y.y = vl_zero;
}

inline void Mat2f::MakeDiag(float k)
{
    x.x = k;          x.y = vl_zero;
    y.x = vl_zero;    y.y = k;
}

inline void Mat2f::MakeIdentity()
{
    x.x = vl_one;     x.y = vl_zero;
    y.x = vl_zero;    y.y = vl_one;
}

inline void Mat2f::MakeBlock(float k)
{
    x.x = k; x.y = k;
    y.x = k; y.y = k;
}

inline Mat2f::Mat2f(ZeroOrOne k)
{
    MakeDiag(float(k));
}

inline Mat2f::Mat2f(Block k)
{
    MakeBlock(float(k));
}

inline Mat2f& Mat2f::operator = (ZeroOrOne k)
{
    MakeDiag(float(k));

    return *this;
}
      
inline Mat2f& Mat2f::operator = (Block k)
{
    MakeBlock(float(k));

    return *this;
}
      
inline Mat2f& Mat2f::operator = (const Mat2f& m)
{
    x = m.x;
    y = m.y;
    
    return *this;
}
      
inline Mat2f& Mat2f::operator += (const Mat2f& m)
{
    x += m.x;
    y += m.y;
    
    return *this;
}

inline Mat2f& Mat2f::operator -= (const Mat2f& m)
{
    x -= m.x;
    y -= m.y;

    return *this;
}

inline Mat2f& Mat2f::operator *= (const Mat2f& m)
{
    Vec2f t;
 
    t = x.x * m.x + x.y * m.y;
    y = y.x * m.x + y.y * m.y;
    x = t;

    return *this;
}

inline Mat2f& Mat2f::operator *= (float s)
{
    x *= s;
    y *= s;

    return *this;
}

inline Mat2f& Mat2f::operator /= (float s)
{
    x /= s;
    y /= s;

    return *this;
}


inline Mat2f Mat2f::operator + (const Mat2f& m) const
{
    Mat2f result;
    
    result.x = x + m.x;
    result.y = y + m.y;

    return result; 
}

inline Mat2f Mat2f::operator - (const Mat2f& m) const
{
    Mat2f result;
    
    result.x = x - m.x;
    result.y = y - m.y;

    return result; 
}

inline Mat2f Mat2f::operator - () const
{
    Mat2f result;
    
    result.x = -x;
    result.y = -y;

    return result; 
}

inline Mat2f Mat2f::operator * (const Mat2f& m) const
{
    Mat2f result;

    result.x.x = x.x * m.x.x + x.y * m.y.x;
    result.x.y = x.x * m.x.y + x.y * m.y.y;
    result.y.x = y.x * m.x.x + y.y * m.y.x;
    result.y.y = y.x * m.x.y + y.y * m.y.y;
    
    return result;
}

inline Mat2f Mat2f::operator * (float s) const
{
    Mat2f result;
    
    result.x = x * s;
    result.y = y * s;

    return result; 
}

inline Mat2f Mat2f::operator / (float s) const
{
    Mat2f result;
    
    result.x = x / s;
    result.y = y / s;

    return result; 
}

inline Mat2f  operator *  (float s, const Mat2f& m)
{
    return m * s;
}

inline Vec2f operator * (const Mat2f& m, const Vec2f& v)
{
    Vec2f result;

    result.x = m.x.x * v.x + m.x.y * v.y;
    result.y = m.y.x * v.x + m.y.y * v.y;

    return result;
}

inline Vec2f operator * (const Vec2f& v, const Mat2f& m)          
{
    Vec2f result;

    result.x = v.x * m.x.x + v.y * m.y.x;
    result.y = v.x * m.x.y + v.y * m.y.y;

    return result;
}

inline Vec2f& operator *= (Vec2f& v, const Mat2f& m)      
{
    float t;
    
    t   = v.x * m.x.x + v.y * m.y.x;
    v.y = v.x * m.x.y + v.y * m.y.y;
    v.x = t;

    return v;
}


inline Mat2f trans(const Mat2f& m)
{
    Mat2f result;

    result.x.x = m.x.x; result.x.y = m.y.x; 
    result.y.x = m.x.y; result.y.y = m.y.y; 
        
    return result;
}

inline float trace(const Mat2f& m)
{
    return m.x.x + m.y.y;
}
            
inline Mat2f adj(const Mat2f& m)            
{
    Mat2f result;

    result.x = -cross(m.y);
    result.y =  cross(m.x);
            
    return result;
}

#ifdef VL_ROW_ORIENT
inline Vec2f xform(const Mat2f& m, const Vec2f& v)
{ return v * m; }
inline Mat2f xform(const Mat2f& m, const Mat2f& n)
{ return n * m; }
#else
inline Vec2f xform(const Mat2f& m, const Vec2f& v)
{ return m * v; }
inline Mat2f xform(const Mat2f& m, const Mat2f& n)
{ return m * n; }
#endif

#endif

/*
    File:           VLMat3.h

    Function:       Defines a 3 x 3 matrix.
                    
    Copyright:      Andrew Willmott
*/

#ifndef VL_Mat3_H
#define VL_Mat3_H





// --- Mat3 Class -------------------------------------------------------------


class Vec4f;

class Mat3f 
{
public:
    // Constructors
    Mat3f();
    Mat3f(float a, float b, float c,
          float d, float e, float f,
          float g, float h, float i);
    Mat3f(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2);
    Mat3f(const Mat3f& m);
    explicit Mat3f(const Mat2f& m, float scale = float(vl_1));
    Mat3f(ZeroOrOne k);
    Mat3f(Block k);

    // Accessor functions
    int           Elts() const { return 9; };
    int           Rows() const { return 3; };
    int           Cols() const { return 3; };

    Vec3f&       operator [] (int i);   
    const Vec3f& operator [] (int i) const;

    float&       operator () (int i, int j);
    float        operator () (int i, int j) const;
    
    float*       Ref();
    const float* Ref() const;

    // Assignment operators
    Mat3f&        operator =  (const Mat3f& m);      
    Mat3f&        operator =  (ZeroOrOne k);     
    Mat3f&        operator =  (Block k);     
    Mat3f&        operator += (const Mat3f& m);      
    Mat3f&        operator -= (const Mat3f& m);     
    Mat3f&        operator *= (const Mat3f& m);     
    Mat3f&        operator *= (float s);                
    Mat3f&        operator /= (float s);       
    
    // Comparison operators
    bool         operator == (const Mat3f& m) const; // M == N?
    bool         operator != (const Mat3f& m) const; // M != N?
    
    // Arithmetic operators
    Mat3f         operator + (const Mat3f& m) const;  // M + N
    Mat3f         operator - (const Mat3f& m) const;  // M - N
    Mat3f         operator - () const;                // -M
    Mat3f         operator * (const Mat3f& m) const;  // M * N
    Mat3f         operator * (float s) const;        // M * s
    Mat3f         operator / (float s) const;        // M / s

    // Initialisers
    void          MakeZero();                         // Zero matrix
    void          MakeIdentity();                     // Identity matrix
    void          MakeDiag(float k = vl_one);        // I
    void          MakeBlock(float k = vl_one);       // all elts = k

    // Vector Transforms
    Mat3f&        MakeRot(const Vec3f& axis, Real theta);
    Mat3f&        MakeRot(const Vec4f& q);           // Rotate by quaternion
    Mat3f&        MakeScale(const Vec3f& s);

    // Homogeneous Transforms
    Mat3f&        MakeHRot(Real theta);           // Rotate by theta rads
    Mat3f&        MakeHScale(const Vec2f& s);    // Scale by s
    Mat3f&        MakeHTrans(const Vec2f& t);    // Translation by t
    
    // Data
    Vec3f x; 
    Vec3f y; 
    Vec3f z; 
};


// --- Matrix operators -------------------------------------------------------

Vec3f&  operator *= (Vec3f& v, const Mat3f& m);        // v *= m
Vec3f   operator * (const Mat3f& m, const Vec3f& v);   // m * v
Vec3f   operator * (const Vec3f& v, const Mat3f& m);   // v * m
Mat3f    operator * (const float s, const Mat3f& m);    // s * m

Mat3f    trans(const Mat3f& m);                   // Transpose            
float   trace(const Mat3f& m);                   // Trace
Mat3f    adj  (const Mat3f& m);                   // Adjoint
float   det  (const Mat3f& m);                   // Determinant
Mat3f    inv  (const Mat3f& m);                   // Inverse
Mat3f    abs  (const Mat3f& m);                   // abs(m_ij)
Mat3f    oprod(const Vec3f& a, const Vec3f& b); // Outer product

// The xform functions help avoid dependence on whether row or column
// vectors are used to represent points and vectors.
Vec3f    xform(const Mat3f& m, const Vec3f& v); // Transform of v by m
Vec2f    xform(const Mat3f& m, const Vec2f& v); // Hom. xform of v by m
Mat3f    xform(const Mat3f& m, const Mat3f& n); // Xform v -> m(n(v))


// --- Inlines ----------------------------------------------------------------

inline Mat3f::Mat3f()
{
}

inline Mat3f::Mat3f(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2)
{
    x = v0;
    y = v1;
    z = v2;
}

inline Mat3f::Mat3f(const Mat3f& m)
{
    x = m.x;
    y = m.y;
    z = m.z;
}

inline Vec3f& Mat3f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 3, "(Mat3::[i]) index out of range");
    return (&x)[i];
}

inline const Vec3f& Mat3f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 3, "(Mat3::[i]) index out of range");
    return (&x)[i];
}

inline float& Mat3f::operator () (int i, int j)
{
    CL_RANGE_MSG(i, 0, 3, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 3, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float Mat3f::operator () (int i, int j) const
{
    CL_RANGE_MSG(i, 0, 3, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 3, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float* Mat3f::Ref()
{
    return &x.x;
}

inline const float* Mat3f::Ref() const
{
    return &x.x;
}

inline void Mat3f::MakeZero()
{
    x.x = vl_zero; x.y = vl_zero; x.z = vl_zero;
    y.x = vl_zero; y.y = vl_zero; y.z = vl_zero;
    z.x = vl_zero; z.y = vl_zero; z.z = vl_zero;
}

inline void Mat3f::MakeDiag(float k)
{
    x.x = k;          x.y = vl_zero;    x.z = vl_zero;
    y.x = vl_zero;    y.y = k;          y.z = vl_zero;
    z.x = vl_zero;    z.y = vl_zero;    z.z = k;
}

inline void Mat3f::MakeIdentity()
{
    x.x = vl_one;     x.y = vl_zero;    x.z = vl_zero;
    y.x = vl_zero;    y.y = vl_one;     y.z = vl_zero;
    z.x = vl_zero;    z.y = vl_zero;    z.z = vl_one;
}

inline void Mat3f::MakeBlock(float k)
{
    x.x = k; x.y = k; x.z = k;
    y.x = k; y.y = k; y.z = k;
    z.x = k; z.y = k; z.z = k;
}


inline Mat3f::Mat3f(ZeroOrOne k)
{
    MakeDiag(float(k));
}

inline Mat3f::Mat3f(Block k)
{
    MakeBlock(float(k));
}

inline Mat3f& Mat3f::operator = (ZeroOrOne k)
{
    MakeDiag(float(k));

    return *this;
}
    
inline Mat3f& Mat3f::operator = (Block k)
{
    MakeBlock(float(k));

    return *this;
}
    
inline Mat3f operator *  (const float s, const Mat3f& m)
{
    return m * s;
}

inline Vec3f operator * (const Mat3f& m, const Vec3f& v)
{
    Vec3f result;
    
    result.x = v.x * m.x.x + v.y * m.x.y + v.z * m.x.z;
    result.y = v.x * m.y.x + v.y * m.y.y + v.z * m.y.z;
    result.z = v.x * m.z.x + v.y * m.z.y + v.z * m.z.z;
    
    return result;
}

inline Vec3f operator * (const Vec3f& v, const Mat3f& m)  
{
    Vec3f result;
    
    result.x = v.x * m.x.x + v.y * m.y.x + v.z * m.z.x;
    result.y = v.x * m.x.y + v.y * m.y.y + v.z * m.z.y;
    result.z = v.x * m.x.z + v.y * m.y.z + v.z * m.z.z;
    
    return result;
}

inline Vec3f& operator *= (Vec3f& v, const Mat3f& m)          
{   
    float t0, t1;
    
    t0  = v.x * m.x.x + v.y * m.y.x + v.z * m.z.x;
    t1  = v.x * m.x.y + v.y * m.y.y + v.z * m.z.y;
    v.z = v.x * m.x.z + v.y * m.y.z + v.z * m.z.z;
    v.x = t0;
    v.y = t1;

    return v;
}

#ifdef VL_ROW_ORIENT
inline Vec2f xform(const Mat3f& m, const Vec2f& v)
{ return proj(Vec3f(v, 1.0) * m); }
inline Vec3f xform(const Mat3f& m, const Vec3f& v)
{ return v * m; }
inline Mat3f xform(const Mat3f& m, const Mat3f& n)
{ return n * m; }
#else
inline Vec2f xform(const Mat3f& m, const Vec2f& v)
{ return proj(m * Vec3f(v, 1.0)); }
inline Vec3f xform(const Mat3f& m, const Vec3f& v)
{ return m * v; }
inline Mat3f xform(const Mat3f& m, const Mat3f& n)
{ return m * n; }
#endif

#endif


/*
    File:           VLMat4.h

    Function:       Defines a 4 x 4 matrix.
                    
    Copyright:      Andrew Willmott
 */

#ifndef VL_Mat4_H
#define VL_Mat4_H





// --- Mat4 Class -------------------------------------------------------------

class Vec3f;

class Mat4f 
{
public:
    
    // Constructors
    
    Mat4f();
    Mat4f(float a, float b, float c, float d,
          float e, float f, float g, float h,
          float i, float j, float k, float l,
          float m, float n, float o, float p);
    Mat4f(const Vec4f& v0, const Vec4f& v1, const Vec4f& v2, const Vec4f& v3);
    Mat4f(const Mat4f& m);
    explicit Mat4f(const Mat3f& m, float scale = float(vl_1));
    Mat4f(ZeroOrOne k);
    Mat4f(Block k);
                    
    // Accessor functions   
    int           Elts() const { return 16; };
    int           Rows() const { return  4; };
    int           Cols() const { return  4; };
    
    Vec4f&       operator [] (int i);
    const Vec4f& operator [] (int i) const;

    float&       operator () (int i, int j);
    float        operator () (int i, int j) const;
    
    float*       Ref();
    const float* Ref() const;

    // Assignment operators 
    Mat4f&       operator =  (const Mat4f& m);
    Mat4f&       operator =  (ZeroOrOne k);
    Mat4f&       operator =  (Block k);
    Mat4f&       operator += (const Mat4f& m);
    Mat4f&       operator -= (const Mat4f& m);
    Mat4f&       operator *= (const Mat4f& m);
    Mat4f&       operator *= (float s);
    Mat4f&       operator /= (float s);

    // Comparison operators 
    bool         operator == (const Mat4f& m) const; // M == N?
    bool         operator != (const Mat4f& m) const; // M != N?
    
    // Arithmetic operators 
    Mat4f        operator + (const Mat4f& m) const;  // M + N
    Mat4f        operator - (const Mat4f& m) const;  // M - N
    Mat4f        operator - () const;                // -M
    Mat4f        operator * (const Mat4f& m) const;  // M * N
    Mat4f        operator * (float s) const;        // M * s
    Mat4f        operator / (float s) const;        // M / s

    // Initialisers 
    void        MakeZero();                         // Zero matrix
    void        MakeIdentity();                     // Identity matrix
    void        MakeDiag(float k = vl_one);        // I
    void        MakeBlock(float k = vl_one);       // All elts = k

    // Homogeneous Transforms   
    Mat4f&      MakeHRot(const Vec3f& axis, Real theta); // Rotate by theta radians about axis   
    Mat4f&      MakeHRot(const Vec4f& q); // Rotate by quaternion   
    Mat4f&      MakeHScale(const Vec3f& s);    // Scale by components of s

    Mat4f&      MakeHTrans(const Vec3f& t);    // Translation by t

    Mat4f&      Transpose();                    // transpose in place
    Mat4f&      AddShift(const Vec3f& t);      // concat shift
    
    // Data
    Vec4f x;
    Vec4f y;
    Vec4f z;
    Vec4f w;
};


// --- Matrix operators -------------------------------------------------------

Vec4f   operator *  (const Mat4f& m, const Vec4f& v); // m * v
Vec4f   operator *  (const Vec4f& v, const Mat4f& m); // v * m
Vec4f&  operator *= (Vec4f& a, const Mat4f& m);       // v *= m
Mat4f    operator *  (float s, const Mat4f& m);        // s * m

Mat4f    trans(const Mat4f& m);              // Transpose            
float   trace(const Mat4f& m);              // Trace
Mat4f    adj  (const Mat4f& m);              // Adjoint
float   det  (const Mat4f& m);              // Determinant
Mat4f    inv  (const Mat4f& m);              // Inverse
Mat4f    oprod(const Vec4f& a, const Vec4f& b); // Outer product
Mat4f    abs  (const Mat4f& m);              // abs(m_ij)

// The xform functions help avoid dependence on whether row or column
// vectors are used to represent points and vectors.
Vec4f    xform(const Mat4f& m, const Vec4f& v); // Transform of v by m
Vec3f    xform(const Mat4f& m, const Vec3f& v); // Hom. xform of v by m
Mat4f    xform(const Mat4f& m, const Mat4f& n); // Xform v -> m(n(v))
    

// --- Inlines ----------------------------------------------------------------

inline Mat4f::Mat4f()
{
}

inline Mat4f::Mat4f(const Vec4f& v0, const Vec4f& v1, const Vec4f& v2, const Vec4f& v3)
{
    x = v0;
    y = v1;
    z = v2;
    w = v3;
}

inline Mat4f::Mat4f(const Mat4f& m)
{
    x = m.x;
    y = m.y;
    z = m.z;
    w = m.w;
}

inline Vec4f& Mat4f::operator [] (int i)
{
    CL_RANGE_MSG(i, 0, 4, "(Mat4::[i]) index out of range");
    return (&x)[i];
}

inline const Vec4f& Mat4f::operator [] (int i) const
{
    CL_RANGE_MSG(i, 0, 4, "(Mat4::[i]) index out of range");
    return (&x)[i];
}

inline float& Mat4f::operator () (int i, int j)
{
    CL_RANGE_MSG(i, 0, 4, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 4, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float Mat4f::operator () (int i, int j) const
{
    CL_RANGE_MSG(i, 0, 4, "(Mat2::(i,j)) index out of range");
    CL_RANGE_MSG(j, 0, 4, "(Mat2::(i,j)) index out of range");
    
    return (&x)[i][j];
}

inline float* Mat4f::Ref()
{
    return &x.x;
}

inline const float* Mat4f::Ref() const
{
    return &x.x;
}

inline Mat4f::Mat4f(ZeroOrOne k)
{
    MakeDiag(float(k));
}

inline Mat4f::Mat4f(Block k)
{
    MakeBlock(float(k));
}

inline Mat4f& Mat4f::operator = (ZeroOrOne k)
{
    MakeDiag(float(k));

    return *this;
}

inline Mat4f& Mat4f::operator = (Block k)
{
    MakeBlock(float(k));

    return *this;
}

inline Mat4f operator * (float s, const Mat4f& m)
{
    return m * s;
}

#ifdef VL_ROW_ORIENT
inline Vec3f xform(const Mat4f& m, const Vec3f& v)
{ return proj(Vec4f(v, 1.0) * m); }
inline Vec4f xform(const Mat4f& m, const Vec4f& v)
{ return v * m; }
inline Mat4f xform(const Mat4f& m, const Mat4f& n)
{ return n * m; }
#else
inline Vec3f xform(const Mat4f& m, const Vec3f& v)
{ return proj(m * Vec4f(v, 1.0)); }
inline Vec4f xform(const Mat4f& m, const Vec4f& v)
{ return m * v; }
inline Mat4f xform(const Mat4f& m, const Mat4f& n)
{ return m * n; }
#endif

#endif
