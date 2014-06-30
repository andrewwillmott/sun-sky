#include <stdint.h>

typedef float Real;

    #define CL_ASSERT(M_X)
    #define CL_ASSERT_MSG(M_X, M_MSG)
    #define CL_RANGE_MSG(M_X, M_0, M_1, M_MSG)
    #define CL_ERROR(M_X)
    #define CL_WARNING(M_X)
#include <stdlib.h>
#include <math.h>



typedef float Real;





    const Real vl_pi = Real(3.14159265358979323846);
    const Real vl_halfPi = Real(vl_pi / 2.0);


    const Real vl_twoPi = Real(2) * vl_pi;





enum ZeroOrOne { vl_zero = 0, vl_0 = 0, vl_one = 1, vl_I = 1, vl_1 = 1 };
enum Block { vl_Z = 0, vl_B = 1, vl_block = 1 };
enum Axis { vl_x, vl_y, vl_z, vl_w };
typedef Axis vl_axis;

const uint32_t VL_REF_FLAG = uint32_t(1) << (sizeof(uint32_t) * 8 - 1);
const uint32_t VL_REF_MASK = (~VL_REF_FLAG);
enum { VL_SV_END = -1 };
const uint32_t VL_SV_MAX_INDEX = (1 << 30);
const int VL_MAX_STEPS = 10000;
    inline double vl_rand()
    { return(rand() / (RAND_MAX + 1.0)); }
    inline void vl_srand(uint32_t seed)
    { srand(seed); }
inline float len(float x)
{ return (fabsf(x)); }
inline double len(double x)
{ return (fabs(x)); }

template<class T_Value> inline T_Value sqr(T_Value x)
{ return(x * x); }

inline float sqrlen(float r)
{ return(sqr(r)); }
inline double sqrlen(double r)
{ return(sqr(r)); }

inline float lerp(float a, float b, float s)
{ return((1.0f - s) * a + s * b); }
inline double lerp(double a, double b, double s)
{ return((1.0 - s) * a + s * b); }

template<class T_Value>
    inline T_Value lerp(T_Value x, T_Value y, Real s)
    { return(x + (y - x) * s); }

template<>
    inline float lerp<float>(float x, float y, Real s)
    { return(float(x + (y - x) * s)); }

inline double sign(double d)
{
    if (d < 0)
        return(-1.0);
    else
        return(1.0);
}

inline double sinc(double x)
{
  if (fabs(x) < 1.0e-6)
    return 1.0;
  else
    return(sin(x) / x);
}

inline float sincf(float x)
{
  if (fabsf(x) < 1.0e-6f)
    return 1.0f;
  else
    return(sinf(x) / x);
}


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


template<class T_Value> inline T_Value clip(T_Value x, T_Value min, T_Value max)
{
    if (x <= min)
        return(min);
    else if (x >= max)
        return(max);
    else
        return(x);
}



inline void SetReal(float &a, double b)
{ a = float(b); }
inline void SetReal(double &a, double b)
{ a = b; }

template <class S, class T> inline void ConvertVec(const S &u, T &v)
{
    for (int i = 0; i < u.Elts(); i++)
        v[i] = u[i];
}

template <class T> inline void ConvertVec(const T &u, T &v)
{ v = u; }

template <class S, class T> inline void ConvertMat(const S &m, T &n)
{
    for (int i = 0; i < m.Rows(); i++)
        ConvertVec(m[i], n[i]);
}

template <class T> inline void ConvertMat(const T &m, T &n)
{ n = m; }





class Vec2f
{
public:


    Vec2f();
    Vec2f(float x, float y);
    Vec2f(const Vec2f& v);
    Vec2f(ZeroOrOne k);
    Vec2f(Axis k);
    explicit Vec2f(float s);
    explicit Vec2f(const float v[]);


    float& operator [] (int i);
    const float& operator [] (int i) const;

    int Elts() const { return(2); };
    float* Ref() const;


    Vec2f& operator = (const Vec2f& a);
    Vec2f& operator = (ZeroOrOne k);
    Vec2f& operator = (Axis k);

    Vec2f& operator += (const Vec2f& a);
    Vec2f& operator -= (const Vec2f& a);
    Vec2f& operator *= (const Vec2f& a);
    Vec2f& operator *= (float s);
    Vec2f& operator /= (const Vec2f& a);
    Vec2f& operator /= (float s);


    bool operator == (const Vec2f& a) const;
    bool operator != (const Vec2f& a) const;


    Vec2f operator + (const Vec2f& a) const;
    Vec2f operator - (const Vec2f& a) const;
    Vec2f operator - () const;
    Vec2f operator * (const Vec2f& a) const;
    Vec2f operator * (float s) const;
    Vec2f operator / (const Vec2f& a) const;
    Vec2f operator / (float s) const;


    Vec2f& MakeZero();
    Vec2f& MakeUnit(int i, float k = vl_one);
    Vec2f& MakeBlock(float k = vl_one);

    Vec2f& Normalize();


    float elt[2];
};




Vec2f operator * (float s, const Vec2f& v);
float dot(const Vec2f& a, const Vec2f& b);
float len(const Vec2f& v);
float sqrlen(const Vec2f& v);
Vec2f norm(const Vec2f& v);
Vec2f norm_safe(const Vec2f& v);
void normalize(Vec2f& v);
Vec2f cross(const Vec2f& v);





inline float& Vec2f::operator [] (int i)
{
    ;
    return(elt[i]);
}

inline const float& Vec2f::operator [] (int i) const
{
    ;
    return(elt[i]);
}

inline Vec2f::Vec2f()
{
}

inline Vec2f::Vec2f(float x, float y)
{
    elt[0] = x;
    elt[1] = y;
}

inline Vec2f::Vec2f(const Vec2f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];
}

inline Vec2f::Vec2f(float s)
{
    elt[0] = s;
    elt[1] = s;
}

inline Vec2f::Vec2f(const float v[])
{
    elt[0] = v[0];
    elt[1] = v[1];
}

inline float* Vec2f::Ref() const
{
    return((float *) elt);
}

inline Vec2f& Vec2f::operator = (const Vec2f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];

    return(*this);
}

inline Vec2f& Vec2f::operator += (const Vec2f& v)
{
    elt[0] += v[0];
    elt[1] += v[1];

    return(*this);
}

inline Vec2f& Vec2f::operator -= (const Vec2f& v)
{
    elt[0] -= v[0];
    elt[1] -= v[1];

    return(*this);
}

inline Vec2f& Vec2f::operator *= (const Vec2f& v)
{
    elt[0] *= v[0];
    elt[1] *= v[1];

    return(*this);
}

inline Vec2f& Vec2f::operator *= (float s)
{
    elt[0] *= s;
    elt[1] *= s;

    return(*this);
}

inline Vec2f& Vec2f::operator /= (const Vec2f& v)
{
    elt[0] /= v[0];
    elt[1] /= v[1];

    return(*this);
}

inline Vec2f& Vec2f::operator /= (float s)
{
    elt[0] /= s;
    elt[1] /= s;

    return(*this);
}

inline Vec2f Vec2f::operator + (const Vec2f& a) const
{
    Vec2f result;

    result[0] = elt[0] + a[0];
    result[1] = elt[1] + a[1];

    return(result);
}

inline Vec2f Vec2f::operator - (const Vec2f& a) const
{
    Vec2f result;

    result[0] = elt[0] - a[0];
    result[1] = elt[1] - a[1];

    return(result);
}

inline Vec2f Vec2f::operator - () const
{
    Vec2f result;

    result[0] = -elt[0];
    result[1] = -elt[1];

    return(result);
}

inline Vec2f Vec2f::operator * (const Vec2f& a) const
{
    Vec2f result;

    result[0] = elt[0] * a[0];
    result[1] = elt[1] * a[1];

    return(result);
}

inline Vec2f Vec2f::operator * (float s) const
{
    Vec2f result;

    result[0] = elt[0] * s;
    result[1] = elt[1] * s;

    return(result);
}

inline Vec2f operator * (float s, const Vec2f& v)
{
    return(v * s);
}

inline Vec2f Vec2f::operator / (const Vec2f& a) const
{
    Vec2f result;

    result[0] = elt[0] / a[0];
    result[1] = elt[1] / a[1];

    return(result);
}

inline Vec2f Vec2f::operator / (float s) const
{
    Vec2f result;

    result[0] = elt[0] / s;
    result[1] = elt[1] / s;

    return(result);
}

inline float dot(const Vec2f& a, const Vec2f& b)
{
    return(a[0] * b[0] + a[1] * b[1]);
}

inline Vec2f cross(const Vec2f& a)
{
    Vec2f result;

    result[0] = -a[1];
    result[1] = a[0];

    return(result);
}

inline float len(const Vec2f& v)
{
    return(sqrt(dot(v, v)));
}

inline float sqrlen(const Vec2f& v)
{
    return(dot(v, v));
}

inline Vec2f norm(const Vec2f& v)
{
    ;
    return(v / len(v));
}

inline Vec2f norm_safe(const Vec2f& v)
{
    return(v / (len(v) + float(1e-8)));
}

inline void normalize(Vec2f& v)
{
    v /= len(v);
}

inline Vec2f& Vec2f::MakeUnit(int i, float k)
{
    if (i == 0)
    { elt[0] = k; elt[1] = vl_zero; }
    else if (i == 1)
    { elt[0] = vl_zero; elt[1] = k; }
    else
        ;
    return(*this);
}

inline Vec2f& Vec2f::MakeZero()
{
    elt[0] = vl_zero; elt[1] = vl_zero;
    return(*this);
}

inline Vec2f& Vec2f::MakeBlock(float k)
{
    elt[0] = k; elt[1] = k;
    return(*this);
}

inline Vec2f& Vec2f::Normalize()
{
    ;
    *this /= len(*this);
    return(*this);
}


inline Vec2f::Vec2f(ZeroOrOne k)
{
    elt[0] = float(k);
    elt[1] = float(k);
}

inline Vec2f::Vec2f(Axis k)
{
    MakeUnit(k, vl_one);
}

inline Vec2f& Vec2f::operator = (ZeroOrOne k)
{
    elt[0] = float(k); elt[1] = float(k);

    return(*this);
}

inline Vec2f& Vec2f::operator = (Axis k)
{
    MakeUnit(k, vl_1);

    return(*this);
}

inline bool Vec2f::operator == (const Vec2f& a) const
{
    return(elt[0] == a[0] && elt[1] == a[1]);
}

inline bool Vec2f::operator != (const Vec2f& a) const
{
    return(elt[0] != a[0] || elt[1] != a[1]);
}
class Vec2f;

class Vec3f
{
public:

    Vec3f();
    Vec3f(float x, float y, float z);
    Vec3f(const Vec3f& v);
    Vec3f(const Vec2f& v, float w);
    Vec3f(ZeroOrOne k);
    Vec3f(Axis a);
    explicit Vec3f(float s);
    explicit Vec3f(const float v[]);


    float& operator [] (int i);
    const float& operator [] (int i) const;

    int Elts() const { return(3); };
    float* Ref() const;


    Vec3f& operator = (const Vec3f& a);
    Vec3f& operator = (ZeroOrOne k);
    Vec3f& operator += (const Vec3f& a);
    Vec3f& operator -= (const Vec3f& a);
    Vec3f& operator *= (const Vec3f& a);
    Vec3f& operator *= (float s);
    Vec3f& operator /= (const Vec3f& a);
    Vec3f& operator /= (float s);


    bool operator == (const Vec3f& a) const;
    bool operator != (const Vec3f& a) const;
    bool operator < (const Vec3f& a) const;
    bool operator >= (const Vec3f& a) const;


    Vec3f operator + (const Vec3f& a) const;
    Vec3f operator - (const Vec3f& a) const;
    Vec3f operator - () const;
    Vec3f operator * (const Vec3f& a) const;
    Vec3f operator * (float s) const;
    Vec3f operator / (const Vec3f& a) const;
    Vec3f operator / (float s) const;


    Vec3f& MakeZero();
    Vec3f& MakeUnit(int i, float k = vl_one);
    Vec3f& MakeBlock(float k = vl_one);

    Vec3f& Normalize();


    float elt[3];
};




Vec3f operator * (float s, const Vec3f& v);
float dot(const Vec3f& a, const Vec3f& b);

float len (const Vec3f& v);
float sqrlen (const Vec3f& v);
Vec3f norm (const Vec3f& v);
Vec3f norm_safe(const Vec3f& v);
Vec3f cross (const Vec3f& a, const Vec3f& b);
Vec3f cross_x (const Vec3f& v);
Vec3f cross_y (const Vec3f& v);
Vec3f cross_z (const Vec3f& v);
Vec2f proj (const Vec3f& v);

void normalize(Vec3f& v);






inline float& Vec3f::operator [] (int i)
{
    ;
    return(elt[i]);
}

inline const float& Vec3f::operator [] (int i) const
{
    ;
    return(elt[i]);
}

inline Vec3f::Vec3f()
{
}

inline Vec3f::Vec3f(float x, float y, float z)
{
    elt[0] = x;
    elt[1] = y;
    elt[2] = z;
}

inline Vec3f::Vec3f(const Vec3f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
}

inline Vec3f::Vec3f(const Vec2f& v, float w)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = w;
}

inline Vec3f::Vec3f(float s)
{
    elt[0] = s;
    elt[1] = s;
    elt[2] = s;
}

inline Vec3f::Vec3f(const float v[])
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
}

inline float* Vec3f::Ref() const
{
    return((float* ) elt);
}

inline Vec3f& Vec3f::operator = (const Vec3f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];

    return(*this);
}

inline Vec3f& Vec3f::operator += (const Vec3f& v)
{
    elt[0] += v[0];
    elt[1] += v[1];
    elt[2] += v[2];

    return(*this);
}

inline Vec3f& Vec3f::operator -= (const Vec3f& v)
{
    elt[0] -= v[0];
    elt[1] -= v[1];
    elt[2] -= v[2];

    return(*this);
}

inline Vec3f& Vec3f::operator *= (const Vec3f& a)
{
    elt[0] *= a[0];
    elt[1] *= a[1];
    elt[2] *= a[2];

    return(*this);
}

inline Vec3f& Vec3f::operator *= (float s)
{
    elt[0] *= s;
    elt[1] *= s;
    elt[2] *= s;

    return(*this);
}

inline Vec3f& Vec3f::operator /= (const Vec3f& a)
{
    elt[0] /= a[0];
    elt[1] /= a[1];
    elt[2] /= a[2];

    return(*this);
}

inline Vec3f& Vec3f::operator /= (float s)
{
    elt[0] /= s;
    elt[1] /= s;
    elt[2] /= s;

    return(*this);
}

inline Vec3f Vec3f::operator + (const Vec3f& a) const
{
    Vec3f result;

    result[0] = elt[0] + a[0];
    result[1] = elt[1] + a[1];
    result[2] = elt[2] + a[2];

    return(result);
}

inline Vec3f Vec3f::operator - (const Vec3f& a) const
{
    Vec3f result;

    result[0] = elt[0] - a[0];
    result[1] = elt[1] - a[1];
    result[2] = elt[2] - a[2];

    return(result);
}

inline Vec3f Vec3f::operator - () const
{
    Vec3f result;

    result[0] = -elt[0];
    result[1] = -elt[1];
    result[2] = -elt[2];

    return(result);
}

inline Vec3f Vec3f::operator * (const Vec3f& a) const
{
    Vec3f result;

    result[0] = elt[0] * a[0];
    result[1] = elt[1] * a[1];
    result[2] = elt[2] * a[2];

    return(result);
}

inline Vec3f Vec3f::operator * (float s) const
{
    Vec3f result;

    result[0] = elt[0] * s;
    result[1] = elt[1] * s;
    result[2] = elt[2] * s;

    return(result);
}

inline Vec3f Vec3f::operator / (const Vec3f& a) const
{
    Vec3f result;

    result[0] = elt[0] / a[0];
    result[1] = elt[1] / a[1];
    result[2] = elt[2] / a[2];

    return(result);
}

inline Vec3f Vec3f::operator / (float s) const
{
    Vec3f result;

    result[0] = elt[0] / s;
    result[1] = elt[1] / s;
    result[2] = elt[2] / s;

    return(result);
}

inline Vec3f operator * (float s, const Vec3f& v)
{
    return(v * s);
}

inline Vec3f& Vec3f::MakeUnit(int n, float k)
{
    switch (n)
    {
    case 0:
        { elt[0] = k; elt[1] = vl_zero; elt[2] = vl_zero; } break;
    case 1:
        { elt[0] = vl_zero; elt[1] = k; elt[2] = vl_zero; } break;
    case 2:
        { elt[0] = vl_zero; elt[1] = vl_zero; elt[2] = k; } break;
    default:
        ;
    }
    return(*this);
}

inline Vec3f& Vec3f::MakeZero()
{
    elt[0] = vl_zero; elt[1] = vl_zero; elt[2] = vl_zero;
    return(*this);
}

inline Vec3f& Vec3f::MakeBlock(float k)
{
    elt[0] = k; elt[1] = k; elt[2] = k;
    return(*this);
}

inline Vec3f& Vec3f::Normalize()
{
    ;
    *this /= len(*this);
    return(*this);
}


inline Vec3f::Vec3f(ZeroOrOne k)
{
    elt[0] = float(k); elt[1] = float(k); elt[2] = float(k);
}

inline Vec3f& Vec3f::operator = (ZeroOrOne k)
{
    elt[0] = float(k); elt[1] = float(k); elt[2] = float(k);

    return(*this);
}

inline Vec3f::Vec3f(Axis a)
{
    MakeUnit(a, vl_one);
}


inline bool Vec3f::operator == (const Vec3f& a) const
{
    return(elt[0] == a[0] && elt[1] == a[1] && elt[2] == a[2]);
}

inline bool Vec3f::operator != (const Vec3f& a) const
{
    return(elt[0] != a[0] || elt[1] != a[1] || elt[2] != a[2]);
}

inline bool Vec3f::operator < (const Vec3f& a) const
{
    return(elt[0] < a[0] && elt[1] < a[1] && elt[2] < a[2]);
}

inline bool Vec3f::operator >= (const Vec3f& a) const
{
    return(elt[0] >= a[0] && elt[1] >= a[1] && elt[2] >= a[2]);
}


inline float dot(const Vec3f& a, const Vec3f& b)
{
    return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline float len(const Vec3f& v)
{
    return(sqrt(dot(v, v)));
}

inline float sqrlen(const Vec3f& v)
{
    return(dot(v, v));
}

inline Vec3f norm(const Vec3f& v)
{
    ;
    return(v / len(v));
}

inline Vec3f norm_safe(const Vec3f& v)
{
    return(v / (len(v) + float(1e-8)));
}

inline void normalize(Vec3f& v)
{
    v /= len(v);
}

inline Vec3f cross(const Vec3f& a, const Vec3f& b)
{
    Vec3f result;

    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];

    return(result);
}

inline Vec3f cross_x(const Vec3f& v)
{ return Vec3f(0, v[2], -v[1]); }

inline Vec3f cross_y(const Vec3f& v)
{ return Vec3f(-v[2], 0, v[0]); }

inline Vec3f cross_z(const Vec3f& v)
{ return Vec3f(v[1], -v[0], 0); }

inline Vec2f proj(const Vec3f& v)
{
    Vec2f result;

    ;

    result[0] = v[0] / v[2];
    result[1] = v[1] / v[2];

    return(result);
}
class Vec3f;

class Vec4f
{
public:

    Vec4f();
    Vec4f(float x, float y, float z, float w);
    Vec4f(const Vec4f& v);
    Vec4f(const Vec3f& v, float w);
    Vec4f(ZeroOrOne k);
    Vec4f(Axis k);
    explicit Vec4f(float s);
    explicit Vec4f(const float v[]);


    float& operator [] (int i);
    const float& operator [] (int i) const;

    int Elts() const { return(4); };
    float* Ref() const;


    Vec4f& operator = (const Vec4f& a);
    Vec4f& operator = (ZeroOrOne k);
    Vec4f& operator = (Axis k);
    Vec4f& operator += (const Vec4f& a);
    Vec4f& operator -= (const Vec4f& a);
    Vec4f& operator *= (const Vec4f& a);
    Vec4f& operator *= (float s);
    Vec4f& operator /= (const Vec4f& a);
    Vec4f& operator /= (float s);


    bool operator == (const Vec4f& a) const;
    bool operator != (const Vec4f& a) const;


    Vec4f operator + (const Vec4f& a) const;
    Vec4f operator - (const Vec4f& a) const;
    Vec4f operator - () const;
    Vec4f operator * (const Vec4f& a) const;
    Vec4f operator * (float s) const;
    Vec4f operator / (const Vec4f& a) const;
    Vec4f operator / (float s) const;


    Vec4f& MakeZero();
    Vec4f& MakeUnit(int i, float k = vl_one);
    Vec4f& MakeBlock(float k = vl_one);

    Vec4f& Normalize();


    float elt[4];
};




Vec4f operator * (float s, const Vec4f& v);
float dot (const Vec4f& a, const Vec4f& b);
Vec4f cross(const Vec4f& a, const Vec4f& b, const Vec4f& c);

float len (const Vec4f& v);
float sqrlen (const Vec4f& v);
Vec4f norm (const Vec4f& v);
Vec4f norm_safe(const Vec4f& v);
Vec3f proj (const Vec4f& v);

void normalize(Vec4f& v);






inline float& Vec4f::operator [] (int i)
{
    ;
    return(elt[i]);
}

inline const float& Vec4f::operator [] (int i) const
{
    ;
    return(elt[i]);
}


inline Vec4f::Vec4f()
{
}

inline Vec4f::Vec4f(float x, float y, float z, float w)
{
    elt[0] = x;
    elt[1] = y;
    elt[2] = z;
    elt[3] = w;
}

inline Vec4f::Vec4f(const Vec4f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
    elt[3] = v[3];
}

inline Vec4f::Vec4f(const Vec3f& v, float w)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
    elt[3] = w;
}

inline Vec4f::Vec4f(float s)
{
    elt[0] = s;
    elt[1] = s;
    elt[2] = s;
    elt[3] = s;
}

inline Vec4f::Vec4f(const float v[])
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
    elt[3] = v[3];
}

inline float* Vec4f::Ref() const
{
    return (float*) elt;
}

inline Vec4f& Vec4f::operator = (const Vec4f& v)
{
    elt[0] = v[0];
    elt[1] = v[1];
    elt[2] = v[2];
    elt[3] = v[3];

    return(*this);
}

inline Vec4f& Vec4f::operator += (const Vec4f& v)
{
    elt[0] += v[0];
    elt[1] += v[1];
    elt[2] += v[2];
    elt[3] += v[3];

    return(*this);
}

inline Vec4f& Vec4f::operator -= (const Vec4f& v)
{
    elt[0] -= v[0];
    elt[1] -= v[1];
    elt[2] -= v[2];
    elt[3] -= v[3];

    return(*this);
}

inline Vec4f& Vec4f::operator *= (const Vec4f& v)
{
    elt[0] *= v[0];
    elt[1] *= v[1];
    elt[2] *= v[2];
    elt[3] *= v[3];

    return(*this);
}

inline Vec4f& Vec4f::operator *= (float s)
{
    elt[0] *= s;
    elt[1] *= s;
    elt[2] *= s;
    elt[3] *= s;

    return(*this);
}

inline Vec4f& Vec4f::operator /= (const Vec4f& v)
{
    elt[0] /= v[0];
    elt[1] /= v[1];
    elt[2] /= v[2];
    elt[3] /= v[3];

    return(*this);
}

inline Vec4f& Vec4f::operator /= (float s)
{
    elt[0] /= s;
    elt[1] /= s;
    elt[2] /= s;
    elt[3] /= s;

    return(*this);
}


inline Vec4f Vec4f::operator + (const Vec4f& a) const
{
    Vec4f result;

    result[0] = elt[0] + a[0];
    result[1] = elt[1] + a[1];
    result[2] = elt[2] + a[2];
    result[3] = elt[3] + a[3];

    return(result);
}

inline Vec4f Vec4f::operator - (const Vec4f& a) const
{
    Vec4f result;

    result[0] = elt[0] - a[0];
    result[1] = elt[1] - a[1];
    result[2] = elt[2] - a[2];
    result[3] = elt[3] - a[3];

    return(result);
}

inline Vec4f Vec4f::operator - () const
{
    Vec4f result;

    result[0] = -elt[0];
    result[1] = -elt[1];
    result[2] = -elt[2];
    result[3] = -elt[3];

    return(result);
}

inline Vec4f Vec4f::operator * (const Vec4f& a) const
{
    Vec4f result;

    result[0] = elt[0] * a[0];
    result[1] = elt[1] * a[1];
    result[2] = elt[2] * a[2];
    result[3] = elt[3] * a[3];

    return(result);
}

inline Vec4f Vec4f::operator * (float s) const
{
    Vec4f result;

    result[0] = elt[0] * s;
    result[1] = elt[1] * s;
    result[2] = elt[2] * s;
    result[3] = elt[3] * s;

    return(result);
}

inline Vec4f Vec4f::operator / (const Vec4f& a) const
{
    Vec4f result;

    result[0] = elt[0] / a[0];
    result[1] = elt[1] / a[1];
    result[2] = elt[2] / a[2];
    result[3] = elt[3] / a[3];

    return(result);
}

inline Vec4f Vec4f::operator / (float s) const
{
    Vec4f result;

    result[0] = elt[0] / s;
    result[1] = elt[1] / s;
    result[2] = elt[2] / s;
    result[3] = elt[3] / s;

    return(result);
}

inline Vec4f operator * (float s, const Vec4f& v)
{
    return(v * s);
}

inline Vec4f& Vec4f::MakeZero()
{
    elt[0] = vl_zero; elt[1] = vl_zero; elt[2] = vl_zero; elt[3] = vl_zero;
    return(*this);
}

inline Vec4f& Vec4f::MakeBlock(float k)
{
    elt[0] = k; elt[1] = k; elt[2] = k; elt[3] = k;
    return(*this);
}

inline Vec4f& Vec4f::Normalize()
{
    ;
    *this /= len(*this);
    return(*this);
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

    return(*this);
}

inline Vec4f& Vec4f::operator = (Axis k)
{
    MakeUnit(k, vl_1);

    return(*this);
}


inline float dot(const Vec4f& a, const Vec4f& b)
{
    return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
}

inline float len(const Vec4f& v)
{
    return(sqrt(dot(v, v)));
}

inline float sqrlen(const Vec4f& v)
{
    return(dot(v, v));
}

inline Vec4f norm(const Vec4f& v)
{
    ;
    return(v / len(v));
}

inline Vec4f norm_safe(const Vec4f& v)
{
    return(v / (len(v) + float(1e-8)));
}

inline void normalize(Vec4f& v)
{
    v /= len(v);
}

class Mat2f
{
public:

    Mat2f();
    Mat2f(float a, float b, float c, float d);
    Mat2f(const Mat2f& m);
    Mat2f(ZeroOrOne k);
    Mat2f(Block k);


    int Rows() const { return(2); };
    int Cols() const { return(2); };

    Vec2f& operator [] (int i);
    const Vec2f& operator [] (int i) const;

    float& operator () (int i, int j);
    float operator () (int i, int j) const;

    float* Ref() const;


    Mat2f& operator = (const Mat2f& m);
    Mat2f& operator = (ZeroOrOne k);
    Mat2f& operator = (Block k);
    Mat2f& operator += (const Mat2f& m);
    Mat2f& operator -= (const Mat2f& m);
    Mat2f& operator *= (const Mat2f& m);
    Mat2f& operator *= (float s);
    Mat2f& operator /= (float s);


    bool operator == (const Mat2f& m) const;
    bool operator != (const Mat2f& m) const;


    Mat2f operator + (const Mat2f& m) const;
    Mat2f operator - (const Mat2f& m) const;
    Mat2f operator - () const;
    Mat2f operator * (const Mat2f& m) const;
    Mat2f operator * (float s) const;
    Mat2f operator / (float s) const;


    void MakeZero();
    void MakeIdentity();
    void MakeDiag(float k = vl_one);
    void MakeBlock(float k = vl_one);


    Mat2f& MakeRot(Real theta);
    Mat2f& MakeScale(const Vec2f& s);


    Vec2f row[2];
};





Vec2f& operator *= (Vec2f& v, const Mat2f& m);
Vec2f operator * (const Mat2f& m, const Vec2f& v);
Vec2f operator * (const Vec2f& v, const Mat2f& m);
Mat2f operator * (float s, const Mat2f& m);

Mat2f trans(const Mat2f& m);
float trace(const Mat2f& m);
Mat2f adj(const Mat2f& m);
float det(const Mat2f& m);
Mat2f inv(const Mat2f& m);
Mat2f oprod(const Vec2f& a, const Vec2f& b);



Vec2f xform(const Mat2f& m, const Vec2f& v);
Mat2f xform(const Mat2f& m, const Mat2f& n);





inline Vec2f& Mat2f::operator [] (int i)
{
    ;
    return(row[i]);
}

inline const Vec2f& Mat2f::operator [] (int i) const
{
    ;
    return(row[i]);
}

inline float& Mat2f::operator () (int i, int j)
{
    ;
    ;

    return(row[i][j]);
}

inline float Mat2f::operator () (int i, int j) const
{
    ;
    ;

    return(row[i][j]);
}

inline float* Mat2f::Ref() const
{
    return((float*) row);
}

inline Mat2f::Mat2f()
{
}

inline Mat2f::Mat2f(float a, float b, float c, float d)
{
    row[0][0] = a; row[0][1] = b;
    row[1][0] = c; row[1][1] = d;
}

inline Mat2f::Mat2f(const Mat2f& m)
{
    row[0] = m[0];
    row[1] = m[1];
}


inline void Mat2f::MakeZero()
{
    row[0][0] = vl_zero; row[0][1] = vl_zero;
    row[1][0] = vl_zero; row[1][1] = vl_zero;
}

inline void Mat2f::MakeDiag(float k)
{
    row[0][0] = k; row[0][1] = vl_zero;
    row[1][0] = vl_zero; row[1][1] = k;
}

inline void Mat2f::MakeIdentity()
{
    row[0][0] = vl_one; row[0][1] = vl_zero;
    row[1][0] = vl_zero; row[1][1] = vl_one;
}

inline void Mat2f::MakeBlock(float k)
{
    row[0][0] = k; row[0][1] = k;
    row[1][0] = k; row[1][1] = k;
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

    return(*this);
}

inline Mat2f& Mat2f::operator = (Block k)
{
    MakeBlock(float(k));

    return(*this);
}

inline Mat2f& Mat2f::operator = (const Mat2f& m)
{
    row[0] = m[0];
    row[1] = m[1];

    return(*this);
}

inline Mat2f& Mat2f::operator += (const Mat2f& m)
{
    row[0] += m[0];
    row[1] += m[1];

    return(*this);
}

inline Mat2f& Mat2f::operator -= (const Mat2f& m)
{
    row[0] -= m[0];
    row[1] -= m[1];

    return(*this);
}

inline Mat2f& Mat2f::operator *= (const Mat2f& m)
{
    *this = *this * m;

    return(*this);
}

inline Mat2f& Mat2f::operator *= (float s)
{
    row[0] *= s;
    row[1] *= s;

    return(*this);
}

inline Mat2f& Mat2f::operator /= (float s)
{
    row[0] /= s;
    row[1] /= s;

    return(*this);
}


inline Mat2f Mat2f::operator + (const Mat2f& m) const
{
    Mat2f result;

    result[0] = row[0] + m[0];
    result[1] = row[1] + m[1];

    return(result);
}

inline Mat2f Mat2f::operator - (const Mat2f& m) const
{
    Mat2f result;

    result[0] = row[0] - m[0];
    result[1] = row[1] - m[1];

    return(result);
}

inline Mat2f Mat2f::operator - () const
{
    Mat2f result;

    result[0] = -row[0];
    result[1] = -row[1];

    return(result);
}

inline Mat2f Mat2f::operator * (const Mat2f& m) const
{




    Mat2f result;

    result[0][0] = row[0][0] * m.row[0][0] + row[0][1] * m.row[1][0];
    result[0][1] = row[0][0] * m.row[0][1] + row[0][1] * m.row[1][1];
    result[1][0] = row[1][0] * m.row[0][0] + row[1][1] * m.row[1][0];
    result[1][1] = row[1][0] * m.row[0][1] + row[1][1] * m.row[1][1];

    return(result);




}

inline Mat2f Mat2f::operator * (float s) const
{
    Mat2f result;

    result[0] = row[0] * s;
    result[1] = row[1] * s;

    return(result);
}

inline Mat2f Mat2f::operator / (float s) const
{
    Mat2f result;

    result[0] = row[0] / s;
    result[1] = row[1] / s;

    return(result);
}

inline Mat2f operator * (float s, const Mat2f& m)
{
    return(m * s);
}

inline Vec2f operator * (const Mat2f& m, const Vec2f& v)
{
    Vec2f result;

    result[0] = m[0][0] * v[0] + m[0][1] * v[1];
    result[1] = m[1][0] * v[0] + m[1][1] * v[1];

    return(result);
}

inline Vec2f operator * (const Vec2f& v, const Mat2f& m)
{
    Vec2f result;

    result[0] = v[0] * m[0][0] + v[1] * m[1][0];
    result[1] = v[0] * m[0][1] + v[1] * m[1][1];

    return(result);
}

inline Vec2f& operator *= (Vec2f& v, const Mat2f& m)
{
    float t;

    t = v[0] * m[0][0] + v[1] * m[1][0];
    v[1] = v[0] * m[0][1] + v[1] * m[1][1];
    v[0] = t;

    return(v);
}


inline Mat2f trans(const Mat2f& m)
{
    Mat2f result;

    result[0][0] = m[0][0]; result[0][1] = m[1][0];
    result[1][0] = m[0][1]; result[1][1] = m[1][1];

    return(result);
}

inline float trace(const Mat2f& m)
{
    return(m[0][0] + m[1][1]);
}

inline Mat2f adj(const Mat2f& m)
{
    Mat2f result;

    result[0] = -cross(m[1]);
    result[1] = cross(m[0]);

    return(result);
}







inline Vec2f xform(const Mat2f& m, const Vec2f& v)
{ return(m * v); }
inline Mat2f xform(const Mat2f& m, const Mat2f& n)
{ return(m * n); }
class Vec4f;

class Mat3f
{
public:

    Mat3f();
    Mat3f(float a, float b, float c,
          float d, float e, float f,
          float g, float h, float i);
    Mat3f(const Mat3f& m);
    Mat3f(ZeroOrOne k);
    Mat3f(Block k);


    int Rows() const { return(3); };
    int Cols() const { return(3); };

    Vec3f& operator [] (int i);
    const Vec3f& operator [] (int i) const;

    float& operator () (int i, int j);
    float operator () (int i, int j) const;

    float* Ref() const;


    Mat3f& operator = (const Mat3f& m);
    Mat3f& operator = (ZeroOrOne k);
    Mat3f& operator = (Block k);
    Mat3f& operator += (const Mat3f& m);
    Mat3f& operator -= (const Mat3f& m);
    Mat3f& operator *= (const Mat3f& m);
    Mat3f& operator *= (float s);
    Mat3f& operator /= (float s);


    bool operator == (const Mat3f& m) const;
    bool operator != (const Mat3f& m) const;


    Mat3f operator + (const Mat3f& m) const;
    Mat3f operator - (const Mat3f& m) const;
    Mat3f operator - () const;
    Mat3f operator * (const Mat3f& m) const;
    Mat3f operator * (float s) const;
    Mat3f operator / (float s) const;


    void MakeZero();
    void MakeIdentity();
    void MakeDiag(float k = vl_one);
    void MakeBlock(float k = vl_one);


    Mat3f& MakeRot(const Vec3f& axis, Real theta);
    Mat3f& MakeRot(const Vec4f& q);
    Mat3f& MakeScale(const Vec3f& s);


    Mat3f& MakeHRot(Real theta);
    Mat3f& MakeHScale(const Vec2f& s);
    Mat3f& MakeHTrans(const Vec2f& t);


    Vec3f row[3];
};




Vec3f& operator *= (Vec3f& v, const Mat3f& m);
Vec3f operator * (const Mat3f& m, const Vec3f& v);
Vec3f operator * (const Vec3f& v, const Mat3f& m);
Mat3f operator * (const float s, const Mat3f& m);

Mat3f trans(const Mat3f& m);
float trace(const Mat3f& m);
Mat3f adj(const Mat3f& m);
float det(const Mat3f& m);
Mat3f inv(const Mat3f& m);
Mat3f oprod(const Vec3f& a, const Vec3f& b);



Vec3f xform(const Mat3f& m, const Vec3f& v);
Vec2f xform(const Mat3f& m, const Vec2f& v);
Mat3f xform(const Mat3f& m, const Mat3f& n);




inline Mat3f::Mat3f()
{
}

inline Vec3f& Mat3f::operator [] (int i)
{
    ;
    return(row[i]);
}

inline const Vec3f& Mat3f::operator [] (int i) const
{
    ;
    return(row[i]);
}

inline float& Mat3f::operator () (int i, int j)
{
    ;
    ;

    return(row[i][j]);
}

inline float Mat3f::operator () (int i, int j) const
{
    ;
    ;

    return(row[i][j]);
}

inline float* Mat3f::Ref() const
{
    return((float *) row);
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

    return(*this);
}

inline Mat3f& Mat3f::operator = (Block k)
{
    MakeBlock(float(k));

    return(*this);
}

inline Mat3f operator * (const float s, const Mat3f& m)
{
    return(m * s);
}

inline Vec3f operator * (const Mat3f& m, const Vec3f& v)
{
    Vec3f result;

    result[0] = v[0] * m[0][0] + v[1] * m[0][1] + v[2] * m[0][2];
    result[1] = v[0] * m[1][0] + v[1] * m[1][1] + v[2] * m[1][2];
    result[2] = v[0] * m[2][0] + v[1] * m[2][1] + v[2] * m[2][2];

    return(result);
}

inline Vec3f operator * (const Vec3f& v, const Mat3f& m)
{
    Vec3f result;

    result[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
    result[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
    result[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];

    return(result);
}

inline Vec3f& operator *= (Vec3f& v, const Mat3f& m)
{
    float t0, t1;

    t0 = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
    t1 = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
    v[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
    v[0] = t0;
    v[1] = t1;

    return(v);
}
inline Vec2f xform(const Mat3f& m, const Vec2f& v)
{ return(proj(m * Vec3f(v, 1.0))); }
inline Vec3f xform(const Mat3f& m, const Vec3f& v)
{ return(m * v); }
inline Mat3f xform(const Mat3f& m, const Mat3f& n)
{ return(m * n); }
class Vec3f;

class Mat4f
{
public:



    Mat4f();
    Mat4f(float a, float b, float c, float d,
          float e, float f, float g, float h,
          float i, float j, float k, float l,
          float m, float n, float o, float p);
    Mat4f(const Mat4f& m);
    Mat4f(ZeroOrOne k);
    Mat4f(Block k);


    int Rows() const { return(4); };
    int Cols() const { return(4); };

    Vec4f& operator [] (int i);
    const Vec4f& operator [] (int i) const;

    float& operator () (int i, int j);
    float operator () (int i, int j) const;

    float* Ref() const;


    Mat4f& operator = (const Mat4f& m);
    Mat4f& operator = (ZeroOrOne k);
    Mat4f& operator = (Block k);
    Mat4f& operator += (const Mat4f& m);
    Mat4f& operator -= (const Mat4f& m);
    Mat4f& operator *= (const Mat4f& m);
    Mat4f& operator *= (float s);
    Mat4f& operator /= (float s);


    bool operator == (const Mat4f& m) const;
    bool operator != (const Mat4f& m) const;


    Mat4f operator + (const Mat4f& m) const;
    Mat4f operator - (const Mat4f& m) const;
    Mat4f operator - () const;
    Mat4f operator * (const Mat4f& m) const;
    Mat4f operator * (float s) const;
    Mat4f operator / (float s) const;


    void MakeZero();
    void MakeIdentity();
    void MakeDiag(float k = vl_one);
    void MakeBlock(float k = vl_one);


    Mat4f& MakeHRot(const Vec3f& axis, Real theta);
    Mat4f& MakeHRot(const Vec4f& q);
    Mat4f& MakeHScale(const Vec3f& s);

    Mat4f& MakeHTrans(const Vec3f& t);

    Mat4f& Transpose();
    Mat4f& AddShift(const Vec3f& t);


    Vec4f row[4];
};




Vec4f operator * (const Mat4f& m, const Vec4f& v);
Vec4f operator * (const Vec4f& v, const Mat4f& m);
Vec4f& operator *= (Vec4f& a, const Mat4f& m);
Mat4f operator * (float s, const Mat4f& m);

Mat4f trans(const Mat4f& m);
float trace(const Mat4f& m);
Mat4f adj(const Mat4f& m);
float det(const Mat4f& m);
Mat4f inv(const Mat4f& m);
Mat4f oprod(const Vec4f& a, const Vec4f& b);



Vec4f xform(const Mat4f& m, const Vec4f& v);
Vec3f xform(const Mat4f& m, const Vec3f& v);
Mat4f xform(const Mat4f& m, const Mat4f& n);




inline Mat4f::Mat4f()
{
}

inline Vec4f& Mat4f::operator [] (int i)
{
    ;
    return(row[i]);
}

inline const Vec4f& Mat4f::operator [] (int i) const
{
    ;
    return(row[i]);
}

inline float& Mat4f::operator () (int i, int j)
{
    ;
    ;

    return(row[i][j]);
}

inline float Mat4f::operator () (int i, int j) const
{
    ;
    ;

    return(row[i][j]);
}

inline float* Mat4f::Ref() const
{
    return (float*) row[0].elt;
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

    return(*this);
}

inline Mat4f& Mat4f::operator = (Block k)
{
    MakeBlock(float(k));

    return(*this);
}

inline Mat4f operator * (float s, const Mat4f& m)
{
    return(m * s);
}
inline Vec3f xform(const Mat4f& m, const Vec3f& v)
{ return(proj(m * Vec4f(v, 1.0))); }
inline Vec4f xform(const Mat4f& m, const Vec4f& v)
{ return(m * v); }
inline Mat4f xform(const Mat4f& m, const Mat4f& n)
{ return(m * n); }
inline Mat2f Rot2f(Real theta)
                { Mat2f result; result.MakeRot(theta); return(result); }
inline Mat2f Scale2f(const Vec2f &s)
                { Mat2f result; result.MakeScale(s); return(result); }

inline Mat3f Rot3f(const Vec3f &axis, Real theta)
                { Mat3f result; result.MakeRot(axis, theta); return(result); }
inline Mat3f Rot3f(const Vec4f &q)
                { Mat3f result; result.MakeRot(q); return(result); }
inline Mat3f Scale3f(const Vec3f &s)
                { Mat3f result; result.MakeScale(s); return(result); }
inline Mat3f HRot3f(Real theta)
                { Mat3f result; result.MakeHRot(theta); return(result); }
inline Mat3f HScale3f(const Vec2f &s)
                { Mat3f result; result.MakeHScale(s); return(result); }
inline Mat3f HTrans3f(const Vec2f &t)
                { Mat3f result; result.MakeHTrans(t); return(result); }

inline Mat4f HRot4f(const Vec3f &axis, Real theta)
                { Mat4f result; result.MakeHRot(axis, theta); return(result); }
inline Mat4f HRot4f(const Vec4f &q)
                { Mat4f result; result.MakeHRot(q); return(result); }
inline Mat4f HScale4f(const Vec3f &s)
                { Mat4f result; result.MakeHScale(s); return(result); }
inline Mat4f HTrans4f(const Vec3f &t)
                { Mat4f result; result.MakeHTrans(t); return(result); }
