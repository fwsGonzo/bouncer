#include <library/math/vector.hpp>

#include <library/log.hpp>
#include <cmath>

namespace library
{
	static const double MIN_V2 = 1e-7;
	
	// vec2 utility
	
	vec2& vec2::rotate(vector_t angle)
	{
		vector_t x = this->x;
		this->x = cos(angle) * x - sin(angle) * this->y;
		this->y = sin(angle) * x + cos(angle) * this->y;
		return *this;
	}
	vec2 vec2::rotated(vector_t angle) const
	{
		return vec2(
			cos(angle) * this->x - sin(angle) * this->y,
			sin(angle) * this->x + cos(angle) * this->y);
	}
	
	vec2::vector_t vec2::length() const
	{
		return sqrtf(length_squared());
	}
	
	vec2 vec2::normalized() const
	{
		return vec2(*this).normalize();
	}
	
	vec2& vec2::normalize()
	{
		vector_t L = length();
		if (L == 0) return *this;
		// normalize to unit length vector
		L = 1.0 / L;
		x *= L; y *= L;
		return *this;
	}
	
	vec2::vector_t vec2::dot(const vec2& v) const
	{
		return x * v.x + y * v.y;
	}
	vec2::vector_t vec2::determinant(const vec2& v) const
	{
		return x * v.y - y * v.x;
	}
	vec2 vec2::frac() const
	{
		// FIXME: does not work for negative numbers
		return vec2(x - (int)x, y - (int)y);
	}
	vec2 vec2::reflect(const vec2& normal) const
	{
		// http://stackoverflow.com/questions/5454661/reflection-how-do-i-do-it
		return *this - 2 * dot(normal) * normal;
	}
	
	// exponentiation
	vec2& vec2::pow(double e)
	{
		this->x = powf(this->x, e);
		this->y = powf(this->y, e);
		return *this;
	}
	vec2& vec2::pow(const vec2& v)
	{
		this->x = powf(this->x, v.x);
		this->y = powf(this->y, v.y);
		return *this;
	}
	
	// linear interpolation
	vec2 vec2::mix(const vec2& v, float mixlevel) const
	{
		return vec2(
			this->x * (1.0 - mixlevel) + v.x * mixlevel,
			this->y * (1.0 - mixlevel) + v.y * mixlevel
		);
	}
	
	// unary - (negate)
	vec2 vec2::operator - () const
	{
		return vec2(-this->x, -this->y);
	}
	
	// arithmetic operators
	vec2& vec2::operator += (const vec2& v)
	{
		this->x += v.x; this->y += v.y;
		return *this;
	}
	vec2& vec2::operator -= (const vec2& v)
	{
		this->x -= v.x; this->y -= v.y;
		return *this;
	}
	vec2& vec2::operator *= (const vec2& v)
	{
		this->x *= v.x; this->y *= v.y;
		return *this;
	}
	vec2& vec2::operator /= (const vec2& v)
	{
		this->x /= v.x; this->y /= v.y;
		return *this;
	}
	
	vec2 vec2::operator + (const vec2& v) const
	{
		return vec2(this->x + v.x, this->y + v.y);
	}
	vec2 vec2::operator - (const vec2& v) const
	{
		return vec2(this->x - v.x, this->y - v.y);
	}
	vec2 vec2::operator * (const vec2& v) const
	{
		return vec2(this->x * v.x, this->y * v.y);
	}
	vec2 vec2::operator / (const vec2& v) const
	{
		return vec2(this->x / v.x, this->y / v.y);
	}
	
	vec2& vec2::operator +=(const vector_t f)
	{
		this->x += f; this->y += f;
		return *this;
	}
	vec2& vec2::operator -=(const vector_t f)
	{
		this->x -= f; this->y -= f;
		return *this;
	}
	vec2& vec2::operator *=(const vector_t f)
	{
		this->x *= f; this->y *= f;
		return *this;
	}
	vec2& vec2::operator /=(const vector_t f)
	{
		this->x /= f; this->y /= f;
		return *this;
	}
	
	vec2  vec2::operator + (const vector_t f) const
	{
		return vec2(this->x + f, this->y + f);
	}
	vec2  vec2::operator - (const vector_t f) const
	{
		return vec2(this->x - f, this->y - f);
	}
	vec2  vec2::operator * (const vector_t f) const
	{
		return vec2(this->x * f, this->y * f);
	}
	vec2  vec2::operator / (const vector_t f) const
	{
		return vec2(this->x / f, this->y / f);
	}
	
	// boolean equality operators
	bool vec2::operator == (const vec2& v) const
	{
		return fabs(x - v.x) < MIN_V2 && fabs(y - v.y) < MIN_V2;
	}
	bool vec2::operator != (const vec2& v) const
	{
		return !(*this == v);
	}
	
	// a and b should be unit vectors
	// t is a value in the range [0, 1]
	vec2 slerp(vec2 va, vec2 vb, float t)
	{
		float theta = acos(va.dot(vb));
		return (va * sin((1 - t) * theta) + vb * sin(t * theta)) / sin(theta);
	}
	
	// rotates <from> towards <to> by dTheta radians
	vec2 rotateTowards(vec2 from, vec2 to, vec2::vector_t dTheta)
	{
		float theta = acos(from.dot(to));
		
		if (theta <= dTheta) return to;
		return slerp(from, to, theta / dTheta);
	}
	
	// log output functions
	
	// write vector-2 to log using format (x, y)
	Log& operator<< (Log& out, const vec2& v)
	{
		return out << "(" << v.x << ", " << v.y << ")";
	}
	
	// write vector-2 to cout using format (x, y)
	std::ostream& operator<< (std::ostream& out, const vec2& v)
	{
		return out << "(" << v.x << ", " << v.y << ")";
	}
	
}
