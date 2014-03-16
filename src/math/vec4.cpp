#include <library/math/vector.hpp>

#include <library/log.hpp>
#include <cmath>

namespace library
{
	// vec4 constructors
	vec4::vec4() : vec3(0.0)
	{
		w = 0.0;
	}
	vec4::vec4(vector_t v) : vec3(v)
	{
		w = v;
	}
	vec4::vec4(vector_t v, vector_t W) : vec3(v, v, v)
	{
		w = W;
	}
	vec4::vec4(vector_t X, vector_t Y, vector_t Z, vector_t W) : vec3(X, Y, Z)
	{
		w = W;
	}
	vec4::vec4(const vec3& v, vector_t W) : vec3(v)
	{
		w = W;
	}
	
	// 4-vector swizzle functions
	vec3 vec4::xyz() const
	{
		return vec3(x, y, z);
	}
	
	// 4-vector operators
	
	// unary - (negate)
	const vec4 vec4::operator - () const
	{
		return vec4(-x, -y, -z, w);
	}
	
	// log output functions
	
	// write vector-4 to log using format (x, y, z, w)
	Log& operator<< (Log& out, const vec4& v)
	{
		return out << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
	}
	
	// write vector-4 to cout using format (x, y, z, w)
	std::ostream& operator<< (std::ostream& out, const vec4& v)
	{
		return out << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
	}
	
}
