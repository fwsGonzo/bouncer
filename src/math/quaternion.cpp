#include <library/math/quaternion.hpp>

#include <library/math/matrix.hpp>
#include <library/math/vector.hpp>
#include <library/log.hpp>
#include <cmath>
#include <cassert>

namespace library
{
	static const float TOLERANCE = 0.0001;
	
	////////////////////////////////
	////      CONSTRUCTORS      ////
	////////////////////////////////
	
	Quaternion::Quaternion()
	{
		this->real = 1.0f;
		this->imag = vec3(0.0f);
	}
	Quaternion::Quaternion(const Quaternion& q)
	{
		this->real = q.real;
		this->imag = q.imag;
	}
	Quaternion::Quaternion(const float real, const vec3& imaginary)
	{
		this->real = real;
		this->imag = imaginary;
	}
	
	// quaternion from axis-angle
	Quaternion::Quaternion(const vec3& angles)
	{
		float cos_z_2 = cosf(0.5 * angles.z);
		float cos_y_2 = cosf(0.5 * angles.y);
		float cos_x_2 = cosf(0.5 * angles.x);
		
		float sin_z_2 = sinf(0.5 * angles.z);
		float sin_y_2 = sinf(0.5 * angles.y);
		float sin_x_2 = sinf(0.5 * angles.x);
		
		real   = cos_z_2 * cos_y_2 * cos_x_2 + sin_z_2 * sin_y_2 * sin_x_2;
		imag.x = cos_z_2 * cos_y_2 * sin_x_2 - sin_z_2 * sin_y_2 * cos_x_2;
		imag.y = cos_z_2 * sin_y_2 * cos_x_2 + sin_z_2 * cos_y_2 * sin_x_2;
		imag.z = sin_z_2 * cos_y_2 * cos_x_2 - cos_z_2 * sin_y_2 * sin_x_2;		
	}
	
	// quaternion from axis-angle
	Quaternion::Quaternion(const vec3& norm, const float theta)
	{
		// start by setting w = cos(theta / 2)
		real = cosf(theta / 2);
		// then, .xyz part = (x, y, z) * sin(theta / 2)
		imag = norm * sinf(theta / 2);
	}
	
	////////////////////////////////
	////    HELPER FUNCTIONS    ////
	////////////////////////////////
	
	float Quaternion::length_squared() const
	{
		return real*real + imag.x*imag.x + imag.y*imag.y + imag.z*imag.z;
	}
	float Quaternion::length() const
	{
		return sqrtf(this->length_squared());
	}
	
	bool Quaternion::isNormalized() const
	{
		return (fabsf(length_squared() - 1.0f) < TOLERANCE);
	}
	
	Quaternion& Quaternion::normalize()
	{
		return *this /= this->length();
	}
	Quaternion Quaternion::normalized() const
	{
		return *this / this->length();
	}
	
	// invert this quaternion
	Quaternion& Quaternion::invert()
	{
		return this->conjugate() /= length_squared();
	}
	// invert a quaternion
	Quaternion Quaternion::invert() const
	{
		return Quaternion(*this).invert();
	}
	
	// computes conjugate of this quaternion
	Quaternion& Quaternion::conjugate()
	{
		imag = -imag; return *this;
	}
	Quaternion Quaternion::conjugate() const
	{
		return Quaternion(*this).conjugate();
	}
	
	// logarithm of quaternion (q = [cos a, axis sin a])
	Quaternion Quaternion::log() const
	{
		float a    = acosf(real);
		float sina = sinf(a);
		
		vec3 retv(0.0);
		if (sina > 0) retv = imag / sina * a;
		
		return Quaternion(0.0, retv);
	}
	
	// exponential function of e^quaternion = exp(axis * a)
	Quaternion Quaternion::exp() const
	{
		float a    = imag.length();
		float sina = sinf(a);
		float cosa = cosf(a);
		
		vec3 retv(0.0);
		if (a > 0) retv = imag / a * sina;
		
		return Quaternion(cosa, retv);
	}
	
	// dot-product: q1 x q2
	float Quaternion::dot(const Quaternion& q) const
	{
		return real * q.real + imag.x * q.imag.x + imag.y * q.imag.y + imag.z * q.imag.z;
	}
	float Quaternion::dot(const Quaternion& q1, const Quaternion& q2)
	{
		return q1.dot(q2);
	}
	
	///////////////////////////////
	////     INTERPOLATION     ////
	///////////////////////////////
	
	// linear interpolation
	Quaternion Quaternion::lerp(const Quaternion& q1, const Quaternion& q2, const float t)
	{
		return (q1 * (1.0 - t) + q2 * t).normalized();
	}
	
	// spherical linear interpolation
	Quaternion Quaternion::slerp(const Quaternion& q1, const Quaternion& q2, const float t)
	{
		Quaternion q3;
		float dotp = Quaternion::dot(q1, q2);
		
		// dot = cos(theta)
		// if (dot < 0), q1 and q2 are more than 90 degrees apart,
		// so we can invert one to reduce spinning
		
		if (dotp < 0)
		{
			dotp = -dotp;
			q3   = -q2;
		}
		else q3 = q2;
		
		if (dotp < 0.95f)
		{
			float angle = acosf(dotp);
			return (q1 * sinf(angle * (1 - t)) + q3 * sinf(angle * t)) / sinf(angle);
		}
		// if the angle is small, use linear interpolation
		return lerp(q1, q3, t);
	}
	
	// helper for squad()
	static Quaternion slerpNoInvert(const Quaternion &q1, const Quaternion &q2, const float t)
	{
		float dotp = q1.dot(q2);
		
		if (dotp > -0.95f && dotp < 0.95f)
		{
			float angle = acosf(dotp);			
			return (q1 * sinf(angle * (1 - t)) + q2 * sinf(angle * t)) / sinf(angle);
		}
		// if the angle is small, use linear interpolation					
		return Quaternion::lerp(q1, q2, t);
	}
	
	// spherical quadrangle interpolation
	Quaternion Quaternion::squad(const Quaternion& q1, const Quaternion& q2, const Quaternion& a, const Quaternion& b, const float t)
	{
		Quaternion c = slerpNoInvert(q1, q2, t),
			       d = slerpNoInvert(a, b, t);
		
		return slerpNoInvert(c, d, 2 * t * (1 - t));
	}
	
	// Shoemake-Bezier interpolation (De Castlejau algorithm)
	Quaternion Quaternion::bezier(const Quaternion& q1, const Quaternion& q2, const Quaternion& a, const Quaternion& b, const float t)
	{
		// level 1
		Quaternion 	q11 = slerpNoInvert(q1, a,  t),
					q12 = slerpNoInvert(a,  b,  t),
					q13 = slerpNoInvert(b,  q2, t);
		// level 2 and 3
		return slerpNoInvert(slerpNoInvert(q11, q12, t), slerpNoInvert(q12, q13, t), t);
	}
	
	// calculate control point for spline interpolation
	Quaternion Quaternion::spline(const Quaternion& qnm1, const Quaternion& qn, const Quaternion& qnp1)
	{
		Quaternion qni(qn.real, -qn.imag);
		return qn * (( (qni * qnm1).log() + (qni * qnp1).log() ) / -4).exp();
	}
	
	////////////////////////////////
	////    OUTPUT FUNCTIONS    ////
	////////////////////////////////
	
	// rotate vector
	vec3 Quaternion::rotate(const vec3& vector) const
	{
		// Q = this * {0, vector} * this.conjugate
		// return Q.imaginary
		return (*this * Quaternion(0, vector) * this->conjugate()).imag;
	}
	
	// to axis-angle
	void Quaternion::toAxisAngle(vec3& axis, float& angle) const
	{
		// precompute stuff
		angle = acosf(real);
		float sinf_theta_inv = 1.0 / sinf(angle);
		// vector part
		axis = imag * sinf_theta_inv;
		// angle part
		angle *= 2.0;
	}
	
	// to 4x4 matrix
	mat4 Quaternion::toMatrix() const
	{
		assert(length() > 0.9999 && length() < 1.0001);
		
		mat4::matrix_t m[] =
		{
			// 1 - 2(y^2 + z^2), 2(xy - wz), 2(xz + wy), 0
			1 - 2 * imag.y * imag.y - 2 * imag.z * imag.z, 2 * imag.x * imag.y - 2 * real * imag.z, 2 * imag.x * imag.z + 2 * real * imag.y, 0,
			// 2xy + 2wz, 1 - 2x^2 - 2z^2, 2yz + 2wx, 0
			2 * imag.x * imag.y + 2 * real * imag.z, 1 - 2 * imag.x * imag.x - 2 * imag.z * imag.z, 2 * imag.y * imag.z + 2 * real * imag.x, 0,
			// 2xz - 2wy, 2yz - 2wx, 1 - 2x^2 - 2y^2, 0
			2 * imag.x * imag.z - 2 * real * imag.y, 2 * imag.y * imag.z - 2 * real * imag.x, 1 - 2 * imag.x * imag.x - 2 * imag.y * imag.y, 0,
			0.0, 0.0, 0.0, 1.0
		};		
		return mat4(m);
	}
	
	///////////////////////////////
	////       OPERATORS       ////
	///////////////////////////////
	
	// assignment operator
	Quaternion& Quaternion::operator = (const Quaternion& q)
	{
		this->real = q.real;
		this->imag = q.imag;
		return *this;
	}
	
	// index operator (returning quaternion elements)
	float Quaternion::operator [] (int index) const
	{
		switch (index)
		{
			case 0:
				return real;
			case 1:
				return imag.x;
			case 2:
				return imag.y;
			case 3:
				return imag.z;
		}
		return 0.0;
	}
	
	// arithmetic operators
	const Quaternion Quaternion::operator + (const Quaternion& q) const
	{
		return Quaternion(real + q.real, imag + q.imag);
	}
	const Quaternion Quaternion::operator - (const Quaternion& q) const
	{
		return Quaternion(real - q.real, imag - q.imag);
	}
	
	// multiplication q1 * q2
	Quaternion& Quaternion::operator *= (const Quaternion& q)
	{
		// [ w = realA * realB - dot(imagA, imagB), v = realB * imagA + realA * imagB + cross(imagA, imagB) ]
		this->real   = real * q.real - imag.dot(q.imag);
		//this->imag = real * q.imag + q.real * imag + cross(imag, q.imag);
		this->imag.x = real * q.imag.x + imag.x * q.real   + imag.y * q.imag.z - imag.z * q.imag.y;
		this->imag.y = real * q.imag.y - imag.x * q.imag.z + imag.y * q.real   + imag.z * q.imag.x;
		this->imag.z = real * q.imag.z + imag.x * q.imag.y - imag.y * q.imag.x + imag.z * q.real;
		
		return *this;
	}
	
	Quaternion Quaternion::operator * (const Quaternion& q) const
	{
		// return copy * q
		return Quaternion(*this) *= q;
	}
	
	// division
	Quaternion& Quaternion::operator /= (const Quaternion& q)
	{
		return *this *= q.invert();
	}
	Quaternion  Quaternion::operator /  (const Quaternion& q) const
	{
		return *this * q.invert();
	}
	
	// scaling
	const Quaternion  Quaternion::operator *  (const float scale) const
	{
		return Quaternion(real * scale, imag * scale);
	}
	const Quaternion  Quaternion::operator /  (const float scale) const
	{
		return Quaternion(real / scale, imag / scale);
	}
	// compound scaling
	Quaternion& Quaternion::operator *= (const float scale)
	{
		real *= scale;  imag *= scale;
		return *this;
	}
	Quaternion& Quaternion::operator /= (const float scale)
	{
		real /= scale;  imag /= scale;
		return *this;
	}
	
	// negation
	const Quaternion Quaternion::operator - () const
	{
		return Quaternion(-real, -imag);
	}
	
	
	// boolean equality operator
	bool Quaternion::operator == (const Quaternion& q) const
	{
		// either same real and imaginary,
		if (this->real == q.real && this->imag == q.imag) return true;
		// or same -real and -imaginary
		if (this->real == -q.real && this->imag == -q.imag) return true;
		
		return false;
	}
	
	// boolean inequality operator
	bool Quaternion::operator != (const Quaternion& q) const
	{
		return !(*this == q);
	}
	
	
	// log output functions
	
	// write Quaternion to log using format [w, (x, y, z)]
	Log& operator<< (Log& out, const Quaternion& q)
	{
		return out << "[" << q.real << " " << q.imag << "]";
	}
	
}
