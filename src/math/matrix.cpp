/**
 * 4x4 Matrix implementation
 * 
**/
#include <library/math/matrix.hpp>
#include <cmath>

namespace library
{
	mat4::mat4() {}
	
	mat4::mat4(const mat4& matrix)
	{
		for (int i = 0; i < ELEMENTS; i++)
			m[i] = matrix[i];
	}
	mat4::mat4(matrix_t matrix[ELEMENTS])
	{
		for (int i = 0; i < ELEMENTS; i++)
			m[i] = matrix[i];
		
	}
	
	// scale matrices
	mat4::mat4(matrix_t sc)
	{
		m[0] = sc;  m[4] = 0.0; m[ 8] = 0.0; m[12] = 0.0;
		m[1] = 0.0; m[5] = sc;  m[ 9] = 0.0; m[13] = 0.0;
		m[2] = 0.0; m[6] = 0.0; m[10] = sc;  m[14] = 0.0;
		m[3] = 0.0; m[7] = 0.0; m[11] = 0.0; m[15] = 1.0;
	}
	mat4::mat4(matrix_t sx, matrix_t sy, matrix_t sz)
	{
		m[0] = sx;  m[4] = 0.0; m[ 8] = 0.0; m[12] = 0.0;
		m[1] = 0.0; m[5] = sy;  m[ 9] = 0.0; m[13] = 0.0;
		m[2] = 0.0; m[6] = 0.0; m[10] = sz;  m[14] = 0.0;
		m[3] = 0.0; m[7] = 0.0; m[11] = 0.0; m[15] = 1.0;
	}
	// vector based matrix (right, up, forward)
	mat4::mat4(const vec3& xa, const vec3& ya, const vec3& za)
	{
		m[0] = xa.x; m[4] = ya.x; m[ 8] = za.x; m[12] = 0.0;
		m[1] = xa.y; m[5] = ya.y; m[ 9] = za.y; m[13] = 0.0;
		m[2] = xa.z; m[6] = ya.z; m[10] = za.z; m[14] = 0.0;
		m[3] = 0.0;  m[7] = 0.0;  m[11] = 0.0;  m[15] = 1.0;
	}
	// translation matrix constructor
	mat4::mat4(const vec3& translation)
	{
		m[0] = 1.0; m[4] = 0.0; m[ 8] = 0.0; m[12] = translation.x;
		m[1] = 0.0; m[5] = 1.0; m[ 9] = 0.0; m[13] = translation.y;
		m[2] = 0.0; m[6] = 0.0; m[10] = 1.0; m[14] = translation.z;
		m[3] = 0.0; m[7] = 0.0; m[11] = 0.0; m[15] = 1.0;
	}
	
	mat4::matrix_t* mat4::data()
	{
		return this->m;
	}
	
	mat4& mat4::identity()
	{
		m[0] = 1.0; m[4] = 0.0; m[ 8] = 0.0; m[12] = 0.0;
		m[1] = 0.0; m[5] = 1.0; m[ 9] = 0.0; m[13] = 0.0;
		m[2] = 0.0; m[6] = 0.0; m[10] = 1.0; m[14] = 0.0;
		m[3] = 0.0; m[7] = 0.0; m[11] = 0.0; m[15] = 1.0;
		return *this;
	}
	
	mat4& mat4::bias()
	{
		m[0] = 0.5; m[4] = 0.0; m[8 ] = 0.0; m[12] = 0.5;
		m[1] = 0.0; m[5] = 0.5; m[9 ] = 0.0; m[13] = 0.5;
		m[2] = 0.0; m[6] = 0.0; m[10] = 0.5; m[14] = 0.5;
		m[3] = 0.0; m[7] = 0.0; m[11] = 0.0; m[15] = 1.0;
		return *this;
	}
	
	vec4 mat4::vright() const
	{
		return vec4(m[0], m[1], m[2], m[3]);
	}
	vec4 mat4::vup() const
	{
		return vec4(m[4], m[5], m[6], m[7]);
	}
	vec4 mat4::vforward() const
	{
		return vec4(m[8], m[9], m[10], m[11]);
	}
	vec4 mat4::vtranslate() const
	{
		return vec4(m[12], m[13], m[14], m[15]);
	}
	
	// transform 3-vector (w = 1.0)
	vec3 mat4::operator* (const vec3& v) const
	{
		return vec3(
			v.x * m[0] + v.y * m[4] + v.z * m[ 8] + m[12],
			v.x * m[1] + v.y * m[5] + v.z * m[ 9] + m[13],
			v.x * m[2] + v.y * m[6] + v.z * m[10] + m[14]
		);
	}
	// transform 4-vector
	vec4 mat4::operator* (const vec4& v) const
	{
		return vec4(
			v.x * m[0] + v.y * m[4] + v.z * m[ 8] + v.w * m[12],
			v.x * m[1] + v.y * m[5] + v.z * m[ 9] + v.w * m[13],
			v.x * m[2] + v.y * m[6] + v.z * m[10] + v.w * m[14],
			v.x * m[3] + v.y * m[7] + v.z * m[11] + v.w * m[15]
		);
	}
	
	mat4 mat4::operator * (const mat4& matrix) const
	{
		return mat4(*this) *= matrix;
	}
	
	mat4& mat4::operator *= (const mat4& matrix)
	{
		mat4 temp(*this);
		
		m[0] = temp.m[0 ] * matrix.m[0 ]
			 + temp.m[4 ] * matrix.m[1 ]
			 + temp.m[8 ] * matrix.m[2 ]
			 + temp.m[12] * matrix.m[3 ];
		
		m[1] = temp.m[1 ] * matrix.m[0 ]
			 + temp.m[5 ] * matrix.m[1 ]
			 + temp.m[9 ] * matrix.m[2 ]
			 + temp.m[13] * matrix.m[3 ];
		
		m[2] = temp.m[2 ] * matrix.m[0 ]
			 + temp.m[6 ] * matrix.m[1 ]
			 + temp.m[10] * matrix.m[2 ]
			 + temp.m[14] * matrix.m[3 ];
		
		m[3] = temp.m[3 ] * matrix.m[0 ]
			 + temp.m[7 ] * matrix.m[1 ]
			 + temp.m[11] * matrix.m[2 ]
			 + temp.m[15] * matrix.m[3 ];
		
		m[4] = temp.m[0 ] * matrix.m[4 ]
			 + temp.m[4 ] * matrix.m[5 ]
			 + temp.m[8 ] * matrix.m[6 ]
			 + temp.m[12] * matrix.m[7 ];
		
		m[5] = temp.m[1 ] * matrix.m[4 ]
			 + temp.m[5 ] * matrix.m[5 ]
			 + temp.m[9 ] * matrix.m[6 ]
			 + temp.m[13] * matrix.m[7 ];
		
		m[6] = temp.m[2 ] * matrix.m[4 ]
			 + temp.m[6 ] * matrix.m[5 ]
			 + temp.m[10] * matrix.m[6 ]
			 + temp.m[14] * matrix.m[7 ];
		
		m[7] = temp.m[3 ] * matrix.m[4 ]
			 + temp.m[7 ] * matrix.m[5 ]
			 + temp.m[11] * matrix.m[6 ]
			 + temp.m[15] * matrix.m[7 ];
		
		m[8] = temp.m[0 ] * matrix.m[8 ]
			 + temp.m[4 ] * matrix.m[9 ]
			 + temp.m[8 ] * matrix.m[10]
			 + temp.m[12] * matrix.m[11];
		
		m[9] = temp.m[1 ] * matrix.m[8 ]
			 + temp.m[5 ] * matrix.m[9 ]
			 + temp.m[9 ] * matrix.m[10]
			 + temp.m[13] * matrix.m[11];
		
		m[10]= temp.m[2 ] * matrix.m[8 ]
			 + temp.m[6 ] * matrix.m[9 ]
			 + temp.m[10] * matrix.m[10]
			 + temp.m[14] * matrix.m[11];
		
		m[11]= temp.m[3 ] * matrix.m[8 ]
			 + temp.m[7 ] * matrix.m[9 ]
			 + temp.m[11] * matrix.m[10]
			 + temp.m[15] * matrix.m[11];

		m[12]= temp.m[0 ] * matrix.m[12]
			 + temp.m[4 ] * matrix.m[13]
			 + temp.m[8 ] * matrix.m[14]
			 + temp.m[12] * matrix.m[15];

		m[13]= temp.m[1 ] * matrix.m[12]
			 + temp.m[5 ] * matrix.m[13]
			 + temp.m[9 ] * matrix.m[14]
			 + temp.m[13] * matrix.m[15];

		m[14]= temp.m[2 ] * matrix.m[12]
			 + temp.m[6 ] * matrix.m[13]
			 + temp.m[10] * matrix.m[14]
			 + temp.m[14] * matrix.m[15];

		m[15]= temp.m[3 ] * matrix.m[12]
			 + temp.m[7 ] * matrix.m[13]
			 + temp.m[11] * matrix.m[14]
			 + temp.m[15] * matrix.m[15];
		
		return *this;
	}
	
	mat4& mat4::translate(matrix_t x, matrix_t y, matrix_t z)
	{
		// translate with delta +(x, y, z)
		m[12] +=  m[0 ] * x + m[4 ] * y + m[8 ] * z;
		m[13] +=  m[1 ] * x + m[5 ] * y + m[9 ] * z;
		m[14] +=  m[2 ] * x + m[6 ] * y + m[10] * z;
		//m[15] +=  m[3 ] * x + m[7 ] * y + m[11] * z;
		return *this;
	}
	
	mat4& mat4::translate_xy(matrix_t x, matrix_t y)
	{
		// translate with delta +(x, y)
		m[12] +=  m[0 ] * x + m[4] * y;
		m[13] +=  m[1 ] * x + m[5] * y;
		m[14] +=  m[2 ] * x + m[6] * y;
		//m[15] +=  m[3 ] * x + m[7] * y;
		return *this;
	}
	
	mat4& mat4::translate_xz(matrix_t x, matrix_t z)
	{
		// translate with delta +(x, y)
		m[12] +=  m[0 ] * x + m[ 8] * z;
		m[13] +=  m[1 ] * x + m[ 9] * z;
		m[14] +=  m[2 ] * x + m[10] * z;
		//m[15] +=  m[3 ] * x + m[7 ] * y;
		return *this;
	}
	
	mat4& mat4::scale(matrix_t scale)
	{
		return *this *= mat4(scale);
	}
	
	mat4& mat4::scale(matrix_t sx, matrix_t sy, matrix_t sz)
	{
		return *this *= mat4(sx, sy, sz);
	}
	
	mat4& mat4::rotateZYX(matrix_t ax, matrix_t ay, matrix_t az)
	{
		return *this *= rotationMatrix(ax, ay, az);
	}
	
	// specialized 4x4 transpose
	// returns this matrix transposed
	mat4& mat4::transpose()
	{
		matrix_t temp;
		
		// swap y, x
		temp = m[1]; m[1] = m[4]; m[4] = temp;
		temp = m[2]; m[2] = m[8]; m[8] = temp;
		temp = m[6]; m[6] = m[9]; m[9] = temp;
		
		// invert (tx, ty, tz)
	  	m[12] *= -1.0;
	  	m[13] *= -1.0;
	  	m[14] *= -1.0;
		
		return *this;
	}
	// returns new transposed matrix
	mat4 mat4::transposed() const
	{
		mat4 mat(*this);
		return mat.transpose();
	}
	
	// batch transform vertices from memory location
	void mat4::batch(void* first, int stride, int count)
	{
		for (int i = 0; i < count; i++)
		{
			float* v = (float*) ((char*)first + i * stride);
			v[0] = v[0] * m[0] + v[1] * m[4] + v[2] * m[ 8] + m[12];
			v[1] = v[0] * m[1] + v[1] * m[5] + v[2] * m[ 9] + m[13];
			v[2] = v[0] * m[2] + v[1] * m[6] + v[2] * m[10] + m[14];
		}
	}
	
	// utility functions
	
	mat4::matrix_t mat4::determinant() const
	{
		return 
		  m[ 0] * (m[5]*(m[10]*m[15]-m[14]*m[11]) + m[9]*(m[14]*m[ 7]-m[ 6]*m[15]) + m[13]*(m[ 6]*m[11]-m[10]*m[ 7]))
		+ m[ 4] * (m[1]*(m[14]*m[11]-m[10]*m[15]) + m[9]*(m[ 2]*m[15]-m[14]*m[ 3]) + m[13]*(m[10]*m[ 3]-m[ 2]*m[11]))
		+ m[ 8] * (m[1]*(m[ 6]*m[15]-m[14]*m[ 7]) + m[5]*(m[14]*m[ 3]-m[ 2]*m[15]) + m[13]*(m[ 2]*m[ 7]-m[ 6]*m[ 3]))
		+ m[12] * (m[1]*(m[10]*m[ 7]-m[ 6]*m[11]) + m[5]*(m[ 2]*m[11]-m[10]*m[ 3]) + m[ 9]*(m[ 6]*m[ 3]-m[ 2]*m[ 7]));
	}
	
	
	// returns translation (tx, ty, tz) from a view matrix (rotation + translation)
	vec3 mat4::transVector() const
	{
		return vec3(m[12], m[13], m[14]);
	}
	
	// returns rotation (rx, ry, rz) from a view matrix (rotation + translation)
	vec3 mat4::lookVector() const
	{
		return vec3(m[8], m[9], -m[10]);
	}
	
	mat4 mat4::rotation() const
	{
		mat4 rot(*this);
		rot.m[12] = rot.m[13] = rot.m[14] = 0.0;
		rot.m[15] = 1.0;
		return rot;
	}
	
	// "global" constructors
	
	// orthographic projection matrix (center top-left)
	mat4 ortho2dMatrix(mat4::matrix_t width, mat4::matrix_t height, mat4::matrix_t znear, mat4::matrix_t zfar)
	{
		mat4::matrix_t depth = zfar - znear;
		/*
			m[0] = 2/width;  m[4] = 0;         m[ 8] = 0;         m[12] = -1.0;
			m[1] = 0;        m[5] = 2/height;  m[ 9] = 0;         m[13] = -1.0;
			m[2] = 0;        m[6] = 0;         m[10] = -2/depth;  m[14] = -(zfar+znear) / depth;
			m[3] = 0;        m[7] = 0;         m[11] = 0;         m[15] =  1.0;
		*/
		mat4::matrix_t a[] =
		{	2/width, 0, 0, 0,
			0, -2/height, 0, 0,
			0, 0, -2/depth, 0,
			-1.0, 1.0, -(zfar+znear) / depth, 1.0
		};
		return mat4(a);
	}
	
	// orthographic projection matrix (center screen)
	mat4 orthoMatrix(mat4::matrix_t width, mat4::matrix_t height, mat4::matrix_t znear, mat4::matrix_t zfar)
	{
		mat4::matrix_t ndcw =  width;  // 
		mat4::matrix_t ndch = -height; // inverting Y-axis
		mat4::matrix_t depth = zfar - znear;
		/*
			m[0] = 2/ndcw;  m[4] = 0;       m[ 8] = 0;         m[12] = -1.0;
			m[1] = 0;       m[5] = 2/ndch;  m[ 9] = 0;         m[13] =  1.0;
			m[2] = 0;       m[6] = 0;       m[10] = -2/depth;  m[14] = -(zfar+znear) / depth;
			m[3] = 0;       m[7] = 0;       m[11] = 0;         m[15] =  1.0;
		*/
		mat4::matrix_t a[] =
		{
			2/ndcw, 0, 0, 0,
			0, 2/ndch, 0, 0,
			0, 0, -2/depth, 0,
			-1.0, 1.0, -(zfar+znear) / depth, 1.0
		};
		return mat4(a);
	}
	
	// FOV perspective matrix (frustum)
	mat4 perspectiveMatrix(mat4::matrix_t fov, mat4::matrix_t aspect, mat4::matrix_t znear, mat4::matrix_t zfar)
	{
		const mat4::matrix_t pio360 = 4.0 * atan(1.0) / 360.0;
		
		mat4::matrix_t h = 1.0 / tan(fov * pio360);
		mat4::matrix_t negd = znear - zfar;
		
		mat4::matrix_t m[mat4::ELEMENTS];
		m[ 0] = h / aspect; m[ 1] = 0; m[ 2] = 0;                       m[ 3] =  0;
		m[ 4] = 0;          m[ 5] = h; m[ 6] = 0;                       m[ 7] =  0;
		m[ 8] = 0;          m[ 9] = 0; m[10] = (zfar + znear) / negd;   m[11] = -1;
		m[12] = 0;          m[13] = 0; m[14] = 2 * znear * zfar / negd; m[15] =  0;
		
		return mat4(m);
	}
	
	mat4 directionMatrix(const vec3& direction, const vec3& up)
	{
		vec3 xa = cross(up, direction).normalize();
		vec3 ya = cross(direction, xa).normalize();
		
		return mat4(xa, ya, direction);
	}
	
	mat4 rotationMatrix(mat4::matrix_t ax, mat4::matrix_t ay, mat4::matrix_t az)
	{
		//////////////////////////////////////////////////////////////////////////////
		// convert Euler angles(x,y,z) to axes(left, up, forward)					//
		// Each column of the rotation matrix represents left, up and forward axis. //
		// The order of rotation is Roll->Yaw->Pitch (Rx*Ry*Rz)
		// Rx: rotation about X-axis, pitch
		// Ry: rotation about Y-axis, yaw(heading)
		// Rz: rotation about Z-axis, roll
		//    Rx           Ry          Rz
		// |1  0   0| | Cy  0 Sy| |Cz -Sz 0|   | CyCz        -CySz         Sy  |	//
		// |0 Cx -Sx|*|  0  1  0|*|Sz  Cz 0| = | SxSyCz+CxSz -SxSySz+CxCz -SxCy|	//
		// |0 Sx  Cx| |-Sy  0 Cy| | 0   0 1|   |-CxSyCz+SxSz  CxSySz+SxCz  CxCy|	//
		//////////////////////////////////////////////////////////////////////////////
		
		// rotation angle about X-axis (pitch)
		mat4::matrix_t sx = sin(ax);
		mat4::matrix_t cx = cos(ax);
		
		// rotation angle about Y-axis (yaw)
		mat4::matrix_t sy = sin(ay);
		mat4::matrix_t cy = cos(ay);
		
		mat4::matrix_t m[mat4::ELEMENTS];
		
		if (az != 0.0)
		{
			// rotation angle about Z-axis (roll)
			mat4::matrix_t sz = sin(az);
			mat4::matrix_t cz = cos(az);
			
			// left vector
			m[0] =  cy * cz;
			m[1] =  sx * sy * cz + cx * sz;
			m[2] = -cx * sy * cz + sx * sz;
			m[3] = 0.0; // w1
			
			// up vector
			m[4] = -cy * sz;
			m[5] = -sx * sy * sz + cx * cz;
			m[6] =  cx * sy * sz + sx * cz;
			m[7] = 0.0; // w2
		}
		else
		{
			// sz = 0.0 (sin(0) = 0)
			// cz = 1.0 (cos(0) = 1)
			
			// left vector
			m[0] =  cy;
			m[1] =  sx * sy;
			m[2] = -cx * sy;
			m[3] = 0.0; // w1
			
			// up vector
			m[4] = 0.0;
			m[5] = cx;
			m[6] = sx;
			m[7] = 0.0; // w2
		}
		
		// forward vector
		m[ 8] =  sy;
		m[ 9] = -sx * cy;
		m[10] =  cx * cy;
		m[11] = 0.0; // w3
		
		// translation (tx, ty, tz, w)
		m[12] = 0.0; // tx
		m[13] = 0.0; // ty
		m[14] = 0.0; // tz
		m[15] = 1.0; // w4
		
		return mat4(m);
	}
	
}
