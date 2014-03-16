#include <library/math/frustum.hpp>

#include <library/math/matrix.hpp>
#include <cmath>

namespace library
{
	
	typedef enum
	{
		E_RIGHT   = 0,  // The RIGHT side of the frustum
		E_LEFT    = 1,  // The LEFT  side of the frustum
		E_BOTTOM  = 2,  // The BOTTOM side of the frustum
		E_TOP     = 3,  // The TOP side of the frustum
		E_BACK    = 4,  // The BACK side of the frustum
		E_FRONT   = 5   // The FRONT side of the frustum
		
	} E_FRUSTUM_SIDE;
	
	typedef enum
	{
		A = 0,  		// The X value of the plane's normal
		B = 1,    		// The Y value of the plane's normal
		C = 2,     		// The Z value of the plane's normal
		D = 3         	// The distance the plane is from the origin
		
	} E_PLANE_DATA;
	
	// NORMALIZE PLANE
	// This normalizes a plane (A side) from a given frustum.
	void Frustum::normalizePlane(int side)
	{
		// Here we calculate the magnitude of the normal to the plane (point A B C)
		// Remember that (A, B][C] is that same thing as the normal's (X, Y, Z).
		// To calculate magnitude you use the equation:  magnitude = sqrt( x^2 + y^2 + z^2)
		float inv_mag = 1.0 / sqrtf(frustum[side][A] * frustum[side][A] +
									frustum[side][B] * frustum[side][B] +
									frustum[side][C] * frustum[side][C]);
		
		// Then we divide the plane's values by it's magnitude.
		frustum[side][A] *= inv_mag;
		frustum[side][B] *= inv_mag;
		frustum[side][C] *= inv_mag;
		frustum[side][D] *= inv_mag;
		
	}
	
	void Frustum::calculate(const mat4& matproj, const mat4& matview)
	{
		// We have our modelview and projection matrix, if we combine these 2 matrices,
		// it will give us our clipping planes.
		calculate(matproj * matview);
	}
	
	// CALCULATE FRUSTUM 
	// This extracts our frustum from the projection and modelview matrix.
	void Frustum::calculate(const mat4& clip)
	{
		// Now we actually want to get the sides of the frustum.  To do this we take
		// the clipping planes we received above and extract the sides from them.
		
		// This will extract the RIGHT side of the frustum
		frustum[E_RIGHT][A] = clip[ 3] - clip[ 0];
		frustum[E_RIGHT][B] = clip[ 7] - clip[ 4];
		frustum[E_RIGHT][C] = clip[11] - clip[ 8];
		frustum[E_RIGHT][D] = clip[15] - clip[12];
		
		// Now that we have a normal (A,B,C) and a distance (D) to the plane,
		// we want to normalize that normal and distance.
		
		// Normalize the RIGHT side
		normalizePlane(E_RIGHT);
		
		// This will extract the LEFT side of the frustum
		frustum[E_LEFT][A] = clip[ 3] + clip[ 0];
		frustum[E_LEFT][B] = clip[ 7] + clip[ 4];
		frustum[E_LEFT][C] = clip[11] + clip[ 8];
		frustum[E_LEFT][D] = clip[15] + clip[12];
		
		// Normalize the LEFT side
		normalizePlane(E_LEFT);
		
		// This will extract the BOTTOM side of the frustum
		frustum[E_BOTTOM][A] = clip[ 3] + clip[ 1];
		frustum[E_BOTTOM][B] = clip[ 7] + clip[ 5];
		frustum[E_BOTTOM][C] = clip[11] + clip[ 9];
		frustum[E_BOTTOM][D] = clip[15] + clip[13];
		
		// Normalize the BOTTOM side
		normalizePlane(E_BOTTOM);
		
		// This will extract the TOP side of the frustum
		frustum[E_TOP][A] = clip[ 3] - clip[ 1];
		frustum[E_TOP][B] = clip[ 7] - clip[ 5];
		frustum[E_TOP][C] = clip[11] - clip[ 9];
		frustum[E_TOP][D] = clip[15] - clip[13];
		
		// Normalize the TOP side
		normalizePlane(E_TOP);
		
		// This will extract the BACK side of the frustum
		frustum[E_BACK][A] = clip[ 3] - clip[ 2];
		frustum[E_BACK][B] = clip[ 7] - clip[ 6];
		frustum[E_BACK][C] = clip[11] - clip[10];
		frustum[E_BACK][D] = clip[15] - clip[14];
		
		// Normalize the BACK side
		normalizePlane(E_BACK);
		
		// This will extract the FRONT side of the frustum
		frustum[E_FRONT][A] = clip[ 3] + clip[ 2];
		frustum[E_FRONT][B] = clip[ 7] + clip[ 6];
		frustum[E_FRONT][C] = clip[11] + clip[10];
		frustum[E_FRONT][D] = clip[15] + clip[14];
		
		// Normalize the FRONT side
		normalizePlane(E_FRONT);
	}
	
	// The code below will allow us to make checks within the frustum.  For example,
	// if we want to see if a point, a sphere, or a cube lies inside of the frustum.
	// Because all of our planes point INWARDS (The normals are all pointing inside the frustum)
	// we then can assume that if a point is in FRONT of all of the planes, it's inside.
	
	// POINT IN FRUSTUM 
	// This determines if a point is inside of the frustum
	bool Frustum::point(float x, float y, float z) const
	{
		// Go through all the sides of the frustum
		for (int i = 0; i < 6; i++)
		{
			// Calculate the plane equation and check if the point is behind a side of the frustum
			if (frustum[i][A] * x + frustum[i][B] * y + frustum[i][C] * z + frustum[i][D] <= 0)
			{
				// The point was behind a side, so it ISN'T in the frustum
				return false;
			}
		}
		// The point was inside of the frustum (In front of ALL the sides of the frustum)
		return true;
	}
	
	// SPHERE IN FRUSTUM 
	// This determines if a sphere is inside of our frustum by it's center and radius.
	bool Frustum::sphere(frustum_t x, frustum_t y, frustum_t z, frustum_t radius) const
	{
		for (int i = 0; i < 6; i++)
		{
			if (frustum[i][0] * x + frustum[i][1] * y + frustum[i][2] * z + frustum[i][3] <= radius)
				return false;
		}
		return true;
	}
	
	// CUBE IN FRUSTUM 
	//  This determines if a cube is in or around our frustum by it's center and 1/2 it's length
	bool Frustum::cube(float x, float y, float z, float size) const
	{
		float size_y = size * 0.5;
		
		for (int i = 0; i < 6; i++)
		{
			if (frustum[i][A] * (x - size) + frustum[i][B] * (y - size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (y - size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (y + size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (y + size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (y - size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (y - size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (y + size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (y + size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			// If we get here, it isn't in the frustum
			return false;
			
		} // next plane
		
		return true;
	}
	
	// Column in frustum (multiple levels of "Sectors" stacked on top)
	
	bool Frustum::column(float x, float z, int cy, int height, float size) const
	{
		float size_y = cy + height;
		float base_y = cy;
		
		for (int i = 0; i < 6; i++)
		{
			if (frustum[i][A] * (x - size) + frustum[i][B] * (base_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (base_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (size_y) + frustum[i][C] * (z - size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (base_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (base_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x - size) + frustum[i][B] * (size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			if (frustum[i][A] * (x + size) + frustum[i][B] * (size_y) + frustum[i][C] * (z + size) + frustum[i][D] > 0)
				continue;
			
			// not in frustum
			return false;
		}
		
		return true;
	}
	
	// create perspective frustum
	// field of view, aspect (width / height), znear plane, zfar plane
	void Frustum::createPerspective(float fovy, float aspect, float znear, float zfar)
	{
		float an = fovy * (3.141592653589f / 180.0f);
		float si = sin(an);
		float co = cos(an);
		
		frustum[0][0] =  0.0f;  frustum[0][1] = -co;    frustum[0][2] =  si;		frustum[0][3] =  0.0f;
		frustum[1][0] =  0.0f;  frustum[1][1] =  co;    frustum[1][2] =  si;		frustum[1][3] =  0.0f;
		frustum[2][0] =  co;    frustum[2][1] =  0.0f;  frustum[2][2] =  si*aspect; frustum[2][3] =  0.0f;
		frustum[3][0] = -co;    frustum[3][1] =  0.0f;  frustum[3][2] =  si*aspect; frustum[3][3] =  0.0f;
		frustum[4][0] =  0.0f;  frustum[4][1] =  0.0f;  frustum[4][2] =  1.0f;		frustum[4][3] =  zfar;
		frustum[5][0] =  0.0f;  frustum[5][1] =  0.0f;  frustum[5][2] = -1.0f;		frustum[5][3] = -znear;
	}
	
}
