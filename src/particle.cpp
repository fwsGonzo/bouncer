#include "particle.hpp"

#include <library/math/toolbox.hpp>

using namespace library;

static const int TAILSIZE  = 4;
static const int PARTICLE_SIZE_MIN = 16;
static const int PARTICLE_SIZE_MAX = 32;

static const float GRAVITY = 1.0;

Particle::Particle(int w, int h)
{
	// assign random position
	vec2 pos;
	pos.x = toolbox::rnd((float) w);
	pos.y = toolbox::rnd((float) h);
	
	// assign random direction
	direction.x = toolbox::rndNorm(1.0);
	direction.y = toolbox::rndNorm(1.0);
	direction.normalize();
	direction *= 0.5;
	
	// assign random size
	size = toolbox::rnd(PARTICLE_SIZE_MAX - PARTICLE_SIZE_MIN)
		+ PARTICLE_SIZE_MIN;
	
	// assign random color
	r = toolbox::rnd(256);
	g = toolbox::rnd(256);
	b = toolbox::rnd(256);
	
	// create tail
	for (int i = 0; i < TAILSIZE; i++)
		tail.emplace_back(pos);
}

void Particle::translate(std::vector<Particle*>& particles, vec2& bounds)
{
	float rad1 = getSize() / 2.0;
	
	/// determine gravity, and un-stuck ///
	for (Particle* p : particles)
	{
		// avoid applying physics from itself
		if (p == this) continue;
		
		// get distance (and direction) from particle A to B
		vec2  delta = p->getPosition(0) - getPosition(0);
		float dist  = delta.length();
		float rad2  = p->getSize() / 2.0;
		
		if (dist > rad1 + rad2)
		{
			// fast normalize
			delta /= dist;
			// apply gravity
			this->direction += delta * GRAVITY * (rad1 + rad2) / (dist * dist);
		}
		else
		{
			// minimum translation distance
			float mtd = (rad1 + rad2 - dist) * 0.5;
			
			// get new normalized direction
			delta /= dist;
			
			// simplest un-stuck
			getPosition(0) -= delta * mtd;
		}
	}
	
	vec2 pos = getPosition(0) + direction;
	
	/// collisions & mass ///
	for (Particle* p : particles)
	{
		// avoid applying physics from itself
		if (p == this) continue;
		
		// get distance (and direction) from particle A to B
		vec2  delta = p->getPosition(0) - pos;
		float dist  = delta.length();
		float rad2  = p->getSize() / 2;
		
		if (dist <= rad1 + rad2)
		{
			// get new normalized direction
			delta /= dist;
			
			float aci = direction.dot(delta);
			float bci = p->direction.dot(delta);
			
			// perpendicular components from collision
			vec2 pa = (bci - aci) * delta;
			vec2 pb = (aci - bci) * delta;
			
			// accelerate apart
			direction += pa;
			p->direction += pb;
		}
	}
	
	pos = getPosition(0) + direction;
	
	// negate directions if outside of the boundries
	if (pos.x < rad1 || pos.x + rad1 > bounds.x)
		direction.x = -direction.x;
	if (pos.y < rad1 || pos.y + rad1 > bounds.y)
		direction.y = -direction.y;
	
	tail[0] += direction;
	
	// un-stuck from boundries
	// X-axis
	vec2& rpos = getPosition(0);
	if (rpos.x < rad1)
		rpos.x = rpos.x + rad1;
	if (rpos.x + rad1 > bounds.x)
		rpos.x = bounds.x - rad1;
	// Y-axis
	if (rpos.y < rad1)
		rpos.y = rad1;
	if (pos.y + rad1 > bounds.y)
		rpos.y = bounds.y - rad1;
	
	// interpolate tail positions
	int TS = tail.size()-1;
	for (int i = 1; i <= TS; i++)
		getPosition(i) = 
		getPosition(i).mix(getPosition(0), 0.2 + 0.3 * (1.0 - (float)i / TS));
}

Fl_Color Particle::getColor(int index)
{
	float age = 1.0 - index / (float) tail.size();
	
	int r = this->r * age;
	int g = this->g * age;
	int b = this->b * age;
	return (r << 24) + (g << 16) + (b << 8);
}

