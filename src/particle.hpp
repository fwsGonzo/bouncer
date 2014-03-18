#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <library/math/vector.hpp>
#include <vector>

typedef unsigned int  Fl_Color;
typedef library::vec2 Tail;

struct Particle
{
	// create a particle inside (0, 0, w, h)
	Particle(int w, int h);
	
	// perform magic/physics
	void translate(std::vector<Particle*>& particles, library::vec2& bounds);
	
	inline int getTailCount() const
	{
		return tail.size();
	}
	
	inline int getSize() const
	{
		return this->size;
	}
	inline library::vec2& getPosition(int t)
	{
		return tail[t];
	}
	Fl_Color getColor(int tailIndex);
	
private:
	// first element in tail is the main ball
	std::vector<Tail> tail;
	// direction/velocity vector
	library::vec2 direction;
	int size;
	int r, g, b;
};

#endif
