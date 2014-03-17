#include "canvas.hpp"

#include <library/math/vector.hpp>
#include "particle.hpp"

using namespace library;

Canvas::Canvas(int X, int Y, int W, int H, int particleCount, int FPS)
	: Fl_Widget(X,Y,W,H,nullptr)
{
	framerate = 1.0 / FPS;
	
	while (particleCount--)
		particles.push_back(new Particle(W, H));
	
	Fl::add_timeout(framerate, cb_timer, (void*) this);	
}
Canvas::~Canvas()
{
	for (Particle* p : particles)
		delete p;
}

void Canvas::draw()
{
	// clear screen with black color
	fl_color(0);
	fl_rectf(0, 0, w(), h());
	
	vec2 bounds(w(), h());
	
	for (Particle* p : particles)
	{
		/// physics time-step ///
		p->translate(particles, bounds);
		
		int size = p->getSize();
		
		/// render particle ///
		for (int i = p->getTailCount()-1; i >= 0; i--)
		{
			// set color
			fl_color(p->getColor(i));
			// get position
			vec2& pos = p->getPosition(i);
			// draw part of tail
			fl_pie(pos.x - size / 2, pos.y - size / 2, size, size, 0, 360);
		}
	}
}
