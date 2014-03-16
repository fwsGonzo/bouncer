#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <FL/fl_draw.H>
#include <library/math/toolbox.hpp>
#include <library/math/vector.hpp>
#include <vector>
#include <iostream>

using namespace library;

static const int PARTICLES = 80;

static const double    FPS = 60.0;
static const int WINDOW_W  = 640;
static const int WINDOW_H  = 480;
static const int TAILSIZE  = 8;
static const int PARTICLE_SIZE_MIN = 10;
static const int PARTICLE_SIZE_MAX = 32;

static const float GRAVITY = 8.0;
static const float GRAVITY_CURVE = 2.0;

struct Tail
{
	library::vec2 position;
	
	Tail(vec2& pos) : position(pos) {}
};

struct Particle
{
	// first element in tail is the main ball
	std::vector<Tail> tail;
	// direction/velocity vector
	library::vec2 direction;
	int size;
	int r, g, b;
	
	Particle(int w, int h)
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
	
	inline vec2& getPosition(int t)
	{
		return tail[t].position;
	}
	inline int getSize() const
	{
		return size;
	}
	
	void translate(std::vector<Particle>& particles, vec2& bounds)
	{
		float rad1 = getSize() / 2;
		
		for (Particle& p : particles)
		{
			if (&p == this) continue;
			
			vec2  delta = p.getPosition(0) - getPosition(0);
			float dist  = delta.length();
			float rad2 = p.getSize() / 2;
			
			if (dist > rad1 + rad2)
			{
				// fast normalize
				delta /= dist;
				// apply gravity
				this->direction += delta * GRAVITY / powf(dist, GRAVITY_CURVE);
			}
			else if (dist < 0.1)
			{
				// minimum translation distance
				vec2 mtd = delta * (rad1 + rad2);
				
				// push apart
				getPosition(0) += mtd * 0.5;
				p.getPosition(0) += mtd * 0.5;
			}
			else
			{
				// minimum translation distance
				float mtd = (rad1 + rad2 - dist) / dist;
				
				// get new normalized direction
				delta /= dist;
				
				float aci = direction.dot(delta);
				float bci = p.direction.dot(delta);
				
				// perpendicular components
				vec2 pa = (bci - aci) * delta;
				vec2 pb = (aci - bci) * delta;
				
				// un-stuck from particle
				getPosition(0) += mtd * pa;
				p.getPosition(0) += mtd * pb;
				// accelerate apart
				direction += pa;
				p.direction += pb;
			}
		}
		
		vec2  pos = getPosition(0) + direction;
		
		// bounce away from boundries
		if (pos.x < rad1 || pos.x + rad1 > bounds.x)
			direction.x = -direction.x;
		if (pos.y < rad1 || pos.y + rad1 > bounds.y)
			direction.y = -direction.y;
		
		tail[0].position += direction;
		
		// un-stuck from boundries
		vec2& rpos = getPosition(0);
		if (rpos.x < rad1)
			rpos.x = rpos.x + rad1;
		if (rpos.x + rad1 > bounds.x)
			rpos.x = bounds.x - rad1;
		
		if (rpos.y < rad1)
			rpos.y = rad1;
		if (pos.y + rad1 > bounds.y)
			rpos.y = bounds.y - rad1;
		
		int TS = tail.size()-1;
		for (int i = 1; i <= TS; i++)
			getPosition(i) = 
			getPosition(i).mix(getPosition(0), 0.2 + 0.3 * (1.0 - (float)i / TS));
		
	}
	
	void renderTail(int index)
	{
		fl_color(getColor(index));
		vec2& pos  = getPosition(index);
		int   size = getSize();
		fl_pie(pos.x - size / 2, pos.y - size / 2, size, size, 0, 360);
	}
	void render()
	{
		// render from lowest to highest Z-value
		for (int i = tail.size()-1; i >= 0; i--)
		{
			renderTail(i);
		}
	}
	
	Fl_Color getColor(int tail)
	{
		float age = 1.0 - tail / (float) TAILSIZE;
		
		int r = this->r * age;
		int g = this->g * age;
		int b = this->b * age;
		return (r << 24) + (g << 16) + (b << 8);
	}
};

class Canvas : public Fl_Widget
{
public:
	Canvas(int X, int Y, int W, int H, int particleCount)
		: Fl_Widget(X,Y,W,H,nullptr)
	{
		Fl::add_timeout(1.0 / FPS, cb_timer, (void*) this);	
		frame = 0;
		
		while (particleCount--)
			particles.emplace_back(W, H);
	}
	
	void draw()
	{
		frame ++;
		//std::cout << "Frame: " << frame << std::endl;
		fl_color(0);
		fl_rectf(0, 0, w(), h());
		
		vec2 bounds(w(), h());
		
		for (Particle& p : particles)
		{
			p.translate(particles, bounds);
			p.render();
		}
	}
	
	// FPS TIMER CALLBACK
	//
	static void cb_timer(void* userdata)
	{
		Canvas* canvas = (Canvas*) userdata;
		canvas->redraw();
		Fl::repeat_timeout(1.0 / FPS, cb_timer, userdata);
    }
	
private:
	int frame;
	std::vector<Particle> particles;
};

int main()
{
	srand(time(0));
	
	Fl_Double_Window win(WINDOW_W, WINDOW_H);
	Canvas canvas(0, 0, WINDOW_W, WINDOW_H, PARTICLES);
	
	win.end();
	win.resizable(canvas);
	win.show();
	return Fl::run();
}
