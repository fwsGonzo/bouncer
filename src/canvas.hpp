#ifndef CANVAS_HPP
#define CANVAS_HPP

#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <vector>

class Particle;

class Canvas : public Fl_Widget
{
public:
	Canvas(int X, int Y, int W, int H, int particleCount, int FPS);
	~Canvas();
	
	void draw();
	
	// redraw timer (based on framerate)
	static void cb_timer(void* userdata)
	{
		Canvas* canvas = (Canvas*) userdata;
		canvas->redraw();
		Fl::repeat_timeout(canvas->framerate, cb_timer, userdata);
    }
	
private:
	double framerate;
	std::vector<Particle*> particles;
};

#endif
