#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <time.h>

static const int FPS = 60;
static const int PARTICLES = 150;

static const int WINDOW_W  = 800;
static const int WINDOW_H  = 600;

#include "canvas.hpp"

int main()
{
	srand(time(0));
	
	Fl_Double_Window win(WINDOW_W, WINDOW_H);
	Canvas canvas(0, 0, WINDOW_W, WINDOW_H, PARTICLES, FPS);
	
	win.end();
	win.resizable(canvas);
	win.show();
	return Fl::run();
}
