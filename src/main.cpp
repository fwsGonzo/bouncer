#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>
#include <time.h>

static const int FPS = 60;
static const int PARTICLES = 10;

#include "canvas.hpp"

int main()
{
	srand(time(0));
	const int WINDOW_W  = Fl::w() * 0.75;
	const int WINDOW_H  = Fl::h() * 0.75;
	const int WINDOW_X  = (Fl::w() / 2) - (WINDOW_W / 2);
	const int WINDOW_Y  = (Fl::h() / 2) - (WINDOW_H / 2);
	
	Fl_Double_Window win(WINDOW_X, WINDOW_Y, WINDOW_W, WINDOW_H);
	Canvas canvas(0, 0, WINDOW_W, WINDOW_H, PARTICLES, FPS);
	
	win.end();
	win.resizable(canvas);
	win.show();
	return Fl::run();
}
