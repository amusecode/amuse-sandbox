#include <stdio.h>
#include <SDL/SDL.h>
#include "../include/draw.h"
#include "../include/kepler.h"

#define Xres 640
#define Yres 480

void plot(SDL_Surface *screen, double r[3], short R, short G, short B) {
  DrawPixel(screen, 200.0*r[0]+320+r[2], 200.0*r[1]+240+r[2], R, G, B);
}

int main (int argc, char *argv[]) {

  SDL_Surface *screen = NULL;
  double t;
  double r[3];
  double v[3];
  char X[600];
  char Y[600];
  char Z[600];
  char m[300];

  SDL_Init( SDL_INIT_EVERYTHING );
  screen = SDL_SetVideoMode( Xres, Yres, 32, SDL_SWSURFACE );
  
  initialize(1024);
  r[0] = 1.0; r[1] = 0.1; r[2] = -0.1;
  v[0] = -0.1; v[1] = 0.1; v[2] = -0.2;
  //r[0] = 1.0; r[1] = 0.0; r[2] = 0.0;
  //v[0] = 0.0; v[1] = 1; v[2] = 0.0;
  set_position(r);
  set_velocity(v);
  set_mu(1.0);
  for (t = 0.01; t < 20; t+=0.005) {
    if (evolve_d(t)==-1) {
      printf("WARNING: Newton root finding failed @t=%2.3e\n", t);
    }
    else {
      plot(screen, r, 50,50,50);
      
      get_position_s(X, Y, Z, 60);
      sprintf(m, "X = %s", X);
      DrawText(screen, m, 1, 1, 's', 0,0,255);
      sprintf(m, "Y = %s", Y);
      DrawText(screen, m, 3, 1, 's', 0,0,255);
      sprintf(m, "Z = %s", Z);
      DrawText(screen, m, 5, 1, 's', 0,0,255);
      get_position(r);
      //get_position_s(X, Y, Z, 256);
      plot(screen, r, 255, 255, 255);

      SDL_Flip(screen);
    }
  }
 
  SDL_Quit();
  return 0;
}
