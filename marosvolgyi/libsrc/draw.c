/*

libdraw, draw library for Simple Direct Media layer.

    Copyright (C) 2010 Marcell Marosvolgyi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#include "SDL/SDL.h"
#include <math.h>
#include <stdio.h>
#include "./charset.h"

void DrawPixel(SDL_Surface *screen, int x, int y, Uint8 R, Uint8 G, Uint8 B)
{
  //source of this procedure: http://www.libsdl.org/intro.de/usingvideo.html
  Uint32 color = SDL_MapRGB(screen->format, R, G, B);
  
  if ( SDL_MUSTLOCK(screen) ) {
    if ( SDL_LockSurface(screen) < 0 ) {
      return;
    }
  }
  switch (screen->format->BytesPerPixel) {
  case 1: { /* Assuming 8-bpp */
    Uint8 *bufp;
    
    bufp = (Uint8 *)screen->pixels + y*screen->pitch + x;
    *bufp = color;
  }
    break;
    
  case 2: { /* Probably 15-bpp or 16-bpp */
    Uint16 *bufp;
    
    bufp = (Uint16 *)screen->pixels + y*screen->pitch/2 + x;
    *bufp = color;
  }
    break;
    
  case 3: { /* Slow 24-bpp mode, usually not used */
    Uint8 *bufp;
    
    bufp = (Uint8 *)screen->pixels + y*screen->pitch + x;
    *(bufp+screen->format->Rshift/8) = R;
    *(bufp+screen->format->Gshift/8) = G;
    *(bufp+screen->format->Bshift/8) = B;
  }
    break;
    
  case 4: { /* Probably 32-bpp */
    Uint32 *bufp;
    
    bufp = (Uint32 *)screen->pixels + y*screen->pitch/4 + x;
    *bufp = color;
  }
    break;
  }
  if ( SDL_MUSTLOCK(screen) ) {
    SDL_UnlockSurface(screen);
  }
  //SDL_UpdateRect(screen, x, y, 1, 1);
}

void swap(int *x, int *y) {
  int foo;
  foo = *y;
  *y = *x;
  *x = foo;
}

void DrawLine(SDL_Surface *screen,
	      int x0, int y0, int x1, int y1,
	      int R, int G, int B) {
  //Bresenham algorithm
  int steep;
  int deltax, deltay;
  int error, ystep, y, x;

  steep = abs(y1 - y0) > abs(x1 - x0);
  if (steep) {
    swap(&x0, &y0);
    swap(&x1, &y1);
  }
  if (x0 > x1) {
    swap(&x0, &x1);
    swap(&y0, &y1);
  }
  deltax = x1 - x0;
  deltay = abs(y1 -y0);
  error = (int)1.0*deltax/2;
  y = y0;
  if (y0 < y1) {
    ystep = 1;
  }
  else {
    ystep = -1;
  }
  x=x0;
  for (; x<x1; x++) {
    if (steep) {
      DrawPixel(screen, y, x, R, G, B);
    }
    else {
      DrawPixel(screen, x, y, R, G, B);
    }
    error -= deltay;
    if (error < 0) {
      y += ystep;
      error += deltax;
    }
  }
}

void DrawBox(SDL_Surface *screen,
		 int x, int y,
		 int w, int h,
		 Uint8 R, Uint8 G, Uint8 B) {

  int x0, y0, x1, y1;
  int x2, y2, x3, y3;

  x0 = x;
  y0 = y;

  x1 = x+w;
  y1 = y;

  x3 = x;
  y3 = y+h;

  x2 = x+w;
  y2 = y+h;

  DrawLine(screen, x0, y0, x1, y1, R, G, B);
  DrawLine(screen, x1, y1, x2, y2, R, G, B);
  DrawLine(screen, x2, y2, x3, y3, R, G, B);
  DrawLine(screen, x3, y3, x0, y0, R, G, B);

}
void DrawEllipse(SDL_Surface *screen,
		 int Xc, int Yc,
		 int a, int b,
		 double p,
		 int n,
		 Uint8 R, Uint8 G, Uint8 B) {
  int i;
  double x_old, y_old, x, y;
  double cosp, sinp;
  double t;

  cosp = cos(p);
  sinp = sin(p);

  x_old = a * cosp + Xc;
  y_old = a * sinp + Yc;

  for (i = 0; i <= n; i++) {
    t = 1.0 * i/n * 2 * 3.1415926535;
    x =  a * cos(t)*cosp - b * sin(t)*sinp + Xc;
    y =  a * cos(t)*sinp + b * sin(t)*cosp + Yc;
    DrawLine (screen,
	      (int)x_old, (int)y_old,
	      (int)x, (int)y,
	      R, G, B);
    x_old = x;
    y_old = y;
  }
}

void DrawCircle(SDL_Surface *screen,
		int x0, int y0, int radius, int n,
		int R, int G, int B) {
  DrawEllipse(screen, x0, y0, radius, radius, 0.0, n, R, G, B);
}

void Bezier(SDL_Surface *screen,
	    int x0, int y0,
	    int x1, int y1,
	    int x2, int y2) {

  double t;
  double x, y;
  double x_old, y_old;
  double dx, dy, l, inc;
  double f1, f2, f3;

  //B(t) = (1-t^2)P0 _ 2(1-t)tP1+t^2P2 t in [0,1]
  x_old = (int)x0;
  y_old = (int)y0;

  dx = abs(x0-x2);
  dy = abs(y0-y2);

  l = pow(dx*dx+dy*dy, 0.5);
  inc = 10.0/(l==0?1:l);

  for (t = 0; t<1; t+=inc) {
    f1 = (1-t*t);
    f2 = 2*(1-t)*t;
    f3 = t*t;
    x = f1 * x0 + f2 * x1 + f3 * x2;
    y = f1 * y0 + f2 * y1 + f3 * y2;
    DrawLine(screen,
	     (int)x_old, (int)y_old,
	     (int)x, (int)y,
	     255, 255, 255);
    x_old = x;
    y_old = y;
  }
}

void DrawChar(SDL_Surface *screen, int charnr, int row, int col,
	       int R, int G, int B) {
  int i, offset;
  int x, y;
  offset = charnr*64;

  for (i=0; i<64; i++) {
    x = i%8 + col;
    y = (int)i/8 + row;
    if (A[i+offset]) {
      DrawPixel(screen, x, y, R, G, B);
    }
    else {
      DrawPixel(screen, x, y, 0, 0, 0);
    }
  }
}

void DrawChar2 (SDL_Surface *screen, int charnr, int row, int col,
	       int R, int G, int B) {
  int i, offset;
  int x, y;
  offset = charnr*64;

  for (i=0; i<64; i++) {
    x = (int)(2 * (i%8) + col);
    y = (int)(2 * (i/8) + row);
    if (A[i+offset]) {
      DrawPixel(screen, x,   y,   R, G, B);
      DrawPixel(screen, x+1, y,   R, G, B);
      DrawPixel(screen, x+1, y+1, R, G, B);
      DrawPixel(screen, x,   y+1, R, G, B);
    }
    else {
      DrawPixel(screen, x, y, 0, 0, 0);
    }
  }
}
/*
printing char 33!10
printing char 64@11
printing char 35#12
printing char 36$13
printing char 37%14
printing char 94^15
printing char 38&16
printing char 42*17
printing char 40(18
printing char 41)19
printing char 43+20
printing char 45-21
printing char 61=
printing char 63?27
*/
int mapascii(char i) {
  if ((i>64)&(i<65+26)) return i-24;
  if ((i>48)&(i<59)) return i-48;
  if (i==32) return 32;//space
  if (i==33) return 10;//!
  if (i==64) return 11;//@
  if (i==35) return 12;//#
  if (i==36) return 13;//$
  if (i==37) return 14;//%
  if (i==94) return 15;//^
  if (i==38) return 16;//&
  if (i==42) return 17;//*
  if (i==40) return 18;//(
  if (i==41) return 19;//)
  if (i==43) return 20;//+
  if (i==45) return 21;//-
  if (i==61) return 10;//
  if (i==63) return 27;//?
}

void DrawText (SDL_Surface *screen, char charnr[], int row, int col,
	       char size,
	       int R, int G, int B) {
  int go = 1;
  int i = 0;

  while (charnr[i]) {
    //printf("printing char %d\n", charnr[i]);
    if (size == 's') {
      DrawChar(screen,
	       mapascii(charnr[i]),
	       row*8, (col+i)*8,
	       R, G, B);
    }
    else {
      DrawChar2(screen,
		mapascii(charnr[i]),
		row*16, (col+i)*16,
		R, G, B);
    }

    if (++i>60) break;
  }
}

void ScrollText(SDL_Surface *screen, char charnr[], int row, int col,
		int R, int G, int B) {

  int i =0;
  while (charnr[i]) {
    DrawChar(screen,
	     mapascii(charnr[i]),
	     row,
	     col+8*i,
	     R, G, B);
    i++;
  }
}
