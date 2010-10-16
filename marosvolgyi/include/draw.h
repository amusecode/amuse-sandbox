void DrawPixel(SDL_Surface *screen, int x, int y,Uint8 R, Uint8 G,Uint8 B);
void swap(int *x, int *y);
void DrawLine(SDL_Surface *screen, 
	      int x0, int y0, 
	      int x1, int y1, Uint8 R, 
	      Uint8 G,Uint8 B);
void DrawCircle(SDL_Surface *screen, 
		int x0, int y0, int radius, int n, 
		int R, int G, int B);
void DrawEllipse(SDL_Surface *screen,
		 int Xc, int Yx, 
		 int a, int b,
		 double phi,
		 int n,
		 Uint8 R, Uint8 G, Uint8 B);
void DrawText(SDL_Surface *screen, 
	      char character[],
	      int row, int col, 
	      char size,
	      Uint8 R, Uint8 G, Uint8 B 
	      );
