#ifndef _graphics_h
#define _graphics_h

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>

/* Graphics-related variables */
extern Display *display;
extern Window win;
extern Pixmap pixmap;
extern XEvent report;
extern GC gc;
extern unsigned width, height;
extern int quit_flag;

/* Corrected function declarations */
void InitializeGraphics(char *windowName);
void CloseGraphics();
void ClearScreen();
void DrawCircle(float x, float y);
void RefreshDisplay();
void CheckForQuit();

#endif