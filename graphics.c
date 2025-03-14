#include "graphics.h"

Display *display;
Window win;
Pixmap pixmap;
XEvent report;
GC gc;
unsigned width = 800, height = 800;
int quit_flag = 0;

void InitializeGraphics(char *command) {
    char *display_name = getenv("DISPLAY");
    display = XOpenDisplay(display_name);
    if (!display) {
        fprintf(stderr, "%s: cannot connect to X server\n", command);
        exit(1);
    }

    int screen_num = DefaultScreen(display);
    win = XCreateSimpleWindow(display, RootWindow(display, screen_num),
                              0, 0, width, height, 2,
                              BlackPixel(display, screen_num),
                              WhitePixel(display, screen_num));
    pixmap = XCreatePixmap(display, win, width, height, DefaultDepth(display, screen_num));
    gc = XCreateGC(display, win, 0, NULL);

    XMapWindow(display, win);
    XFlush(display);
    XSelectInput(display, win, ExposureMask | KeyPressMask);
}

void ClearScreen() {
    XSetForeground(display, gc, BlackPixel(display, DefaultScreen(display)));
    XFillRectangle(display, pixmap, gc, 0, 0, width, height);
}

void DrawCircle(float x, float y) {
    int i = (int)(x * width);
    int j = height - (int)(y * height);
    int radius = 2;
    XSetForeground(display, gc, WhitePixel(display, DefaultScreen(display)));
    XFillArc(display, pixmap, gc, i - radius, j - radius, 2 * radius, 2 * radius, 0, 360 * 64);
}

void RefreshDisplay() {
    XCopyArea(display, pixmap, win, gc, 0, 0, width, height, 0, 0);
    XFlush(display);
}

void CheckForQuit() {
    while (XPending(display)) {
        XNextEvent(display, &report);
        if (report.type == KeyPress && XLookupKeysym(&report.xkey, 0) == XK_q) {
            quit_flag = 1;
        }
    }
}

void CloseGraphics() {
    XFreeGC(display, gc);
    XFreePixmap(display, pixmap);
    XDestroyWindow(display, win);
    XCloseDisplay(display);
}