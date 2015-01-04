#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef __APPLE__
#include <SDL/SDL.h>
#else
#include <SDL.h>
#endif
#include <vector>
#include <stdint.h>
using namespace std;

typedef Uint32 color_t;
#define RGB(r, g, b) ((color_t)((((unsigned)(r) & 0xFF) << 16) | (((unsigned)(g) & 0xFF) << 8) | (((unsigned)(b) & 0xFF))))
#define GetRValue(c) ((unsigned)((color_t)(c) & 0xFF0000) >> 16)
#define GetGValue(c) ((unsigned)((color_t)(c) & 0xFF00) >> 8)
#define GetBValue(c) ((unsigned)((color_t)(c) & 0xFF))

static Uint8 gammaCorrectionTable[256];

static void initGammaCorrectionTable()
{
    for(size_t i = 0; i < sizeof(gammaCorrectionTable) / sizeof(gammaCorrectionTable[0]); i++)
    {
        double v = (double)i / 0xFF;
        const double a = 0.055;
        if(v <= 0.0031308)
            v *= 12.92;
        else
            v = (1 + a) * pow(v, 1 / 2.4) - a;
        v *= 0xFF;
        v = floor(v + 0.5);
        if(v > 0xFF)
            v = 0xFF;
        else if(v < 0)
            v = 0;
        gammaCorrectionTable[i] = (int)v;
    }
}

static inline color_t correctGamma(color_t c)
{
    return RGB(gammaCorrectionTable[GetRValue(c)], gammaCorrectionTable[GetGValue(c)], gammaCorrectionTable[GetBValue(c)]);
}

const float eps = 5e-6, NearDist = 0.01;

struct polygon_t
{
    float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    color_t color;
};

typedef int_fast32_t fast_int;

fast_int ScreenXRes = 640, ScreenYRes = 480;
fast_int SuperSampleFactor = 1;
const fast_int size = 8;
#define YScaleFactor ((float)ScreenXRes / ScreenYRes)

struct span_t // span [start, end]
{
    color_t value;
    fast_int start, end;
    float startz, deltaz; // z values are really 1 / z so value for x == start + dx is really 1 / (startz + dx * deltaz)
    span_t * next;
};

fast_int sbufferwidth = ScreenXRes / SuperSampleFactor + 1;
span_t * sbuffer = NULL;
span_t ** freeliststart = NULL;
fast_int * freeindexstart = NULL;
span_t ** drawliststart = NULL;

void onresize()
{
    sbufferwidth = ScreenXRes / SuperSampleFactor + 1;
    delete []sbuffer;
    delete []freeliststart;
    delete []freeindexstart;
    delete []drawliststart;
    sbuffer = new span_t[sbufferwidth * ScreenYRes];
    freeliststart = new span_t *[ScreenYRes];
    freeindexstart = new fast_int[ScreenYRes];
    drawliststart = new span_t *[ScreenYRes];
}

void oninit()
{
    initGammaCorrectionTable();
    onresize();
}

void initsbuffer()
{
    for(fast_int y = 0; y < ScreenYRes; y++)
    {
        drawliststart[y] = &sbuffer[0 + y * sbufferwidth];
        sbuffer[0 + y * sbufferwidth].next = NULL;
        sbuffer[0 + y * sbufferwidth].start = 0;
        sbuffer[0 + y * sbufferwidth].end = ScreenXRes - 1;
        sbuffer[0 + y * sbufferwidth].startz = 0;
        sbuffer[0 + y * sbufferwidth].deltaz = 0;
        sbuffer[0 + y * sbufferwidth].value = RGB(0, 0, 0); // background
        freeliststart[y] = NULL;
        freeindexstart[y] = 1;
    }
}

inline span_t * allocspan(fast_int y)
{
    span_t * retval = freeliststart[y];
    if(retval)
    {
        freeliststart[y] = retval->next;
    }
    else if(freeindexstart[y] < sbufferwidth)
    {
        retval = &sbuffer[freeindexstart[y]++ + y * sbufferwidth];
    }
    return retval;
}

inline void freespan(fast_int y, span_t * span)
{
    span->next = freeliststart[y];
    freeliststart[y] = span;
}

inline void insertspan(fast_int start, fast_int end, fast_int y, float startz, float deltaz, color_t value) // all z values are really 1 / z so all compares are inverted
{
    span_t * curspan = drawliststart[y];
    span_t ** curspanp = &drawliststart[y];
    span_t * nextnewspanlist = NULL;
    bool donewithnewspan = false;
    while(curspan)
    {
        if(nextnewspanlist && donewithnewspan)
        {
            donewithnewspan = false;
            span_t * curspan = nextnewspanlist;
            nextnewspanlist = curspan->next;
            start = curspan->start;
            end = curspan->end;
            startz = curspan->startz;
            freespan(y, curspan);
        }
        if(donewithnewspan)
            return;
        if(curspan->end < start) // curspan is to left
        {
            // combine spans while we're at it
            span_t * nextspan = curspan->next;
#if 1
            if(nextspan)
            {
                if(nextspan->end < start && curspan->deltaz == nextspan->deltaz && curspan->value == nextspan->value && fabs(curspan->startz + curspan->deltaz * (nextspan->start - curspan->start) - nextspan->startz) < eps)
                {
                    curspan->end = nextspan->end;
                    curspan->next = nextspan->next;
                    freespan(y, nextspan);
                    continue;
                }
            }
#endif
            curspanp = &curspan->next;
            curspan = curspan->next;
            continue;
        }
        if(curspan->start > end) // curspan is to right
        {
            span_t * newspan = allocspan(y);
            newspan->next = curspan;
            *curspanp = newspan;
            newspan->start = start;
            newspan->end = end;
            newspan->startz = startz;
            newspan->deltaz = deltaz;
            newspan->value = value;
            donewithnewspan = true;
            continue;
        }

        if(curspan->start < start)
        {
            // split curspan in between newspan->start - 1 and newspan->start
            // with curspan becoming right side
            span_t * leftspan = allocspan(y);
            *curspanp = leftspan;
            curspanp = &leftspan->next;
            leftspan->next = curspan;
            leftspan->start = curspan->start;
            leftspan->end = start - 1;
            leftspan->startz = curspan->startz;
            leftspan->deltaz = curspan->deltaz;
            leftspan->value = curspan->value;
            curspan->startz += curspan->deltaz * (start - curspan->start);
            curspan->start = start;
        }
        else if(curspan->start > start)
        {
            // split newspan in between curspan->start - 1 and curspan->start
            // with newspan becoming left side
            span_t * rightspan = allocspan(y);
            rightspan->next = nextnewspanlist;
            nextnewspanlist = rightspan;
            rightspan->start = curspan->start;
            rightspan->end = end;
            rightspan->startz += deltaz * (curspan->start - start);
            end = curspan->start - 1;
        }
        // now curspan->start == start

        if(curspan->end > end)
        {
            // split curspan in between newspan->end and newspan->end + 1
            // with curspan becoming left side
            span_t * rightspan = allocspan(y);
            rightspan->next = curspan->next;
            curspan->next = rightspan;
            rightspan->start = end + 1;
            rightspan->end = curspan->end;
            rightspan->startz = curspan->startz + curspan->deltaz * (rightspan->start - curspan->start); // for end + 1
            rightspan->deltaz = curspan->deltaz;
            rightspan->value = curspan->value;
            curspan->end = end;
            continue;
        }
        else if(curspan->end < end)
        {
            // split newspan in between curspan->end and curspan->end + 1
            span_t * rightspan = allocspan(y);
            rightspan->next = nextnewspanlist;
            nextnewspanlist = rightspan;
            rightspan->start = curspan->end + 1;
            rightspan->end = end;
            rightspan->startz = startz + deltaz * (rightspan->start - start);
            end = curspan->end;
            continue;
        }
        // now curspan->end == end

        float ns_endz = startz + deltaz * (end - start);
        float cs_endz = curspan->startz + curspan->deltaz * (curspan->end - curspan->start);
        if(curspan->startz >= startz && cs_endz >= ns_endz) // new span is behind curspan
        {
            donewithnewspan = true;
            continue;
        }
        if(curspan->startz < startz && cs_endz < ns_endz) // new span is in front of curspan
        {
            // replace curspan with newspan
            curspan->start = start;
            curspan->end = end;
            curspan->startz = startz;
            curspan->deltaz = deltaz;
            curspan->value = value;
            donewithnewspan = true;
            continue;
        }
        // curspan and newspan intersect each other
        float divisor = curspan->deltaz - deltaz;
        // divisor can't be 0 because it would make *_endz equal
        float numerator = startz - deltaz * start - curspan->startz + curspan->deltaz * curspan->start;
        float intx = numerator / divisor;
        if(curspan->startz < startz) // newspan start is in front of curspan start
        {
            // left side is newspan, right side is curspan
            fast_int endx = (fast_int)floor(intx);
            fast_int startx = endx + 1;
            if(endx < start) // left side is null
            {
                donewithnewspan = true;
                continue;
            }
            if(startx > curspan->end) // right side is null
            {
                // replace curspan with newspan
                curspan->start = start;
                curspan->end = end;
                curspan->startz = startz;
                curspan->deltaz = deltaz;
                curspan->value = value;
                donewithnewspan = true;
                continue;
            }
            span_t * newspan = allocspan(y);
            *curspanp = newspan;
            newspan->next = curspan;
            newspan->start = start;
            newspan->end = endx;
            newspan->startz = startz;
            newspan->deltaz = deltaz;
            newspan->value = value;
            curspan->startz += curspan->deltaz * (startx - curspan->start);
            curspan->start = startx;
            donewithnewspan = true;
            continue;
        }
        else // new span start is behind curspan start
        {
            // left side is curspan, right side is newspan
            fast_int endx = (fast_int)floor(intx);
            fast_int startx = endx + 1;
            if(endx < curspan->start) // left side is null
            {
                // replace curspan with newspan
                curspan->start = start;
                curspan->end = end;
                curspan->startz = startz;
                curspan->deltaz = deltaz;
                curspan->value = value;
                donewithnewspan = true;
                continue;
            }
            if(startx > end) // right side is null
            {
                donewithnewspan = true;
                continue;
            }
            span_t * newspan = allocspan(y);
            newspan->next = curspan->next;
            curspan->next = newspan;
            newspan->start = startx;
            newspan->end = end;
            newspan->startz = startz + deltaz * (startx - start);
            newspan->deltaz = deltaz;
            newspan->value = value;
            curspan->end = endx;
            donewithnewspan = true;
            continue;
        }
    }

    while(!donewithnewspan)
    {
        if(nextnewspanlist && donewithnewspan)
        {
            donewithnewspan = false;
            span_t * curspan = nextnewspanlist;
            nextnewspanlist = curspan->next;
            start = curspan->start;
            end = curspan->end;
            startz = curspan->startz;
            freespan(y, curspan);
        }
        *curspanp = curspan = allocspan(y);
        curspan->start = start;
        curspan->end = end;
        curspan->startz = startz;
        curspan->deltaz = deltaz;
        curspan->value = value;
        curspan->next = NULL;
        donewithnewspan = true;
    }
}

inline void drawscanline(vector<span_t *> &curspans, SDL_Surface * screen, fast_int sy) // screen must be locked
{
    Uint8 * startaddr = (Uint8 *)screen->pixels;
    startaddr += screen->pitch * sy;
    fast_int bpp = screen->format->BytesPerPixel;

    for(fast_int dy = 0; dy < SuperSampleFactor; dy++)
        curspans[dy] = drawliststart[sy * SuperSampleFactor + dy];
    for(fast_int sx = 0; sx < screen->w; sx++, startaddr += bpp)
    {
        fast_int r = 0, g = 0, b = 0;
        fast_int startx = sx * SuperSampleFactor, endx = startx + SuperSampleFactor - 1;
        for(fast_int dy = 0; dy < SuperSampleFactor; dy++)
        {
            span_t *& curspan = curspans[dy];
            if(curspan && curspan->end < startx)
                curspan = curspan->next;
            while(curspan)
            {
                fast_int curR = GetRValue(curspan->value);
                fast_int curG = GetGValue(curspan->value);
                fast_int curB = GetBValue(curspan->value);
                fast_int spanStart = curspan->start;
                fast_int spanEnd = curspan->end;
                if(spanStart < startx)
                    spanStart = startx;
                if(spanEnd > endx)
                    spanEnd = endx;
                fast_int length = spanEnd - spanStart + 1;
                r += curR * length;
                g += curG * length;
                b += curB * length;
                if(curspan->end < endx)
                {
                    curspan = curspan->next;
                }
                else
                {
                    break;
                }
            }
        }
        r /= SuperSampleFactor * SuperSampleFactor;
        g /= SuperSampleFactor * SuperSampleFactor;
        b /= SuperSampleFactor * SuperSampleFactor;
        *(Uint32 *)startaddr = correctGamma(RGB(r, g, b));
    }
}

void lockSurface(SDL_Surface *s)
{
    while(0 != SDL_LockSurface(s))
    {
        timespec ts;
        ts.tv_nsec = 1000000;
        ts.tv_sec = 0;
        nanosleep(&ts, NULL);
    }
}

void drawsbuffer(vector<span_t *> &curspans, SDL_Surface * screen)
{
    lockSurface(screen);
    for(fast_int sy=0;sy<screen->h;sy++)
    {
        drawscanline(curspans, screen, sy);
    }
    SDL_UnlockSurface(screen);
}

static inline fast_int getabcd(float * a, float * b, float * c, float * d, float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3)
{
    *a = y3 * (z2 - z1) + y2 * (z1 - z3) + y1 * (z3 - z2);
    *b = x3 * (z1 - z2) + x1 * (z2 - z3) + x2 * (z3 - z1);
    *c = x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2);
    *d = -x3 * y2 * z1 + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 + x1 * y2 * z3;
    if(fabs(*d) < eps) return 0;
    return 1;
}

static inline void GetLineEq(float * a, float * b, float * c, float x1, float y1, float x2, float y2, float insidex, float insidey)
{
    *a = y1 - y2;
    *b = x2 - x1;
    *c = x1 * y2 - x2 * y1;
    if(*a * insidex + *b * insidey + *c < 0.0)
    {
        *a = -*a;
        *b = -*b;
        *c = -*c;
    }
}

void drawtri(const polygon_t * p)
{
    float a, b, c, d;
    float a1, b1, c1, a2, b2, c2, a3, b3, c3;
    float fminx, fmaxx, fminy, fmaxy;
    fast_int minx, miny, maxx, maxy;
    float b1_y, b2_y, b3_y, ia1_x, ia2_x, ia3_x;
    float dx, ix, x, dy, iy, y;
    fast_int xi, yi;
    float d1, d2, d3;
    float curz;
    fast_int finishx;
    float pos;
    float x1, y1, x2, y2, x3, y3;

    if(fabs(p->x1 - p->x2) < eps && fabs(p->y1 - p->y2) < eps && fabs(p->z1 - p->z2) < eps)
        return;

    if(fabs(p->x1 - p->x3) < eps && fabs(p->y1 - p->y3) < eps && fabs(p->z1 == p->z3) < eps)
        return;

    if(fabs(p->x3 - p->x2) < eps && fabs(p->y3 - p->y2) < eps && fabs(p->z3 - p->z2) < eps)
        return;

    if(p->z1 < NearDist || p->z2 < NearDist || p->z3 < NearDist)
        return;

    if(!getabcd(&a, &b, &c, &d, p->x1, p->y1, p->z1, p->x2, p->y2, p->z2, p->x3, p->y3, p->z3))
        return;

    x1 = ScreenXRes / 2 + ScreenXRes / 2 * (p->x1 / p->z1);
    y1 = ScreenYRes / 2 - ScreenYRes / 2 * (p->y1 / p->z1) * YScaleFactor;

    x2 = ScreenXRes / 2 + ScreenXRes / 2 * (p->x2 / p->z2);
    y2 = ScreenYRes / 2 - ScreenYRes / 2 * (p->y2 / p->z2) * YScaleFactor;

    x3 = ScreenXRes / 2 + ScreenXRes / 2 * (p->x3 / p->z3);
    y3 = ScreenYRes / 2 - ScreenYRes / 2 * (p->y3 / p->z3) * YScaleFactor;

    GetLineEq(&a1, &b1, &c1, x1, y1, x2, y2, x3, y3);
    GetLineEq(&a2, &b2, &c2, x2, y2, x3, y3, x1, y1);
    GetLineEq(&a3, &b3, &c3, x3, y3, x1, y1, x2, y2);

    fminx = x1;
    fmaxx = x1;
    fminy = y1;
    fmaxy = y1;
    if(x2 < fminx)
        fminx = x2;
    else if(x2 > fmaxx)
        fmaxx = x2;

    if(y2 < fminy)
        fminy = y2;
    else if(y2 > fmaxy)
        fmaxy = y2;

    if(x3 < fminx)
        fminx = x3;
    else if(x3 > fmaxx)
        fmaxx = x3;

    if(y3 < fminy)
        fminy = y3;
    else if(y3 > fmaxy)
        fmaxy = y3;

    if(fminx < 0)
        fminx = 0;

    if(fmaxx > ScreenXRes - 1.0)
        fmaxx = ScreenXRes - 1.0;

    if(fminy < 0)
        fminy = 0;

    if(fmaxy > ScreenYRes - 1.0)
        fmaxy = ScreenYRes - 1.0;

    minx = (fast_int)floor(fminx);
    maxx = (fast_int)ceil(fmaxx);
    miny = (fast_int)floor(fminy);
    maxy = (fast_int)ceil(fmaxy);

    b1_y = miny * b1;
    b2_y = miny * b2;
    b3_y = miny * b3;

    ia1_x = minx * a1;
    ia2_x = minx * a2;
    ia3_x = minx * a3;

    a /= -d;
    b /= -d;
    c /= -d;
    d = -1.0;

    dx = 2.0 * a / ScreenXRes;
    ix = 2.0 * a * (minx - ScreenXRes / 2.0) / ScreenXRes;
    dy = -2.0 * b / ScreenYRes / YScaleFactor;
    iy = -2.0 * b * (miny - ScreenYRes / 2.0) / ScreenYRes / YScaleFactor;
    y = iy + c;

    for(yi = miny; yi <= maxy; yi++)
    {
        d1 = ia1_x + b1_y + c1;
        d2 = ia2_x + b2_y + c2;
        d3 = ia3_x + b3_y + c3;
        x = ix + y;
        curz = x;
        xi = minx;
        if(a1 > 0 && d1 < 0)
        {
            pos = ceil(d1 / -a1) + minx;
            if(xi < pos)
            {
                if(pos < ScreenXRes)
                    xi = (fast_int)pos;
                else
                    xi = ScreenXRes;
            }
        }
        if(a2 > 0 && d2 < 0)
        {
            pos = ceil(d2 / -a2) + minx;
            if(xi < pos)
            {
                if(pos < ScreenXRes)
                    xi = (fast_int)pos;
                else
                    xi = ScreenXRes;
            }
        }
        if(a3 > 0 && d3 < 0)
        {
            pos = ceil(d3 / -a3) + minx;
            if(xi < pos)
            {
                if(pos < ScreenXRes)
                    xi = (fast_int)pos;
                else
                    xi = ScreenXRes;
            }
        }
        d1 += a1 * (xi - minx);
        d2 += a2 * (xi - minx);
        d3 += a3 * (xi - minx);
        for(; xi <= maxx; xi++)
        {
            if(d1 < 0 && a1 <= 0)
                break;
            if(d2 < 0 && a2 <= 0)
                break;
            if(d3 < 0 && a3 <= 0)
                break;
            if(d1 >= 0 && d2 >= 0 && d3 >= 0)
                break;
            d1 += a1;
            d2 += a2;
            d3 += a3;
        }
        if((d1 > 0 || a1 > 0) && (d2 > 0 || a2 > 0) && (d3 > 0 || a3 > 0))
        {
            curz += dx * (xi - minx);

            finishx = maxx;
            if(a1 < 0 && d1 > 0)
            {
                pos = floor(d1 / -a1) + xi;
                if(finishx > pos)
                {
                    if(pos > xi)
                        finishx = (fast_int)pos;
                    else
                        finishx = xi;
                }
            }
            if(a2 < 0 && d2 > 0)
            {
                pos = floor(d2 / -a2) + xi;
                if(finishx > pos)
                {
                    if(pos > xi)
                        finishx = (fast_int)pos;
                    else
                        finishx = xi;
                }
            }
            if(a3 < 0 && d3 > 0)
            {
                pos = floor(d3 / -a3) + xi;
                if(finishx > pos)
                {
                    if(pos > xi)
                        finishx = (fast_int)pos;
                    else
                        finishx = xi;
                }
            }

            insertspan(xi, finishx, yi, curz, dx, p->color);
        }

        b1_y += b1;
        b2_y += b2;
        b3_y += b3;
        y += dy;
    }
}

float frand()
{
    return (float)rand() / RAND_MAX * 2 - 1;
}

void rotatex(float angle, float * x, float * y, float * z)
{
    float s = sin(angle);
    float c = cos(angle);
    float t = *y * c - *z * s;
    *z = *z * c + *y * s;
    *y = t;
}

void rotatey(float angle, float * x, float * y, float * z)
{
    float s = sin(angle);
    float c = cos(angle);
    float t = *x * c - *z * s;
    *z = *z * c + *x * s;
    *x = t;
}

float angle = 0;

void orientpoint(float * x, float * y, float * z)
{
    rotatex(angle * 2, x, y, z);
    rotatey(angle, x, y, z);
    *z += size * 5;
}

inline void CrossProduct(float * dx, float * dy, float * dz, float ax, float ay, float az, float bx, float by, float bz)
{
    *dx = ay * bz - az * by;
    *dy = az * bx - ax * bz;
    *dz = ax * by - ay * bx;
}

void drawface(float x00, float y00, float z00, float x01, float y01, float z01, float x10, float y10, float z10, color_t color)
{
    orientpoint(&x00, &y00, &z00);
    orientpoint(&x01, &y01, &z01);
    orientpoint(&x10, &y10, &z10);
    float avgx = (x01 + x10) * 0.5, avgy = (y01 + y10) * 0.5, avgz = (z01 + z10) * 0.5;
    float nx, ny, nz;
    CrossProduct(&nx, &ny, &nz, x01 - x00, y01 - y00, z01 - z00, x10 - x00, y10 - y00, z10 - z00);
    if(nx * avgx + ny * avgy + nz * avgz <= 0) return;
    const float lx = -1, ly = 1, lz = -2;
    float divisor = -sqrt((nx * nx + ny * ny + nz * nz) * (lx * lx + ly * ly + lz * lz));
    float d = (nx * lx + ny * ly + nz * lz) / divisor;
    if(d < 0) d = 0;
    d *= 0.6;
    d += 0.4;
    if(d > 1) d = 1;
    fast_int r = (fast_int)(GetRValue(color) * d), g = (fast_int)(GetGValue(color) * d), b = (fast_int)(GetBValue(color) * d);
    color = RGB(r, g, b);
    polygon_t p;
    p.x1 = x00;
    p.y1 = y00;
    p.z1 = z00;
    p.x2 = x01;
    p.y2 = y01;
    p.z2 = z01;
    p.x3 = x10;
    p.y3 = y10;
    p.z3 = z10;
    p.color = color;
    drawtri(&p);
    float x11 = x10 + (x01 - x00);
    float y11 = y10 + (y01 - y00);
    float z11 = z10 + (z01 - z00);
    p.x1 = x10;
    p.y1 = y10;
    p.z1 = z10;
    p.x2 = x01;
    p.y2 = y01;
    p.z2 = z01;
    p.x3 = x11;
    p.y3 = y11;
    p.z3 = z11;
    drawtri(&p);
}

void drawcube(fast_int facemask, float minx, float maxx, float miny, float maxy, float minz, float maxz, color_t color)
{
    if(facemask & 0x01) drawface(minx, miny, minz, minx, maxy, minz, minx, miny, maxz, color);
    if(facemask & 0x02) drawface(maxx, miny, minz, maxx, miny, maxz, maxx, maxy, minz, color);
    if(facemask & 0x04) drawface(minx, miny, minz, minx, miny, maxz, maxx, miny, minz, color);
    if(facemask & 0x08) drawface(minx, maxy, minz, maxx, maxy, minz, minx, maxy, maxz, color);
    if(facemask & 0x10) drawface(minx, miny, minz, maxx, miny, minz, minx, maxy, minz, color);
    if(facemask & 0x20) drawface(minx, miny, maxz, minx, maxy, maxz, maxx, miny, maxz, color);
}

double Timer()
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char ** argv)
{
    // initialize SDL video
    if(SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        printf("Unable to init SDL: %s\n", SDL_GetError());
        return 1;
    }

    // make sure SDL cleans up before exit
    atexit(SDL_Quit);

    // create a new window
#ifdef DEBUG
    SDL_Surface * screen = SDL_SetVideoMode(ScreenXRes, ScreenYRes, 32, /*SDL_FULLSCREEN |*/ SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_RESIZABLE);
#else
    SDL_Surface * screen = SDL_SetVideoMode(ScreenXRes, ScreenYRes, 32, /*SDL_FULLSCREEN |*/ SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_RESIZABLE);
#endif
    if(!screen)
    {
        printf("Unable to set %ix%i video: %s\n", (int)ScreenXRes, (int)ScreenYRes, SDL_GetError());
        return 1;
    }
    ScreenXRes *= SuperSampleFactor;
    ScreenYRes *= SuperSampleFactor;
    vector<span_t *> curspans;
    curspans.resize(SuperSampleFactor);

    oninit();

    double starttime = Timer();

    // program main loop
    bool done = false;
    fast_int frames = 0;
    while(!done)
    {
        frames++;
        double curtime = Timer();
        if(curtime - 0.1 > starttime)
        {
            float fps = frames / (curtime - starttime);
            frames = 0;
            starttime = curtime;
#if defined(DEBUG) || 1
            printf("%ix%i : %ix : FPS : %g          \r", screen->w, screen->h, (int)SuperSampleFactor, fps);
            fflush(stdout);
#endif
        }
        // message processing loop
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            // check for messages
            switch(event.type)
            {
                // exit if the window is closed
            case SDL_QUIT:
                done = true;
                break;

                // check for keypresses
            case SDL_KEYDOWN:
            {
                // exit if ESCAPE is pressed
                if(event.key.keysym.sym == SDLK_ESCAPE)
                    done = true;
                if(event.key.keysym.sym >= SDLK_0 && event.key.keysym.sym <= SDLK_9)
                {
                    SuperSampleFactor = event.key.keysym.sym - SDLK_0;
                    if(SuperSampleFactor == 0)
                        SuperSampleFactor = 16;
                    ScreenXRes = screen->w;
                    ScreenYRes = screen->h;
                    ScreenXRes *= SuperSampleFactor;
                    ScreenYRes *= SuperSampleFactor;
                    curspans.resize(SuperSampleFactor);
                    onresize();
                }
                break;
            }

            case SDL_VIDEORESIZE:
            {
                ScreenXRes = event.resize.w;
                ScreenYRes = event.resize.h;
                screen = SDL_SetVideoMode(ScreenXRes, ScreenYRes, 32, SDL_HWSURFACE | SDL_DOUBLEBUF | SDL_RESIZABLE);
                ScreenXRes *= SuperSampleFactor;
                ScreenYRes *= SuperSampleFactor;
                curspans.resize(SuperSampleFactor);
                onresize();
                break;
            }
            } // end switch
        } // end of message processing

        initsbuffer();
        srand(2);
        angle = fmod(Timer() / 16, 2 * M_PI);
        for(fast_int x=-size;x<size;x++)
        {
            for(fast_int y=-size;y<size;y++)
            {
                for(fast_int z=-size;z<size;z++)
                {
                    fast_int facemask = 0;
                    if(x == -size) facemask |= 0x01;
                    if(x == size - 1) facemask |= 0x02;
                    if(y == -size) facemask |= 0x04;
                    if(y == size - 1) facemask |= 0x08;
                    if(z == -size) facemask |= 0x10;
                    if(z == size - 1) facemask |= 0x20;
                    //facemask |= 0x3F;
                    if(facemask == 0) {z = size - 1; z--; continue;} // skip invisible space
                    drawcube(facemask, x, x + 1, y, y + 1, z, z + 1, RGB((x * 127 / size + 128) & 0xFF, (y * 127 / size + 128) & 0xFF, (z * 127 / size + 128) & 0xFF));
                }
            }
        }

        drawsbuffer(curspans, screen);
        SDL_Flip(screen);
    } // end main loop

    printf("Exited cleanly\n");
    return 0;
}
