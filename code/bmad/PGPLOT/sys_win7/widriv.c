/*
 * PGPLOT driver Version 0.1 - Mar 13, 2001
 *     for native Win32
 *     by  "Thomas Ott" <ott@mpe.mpg.de>
 *
 *     WWW: http://www.mpe.mpg.de/~tho
 *
 * Version 0.2 - Apr 14, 2002
 *   Fixed the following memory leaks:
 *    - Always call EndPaint in reply to WM_PAINT message
 *    - Add call to ReleaseDC in ifunc = 16
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <windows.h>

#define WI_IDENT "PGPLOT /windows"       /* Name to prefix messages to user */

#define MAXCOLOR	256
#define MAXCOLOR1	(MAXCOLOR-1)

/* Default size of plot window */
#define DEFX 640
#define DEFY 480

/* globals */
int firsttime = 1;
int painting = 0;
int wantcursor = 0;
int rubberband = 0;
int cursorx, cursory, cursormode = 0, xmouseref, ymouseref;
int rubberx1, rubberx2, rubbery1, rubbery2, rubbermode;
int cursorkey;
int WINDOW_ID = 0;
unsigned char curcol = 0;
long npnt, nVertices = 0;
int dw, dh;
HWND windows[10];
HDC memdc[10], windowdc;
RECT windowsize[10];
HGDIOBJ defpen;
HANDLE Creation, Mouse, NotRepainting;
HBITMAP bitmaps[10];
COLORREF colors[MAXCOLOR];
POINT *polypoints;
HANDLE globMutex;

void drawRubberBand(HDC hdc, int fromx, int fromy, int tox, int toy, int mode) {
	switch (mode) {
		case 1:
			MoveToEx(hdc, fromx, fromy, NULL);
			LineTo(hdc, tox, toy);
			break;
		case 2:
			MoveToEx(hdc, fromx, fromy, NULL);
			LineTo(hdc, fromx, toy);
			LineTo(hdc, tox, toy);
			LineTo(hdc, tox, fromy);
			LineTo(hdc, fromx, fromy);
			break;
		case 3:
			break;
		case 4:
			break;
		case 5:
			MoveToEx(hdc, windowsize[WINDOW_ID].left, toy, NULL);
			LineTo(hdc, windowsize[WINDOW_ID].right, toy);
			break;
		case 6:
			MoveToEx(hdc, tox, windowsize[WINDOW_ID].top, NULL);
			LineTo(hdc, tox, windowsize[WINDOW_ID].bottom);
			break;
		case 7:
			MoveToEx(hdc, windowsize[WINDOW_ID].left, toy, NULL);
			LineTo(hdc, windowsize[WINDOW_ID].right, toy);
			MoveToEx(hdc, tox, windowsize[WINDOW_ID].top, NULL);
			LineTo(hdc, tox, windowsize[WINDOW_ID].bottom);
			break;
		default: break;
	}
	if (mode > 0) {
		if (rubberband == 1) rubberband = 0;
		else rubberband = 1;
		rubberx1 = fromx;
		rubberx2 = tox;
		rubbery1 = fromy;
		rubbery2 = toy;
		rubbermode = mode;
	}
}

void deleteRubberBand(HDC hdc) {
	if ((rubberx1 == -1) && (rubberx2 == -1) && (rubbery1 == -1) && (rubbery2 == -1)) return;
	switch (rubbermode) {
		case 1:
			MoveToEx(hdc, rubberx1, rubbery1, NULL);
			LineTo(hdc, rubberx2, rubbery2);
			break;
		case 2:
			MoveToEx(hdc, rubberx1, rubbery1, NULL);
			LineTo(hdc, rubberx1, rubbery2);
			LineTo(hdc, rubberx2, rubbery2);
			LineTo(hdc, rubberx2, rubbery1);
			LineTo(hdc, rubberx1, rubbery1);
			break;
		case 3:
			break;
		case 4:
			break;
		case 5:
			MoveToEx(hdc, windowsize[WINDOW_ID].left, rubbery2, NULL);
			LineTo(hdc, windowsize[WINDOW_ID].right, rubbery2);
			break;
		case 6:
			MoveToEx(hdc, rubberx2, windowsize[WINDOW_ID].top, NULL);
			LineTo(hdc, rubberx2, windowsize[WINDOW_ID].bottom);
			break;
		case 7:
			MoveToEx(hdc, windowsize[WINDOW_ID].left, rubbery2, NULL);
			LineTo(hdc, windowsize[WINDOW_ID].right, rubbery2);
			MoveToEx(hdc, rubberx2, windowsize[WINDOW_ID].top, NULL);
			LineTo(hdc, rubberx2, windowsize[WINDOW_ID].bottom);
			break;
		default: break;
	}
	rubberx1 = rubberx2 = rubbery1 = rubbery2 = -1;
}

LRESULT CALLBACK PgplotGraphProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	HDC hdc;
	PAINTSTRUCT ps;
	int i, id = -1;

	for (i = 1; i < 10; i++) if (hwnd == windows[i]) id = i;
	if (id == -1) return DefWindowProc(hwnd, message, wParam, lParam);

	switch(message) {
		case WM_PAINT:
			if (WaitForSingleObject(globMutex, 500) != WAIT_OBJECT_0) return 0;
			if (GetUpdateRect(hwnd, NULL, FALSE)) {
				if (memdc[id] != NULL) {
					if ((hdc = BeginPaint(hwnd, &ps)) != NULL) {
						BitBlt(hdc, 0, 0, windowsize[id].right, windowsize[id].bottom, memdc[id], 0, 0, SRCCOPY);
					}
					EndPaint(hwnd, &ps);
				}
			}
			ReleaseMutex(globMutex);
			return 0;
		case WM_MOUSEMOVE: if (wantcursor) {
			deleteRubberBand(windowdc);
			cursorx = LOWORD(lParam);
			cursory = HIWORD(lParam);
			drawRubberBand(windowdc, xmouseref, ymouseref, cursorx, cursory, cursormode);
			return 0;
		}
			break;
		case WM_CHAR: if (wantcursor) {
			cursorkey = (TCHAR)wParam;
			cursorx = LOWORD(lParam);
			cursory = HIWORD(lParam);
			SetEvent(Mouse);
			return 0;
		}
			break;
		case WM_LBUTTONDOWN: if (wantcursor) {
			cursorkey = 'A';
			cursorx = LOWORD(lParam);
			cursory = HIWORD(lParam);
			SetEvent(Mouse);
			return 0;
		}
			break;
		case WM_MBUTTONDOWN: if (wantcursor) {
			cursorkey = 'D';
			cursorx = LOWORD(lParam);
			cursory = HIWORD(lParam);
			SetEvent(Mouse);
			return 0;
		}
			break;
		case WM_RBUTTONDOWN: if (wantcursor) {
			cursorkey = 'X';
			cursorx = LOWORD(lParam);
			cursory = HIWORD(lParam);
			SetEvent(Mouse);
			return 0;
		}
			break;
		case WM_DESTROY: if (wantcursor) {
			cursorkey = 'X';
			SetEvent(Mouse);
		}
			windows[id] = NULL;
			DeleteObject(bitmaps[id]);
			bitmaps[id] = NULL;
			if (memdc != NULL) DeleteDC(memdc[id]);
			memdc[id] = NULL;
			return 0;
		default: break;
	}
	return DefWindowProc(hwnd, message, wParam, lParam);
}

DWORD WINAPI GraphProc(LPVOID params) {
	MSG msg;
	char title[30];
	RECT r1, r2;

	sprintf(title, "PGPLOT Window %i", WINDOW_ID);
	windows[WINDOW_ID] = CreateWindow("PGPLOT", title, WS_OVERLAPPEDWINDOW|WS_VISIBLE, 1, 1, DEFX, DEFY, NULL, NULL, NULL, NULL);
	ShowWindow(windows[WINDOW_ID], SW_SHOWNORMAL);
//	GetClientRect(windows[WINDOW_ID], &windowsize[WINDOW_ID]);
	GetClientRect(windows[WINDOW_ID], &r1);
	GetWindowRect(windows[WINDOW_ID], &r2);
	dw = r2.right - r2.left - (r1.right - r1.left);
	dh = r2.bottom - r2.top - (r1.bottom - r1.top);

	SetEvent(Creation);
// 	waiting = 0;

	while( GetMessage( &msg, NULL, 0, 0 ) ) {
		TranslateMessage( &msg );
		DispatchMessage( &msg );
	}

	return 0;
}

static void Initial()
{
	int i;
	WNDCLASS wndclass;

	for (i = 0; i < MAXCOLOR; i++) colors[i] = i + i * 256 + i * 256 * 256;
	colors[0]  = 0x00000000;  // current background color
	colors[1]  = 0x00ffffff;  // current foreground color
	colors[2]  = 0x000000ff;  // Red  
	colors[3]  = 0x0000ff00;  // Green
	colors[4]  = 0x00ff0000;  // Blue 
	colors[5]  = 0x00ffff00;  // Cyan (Green + Blue)
	colors[6]  = 0x00ff00ff;  // Magenta (Red + Blue)
	colors[7]  = 0x0000ffff;  // Yellow (Red + Green)
	colors[8]  = 0x00007fff;  // Orange (Red + Yellow)
	colors[9]  = 0x0000ff7f;  // Green + Yellow
	colors[10] = 0x007fff00;  // Green + Cyan
	colors[11] = 0x00ff7f00;  // Blue + Cyan
	colors[12] = 0x00ff007f;  // Blue + Magenta
	colors[13] = 0x007f00ff;  // Red + Magenta
	colors[14] = 0x00555555;  // Dark Gray
	colors[15] = 0x00aaaaaa;  // Light Gray

	wndclass.style = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc = &PgplotGraphProc;
	wndclass.cbClsExtra = 0;
	wndclass.cbWndExtra = 2 * sizeof(void FAR *);
	wndclass.hInstance = NULL;
	wndclass.hIcon = LoadIcon(NULL, IDI_WINLOGO);
	wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndclass.hbrBackground = (HBRUSH)GetStockObject(BLACK_BRUSH);
	wndclass.lpszMenuName = NULL;
	wndclass.lpszClassName = "PGPLOT";
	RegisterClass(&wndclass);

	NotRepainting = CreateEvent(NULL, TRUE, TRUE, NULL);
	globMutex = CreateMutex(NULL, FALSE, "PGPLOTmutex");

	for (i = 0; i < 10; i++) {
		memdc[i] = NULL;
		bitmaps[i] = NULL;
	}

	firsttime = 0;
}

void widriv_(int *opcode,float  rbuf[], int *nbuf, char *chr, int *lchr,int *Mode,int  len)

{
	int i;
	if((WINDOW_ID < 0) && (*opcode != 1) && (*opcode != 2) && (*opcode != 3) && (*opcode != 4) && (*opcode != 5) && (*opcode != 7) && (*opcode != 8) && (*opcode != 9) && (*opcode != 27)) return;

/*
 * Branch on the specified PGPLOT opcode.
 */
	switch(*opcode) {

/*--- IFUNC=1, Return device name ---------------------------------------*/
	case 1:
		sprintf(chr, "WINDOWS (Win32 native window)");
		*lchr = strlen(chr);
		for(i = *lchr; i < len; i++)
			chr[i] = ' ';
		break;
/*--- IFUNC=2, Return maximum dimensions of view surface, and range of color index*/
	case 2:
		rbuf[0] = 0.0;
		rbuf[1] = -1.0;  /* Report no effective max plot width */
		rbuf[2] = 0.0;
		rbuf[3] = -1.0;  /* Report no effective max plot height */
		rbuf[4] = 0.0;
		rbuf[5] = (float)MAXCOLOR1;
		*nbuf = 6;
		break;
/*--- IFUNC=3, Return device scale ----------------------------------------------*/
	case 3:
		rbuf[0] = 75.0;	/* assuming landscape paper setting */
		rbuf[1] = 75.0;
		rbuf[2] = 1.0;		/* Device coordinates per pixel */
		*nbuf = 3;
		break;
/*--- IFUNC=4, Return device capabilities ---------------------------------------*/
	case 4:
		chr[0] = 'I'; /* Interactive device */
		chr[1] = 'C'; /* Cursor is available */
		chr[2] = 'N'; /* Dashed lines not available */
		chr[3] = 'A'; /* Area fill available */
		chr[4] = 'N'; /* Thick lines not available */
		chr[5] = 'R'; /* Rectangle fill available */
		chr[6] = 'P'; /* Line of pixels is available */
		chr[7] = 'N'; /* PGPLOT window is not deleted when closed */
		chr[8] = 'Y'; /* Can return color representation */
		chr[9] = 'N'; /* Not used */
		chr[10]= 'N'; /* Area-scroll not available */
		*lchr = 11;
		break;
/*--- IFUNC=5, Return default device/file name ----------------------------------*/
	case 5:
		chr[0] = ' ';  /* Default name is "" */
		*lchr = 0;
		break;
/*--- IFUNC=6, Return default size of view surface ------------------------------*/
	case 6:
		if (!IsIconic(windows[WINDOW_ID])) GetClientRect(windows[WINDOW_ID], &windowsize[WINDOW_ID]);
		rbuf[0] = 0.0;
		rbuf[1] = (float)windowsize[WINDOW_ID].right;
		rbuf[2] = 0.0;
		rbuf[3] = (float)windowsize[WINDOW_ID].bottom;
		*nbuf = 4;

		break;
/*--- IFUNC=7, Return miscellaneous defaults ------------------------------------*/
	case 7:
		rbuf[0] = 1.0;
		*nbuf = 1;
		break;
/*--- IFUNC=8, Select device ----------------------------------------------------*/
	case 8:
		WINDOW_ID =(int)(rbuf[1]+0.5);
		break;
/*--- IFUNC=9, Open workstation -------------------------------------------------*/
	case 9:
		if (firsttime) Initial();
/*
 * Assign the returned device unit number and success indicator.
 * Assume failure to open until the workstation is open.
 */
		rbuf[0] = rbuf[1] = 0.0;
		*nbuf = 2;
/*
 * Prepare the display name.
 */
		if(*lchr >= len) {
			fprintf(stderr, "%s: Display name too long.\n", WI_IDENT);
			return;
		} else {
			chr[*lchr] = '\0';
		}
		WINDOW_ID = WINDOW_ID + 1;
		if (WINDOW_ID == 0) WINDOW_ID = 1;
		if (!IsWindow(windows[WINDOW_ID])) {
/*
 * Create the window.
 */
			DWORD id;


			WaitForSingleObject(globMutex, INFINITE);
			Creation = CreateEvent(NULL, TRUE, FALSE, NULL);
			if (CreateThread(NULL, 0, GraphProc, NULL, 0, &id) == NULL) printf("Could not start thread\n");
			else {
				WaitForSingleObject(Creation, INFINITE);
				curcol = 0;
				rbuf[0] = (float)WINDOW_ID;
				rbuf[1] = 1.0;
				*nbuf = 2;
			}
			CloseHandle(Creation);
			ReleaseMutex(globMutex);
		} else {
/*
 * Activate the window
 */
			curcol = 0;
			rbuf[0] = (float)WINDOW_ID;
			rbuf[1] = 1.0;
			*nbuf = 2;
		}
		break;
/*--- IFUNC=10, Close workstation -----------------------------------------------*/
	case 10:
		WINDOW_ID = -1;
		painting = 0;
		break;
/*--- IFUNC=11, Begin picture ---------------------------------------------------*/
	case 11:
	{
		HDC dc;

		WaitForSingleObject(globMutex, INFINITE);
		if (IsIconic(windows[WINDOW_ID])) ShowWindow(windows[WINDOW_ID], SW_RESTORE);
//		InvalidateRect(windows[WINDOW_ID], NULL, TRUE);
		painting = 1;
		SetWindowPos(windows[WINDOW_ID], HWND_NOTOPMOST, 0, 0, (int)(rbuf[0] + 0.5) + dw, (int)(rbuf[1] + 0.5) + dh, SWP_SHOWWINDOW|SWP_NOMOVE);
		GetClientRect(windows[WINDOW_ID], &windowsize[WINDOW_ID]);
		dc = GetDC(windows[WINDOW_ID]);
		if (memdc[WINDOW_ID] == NULL) {
			memdc[WINDOW_ID] = CreateCompatibleDC(dc);
		}
		defpen = SelectObject(memdc[WINDOW_ID], CreatePen(PS_SOLID, 0, colors[1]));
		if (bitmaps[WINDOW_ID] != NULL) DeleteObject(bitmaps[WINDOW_ID]);
		bitmaps[WINDOW_ID] = CreateCompatibleBitmap(dc, windowsize[WINDOW_ID].right, windowsize[WINDOW_ID].bottom);
		ReleaseDC(windows[WINDOW_ID], dc);
		SelectObject(memdc[WINDOW_ID], bitmaps[WINDOW_ID]);
		FillRect(memdc[WINDOW_ID], &windowsize[WINDOW_ID], GetStockObject(BLACK_BRUSH));
		ReleaseMutex(globMutex);
		break;
	}
/*--- IFUNC=12, Draw line -------------------------------------------------------*/
	case 12:
		WaitForSingleObject(globMutex, INFINITE);
//		if (MoveToEx(memdc[WINDOW_ID], (int)(rbuf[0] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[1] + 0.5), NULL) == 0) printf("MoveToEx failed\n");
//		if (LineTo(memdc[WINDOW_ID], (int)(rbuf[2] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[3] + 0.5)) == 0) printf("LineTo failed from %i to %i\n", (int)(rbuf[2] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[3] + 0.5));
		MoveToEx(memdc[WINDOW_ID], (int)(rbuf[0] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[1] + 0.5), NULL);
		LineTo(memdc[WINDOW_ID], (int)(rbuf[2] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[3] + 0.5));
		ReleaseMutex(globMutex);
		break;
/*--- IFUNC=13, Draw dot --------------------------------------------------------*/
	case 13:
		WaitForSingleObject(globMutex, INFINITE);
		SetPixel(memdc[WINDOW_ID], (int)(rbuf[0] + 0.5), windowsize[WINDOW_ID].bottom-(int)(rbuf[1] + 0.5), colors[curcol]);
		ReleaseMutex(globMutex);
		break;
/*--- IFUNC=14, End picture -----------------------------------------------------*/
	case 14:
		WaitForSingleObject(globMutex, INFINITE);
		DeleteObject(SelectObject(memdc[WINDOW_ID], defpen));
		ReleaseMutex(globMutex);
		painting = 0;
		
	    *nbuf = -1;
		break;
/*--- IFUNC=15, Set color index -------------------------------------------------*/
	case 15: {
		HPEN pen;
		if( ((int)(rbuf[0] + 0.5) < 0)
		 || ((int)(rbuf[0] + 0.5) >= MAXCOLOR) ) break;
		WaitForSingleObject(globMutex, INFINITE);
		pen = CreatePen(PS_SOLID, 0, colors[(int)(rbuf[0] + 0.5)]);
		DeleteObject(SelectObject(memdc[WINDOW_ID], pen));
		curcol = (int)(rbuf[0] + 0.5);
		ReleaseMutex(globMutex);
	}
		break;
/*--- IFUNC=16, Flush buffer ----------------------------------------------------*/
	case 16: {
		HDC DC;

		WaitForSingleObject(globMutex, INFINITE);
		DC = GetDC(windows[WINDOW_ID]);
		BitBlt(DC, 0, 0, windowsize[WINDOW_ID].right, windowsize[WINDOW_ID].bottom, memdc[WINDOW_ID], 0, 0, SRCCOPY);
		ReleaseDC(windows[WINDOW_ID], DC);
		ReleaseMutex(globMutex);
	}
		break;
/*--- IFUNC=17, Read cursor -----------------------------------------------------*/
	case 17: {
		HPEN pen;

		windowdc = GetDC(windows[WINDOW_ID]);
		pen = SelectObject(windowdc, CreatePen(PS_SOLID, 0, colors[curcol]));
		cursorkey = -1;
		cursormode = (int)(rbuf[4] + 0.5);
		if (cursormode > 0) SetROP2(windowdc, R2_XORPEN);
		xmouseref = (int)(rbuf[2] + 0.5);
		ymouseref = windowsize[WINDOW_ID].bottom-(int)(rbuf[3] + 0.5);
		cursorx = (int)(rbuf[0] + 0.5);
		cursory = windowsize[WINDOW_ID].bottom-(int)(rbuf[1] + 0.5);
		wantcursor = 1;
		drawRubberBand(windowdc, xmouseref, ymouseref, cursorx, cursory, cursormode);
		Mouse = CreateEvent(NULL, TRUE, FALSE, NULL);
		WaitForSingleObject(Mouse, INFINITE);
		wantcursor = 0;
		deleteRubberBand(windowdc);
		CloseHandle(Mouse);
		if (cursormode > 0) SetROP2(windowdc, R2_COPYPEN);
		*lchr = 1;
		chr[0] = (char)cursorkey;
		rbuf[0] = (float)cursorx;
		rbuf[1] = (float)(windowsize[WINDOW_ID].bottom - cursory);
		*nbuf = 2;
		cursormode = 0;
		DeleteObject(SelectObject(windowdc, pen));
		ReleaseDC(windows[WINDOW_ID], windowdc);
	}
		break;
/*--- IFUNC=18, Erase alpha screen ----------------------------------------------*/
	case 18:
	    *nbuf = -1;
		break;
/*--- IFUNC=19, Set line style --------------------------------------------------*/
	case 19:
		*nbuf = -1;
		break;
/*--- IFUNC=20, Polygon fill ----------------------------------------------------*/
	case 20:
		WaitForSingleObject(globMutex, INFINITE);
		if (nVertices == 0) {
			nVertices = (int)(rbuf[0] + 0.5);
			npnt = 0;
			polypoints = (POINT *)malloc(nVertices * sizeof(POINT));
		} else {
			polypoints[npnt].x = (int)(rbuf[0] + 0.5);
			polypoints[npnt].y = windowsize[WINDOW_ID].bottom - (int)(rbuf[1] + 0.5);
			npnt++;
		}
		if (npnt == nVertices) {
			HBRUSH b;

			b = SelectObject(memdc[WINDOW_ID], CreateSolidBrush(colors[curcol]));
			Polygon(memdc[WINDOW_ID], polypoints, npnt);
			DeleteObject(SelectObject(memdc[WINDOW_ID], b));
			free(polypoints);
			nVertices = 0;
			npnt = 0;
		}
		ReleaseMutex(globMutex);
		*nbuf = -1;
		break;
/*--- IFUNC=21, Set color representation ----------------------------------------*/
	case 21: {
		int c;
		COLORREF oldcolor, newcolor;

		c = (int)(rbuf[0] + 0.5);
		if ((c < 0) || (c >= MAXCOLOR)) break;
		oldcolor = colors[c];
		newcolor = (int)(rbuf[1]*MAXCOLOR1 + 0.5) + 
					(int)(rbuf[2]*MAXCOLOR1 + 0.5) * 256 +
					(int)(rbuf[3]*MAXCOLOR1 + 0.5) * 256 * 256;
		colors[c] = newcolor;
	}
		break;
/*--- IFUNC=22, Set line width --------------------------------------------------*/
	case 22:
		*nbuf = -1;
		break;
/*--- IFUNC=23, Escape function -------------------------------------------------*/
	case 23:
	    *nbuf = -1;
		break;
/*--- IFUNC=24, Rectangle fill --------------------------------------------------*/
	case 24: {
		RECT r;
		HBRUSH b;

		WaitForSingleObject(globMutex, INFINITE);
		r.left = (int)(rbuf[0] + 0.5);
		r.right = (int)(rbuf[2] + 0.5) + 1;
		r.top = windowsize[WINDOW_ID].bottom - (int)(rbuf[1] + 0.5);
		r.bottom = windowsize[WINDOW_ID].bottom - (int)(rbuf[3] + 0.5) + 1;
		b = CreateSolidBrush(colors[curcol]);
		FillRect(memdc[WINDOW_ID], &r, b);
		DeleteObject(b);
		ReleaseMutex(globMutex);
	}
		break;
/*--- IFUNC=25, Set fill pattern ------------------------------------------------*/
	case 25:
		*nbuf = -1;
		break;
/*--- IFUNC=26, Line of pixels --------------------------------------------------*/
	case 26: {
		int x, y, i;
		unsigned char *pixels;

		WaitForSingleObject(globMutex, INFINITE);
		x = (int)(rbuf[0] + 0.5);
		y = windowsize[WINDOW_ID].bottom - (int)(rbuf[1] + 0.5);
		if ((pixels = (unsigned char *)malloc(*nbuf * sizeof(unsigned char))) == NULL) break;
		for (i = 2; i < *nbuf; i++) {
			pixels[i - 2] = (unsigned char)(rbuf[i] + 0.5);
			SetPixel(memdc[WINDOW_ID], x, y, colors[pixels[i - 2]]);
			x++;
		}
		free(pixels);
		ReleaseMutex(globMutex);
	}
		break;
/*--- IFUNC=27, Scaling information ---------------------------------------------*/
	case 27:
		rbuf[0] = 0.0;
		rbuf[1] = 1.0;
		rbuf[2] = 0.0;
		rbuf[3] = 1.0;
	    *nbuf = 4;
		break;
/*--- IFUNC=28, Draw marker -----------------------------------------------------*/
	case 28:
		*nbuf = -1;
		break;
/*--- IFUNC=29, Query color representation --------------------------------------*/
	case 29:
		if( ((int)(rbuf[0] + 0.5) < 0)
		 || ((int)(rbuf[0] + 0.5) >= MAXCOLOR) ) break;
		rbuf[1] = (colors[(int)(rbuf[0] + 0.5)] & 0x000000ff) / (float)MAXCOLOR1;
		rbuf[2] = ((colors[(int)(rbuf[0] + 0.5)] >> 8) & 0x000000ff) / (float)MAXCOLOR1;
		rbuf[3] = ((colors[(int)(rbuf[0] + 0.5)] >> 16) & 0x000000ff) / (float)MAXCOLOR1;
		*nbuf = 4;
		break;
/*--- IFUNC=30, Scroll rectangle ------------------------------------------------*/
	case 30:
		*nbuf = -1;
		break;
 	default:
		fprintf(stderr, "/WINDOWS: Unexpected opcode=%d in driver.\n", *opcode);
		*nbuf = -1;
		break;
	}
	return;
}
