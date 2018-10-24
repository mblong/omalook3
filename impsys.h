// Definitions for compile options
#define Mac
#define MacOS 1
#define Bundle 1	// set if the application is in an OSX bundle


#define DOTIFF 1

#define SWAP_BYTES 0

#define TWOBYTE short
#define SERIALOUT ".aout"
#define SERIALIN ".ain"
/*
// In versions before 2.1, DATAWORD was short
// like this:

#define DATAWORD short
#define DATABYTES 2
#define MAXDATAVALUE 32000		// the maximum value of a data word -- depends on the dataword size 


// if you want int variables, use this
#define DATAWORD int
#define DATABYTES 4
#define MAXDATAVALUE 0x7FFFFFFF
// for either case, these are the proper format specifiers

#define MAXMSG " Maximum %d at Row %d and Column %d.\n"
#define MINMSG " Minimum %d at Row %d and Column %d.\n"
#define DATAMSG " %d\n"
#define DATAFMT "%d"
*/

// starting with version 2.1, I'm liking float variabes. slower, but it gets rid of all that scale factor stuff

#define FLOAT 1
#define DATAWORD float
#define DATABYTES 4
#define MAXDATAVALUE 0x7FFFFFFF

#define MAXMSG " Maximum %g at Row %d and Column %d.\n"
#define MINMSG " Minimum %g at Row %d and Column %d.\n"
#define DATAMSG " %g\n"
#define DATAFMT "%g"



#define DOFFSET 80	/* -1 for old oma data; 80 for CCD data */
#define SFACTRMSG " Scale factor is %d.\n"

#define __ALLNU__	1

#define TBTHICK			20		/* title bar thickness */
#define NCOLOR 254				/* the number of colors assumed on MAC II*/
#define MAXNOCTRS 10			/* the max number of contours for contour plots */

#define SIZEOFEXTRA	50000		/* amount of extra memory needed before a window will
									open */
#ifndef DO_MACH_O
	#include <MacHeadersCarbon.h>
	#include <fcntl.h>
	#include <Strings.h>
	#include <StdLib.h>
	#include <StdIO.h>
	#include <CType.h>
	#include <unistd.h>
	#include <Carbon.h>
	#include <MacWindows.h>
	
	
	int pprintf(const char *,...);
#else
	#include <Carbon/Carbon.h>
	#include <fcntl.h>
	#include <unistd.h>
	//#define O_RDONLY		0x2			// open the file in read only mode 	//- mm 980420 -//
	#define printf omaprintf
	#define gets omagets
#endif

#define update_pause update_status
#define printmacronum update_status
#define drawtools update_status
#define printparams update_status
#define printxyzstuff printxyzstuff_nib
#define printcmin_cmax update_status

#define LONGBYTES 4				/* the length of long integers */

#define PICTHEADERSIZE	512		/* the number of bytes in a MacDraw PICT file header */


#define BITMAP		1		/* specifies that a window contains a bitmap */
#define LINEGRAPH	2		/* specifies that a window contains linegraphics */
#define ZOOM		3		/* specifies that a window containes a zoomed image */
#define FFT			4		/* specifies that a window contains the FFT of a subregion */
#define AUTO		5		/* Autocorrelation Window */
#define RGBMAP		6		/* RGB window */

#define QBITMAP		-1		/* specifies that a window contains a bitmap */
#define QLINEGRAPH	-2		/* specifies that a window contains linegraphics */
#define QZOOM		-3		/* specifies that a window containes a zoomed image */
#define QFFT		-4		/* specifies that a window contains the FFT of a subregion */
#define QAUTO		-5		/* Autocorrelation Window */
#define QRGBMAP		-6		/* RGB window */
#define QPIVGRAPH	-7		/* PIV window */
#define QROWCOL		2


int 	assigncolor1(),beep(),block(),checkpar(),calc(),
		drawboxes(),
		drawquad(),displa(),dochanplot(),dvect(float,float),
		drawtemplate(),file_error(),
		fileflush(),
		iset(float,float),maxx(),nomemory(),newpen(),openwindow(short ,short ,char* ,short),ploton(),
		pset(float,float),
		reset(),
		setincrementdialog(),setcheckbox(),setmaxbuttons(),setsizebuttons(),setctrvaltxt(),
		setbuttons(),setscreensize(),setviewdialog(),setvalues(),setarrow(),
		setccd(),
		scaly(float,float),scalx(float,float),settext(),set(float,float),
		setup_status_fonts(),savetemplate(),set_lr_point(),set_ul_point(),update_status(),
		vect(float,float);
		
int		getdatafile(),savedatafile(),getpictfile(),
		saveprefs(),loadprefs(),
		saverectangle(),getpalettefile(),activegwnum(),contour(),
		plotchans(),surface(),histogram(),halftone(),do_line_plot(),
		plotone(),pdvect(float,float),setpref(),setfonts();


		
int 	setup_status(),closewindow(),checkscroll(),zoom_image(),erasebox(),stopmacro(),
		find_minmax(),setcolors(),lowres(),printparams(),getwindowdata(),
		display(),try(),savefile(),setdata(),getlin(),integ(),loadmatch(),rotate(),
		getfile(),maknew(),fill_in_command(),checkevents(),nosuch(),dinit(),comdec(),helpdone(),
		clear_macro_to_end(),concatfile(),illval(),arithfile(),setcolordialog();

int 	loadtemplate();

void 	setcontour(),setsurface(),setplotone(),sethist(),template(),
		setprefix();

DATAWORD idat(int,int);



int		keylimit(),getpairs(float,float,float),data_histogram(),setbounds(DATAWORD*,DATAWORD, DATAWORD,int),
		get_vect_dir(float,float,float *,float *),
		piv_reorder(),getpivdata(),swap_peaks(),peaks(TWOBYTE nx,TWOBYTE ny,int m,int rad,DATAWORD* dpt),
		get_points(),checkstoprequest(),
		printvel(float,float,float,float),pivplot(int,float,float,float,float,int),fft_image(),change_vector(Point, short, short);


