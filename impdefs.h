/*
OMAX -- Photometric Image Processing and Display
Copyright (C) 2006  by the Developers of OMA

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

//#include "impsys.h"	/* system dependent parameters */ 
#define TWOBYTE short
#define FLOAT 1
#define DATAWORD float
#define DATABYTES 4

#ifdef DOTIFF
#include "tiffio.h"
#endif

#define PI 3.14159265358979323846
#define Pi 3.14159265358979323846

#define NGWINDS		30		/* the maximum number of graphics windows */

/* Things for loops in macros */
#define NESTDEPTH 20		// Should add checking for overflow; just make big for now
// depth of execute commands
#define EX_NEST_DEPTH 40	// Should add checking for overflow; just make big for now


#define NUMWATCHES	8		/* the number of watch cursors */
#define	MYWATCHNUM	257-1	/* 257 is the id # of my first watch cursor */


#define EOL     0        	/* this marks the end of the command list */
#define CHPERLN 128      	/* maximum number of characters per line */
#define HEADLEN 30       	/* number of bytes in header */
#define TRAILEN 62       	/* number of bytes in trailer */
#define COMLEN 512-HEADLEN-TRAILEN  /* number of bytes in comment buffer */
#define MBUFLEN 10240     	/* number of bytes in macro buffer */
#define VBUFLEN	1024		/* the storage for variable names in macros */
#define MATHLEN 10240		/* number of elements in the file math buffer */
#define DBUFLEN 490000		/* number of elements in the data buffer */
#define MAXCHAN	1248		/* maximum number of channels */
#define MAXTRAK	1248		/* maximum number of tracks */
#define LF '\012'
#define CarRet '\015'
#define BEL '\007'
#define SIT	0				/* the detector is a SIT */
#define CCD 1				/* the detector is a CCD */

#define GRAYMAP 1			/* specify a gray color map */
#define BGRBOW	0			/* a blue to green to red rainbow color map */
#define REDMAP	2			/* red shades color map */
#define GREENMAP 3			/* green shades color map */
#define BLUEMAP  4			/* blue shades color map */
#define FROMAFILE 5                 /* palette from a file */
#define DEFAULTMAP 5
#define FROMAFILE2 6
#define FROMAFILE3 7

#define NUMPAL 8

#define MAXDOFFSET 80		/* the maximum databuffer offset */

#define NUM_TEMP_IMAGES 10	// number of temporary images
			   
/* Define the indices to important locations in the header */

#define NMAX    1
#define LMAX    2
#define NMIN    3
#define LMIN    4
#define NFRAM   5
#define NPTS    6
#define NCHAN   7
#define NTRAK   8
#define NTIME   9
#define NPREP   10
#define NX0     11
#define NY0     12
#define NDX     13
#define NDY     14

/* Define the indices to important locations in the trailer */

#define FREESP	0
#define IDWRDS	1			// use this to indicate byte ordering
#define RULER_CODE 2		/* if trailer[RULER_CODE] == MAGIC_NUMBER, assume a ruler */
#define MAGIC_NUMBER 12345  /*   has been defined. */
#define RULER_SCALE 3		/* A floating point number occupying trailer[3] & [4]. Pixels/Unit. */
#define RULER_UNITS 5		/* The starting location of a text string specifying the 
								name of the units. Occupies trailer[5] to trailer[12] */
#define RUNNUM	13
#define	TOMA2	14
#define	ROMA1	15
#define SFACTR	17
#define NDATE	18
#define DMODE	21
#define NDATW	22
#define SAMT	23
#define SUBFC	24
#define NREAD	25
#define LSYFG	26
#define COLFG	27
#define NDISF	28
#define NDELY	29
#define ACSTAT	30

// use these constants in both bytes of trailer[IDWRDS] to specify the byte ordering of files
// big endian is PowerPC et al
// little endian is intel et al

#define LITTLE_ENDIAN_CODE 127
#define BIG_ENDIAN_CODE 0



#define INIT		0
#define RUN			1
#define INFO		2
#define SEND		3
#define TRANS		4
#define BYE			5
#define FLUSH		6
/*#define TAKE		7 */
#define RECEIVE 	8
#define ASK			9
#define FORCE_INIT	10

#define PHOTOMETRICS_CC200		0
#define PRINCETON_INSTRUMENTS_1	1
#define STAR_1					2
#define SCOPE					3

/* Math Routine Definitions */

#define ADD	1
#define SUB	2
#define MUL	3
#define DIV	4

/* Prefix and Suffix Types */

#define SAVE_DATA 	1
#define GET_DATA	2
#define MACROS_DATA		3
//#define GRAPHICS_DATA	4
#define SETTINGS_DATA	4

#define TIFF_DATA		8
#define TIF_DATA		9
#define CSV_DATA        10
#define FTS_DATA        11
#define RAW_DATA        12
#define PDF_DATA        13

#ifdef DO_MACH_O
	#ifndef SETTINGSFILE
		#define SETTINGSFILE "oma.app/Contents/Resources/OMA Settings"
		#define PALETTEFILE	"oma.app/Contents/Resources/OMA palette.pa1"
		#define PALETTEFILE2 "oma.app/Contents/Resources/OMA palette2.pa1"
		#define PALETTEFILE3 "oma.app/Contents/Resources/OMA palette3.pa1"
	#endif
#else
    #define SETTINGSFILE "OMA Settings"
    #define PALETTEFILE	"OMA Palette"
    #define PALETTEFILE2 "OMA Palette2"
    #define PALETTEFILE3 "OMA Palette3"
#endif

#define NUMBER_LINEOUTS	30

#define OMA_OK 0
#define OMA_FILE -1

#define LOLAB_LENGTH 15	/* Line Out Labels can have 15 characters */
typedef struct {
	char name[LOLAB_LENGTH+1];
} LOLabel;				

typedef struct {
	char name[16];
} Cname;

typedef struct {
	Cname text;
	int (*fnc)();
} ComDef;

#define MAX_VAR_LENGTH 32
#define MAX_VAR 100

typedef struct {
	char vname[MAX_VAR_LENGTH];
	int ivalue;
	float fvalue;
	int is_float;
	char estring[128];
} Variable;

typedef struct {
	char op_char;
	int ivalue;
	float fvalue;
	char estring[128];
} Expression_Element;

typedef struct {
	int window;
	int where;
	char lineongraph;
} RCMarker;


/* Define the Menu numbers */
#define	appleMenu		0				
#define	fileMenu		1
#define	editMenu		2
#define graphMenu		3
#define lineoutMenu		4
#define starMenu		5
#define nu200Menu		6
#define st6Menu			7
#define SensiCamMenu	8

#define	rMenuBar	128	/* menu bar */

#define	mApple		128	/* Apple menu */
#define	iAbout		1

#define	mFile		129	/* File menu */
#define	iNew		1
#define	iClose		4
#define	iQuitSeparator	10
#define	iQuit		11

#define	mEdit		130	/* Edit menu */



#define UNIT_NAME_LENGTH 16		/* length of the unit name for the ruler - 
									don't make this bigger without considering what it does
									to the trailer -- this is saved with each file */


#define NUM_CAMERA_FORMATS  8	/* the max number of headers that can be stored in memory */

#define HISTOSIZE 256			/* the size of the histogram data array used for data (not
								   display) */
								   
#define NPEAK 32

#define PIV	1
