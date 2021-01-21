/*
 *  oma_bits.h
 *  omalook3
 *
 *  Created by Marshall Long on 11/3/07.
 *  Copyright 20018. All rights reserved.
 *
 */


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

#define TWOBYTE short
#define FLOAT 1
#define DATAWORD float
#define DATABYTES 4


#define PALETTEFILE	"/Library/QuickLook/omalook3.qlgenerator/Contents/Resources/OMA palette.pa1"
#define THUMBSIZE 256
#define PREVIEWSIZE 512
#define IMAGEDATA 0
#define PIVDATA 1
#define UNKNOWNTEXT 2
#define UNKNOWN 3
#define OMAPREFS 4
#define HOBJ 5
#define ERROR -1

#define EOL     0        	/* this marks the end of the command list */
#define CHPERLN 128      	/* maximum number of characters per line */
#define HEADLEN 30       	/* number of bytes in header */
#define TRAILEN 62       	/* number of bytes in trailer */
#define COMLEN 512-HEADLEN-TRAILEN  /* number of bytes in comment buffer */

#define MAXWIDE 8000        // largest width image this will try to decode

#define MAXDOFFSET 80		/* the maximum databuffer offset */


			   
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
#define	IS_COLOR	15
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

#define OMA_OK 0
#define OMA_FILE -1

typedef struct {
	char name[16];
} Cname;

typedef struct {
	Cname text;
	int (*fnc)();
} ComDef;





#define UNIT_NAME_LENGTH 16		/* length of the unit name for the ruler - 
									don't make this bigger without considering what it does
									to the trailer -- this is saved with each file */


// oma2 data format

/******************** Constants for Classes ********************/
#define NSPECS  32   // number of integers in an Image specification array
#define NVALUES 16   // number of values associated with an Image (things like min, max, etc.)
#define NRULERCHAR 16   // number of characters in the units of the ruler
#define NUM_IMAGE_PTRS 3    // the number of pointers in the Image Class
// locations within the specs array
enum {ROWS,COLS,X0,Y0,DX,DY,LMAX_,LMIN_,IS_COLOR_,HAVE_MAX,HAS_RULER,
    LRMAX,LRMIN,LGMAX,LGMIN,LBMAX,LBMIN,NFRAMES,SAVE_FORMAT};

// types if integers the data can be saved as
enum {UNSIGNED16=59464,SIGNED16,UNSIGNED8,SIGNED8};

// SAVE_FORMAT is a 2021 addition to the specs array
// This will be used to save data in other than DATAWORD types
// Will assume that the chances of specs[SAVE_FORMAT] accidentally being one of the magic values is small

// locations within the values array
enum {MIN,MAX,RMAX,RMIN,GMAX,GMIN,BMAX,BMIN,RULER_SCALE_};

// Image error codes and command return codes
enum {NO_ERR,SIZE_ERR,FILE_ERR,MEM_ERR,ARG_ERR,CMND_ERR,GET_MACRO_LINE};


#define OMA2_BINARY_DATA_STRING  "OMA2 Binary Data 1.0"

typedef struct {
    int         specs[NSPECS];      // information on Image size, type, etc.
    DATAWORD    values[NVALUES];    // important values (things like min, max, etc.)
    char        unit_text[NRULERCHAR];
    int         error;
    int         is_big_endian;
    int         commentSize;
    int         extraSize;
    char*       comment;
    float*      extra;
    DATAWORD*   data;
    
} oma2data;

OSStatus readHobj( CFURLRef,TWOBYTE **ptrh,TWOBYTE **ptrt, DATAWORD **ptrd);


