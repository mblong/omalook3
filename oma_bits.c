/*
 *  oma_bits.c
 *  omalook3
 *
 *  Created by Marshall Long on 11/3/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


#include <CoreFoundation/CoreFoundation.h>
#include <CoreServices/CoreServices.h>
#include <QuickLook/QuickLook.h>
#include "oma_bits.h"
#include <math.h>

/* ********** */

int get_byte_swap_value(short id,int is_big_endian)
{
	char* cp;
	cp = (char*) &id;
	if( is_big_endian ) {	// running on a PowerPC
		if( *(cp) ==  LITTLE_ENDIAN_CODE && *(cp+1) ==  LITTLE_ENDIAN_CODE)
			return 1;	// must have been saved on an intel machine so have to swap bytes
		else
			return 0;	// must be same kind, leave it alone
	} else {			// running on intel
		if( *(cp) ==  LITTLE_ENDIAN_CODE && *(cp+1) ==  LITTLE_ENDIAN_CODE)
			return 0;	// must have been saved on an intel machine so leave it as is
		else
			return 1;	// must be from a powerPC, have to change it
	}

}

 int two_to_four(DATAWORD* dpt, int num, TWOBYTE scale)
 {
		TWOBYTE* two_byte_ptr;
		int i;
		
		two_byte_ptr = (TWOBYTE*) dpt;
		
		if(scale <= 0) scale = 1;
		for(i=num-1; i>=0; i--){			
			*(dpt+i) = *(two_byte_ptr+i);
			*(dpt+i) *= scale;
		}
		return 0;
}

 void swap_bytes_routine(char* co, int num,int nb)
 {
	int nr;
	char ch;
	if(nb == 2){
		for(nr=0; nr < num; nr += nb) {
			ch = co[nr+1];
			co[nr+1] = co[nr];
			co[nr] = ch;
		}
	} else if (nb == 4){
		for(nr=0; nr < num; nr += nb) {
			ch = co[nr+3];
			co[nr+3] = co[nr];
			co[nr] = ch;
			ch = co[nr+2];
			co[nr+2] = co[nr+1];
			co[nr+1] = ch;
		}
	} else if (nb == 8){
		for(nr=0; nr < num; nr += nb) {
			ch = co[nr+7];
			co[nr+7] = co[nr];
			co[nr] = ch;
			
			ch = co[nr+6];
			co[nr+6] = co[nr+1];
			co[nr+1] = ch;
			
			ch = co[nr+5];
			co[nr+5] = co[nr+2];
			co[nr+2] = ch;
			
			ch = co[nr+4];
			co[nr+4] = co[nr+3];
			co[nr+3] = ch;
		}
	}
	
}

int get_oma_data_type( CFURLRef url)
{
	// turn the url into a file path string
	
	char txt[256];
	FILE *fp;
	int k,n,m;
		
	// turn path into fsref
	FSRef theFSRef;
	FSCatalogInfo myInfo;
	FInfo finfo,*finfo_ptr;

	CFURLGetFSRef(url,&theFSRef);
	FSRefMakePath(&theFSRef,(unsigned char*)txt,255); // this is the path
	FSGetCatalogInfo (&theFSRef,kFSCatInfoFinderInfo,&myInfo,NULL,NULL,NULL);
	finfo_ptr = (FInfo*) &myInfo.finderInfo;
	finfo = *finfo_ptr;
	if( finfo.fdType == 'OPRF' ) return OMAPREFS;	// don't show prefs -- there's no data there
	if( finfo.fdType == 'TEXT' ) {	// this could be PIV data 
		fp = fopen(txt,"r");
		if( fp != NULL) {	
			fscanf( fp,"%d\t%d\t%d\t",&k,&m,&n);
			// read the text label to be sure this is a valid piv file 
			// don't try and compare the \n, since that may be different between os9 and osx
			fread(txt,strlen("x locations, y locations, peaks/location"),1,fp);
			if( memcmp(txt,"x locations, y locations, peaks/location",
						strlen("x locations, y locations, peaks/location")) != 0) {
				// not a valid PIV file
				fclose(fp);
				return(UNKNOWNTEXT);
			} else {
				fclose(fp);
				return(PIVDATA);
			}
		} else {
			return ERROR;
		}
	} 

	if( finfo.fdType == 'ODTA' ) return IMAGEDATA;
	if( strncmp(&txt[strlen(txt)-4],".o2d",4) == 0 ||
        strncmp(&txt[strlen(txt)-4],".O2D",4) == 0)return IMAGEDATA;
    if( strncmp(&txt[strlen(txt)-5],".hobj",5) == 0 ||
       strncmp(&txt[strlen(txt)-5],".HOBJ",5) == 0)return HOBJ;

	// other checks ??
	
	// Going forward, need to have oma set a UTI according to data type so that appropriate previews can be generated
	
	return UNKNOWN;
}

OSStatus get_oma_data( CFURLRef url, TWOBYTE **head, TWOBYTE **trail, DATAWORD **data)
{
	DATAWORD *datpt;
	TWOBYTE	*header,*trailer;   //trailer[TRAILEN/2];
	char	comment[COMLEN] = {0};
	int nbyte,nr;
	int fd;
	int is_big_endian,swap_bytes;
	
	UInt8 txt[512];
    
    int doffset = 80;
	
	//FSRef theFSRef;
	//CFURLGetFSRef(url,&theFSRef);
	//FSRefMakePath(&theFSRef,(unsigned char*)txt,512); // this is the path
    CFURLGetFileSystemRepresentation(url, true, txt, 512);
	
	// get the byte ordering on this machine
	
	CFByteOrder the_order = CFByteOrderGetCurrent();
	if(the_order == CFByteOrderBigEndian){
		// powerPC ordering
		swap_bytes = 0;
		is_big_endian = 1;
	} else {
		// intel ordering
		swap_bytes = 1;
		is_big_endian = 0;
	}
    
	// read in information from the header
	
	if((fd = open((char*)txt,O_RDONLY)) == -1) {
		return -1;
	} else {
		header = (TWOBYTE*) malloc(HEADLEN);
		if(header == NULL)
			return -1;
		trailer = (TWOBYTE*) malloc(TRAILEN);
		if(trailer == NULL)
			return -1;
		read(fd,(char*)header,HEADLEN);
        
        if (strncmp((char*)header, OMA2_BINARY_DATA_STRING,HEADLEN) == 0) { // new save format
            int nspecs,nvalues,nrulerchar;
            oma2data image;
            read(fd,&nspecs,sizeof(int));
            read(fd,&nvalues,sizeof(int));
            read(fd,&nrulerchar,sizeof(int));
            // we'll assume that the number of specs, values, and ruler characters won't decrease in future versions, but allow for them to increase
            read(fd,image.specs,sizeof(int)*nspecs);
            read(fd,image.values,sizeof(DATAWORD)*nvalues);
            read(fd,image.unit_text,nrulerchar);
            read(fd,&image.error,sizeof(int));
            read(fd,&image.is_big_endian,sizeof(int));
            read(fd,&image.commentSize,sizeof(int));
            read(fd,&image.extraSize,sizeof(int));
            if(image.commentSize >0){
                image.comment = malloc(image.commentSize);
                read(fd,image.comment,image.commentSize);
                if(image.comment) free(image.comment);
            }
            if(image.extraSize > 0){
                image.extra = (float*)malloc(image.extraSize*sizeof(float));
                read(fd,image.extra,image.extraSize*sizeof(float));
                if(image.extra) free(image.extra);
            }
            
            header[NCHAN] = image.specs[COLS];
            header[NTRAK] = image.specs[ROWS];
            
            header[LMAX] = image.specs[LMAX_]/image.specs[COLS];
            header[NMAX] = image.specs[LMAX_]%image.specs[COLS];
            header[LMIN] = image.specs[LMIN_]/image.specs[COLS];
            header[NMIN] = image.specs[LMIN_]%image.specs[COLS];
            
            doffset = 0;
            
            trailer[SFACTR] = 1;
            
            trailer[IS_COLOR] = image.specs[IS_COLOR_];
            
            swap_bytes = 0;
            
        } else {
            
            read(fd,comment,COMLEN);
            read(fd,(char*)trailer,TRAILEN);
            
            
            
            swap_bytes = get_byte_swap_value(trailer[IDWRDS],is_big_endian);
            if(swap_bytes) {
                swap_bytes_routine((char*)header,HEADLEN,2);
                swap_bytes_routine((char*)trailer,TRAILEN,2);
            }
		}
		// some checks here for valid OMA data
		// see if the image seems just too big
		if(header[NCHAN] > MAXWIDE || header[NTRAK] > MAXWIDE*3 ){
			free(header);
            free(trailer);
			return -1;
		}
		
		nbyte = (header[NCHAN]*header[NTRAK]) * DATABYTES;
        
        if (nbyte < doffset* DATABYTES) {   // in case of a very small image
            nbyte = doffset* DATABYTES;
        }
		
		datpt = (DATAWORD*) malloc(nbyte);
		if(datpt == NULL){
			free(header);
            free(trailer);
			return -1;
		}
        
        int dataSize = DATABYTES;
        FILE* f=fopen((char*)txt ,"r");
        fseek(f, 0, SEEK_END);
        long size = ftell(f)-HEADLEN+COMLEN+TRAILEN;
        fclose(f);
        if (roundf((float)nbyte/size) == 2.) {
            dataSize = 2;
        }

        
        if(doffset > 0) read(fd,(char*)datpt,doffset* dataSize);   // read and then forget any offset (from old oma files)
        
		nr = read(fd,(char*)datpt, nbyte);
		if(nr == 0) {
			free(header);
            free(trailer);
			free(datpt);
			return -1;
		}
		
		if(dataSize == 2) {
			// this is a 2-byte data file
			// adjust to 4-byte format
			// printf("2-byte input file\n");
			if(swap_bytes)  swap_bytes_routine((char*) datpt, nr,2);
			two_to_four(datpt,nr/2,trailer[SFACTR]);
			trailer[SFACTR]=1;
		} else {
			if(swap_bytes) swap_bytes_routine((char*) datpt, nr, DATABYTES);
		}
        
		close(fd);
        
	}
	*head = header;
	*data = datpt;
    *trail = trailer;
	return noErr;
}

int getpalettefile(char* name, RGBColor *filecolor) 
{

	unsigned short i;
	int fd;
	unsigned char colors[256];

	fd = open(PALETTEFILE,O_RDONLY);
	if(fd == -1) {
		return -1;
	}
	read(fd,colors,256);
	for(i=0; i<256; i++)
		filecolor[i].red = colors[i]<<8;
	read(fd,colors,256);
	for(i=0; i<256; i++)
		filecolor[i].green = colors[i]<<8;
	read(fd,colors,256);
	for(i=0; i<256; i++)
		filecolor[i].blue = colors[i]<<8;

	return 0;
	
}


//-----------------------------------------------
Ptr Get_rgb_from_image_buffer( TWOBYTE *header, TWOBYTE *trailer, DATAWORD *datpt, int pixsiz)
{
	Ptr ptr;
	DATAWORD* point = datpt;
	long k = 0, i,j,n=0;
	int ntrack = header[NTRAK];
	int nchan = header[NCHAN];
	int nth;
	float pix_scale;

	float fpindx;
	DATAWORD crange,cm,ncm1,indx,cmin=0,cmax=0;
	int pindx;
	
	RGBColor color[256];
	DATAWORD *mydatpt;
    
    
	
	if( getpalettefile(PALETTEFILE, color) == -1) {
		return 0;
	}
	cmax = *(datpt+header[NMAX]+header[LMAX]*header[NCHAN]);
	cmin = *(datpt+header[NMIN]+header[LMIN]*header[NCHAN]);
	if(cmin == cmax) {
		mydatpt = datpt;
		cmin = cmax = *mydatpt;
        
		while ( mydatpt < datpt+header[NCHAN]*header[NTRAK]) {
			if ( *mydatpt < cmin )cmin = *mydatpt;
			if ( *mydatpt > cmax )cmax = *mydatpt;
			mydatpt++;
		}
	}
	crange = cmax - cmin;
	ncm1 = 253;
	cm = cmin;
	
	if( pixsiz > 0 ){
		nth = 1;
		pix_scale = 1.0;
	} else {
		nth = abs(pixsiz);
		pix_scale=1.0/nth;
	}
	
	ptr = calloc(header[NCHAN]/nth*header[NTRAK]/nth,4);
	
	if(ptr == NULL){
		return ptr;
	}
	
	if( pixsiz > 0 ) {
		for(i=0; i < ntrack; i++){
			for(j=0; j < nchan; j++){
				indx = *(point+k++) - cm;
				fpindx = (float)indx * ncm1;
				pindx = fpindx/crange;
				if( pindx > ncm1)
					pindx = ncm1;
				if( pindx < 0)
					pindx = 0;
				++pindx;
				*(ptr+n++) =pindx;
				*(ptr+n++) = color[pindx].red/256;
				*(ptr+n++) = color[pindx].green/256;
				*(ptr+n++) = color[pindx].blue/256;
			}
		}
	}else {
		i = 0;
		while(++i < ntrack/nth){
			j = 0;
			while( j++ < nchan/nth){
				indx = *(point+k) - cm;
				k += nth;
				fpindx = (float)indx * ncm1;
				pindx = fpindx/crange;
				if( pindx > ncm1)
					pindx = ncm1;
				if( pindx < 0)
					pindx = 0;
				++pindx;
				//*(ptr+n++) = 0xFF;
				*(ptr+n++) =pindx;
				*(ptr+n++) = color[pindx].red/256;
				*(ptr+n++) = color[pindx].green/256;
				*(ptr+n++) = color[pindx].blue/256;
			}
		k = i * nth * nchan;
		}
	}
	
	return ptr;

}
int scale_pixval(DATAWORD val, DATAWORD cmin, DATAWORD crange)
{
    int pval;
    float fpval;
    
    fpval = (val-cmin) * 255.0;
    pval = fpval/crange;
    if( pval > 255)
        pval = 255;
    if( pval < 0)
        pval = 0;
    return pval;
}


Ptr Get_color_rgb_from_image_buffer( TWOBYTE *header, TWOBYTE *trailer, DATAWORD *datpt, int pixsiz)
{
	Ptr ptr;
	DATAWORD* point = datpt;
	long k = 0, i,j,n=0;
	int ntrack = header[NTRAK];
	int nchan = header[NCHAN];
	int nth,intensity;
	float pix_scale;
	DATAWORD *pt_green,*pt_blue;
    
	//float fpindx;
	DATAWORD crange,cmax,cmin;
    //int pindx;
    
	cmax = *(datpt+header[NMAX]+header[LMAX]*header[NCHAN]);
	cmin = *(datpt+header[NMIN]+header[LMIN]*header[NCHAN]);

    
	crange = cmax - cmin;
	
	if( pixsiz > 0 ){
		nth = 1;
		pix_scale = 1.0;
	} else {
		nth = abs(pixsiz);
		pix_scale=1.0/nth;
	}
	
    ptr = calloc(header[NCHAN]/nth*header[NTRAK]/nth/3,4);
	if(ptr == NULL){
        return ptr;
	}
	pt_green = point + nchan*ntrack/3;
	pt_blue =  pt_green + nchan*ntrack/3;
    
	if( pixsiz > 0 ) {
		for(i=0; i < ntrack/3; i++){
			for(j=0; j < nchan; j++){
				*(ptr+n+1) = scale_pixval(*(point+k), cmin, crange);
				*(ptr+n+2) = scale_pixval(*(pt_green+k), cmin, crange);
				*(ptr+n+3) = scale_pixval(*(pt_blue+k++), cmin, crange);
				intensity = ( (unsigned char) *(ptr+n+1) + (unsigned char) *(ptr+n+2) + (unsigned char) *(ptr+n+3))/3;
				*(ptr+n) = intensity;
				n += 4;
			}
		}
	}else {
		i = 0;
		while(++i < ntrack/nth/3){
			j = 0;
			while( j++ < nchan/nth){
				*(ptr+n+1) = scale_pixval(*(point+k), cmin, crange);
				*(ptr+n+2) = scale_pixval(*(pt_green+k), cmin, crange);
				*(ptr+n+3) = scale_pixval(*(pt_blue+k), cmin, crange);
				intensity = ( (unsigned char) *(ptr+n+1) + (unsigned char) *(ptr+n+2) + (unsigned char) *(ptr+n+3))/3;
				*(ptr+n) = intensity;
				k += nth;
				n += 4;
			}
			k = i * nth * nchan;
		}
	}
	return ptr;
}


//__________________________________________________________________________________

// -----------------------------------------------------------------------------
//	GWorldImageBufferRelease
// -----------------------------------------------------------------------------
//
void
GWorldImageBufferRelease(
	void*						inInfo,
	const void*					inData,
	size_t						inSize )
{
	#pragma unused( inData, inSize )
	//printf("image release\n");
	if(inData!=0) DisposePtr( (Ptr) inData );
}

// -----------------------------------------------------------------------------

int plotpivdata( CGContextRef cgContext, CGRect dstRect, float pen_width, int no_of_velx, int no_of_vely, int fftsize, int boxinc, float* xpeaks, float* ypeaks)
{
	float x0=0,y0=0,x1=0,y1=0,scale=1,dx,dy,len;
	float		arrowx = 5;
	float		arrowy = 2;	
	float sfact;

	int ulx,uly,nc,nt;
	int num=0;

	// Initialize
	CGContextSetRGBFillColor (cgContext, 1, 1, 1, 1);
	CGContextSetRGBStrokeColor (cgContext, 0, 0, 0, 1);
	CGContextSetLineWidth(cgContext, pen_width);
	CGContextSetLineJoin(cgContext,kCGLineJoinRound);
	CGContextSetLineCap(cgContext,kCGLineCapRound);
	CGContextFillRect (cgContext,dstRect);
	CGContextStrokeRect(cgContext,dstRect);
	CGContextBeginPath(cgContext);
	
	nc = (no_of_velx-1)*boxinc+fftsize;
	nt = (no_of_vely-1)*boxinc+fftsize;
	
	x1 = nc;
	y1 = nt;
	// local coordinates run from 0 to x1; 0 to y1
	sfact = dstRect.size.width/x1;
	
	for(uly = 0; uly <= nt - fftsize; uly += boxinc){
		for(ulx = 0; ulx <= nc - fftsize; ulx += boxinc){
			
			x0 = ulx + fftsize/2;		// the center of this box (used in plotting) 
			y0 = uly + fftsize/2;
				
			// plot velocities 
			
			x1 = scale*(xpeaks[num] - fftsize/2);
			y1 = scale*(ypeaks[num] - fftsize/2);
			
			x1 = ulx + xpeaks[num] + x1;
			y1 = uly + ypeaks[num] + y1;
			num += 1;

			CGContextMoveToPoint(cgContext,x0*sfact, dstRect.size.height - y0*sfact);
			CGContextAddLineToPoint(cgContext,x1*sfact, dstRect.size.height - y1*sfact);
			
			dx = x1 - x0;
			dy = y1 - y0;
			len = sqrt(dx*dx +dy*dy);
		
			if(len != 0.0) {
				x0 = x1 - dx/len*arrowx + dy/len*arrowy;
				y0 = y1 - dx/len*arrowy - dy/len*arrowx;
				CGContextMoveToPoint(cgContext,x0*sfact, dstRect.size.height - y0*sfact);
				CGContextAddLineToPoint(cgContext,x1*sfact, dstRect.size.height - y1*sfact);

				x0 = x1 - dx/len*arrowx - dy/len*arrowy;
				y0 = y1 + dx/len*arrowy - dy/len*arrowx;
				CGContextMoveToPoint(cgContext,x0*sfact, dstRect.size.height - y0*sfact);
				CGContextAddLineToPoint(cgContext,x1*sfact, dstRect.size.height - y1*sfact);
			}
		}
	}
	

	//CGContextMoveToPoint(cgContext,100.,0.);
	//CGContextAddLineToPoint(cgContext,0.,100.);
	CGContextStrokePath(cgContext );

	return 0;
}

/*
int pivplt(int n,int index)
{
	extern char scratch[];
    extern char cmnd[];
	
	float x0=0,y0=0,x1=0,y1=0,scale=1;

	int ulx,uly,nc,nt;
	int num=1;
	
	if( sscanf(&cmnd[index],"%f %s",&scale,scratch) < 1)
		scale = n;
		
	if(n == 0)
		scale = pivscale;

	if(scale == 0.0) scale = 1.0;
	pivscale = scale;
	
	scale -= 1.0;	// the scale factor input is a multiplicitive factor 

	nc = (no_of_velx-1)*boxinc+fftsize;
	nt = (no_of_vely-1)*boxinc+fftsize;
	
	x1 = nc;
	y1 = nt;
	
	if( !pivplot(INIT,1.0,0.0,x1,y1,0) ) {
		beep();
		return -2;
	}

	// get the next subarray of the current image 
		
	for(uly = 0; uly <= nt - fftsize; uly += boxinc){
		for(ulx = 0; ulx <= nc - fftsize; ulx += boxinc){
			
			x0 = ulx + fftsize/2;		// the center of this box (used in plotting) 
			y0 = uly + fftsize/2;
				
			// plot velocities 
			
			x1 = scale*(xpeaks[num] - fftsize/2);
			y1 = scale*(ypeaks[num] - fftsize/2);
			
			x1 = ulx + xpeaks[num] + x1;
			y1 = uly + ypeaks[num] + y1;
			
			num += peaks_to_save;
			
			get_vect_dir(x0,y0,&x1,&y1);
			
			pivplot(PLOT,x0,y0,x1,y1,0);
			
		}
	}
	pivplot(CLOSE,x0,y0,x1,y1,0); 
	
	return 0;
		
}
*/
int getpivdata(CFURLRef url, int* no_of_velx, int* no_of_vely, int* fftsize, int* boxinc, float **xpeaks, float **ypeaks)		
{
	int peaks_to_save;		// the number of peaks to saved 
	float  peakval;
	float* xpeak;
	float* ypeak;	
	
	CFStringRef urlpath;
	char txt[256];
	FILE *fp;
	int n,m;
	int	i = 0;
		
	urlpath = CFURLCopyFileSystemPath(url,kCFURLPOSIXPathStyle);
	CFStringGetCString (urlpath,txt,256,kCFStringEncodingMacRoman);

	fp = fopen(txt,"r");
	if( fp != NULL) {	
		fscanf( fp,"%d\t%d\t%d\t",no_of_velx,no_of_vely,&peaks_to_save);
		// read the text label to be sure this is a valid piv file 
		// don't try and compare the \n, since that may be different between os9 and osx
		fread(txt,strlen("x locations, y locations, peaks/location"),1,fp);

		if( memcmp(txt,"x locations, y locations, peaks/location",
					strlen("x locations, y locations, peaks/location")) != 0) {
			//printf("Not a valid PIV file.\n");
			fclose(fp);
			return(-1);
		}
		fscanf( fp,"%d\t%d",fftsize,boxinc);
		fread(txt,strlen("\tfft size, increment\n"),1,fp);			// skip past more text 
		fread(txt,strlen("x coord\ty coord\tvalue\n"),1,fp);		// skip past more text 
		
		// allocate memory 
		xpeak = (float *)malloc(*no_of_velx * *no_of_vely *  4 );	// just get first peak
		ypeak = (float *)malloc(*no_of_velx * *no_of_vely *  4 );

		if( xpeak == 0 || ypeak == 0 ) {
			//nomemory();
			if(xpeak != 0) free(xpeak);
			if(ypeak != 0) free(ypeak);
			return -2;
		}
		m = 0;
        for( i=0; i< *no_of_velx * *no_of_vely; i++) {
			for(n=0; n<peaks_to_save; n++) {
				if( n == 1 ) {
					fscanf( fp,"%f\t%f\t%g",&(xpeak[m]),&(ypeak[m]),&peakval);
					m++;
				} else {
					// for a thumbnail or preview, we don't care about any extra peaks (or the peak value)
					// just read and discard
					fscanf( fp,"%f\t%f\t%g",&peakval,&peakval,&peakval);
				}
			}
        }
    	fclose(fp);
		*xpeaks = xpeak;
		*ypeaks = ypeak;
		return 1;
	}
	else {
		return -3;
	}
}

unsigned long fsize(char* file)
{
    FILE * f = fopen(file, "r");
    if (f == NULL) {
        return 0;
    }
    fseek(f, 0, SEEK_END);
    unsigned long len = (unsigned long)ftell(f);
    fclose(f);
    return len;
}

DATAWORD idat(DATAWORD *datpt, int nt, int nc, int rows, int cols)
{
    DATAWORD *pt;
    int index;
    
    if (datpt == 0) return(0);
    if(nc < 0) nc = 0;
    if(nt < 0) nt = 0;
    if(nc > cols-1) nc = cols-1;
    if(nt > rows-1) nt = rows-1;
    index = nc + nt*cols;
    //if (index >= header[NCHAN]*header[NTRAK]) return(0);    // check for illegal value passed
    pt = datpt + index;
    return(*pt);
}


void decodeHobj(DATAWORD* theImage,int width, int height){
    int nt,nc,tweak;
    unsigned int x0,x1;
    
    for(nt=0; nt<height; nt++) {
        for(nc=0; nc<width;nc++){
            // check for crossed a 256 boundary
            /*
             // This bounds checking is not required because of the way the getpix method works
             int dec = 1;
             int inc = 1;
             if (nc == 0) dec = 0;
             if (nc == width-1) inc = 0;
             x0 = iBuffer.getpix(nt,nc-dec);
             x1 = iBuffer.getpix(nt,nc+inc);
             */
            //x0 = theImage->getpix(nt,nc-1);
            x0 = idat(theImage,nt,nc-1, height,width);
            //x1 = theImage->getpix(nt,nc+1);
            x1 = idat(theImage,nt,nc+1,height,width);
            x1/=256;
            x0/=256;
            if(x1 != x0){
                // we crossed a boundary
                tweak = 256*(x1-x0);
                //theImage->setpix(nt,nc,theImage->getpix(nt,nc)+tweak);
                *(theImage + nt*width + nc) = idat(theImage,nt,nc, height,width)+tweak;
            }
        }
    }
    
}

OSStatus readHobj(CFURLRef url, TWOBYTE **head, TWOBYTE **trail, DATAWORD **data) {
    
    //return -1;
    
    int fd;
    char buffer[256];
    char* pointer;
    UInt8 txt[512];
    DATAWORD *datpt;
    TWOBYTE    *header,*trailer;   //trailer[TRAILEN/2];

    CFURLGetFileSystemRepresentation(url, true, txt, 512);
    
    fd = open((char*)txt,O_RDONLY);
    if(fd == -1) {
        return -1;
    }
    header = (TWOBYTE*) malloc(HEADLEN);
    if(header == NULL)
        return -1;
    trailer = (TWOBYTE*) malloc(TRAILEN);
    if(trailer == NULL)
        return -1;
    
    unsigned long filesize = fsize((char*)txt);
    read(fd,buffer,84);
    pointer = buffer+72;
    swap_bytes_routine(pointer,4,4);   // 4 bytes taken 4 at a time (1 int)
    int imagePixels = *((int*)pointer);
    
    swap_bytes_routine(buffer+80,4,2);   // 4 bytes taken two at a time (2 shorts)
    short rows = *((short*)buffer+40)+1;
    short cols = *((short*)buffer+41)+1;
    
    //unsigned char* ptr2 = new unsigned char[rows*6]; // skip over this part that has row, start pix, end pix for each row
    unsigned char* ptr2 = malloc(rows*6);
    if(ptr2 == 0) {
        close(fd);
        return -1;
    }
    
    read(fd,ptr2,rows*6);    // skip over the header
    //delete[] ptr2;
    free(ptr2);
    
    read(fd,buffer,16);
    swap_bytes_routine(buffer,16,4); // 16 bytes taken 4 at a time (4 ints)
    
    int samplesPerPix = *(int*)buffer;
    /* Notes on color hobj files:
     if there are 3 samples per pixel (an rgb color image), just reading in three times the pixels
     and decoding doesn't work. Assuming the order is RGB, the green in not coded the same way -- it is
     stored as little endian and has none of the crossing 256 boundry stuff. There also seems to be some
     extra information so there is an offset on the on the LHS of the image
     */
    unsigned long bytesPerPix = *((int*)buffer+1);    // not sure what this is
    int binary_file_bytes_per_data_point;
    bytesPerPix = filesize/samplesPerPix/imagePixels;
    
    binary_file_bytes_per_data_point = bytesPerPix;
    
    
    //printf("%d samples per pixel; %d bytes per pixel\n",samplesPerPix,bytesPerPix);
    close(fd);
    
    int bin_header = 84+rows*6 + 16;
    header[NCHAN] = cols;
    header[NTRAK] = rows*samplesPerPix;
    trailer[IS_COLOR]=1;

    // integer arguments to readBinary:
    //  bin_rows, bin_cols, bin_header, binary_file_bytes_per_data_point, swap_bytes_flag, unsigned_flag
    //readBinary(txt,theImage,rows*samplesPerPix,cols,bin_header,(int)bytesPerPix,1, 1);
    int swap_bytes_flag=1;
    
    // this stuff from read binary
    // reopen file and skip over header info
    
    fd = open((char*)txt,O_RDONLY);
    if(fd == -1) {
        return -1;
    }
    if(bin_header > 0) {
        ptr2 = malloc(bin_header);
        if(ptr2 == 0) {
            return -1;
        }
        read(fd,ptr2,bin_header);    // skip over the header
        free(ptr2);
    }
    
    int nbyte = (header[NCHAN]*header[NTRAK]) * binary_file_bytes_per_data_point;
    unsigned int nr,i,npts;
    short *sptr;
    char tc;
    int *iptr;
    
    npts = (header[NCHAN]*header[NTRAK]);
    datpt = (DATAWORD*) malloc(npts*sizeof(DATAWORD));
    
    if( binary_file_bytes_per_data_point == 1) {
        // allocate memory
        ptr2 = malloc(nbyte);
        if(ptr2 == 0) {
            close(fd);
            return -1;
        }
        // Read in the actual data
        nr = read(fd,ptr2, nbyte);
        //printf("%d Bytes read.\n",nr);
        close(fd);
        for(i=0; i<npts; i++){
            *(datpt+i) = *(ptr2+i);
        }
        free(ptr2);
    } else if( binary_file_bytes_per_data_point == sizeof(short)) {
        // allocate memory
        sptr = malloc(nbyte);
        if(sptr == 0) {
            close(fd);
            return -1;
        }
        // Read in the actual data
        nr = read(fd,sptr, nbyte);
        //printf("%d Bytes read.\n",nr);
        close(fd);
        if(swap_bytes_flag){
            // fiddle the byte order
            ptr2 = (unsigned char *)sptr;        // a copy of the data pointer
            for(i=0; i< nr; i+=2){
                tc = *(ptr2);
                *(ptr2) = *(ptr2+1);
                *(++ptr2) = tc;
                ptr2++;
            }
        }
        for(i=0; i<npts; i++){
                *(datpt+i) = *(sptr+i);
        }
        free(sptr);
    }  else if( binary_file_bytes_per_data_point == sizeof(int)){
        // allocate memory
        iptr = malloc(nbyte);
        if(iptr == 0) {
            close(fd);
            return -1;
        }
        // Read in the actual data
        nr = read(fd,iptr, nbyte);
        printf("%d Bytes read.\n",nr);
        close(fd);
        for(i=0; i<npts; i++){
            *(datpt+i) = *(iptr+i);
        }
        free(iptr);
    }
    // decode
    if(bytesPerPix == 2 ){
        decodeHobj(datpt, cols, rows*samplesPerPix);
    }
    // demosaic
    int bayer[2][2] = {0,1,1,2};
    int newrow = rows*samplesPerPix/2;
    int newcol = cols/2;
    int row,col;
    DATAWORD maxvalue,minvalue,value;
    DATAWORD* finalImage = malloc(newrow*3*newcol*sizeof(DATAWORD));
    DATAWORD* green = finalImage + newrow*newcol;
    DATAWORD* blue = finalImage + 2*newrow*newcol;
    maxvalue = minvalue = *datpt;
    for (row=0; row < newrow*2; row++) {
        int r=row/2;
        for (col=0; col < newcol*2; col++){
            int c=col/2;
            switch (bayer[row&1][col&1]){
                case 0:        // red
                    //red.setpix(r,c,theImage->getpix(row,col));   //*datp++;
                    value = idat(datpt,row,col,rows*samplesPerPix,cols);
                    *(finalImage+c+r*newcol)= value;
                    if(value > maxvalue){
                        maxvalue = value;
                        header[NMAX] = c;
                        header[LMAX] = r;
                    } else if (value <minvalue){
                        minvalue = value;
                        header[NMIN] = c;
                        header[LMIN] = r;

                    }
                    break;
                case 1:        // green
                    if(row&1){
                        //green.setpix(r,c,(green.getpix(r,c)+theImage->getpix(row,col))/2.);
                        value = (*(green+c+r*newcol)+idat(datpt,row,col,rows*samplesPerPix,cols))/2.;
                        *(green+c+r*newcol)=value;
                        if(value > maxvalue){
                            maxvalue = value;
                            header[NMAX] = c;
                            header[LMAX] = r + newrow;
                        } else if (value <minvalue){
                            minvalue = value;
                            header[NMIN] = c;
                            header[LMIN] =r + newrow;
                        }
                    } else {
                        //green.setpix(r,c,theImage->getpix(row,col));// =  *datp++;
                        *(green+c+r*newcol)=idat(datpt,row,col,rows*samplesPerPix,cols);
                    }
                    break;
                case 2:        // blue
                    //blue.setpix(r,c,theImage->getpix(row,col));// *pt_blue++ =  *datp++;
                    value = idat(datpt,row,col,rows*samplesPerPix,cols);
                    *(blue+c+r*newcol)= value;
                    if(value > maxvalue){
                        maxvalue = value;
                        header[NMAX] = c;
                        header[LMAX] = r + 2*newrow;
                    } else if (value <minvalue){
                        minvalue = value;
                        header[NMIN] = c;
                        header[LMIN] =r + 2*newrow;
                    }

                    break;
            }
        }
    }
    header[NCHAN] = newcol;
    header[NTRAK] = 3*newrow;
    free(datpt);
    *head = header;
    *data = finalImage;
    *trail = trailer;
    
    return noErr;

}

