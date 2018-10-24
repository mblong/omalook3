#include <CoreFoundation/CoreFoundation.h>
#include <CoreServices/CoreServices.h>
#include <QuickLook/QuickLook.h>
#include "oma_bits.h"

OSStatus get_oma_data( CFURLRef,TWOBYTE **ptrh,TWOBYTE **ptrt, DATAWORD **ptrd);

Ptr Get_rgb_from_image_buffer( TWOBYTE *header, TWOBYTE *trailer, DATAWORD *datpt, int pixsiz);
Ptr Get_color_rgb_from_image_buffer( TWOBYTE *header, TWOBYTE *trailer, DATAWORD *datpt, int pixsiz);
CGImageRef LoadImageFromData(int width, int height, void* imagedata);
void GWorldImageBufferRelease(void*	inInfo, const void*	inData, size_t	inSize );
int get_oma_data_type( CFURLRef url);
int getpivdata(CFURLRef url, int* no_of_velx, int* no_of_vely, int* fftsize, int* boxinc, float **xpeaks, float **ypeaks);
int plotpivdata( CGContextRef cgContext, CGRect dstRect, float pen_width, int no_of_velx, int no_of_vely, int fftsize, int boxinc, float* xpeaks, float* ypeaks);
/* -----------------------------------------------------------------------------
    Generate a thumbnail for file

   This function's job is to create thumbnail for designated file as fast as possible
   ----------------------------------------------------------------------------- */

OSStatus GenerateThumbnailForURL(void *thisInterface, QLThumbnailRequestRef thumbnail, CFURLRef url, CFStringRef contentTypeUTI, CFDictionaryRef options, CGSize maxSize)
{
	CGContextRef cgContext;
	CGSize canvasSize = {400.,400.};
	CGRect dstRect = {0.,0.,399., 399.};
	Ptr rgbdata;
	
	DATAWORD *datpt;
	TWOBYTE *header;
    TWOBYTE *trailer;
	
	int pixsize = 1,maxdimension,width,height,data_type;
	
	CGImageRef image=0;
	
	CGDataProviderRef	provider=0;
	CGColorSpaceRef		colorspace;
	size_t			rowbytes;
	
	datpt = 0;
	header = 0;
	rgbdata = 0;
	
	data_type = get_oma_data_type( url);
	
	if( data_type == IMAGEDATA || data_type == UNKNOWN || data_type == HOBJ){
        if ( data_type == HOBJ) {
            if (readHobj(url, &header, &trailer, &datpt) != noErr) {
                return -1;
            }
        }else {
            if (get_oma_data(url, &header, &trailer, &datpt) != noErr) {
                return -1;
            }
        }
        if(trailer[IS_COLOR]){
            if( header[NCHAN] > header[NTRAK]/3 )
                maxdimension =  header[NCHAN];
            else
                maxdimension =  header[NTRAK]/3;
            
            pixsize = 1 + maxdimension/THUMBSIZE;
            if(pixsize > 1) pixsize = -pixsize;
            
            width = header[NCHAN]/abs(pixsize);
            height = header[NTRAK]/abs(pixsize)/3;

            if( width == 0 || height == 0 ){
                if(header != 0) free(header);
                if(datpt != 0) free(datpt);
                if(trailer != 0) free(trailer);
                return -1;
            }
            
            rgbdata = Get_color_rgb_from_image_buffer(header, trailer, datpt, pixsize);
            
        } else {
		
            if( header[NCHAN] > header[NTRAK] )
                maxdimension =  header[NCHAN];
            else
                maxdimension =  header[NTRAK];
            
            pixsize = 1 + maxdimension/THUMBSIZE;
            if(pixsize > 1) pixsize = -pixsize;
            
            width = header[NCHAN]/abs(pixsize);
            height = header[NTRAK]/abs(pixsize);
            
            if( width == 0 || height == 0 ){
                if(header != 0) free(header);
                if(datpt != 0) free(datpt);
                if(trailer != 0) free(trailer);
                return -1;
            }
                
            rgbdata = Get_rgb_from_image_buffer(header, trailer, datpt, pixsize);
        }
		
		if(header != 0) free(header);
		if(datpt != 0) free(datpt);
        if(trailer != 0) free(trailer);

		if(rgbdata == NULL){
			return -1;
		}
		canvasSize.width = width;
		canvasSize.height = height;
		dstRect.size = canvasSize;

		cgContext = QLThumbnailRequestCreateContext(thumbnail, *(CGSize *)&canvasSize, true, NULL);
		 if(cgContext) {
			
			colorspace = CGColorSpaceCreateDeviceRGB();
			rowbytes = width * 4;
			provider = CGDataProviderCreateWithData( NULL,rgbdata, height * rowbytes,
					GWorldImageBufferRelease );

			if(provider){
				// create an image
				image = CGImageCreate( width, height, 8, 32, rowbytes, colorspace,
						kCGImageAlphaNoneSkipFirst, provider, NULL, false, kCGRenderingIntentDefault );
				
				
				CGContextDrawImage (cgContext,dstRect,image);
			}
			QLThumbnailRequestFlushContext(thumbnail, cgContext);
			CFRelease(cgContext);
		
			if(image !=0) CGImageRelease (image);
			//if(rgbdata != 0) free(rgbdata);
			if(provider !=0) CGDataProviderRelease(provider);
		}	

		return noErr;
	} else if ( data_type == PIVDATA) {
		// plot PIV thumbnail
		int no_of_velx, no_of_vely, fftsize, boxinc;
		float *xpeaks,*ypeaks;
		if(  getpivdata( url, &no_of_velx, &no_of_vely, &fftsize, &boxinc, &xpeaks, &ypeaks) == 1){
			//successfully got PIV dat
			if( no_of_velx > no_of_vely )
				maxdimension =  no_of_velx;
			else
				maxdimension =  no_of_vely;

			
			canvasSize.width =  no_of_velx*THUMBSIZE;
			canvasSize.width /=  maxdimension;
			canvasSize.height = no_of_vely*THUMBSIZE;
			canvasSize.height /=  maxdimension;
			dstRect.size = canvasSize;
			
			cgContext = QLThumbnailRequestCreateContext(thumbnail, *(CGSize *)&canvasSize, false, NULL);
			if(cgContext) {
				plotpivdata(cgContext, dstRect, 0.5, no_of_velx, no_of_vely, fftsize, boxinc, xpeaks, ypeaks);

				QLThumbnailRequestFlushContext(thumbnail, cgContext);
				CFRelease(cgContext);

			}
			free(xpeaks);
			free(ypeaks);
			return noErr;

		}
	}
	return -1;

}

void CancelThumbnailGeneration(void* thisInterface, QLThumbnailRequestRef thumbnail)
{
    // implement only if supported
}

