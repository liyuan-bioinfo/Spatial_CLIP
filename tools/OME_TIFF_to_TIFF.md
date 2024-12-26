# Transform OME-ITFF to the compatible tif image (pyramidal tiff)
- [Introduction](https://andrewjanowczyk.com/converting-an-existing-image-into-an-openslide-compatible-format/)

## Details of installing ImageMagick on Linux system
[Download](https://imagemagick.org/archive/ImageMagick.tar.gz)

[Install](https://github.com/ImageMagick/ImageMagick/blob/main/Install-unix.txt)
    
Unpack

    $ tar xvfz ImageMagick.tar.gz
    
Configure
    
    $ cd ImageMagick-7.1.1
    $ ./configure
  
  In most cases you will simply want to compile ImageMagick with this command:
  
    $ make
    
Install
  
  Now that ImageMagick is configured and built, type:
  
    $ make install
  

  To confirm your installation of the ImageMagick distribution was successful,
  ensure that the installation directory is in your executable search path
  and type:  
    $ display

## Conversion using ImageMagick
    $ convert input -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:output.tif'
    input: source image - can be any kind of image supported by ImageMagick (almost all)
    ptif: specify your image format as pyramid tiff
    output.tif: output image - *MUST* have tif or ptif extension;
    256x256: the tile size    

    $ convert demo/17642_500_f00001_original.jpg -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:demo/output.tiff'
    $ convert 20190715_PDAC_65289_ROI2_4_plex/148Nd_pan-Keratin.ome.tiff -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:20190715_PDAC_65289_ROI2_4_plex_convert/148Nd_pan-Keratin.ome.tiff'
  
