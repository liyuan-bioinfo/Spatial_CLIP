# Transform OME-ITFF to the compatible tif image (pyramidal tiff)
(Introduction)[https://andrewjanowczyk.com/converting-an-existing-image-into-an-openslide-compatible-format/]
(Tool:ImageMagick)[https://github.com/ImageMagick/ImageMagick/blob/main/Install-unix.txt]
Download & Unpack
  
Download
  https://imagemagick.org/archive/ImageMagick.tar.gz.
    
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

# https://andrewjanowczyk.com/converting-an-existing-image-into-an-openslide-compatible-format/
# https://iipimage.sourceforge.io/documentation/images/
# Conversion using ImageMagick
You can also use ImageMagick (version 6.4.7-10 and upwards) to create Tiled Pyramid TIFF. In this case use the convert command. 
For example, to generate a 256x256 pyramid tiled tiff using JPEG compression:

    $ convert input -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:output.tif'
    input: source image - can be any kind of image supported by ImageMagick (almost all)
    ptif: specify your image format as pyramid tiff
    output.tif: output image - *MUST* have tif or ptif extension;
    256x256: the tile size    

    $ convert demo/17642_500_f00001_original.jpg -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:demo/output.tiff'
    $ convert 20190715_PDAC_65289_ROI2_4_plex/148Nd_pan-Keratin.ome.tiff -define tiff:tile-geometry=256x256 -compress jpeg 'ptif:20190715_PDAC_65289_ROI2_4_plex_convert/148Nd_pan-Keratin.ome.tiff'
  