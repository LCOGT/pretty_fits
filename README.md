# colourImager
Python pipeline to combine 3 FITS files into a single colour JPEG/PNG

# Installation
In addition to the `requirement.txt` will need to have:
- [SExtractor](http://www.astromatic.net/software/sextractor) from Astromatics. It finds all the sources (i.e. things that look like stars) in each image.
- [alipy](http://obswww.unige.ch/~tewes/alipy/installation.html) calls SExctractor and does a cross-correlation to match features in the images.
- (optional) [STiff](http://www.astromatic.net/software/stiff) from Astromatics. It produces well balenced greyscale and colour images from FITS files.
- (optional) [ImageMagick](http://www.imagemagick.org) to convert the TIFF file STiff will generate to JPEG/PNG and give us a good level of compression.
- The rest of the requirements can be installed with `pip install -r requirements.txt`.

# Running
- In `image_ali.py` change `DATA_DIR` to where ever you will be putting data
- Put the 3 FITS files you want to combine into a sub-directory of `DATA_DIR`
- Run `python image_ali.py -d <sub-directory name>`
- If you are using STIFF add the `-st` command-line option
- If you are using FPacked data, use the `-z` option
