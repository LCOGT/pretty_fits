# colourImager
Python pipeline to combine 3 FITS files into a single colour JPEG/PNG

# Installation
In addition to the `requirement.txt` will need to have:
- [SExtractor](http://www.astromatic.net/software/sextractor) from Astromatics. It finds all the sources (i.e. things that look like stars) in each image.
- [alipy](http://obswww.unige.ch/~tewes/alipy/installation.html) calls SExctractor and does a cross-correlation to match features in the images.
- The rest of the requirements can be installed with `pip install -r requirements.txt`.

# Running
- Create a directory called `temp` at the same level as `image_ali.py`
- Put the 3 FITS files you want to combine in this directory
- Make the directory `temp` and `image_ali.py` sit in, writable. The process creates temporary files.
- Run `python image_ali.py`
