# Pretty FITS

This package aims to make nice looking Red, Green, Blue JPGs from the corresponding FITS files. FITS data handling is done thanks to [astropy](https://github.com/astropy) and it uses [astroscrappy](https://github.com/astropy/astroscrappy) to remove unsightly cosmic ray hits. It also aligns images using my fork of the [ALIpy wrapper](https://github.py/LCOGT/alipy/) for SExtractor.

## Requirements
- SExtractor (I'm really sorry about this but it is the best way I've found to align the images)
