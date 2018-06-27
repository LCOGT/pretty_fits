import six
from astropy.io import fits
from datetime import datetime
from fits2image.conversions import fits_to_jpg
from numpy import shape
import os
from scipy.ndimage import interpolation

from .quad import calibrate
from .file_utils import make_filename, write_temp_fits, clean_up_temps

def reproject_data(filename, transform):

    inv = transform.inverse()
    (matrix, offset) = inv.matrixform()

    hdu = fits.open(filename)
    data = hdu[1].data
    datashape = shape(data)
    # LCO 0.4m data x and y are flipped
    tel = hdu[1].header['TELID'][0:3]
    # if tel == '0m4':
        # offset = (offset[1], offset[0])
    data = interpolation.affine_transform(data, matrix, offset=offset, output_shape = datashape)
    tempfile = write_temp_fits(data)
    return  tempfile, datashape

def calc_offsets(files, planet=False):
    offsets = []
    ref = files[0]
    imgs = files[1:]
    offsets.append((0,0))
    if planet:
        (refx,refy) = planet_centre_coord(ref)
        for file in files:
            (x,y) = planet_centre_coord(imgs)
            offsetx = int(refx - x)
            offsety = int(refy - y)
            offsets.append((offsetx, offsety))
    else:
        idents = calibrate(ref, imgs)
        for ident in idents:
            (a,b,c,d) = ident.trans.v
            offsets.append((int(c),int(d)))
    return offsets

def planet_centre_coord(filename):
    # Extract photometry info for brightest object i.e. the planet
    hdu = fits.open(filename)
    x = hdu[2].data[0][0]
    y = hdu[2].data[0][1]
    return x,y

def crop_align(filename, offsetx, offsety):
    im = Image.open(filename)
    img = ImageChops.offset(im, offsetx, offsety)
    area = (200, 0, 1300, 1000)
    img.crop(area).save(filename)
    return

def planet_align_process(file_list, out_dir, stamp=False):
    offsets = calc_offsets(file_list)
    for idx, files in enumerate(file_list):
        outfile = make_filename(files[0], out_dir)
        create_jupiter_image(files, outfile)
        (offsetx, offsety) = offsets[idx]
        if stamp:
            label_text = date_obs.strftime('%Y-%m-%d %H:%M:%S')
        crop_align(outfile, offsetx, offsety)
        if stamp:
            add_label(outfile, label_text, label_font='/Library/Fonts/Courier New.ttf')
    return

def align_process(file_list, out_dir, fits_file=False, planet=False, stamp=False, colour=False):
    '''
    Takes a list of FITS files (FZ format) and aligns them based on the 1st file in list
    Saves either a sequence of aligned B&W JPEGs or a single 3 colour JPEG.

    Returns True if files saved successfully

    :param file_list: list of FITS files each with photometry table in a FITS extension.
                        First in the list is the reference image
    :type file_list: list

    :param out_dir: Full path to writable directory on the local files system where processed
                    files will be saved.
    :type out_dir: string

    :param fits: Return FITS files. Default - False will return JPEGs.
                Files will be saved to `out_dir`
    :type fits: boolean

    :param planet: Align on planet
    :type planet: boolean

    :param stamp: Adds the time stamp to the images (JPEG only)
    :type stamp: boolean
    '''
    if planet:
        return planet_align_process(file_list, out_dir, stamp)
    else:
        ref = file_list[0]
        imgs = file_list[1:]
        calibrations = calibrate(ref, imgs)
        temp_files = []
        datashape = shape(fits.getdata(ref))
        for idx, c in enumerate(calibrations):
            if c.ok:
                tmp_file, datashape = reproject_data(c.ukn.filepath, c.trans)
                if not colour:
                    outfile = make_filename(c.ukn.filepath, out_dir)
                    fits_to_jpg(tmp_file, outfile, width=datashape[0], height=datashape[1], median=True, percentile=99.95)
                temp_files.append(tmp_file)
            else:
                temp_files.append(imgs[idx])
        outfile = make_filename(ref, out_dir)
        if colour:
            fits_to_jpg([ref]+temp_files, outfile, width=datashape[0], height=datashape[1], color=True)
        else:
            fits_to_jpg(ref, outfile, width=datashape[0], height=datashape[1], median=True, percentile=99.95)
        clean_up_temps(temp_files)
    return
