'''
ColourImager - Python pipeline to combine 3 FITS files into a single colour JPEG/PNG
Copyright (C) 2015 Edward Gomez, LCOGT

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
from astropy.io import fits
from astroscrappy import detect_cosmics
import PIL.Image as pli
from PIL import ImageFont, ImageDraw
import numpy as np
import alipy
import glob
import logging
import sys, os
import argparse
import tempfile
import shutil
import subprocess

DATA_DIR = '/Users/egomez/Downloads/'

def reshape(data):
    return data.shape

def remove_cr(data):
    '''
    Removes high value pixels which are presumed to be cosmic ray hits.
    '''
    m, imdata = detect_cosmics(data, readnoise=20., gain=1.4, sigclip=5., sigfrac=.5, objlim=6.)
    return imdata

def clean_data(data):
    '''
    - Remove bogus (i.e. negative) pixels
    - Remove Cosmic Rays
    - Find the 99.5% percentile and remove everything above that
    - Subtract the median sky value
    - Scale the images in the range for JPEG
    '''
    # Level out the colour balance in the frames
    logging.warning('--- Begin CR removal ---')
    median = np.median(data)
    data[data<0.]=median
    # Run astroScrappy to remove pesky cosmic rays
    data = remove_cr(data)
    logging.warning('Median=%s' % median)
    logging.warning('Max after median=%s' % data.max())
    return data


def scale_data(data, i):
    # Recalculate the median
    logging.warning('--- Begin Scaling ---')
    data[data<0.]=0.
    median = np.median(data)
    data-= median
    data[data<0.]=0.
    sc_data= data #np.arcsinh(data)
    max_val = np.percentile(sc_data,99.5)
    logging.warning('99.5 =%s' % max_val)
    scaled = sc_data*255./(max_val)
    scaled[scaled>255.]=255.
    logging.warning('Median of scaled=%s' % np.median(scaled))
    logging.warning('Min scaled=%s' % scaled.min())
    return scaled

def select_images(folder='temp', fpacked=False):
    if fpacked:
        filetype = "%s/*.fits.fz"
        images_to_align = sorted(glob.glob(filetype % folder))
        ref_image = images_to_align[0]
        hdr = fits.getheader(ref_image, 1)
    else:
        filetype = "%s/*.fits"
        images_to_align = sorted(glob.glob(filetype % folder))
        ref_image = images_to_align[0]
        hdr = fits.getheader(ref_image)
    objectname = hdr['OBJECT'].strip().replace(' ','_')
    filename = '%s-%s.jpg' % (objectname, hdr['REQNUM'])
    return ref_image, images_to_align, filename

def read_aligned(filelist):
    # Scale the images
    rgb_list =[]
    for i, file_in in enumerate(filelist):
        data, hdrs = fits.getdata(file_in, header=True)
        data = clean_data(data)
        data = scale_data(data, i)
        logging.warning('Shape of %s %s' % (file_in, str(data.shape)))
        rgb_list.append(data)
    return rgb_list

def create_colour_simple(img_list, filename='test.jpg', object_name='Unknown', credit=False, preview=False):
    rgb_list = read_aligned(img_list)
    rgb_cube = np.dstack(rgb_list).astype(np.uint8)[::-1, :, :]  # make cube, flip vertical axis
    colour_img = pli.fromarray(rgb_cube)
    if credit:
        font = ImageFont.truetype("Helvetica.ttf", 40)
        textstamp = 'Las Cumbres Observatory Global Telescope Network - %s' % object_name
        draw = ImageDraw.Draw(colour_img)
        draw.text((50, 10), textstamp, font=font, fill=(255,255,255))
    colour_img.save(filename)
    if preview:
        colour_img.show()

    return True


def read_write_data(filelist):
    '''
    Overwrite FITS files with cleaned and scaled data
    - Data is read into uncompressed FITS file to remove dependency on FPack
    '''
    img_list =[]
    for i, file_in in enumerate(filelist):
        data, hdrs = fits.getdata(file_in, header=True)
        data = clean_data(data)
        # data = scale_data(data, i)
        logging.warning('Shape of %s %s' % (file_in, str(data.shape)))
        new_filename = file_in.split('.')[0] + "_c.fits"
        hdu = fits.PrimaryHDU(data, header=hdrs)
        hdu.writeto(new_filename)
        img_list.append(new_filename)

    return img_list


def create_colour_stiff(img_list, filename='test.jpg'):
    resp = subprocess.call(['stiff']+img_list)
    if resp == 0:
        resp = subprocess.call(['convert', '-quality','70%', '-resize','1500', 'stiff.tif', filename])

    return True


def reproject_files(ref_image, images_to_align, tmpdir='temp/'):
    identifications = alipy.ident.run(ref_image, images_to_align[1:3], visu=False)
    outputshape = alipy.align.shape(ref_image)

    for id in identifications:
        if id.ok:
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True, outdir=tmpdir)

    aligned_images = sorted(glob.glob(tmpdir+"/*_affineremap.fits"))

    img_list = [ref_image]+aligned_images

    return img_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="Comma separated list of 3 files")
    parser.add_argument("-d", "--directory", help="directory of files inside Downloads")
    parser.add_argument("-c", "--credit", help="apply a standard LCOGT credit watermark", action="store_true")
    parser.add_argument("-st", "--stiff", help="use STIFF for combining images", action="store_true")
    parser.add_argument("-p", "--preview", help="show a PIL generated JPEG preview", action="store_true")
    parser.add_argument("-z", "--fpack", help="Are the files rice compressed with fpack", action="store_true")
    args = parser.parse_args()

    folder_name = args.directory
    folder_path = os.path.join(DATA_DIR, folder_name)
    tmpdir = tempfile.mkdtemp()

    ref_image, images_to_align, filename = select_images(folder=folder_path, fpacked=args.fpack)


    if args.stiff:
        img_list = read_write_data(images_to_align)
        img_list = reproject_files(img_list[0], img_list, tmpdir)
        create_colour_stiff(img_list, filename)
    else:
        img_list = reproject_files(ref_image, images_to_align, tmpdir)
        create_colour_simple(img_list, filename, object_name=folder_name, credit=args.credit, preview=args.preview)

    # Remove the temporary files
    shutil.rmtree(tmpdir)
