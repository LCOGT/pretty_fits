'''
ColourImager - Python pipeline to combine 3 FITS files into a single colour JPEG/PNG
Copyright (C) 2017 Edward Gomez, Las Cumbres Observatory

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
from PIL import ImageFont, ImageDraw, Image
from fits2image.conversions import fits_to_jpg
from alipy.ident import run as alirun
from alipy.align import shape

import argparse
import glob
import logging
import numpy as np
import shutil
import subprocess
import sys, os
import tempfile


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

def planet_data_scale(data):
    data[data<0.]=0.
    std = np.sqrt(np.mean(data))
    cutoff = 3*std + np.mean(data)
    return data, cutoff

def scale_for_jupiter(filename):
    hdul = fits.open(filename)
    data = hdul[1].data
    min_val = np.percentile(data,99.9)
    data -= min_val
    data, cutoff = planet_data_scale(data)
    data[data<cutoff] = 0.
    data = data**2.
    scaled_planet = data*256./(data.max())
    return scaled_planet

def scale_for_moons(filename):
    hdul = fits.open(filename)
    data = hdul[1].data
    data[data<0.] = 0.
    data = np.arcsinh(data)
    moon_med = np.median(data)
    moon_max_val = np.percentile(data,99.95)
    data -= moon_med
    scaled_moon = data*256./(moon_max_val)
    return scaled_moon

def create_jupiter_image(infiles, outfile):
    imgs = []
    for infile in infiles:
        jupiter_data = scale_for_jupiter(infile)
        moons_data = scale_for_moons(infile)
        im_j = Image.fromarray(jupiter_data)
        im_m = Image.fromarray(moons_data)
        im = Image.blend(im_j.convert('RGB'), im_m.convert('RGB'), 0.4)
        if im.mode != 'L':
            im = im.convert('L')
        imgs.append(im)
    if len(imgs) == 3:
        # Only make a colour image if we have 3 files
        im_rgb = Image.merge('RGB', imgs)
    else:
        im_rgb = imgs[0]
    im_rgb.save(outfile)
    return

def create_saturn_image(infile, outfile):
    '''
    Makes a single colour Saturn image, rotating and cropping appropriately
    '''
    data = fits.getdata(infile)
    min_val = np.percentile(data,99.9)
    data -= min_val
    scaled_planet = data*256./(data.max())
    tmp_im = Image.fromarray(scaled_planet)
    tmp_im.crop((w/2-dp, h/2-dp, w/2+200, h/2+200))
    tmp_im.transpose(method=Image.ROTATE_90)
    tmp_im.convert('RGB')
    tmp_im.save(outfile)
    return


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
        # data = clean_data(data)
        data = scale_data(data, i)
        logging.warning('Shape of %s %s' % (file_in, str(data.shape)))
        rgb_list.append(data)
    return rgb_list


def read_write_data(filelist):
    '''
    Overwrite FITS files with cleaned and scaled data
    - Data is read into uncompressed FITS file to remove dependency on FPack
    '''
    img_list =[]
    for i, file_in in enumerate(filelist):
        new_filename = file_in.split('.')[0] + "_c.fits"
        if os.path.isfile(new_filename):
            # Reuse previously reduced images
            img_list.append(new_filename)
            continue
        data, hdrs = fits.getdata(file_in, header=True)
        data = clean_data(data)
        # data = scale_data(data, i)
        logging.warning('Shape of %s %s' % (file_in, str(data.shape)))
        hdu = fits.PrimaryHDU(data, header=hdrs)
        hdu.writeto(new_filename)
        img_list.append(new_filename)

    return img_list


def create_colour_stiff(img_list, filename='test.jpg', size=1500, tiff=False):
    resp = subprocess.call(['stiff']+img_list)
    if resp == 0 and not tiff:
        resp = subprocess.call(['convert', '-quality','70%', '-resize',size, 'stiff.tif', filename])
    elif resp == 0 and tiff:
        shutil.copyfile('stiff.tif',filename.replace('jpg','tiff'))

    return True


def reproject_files(ref_image, images_to_align, tmpdir='temp/'):
    identifications = alirun(ref_image, images_to_align[1:3], visu=False)
    outputshape = shape(ref_image)

    for id in identifications:
        if id.ok:
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=False, outdir=tmpdir)

    aligned_images = sorted(glob.glob(tmpdir+"/*_affineremap.fits"))

    img_list = [ref_image]+aligned_images

    return img_list


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_directory", help="directory for input files")
    parser.add_argument("-o", "--out_directory", help="directory for output files")
    parser.add_argument("-st", "--stiff", help="use STIFF for combining images", action="store_true")
    parser.add_argument("-t", "--tiff", help="Create a TIFF file", action="store_true")
    parser.add_argument("-z", "--fpack", help="Are the files rice compressed with fpack", action="store_true")
    parser.add_argument("-s", "--size", help="Size in pixels of x axis", default='1500')
    parser.add_argument("-n", "--name", help="Name of the output file (.jpg will be appended)")
    parser.add_argument("-p", "--planet", help="Planetary processing (no alignment)" , action="store_true")
    args = parser.parse_args()

    folder_path = args.in_directory
    tmpdir = tempfile.mkdtemp()

    ref_image, images_to_align, filename = select_images(folder=folder_path, fpacked=args.fpack)
    if args.planet:
        create_planet_image(images_to_align)
    if args.name:
        name = args.name.replace('.jpg','')
        filename = "{}.jpg".format(name)
    img_list = read_write_data(images_to_align)
    img_list = reproject_files(img_list[0], img_list, tmpdir)
    filename = os.path.join(args.out_directory, filename)
    if args.stiff:
        create_colour_stiff(img_list, filename, size=args.size, tiff=args.tiff)
    else:
        fits_to_jpg(img_list, filename, width=args.size, height=args.size, color=True)

    # Remove the temporary files
    shutil.rmtree(tmpdir)
