from astropy.io import fits
from astroscrappy import detect_cosmics
import PIL.Image as pli
import numpy as np
import alipy
import glob
import logging

def reshape(data):
    return data.shape

def remove_cr(data):
    m, imdata = detect_cosmics(data, readnoise=20., gain=1.4, sigclip=4., sigfrac=.2, objlim=6.)
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
    logging.warning('--- Begin Scaling ---')
    median = np.median(data)
    data[data<0.]=median
    # Run astroScrappy to remove pesky cosmic rays
    data = remove_cr(data)
    logging.warning('Median=%s' % median)
    logging.warning('Max after median=%s' % data.max())
    return data


def scale_data(data):
    # Recalculate the median
    data[data<0.]=0.
    median = np.median(data)
    data-= median
    data[data<0.]=0.
    sc_data= data
    max_val = np.percentile(sc_data,99.5)
    logging.warning('99.5 =%s' % max_val)
    scaled = sc_data*255./(max_val)
    scaled[scaled>255]=255
    logging.warning('Median of scaled=%s' % np.median(scaled))
    logging.warning('Min scaled=%s' % scaled.min())
    return scaled

def select_images(folder='temp/'):
    images_to_align = sorted(glob.glob("%s*.fits" % folder))
    ref_image = images_to_align[0]
    return ref_image, images_to_align

def read_aligned(filelist):
    # Scale the images
    rgb_list =[]
    for file_in in filelist:
        data, hdrs = fits.getdata(file_in, header=True)
        data = clean_data(data)
        data = scale_data(data)
        rgb_list.append(data)
    return rgb_list

def create_colour(rgb_list):
    rgb_cube = np.dstack(rgb_list).astype(np.uint8)[::-1, :, :]  # make cube, flip vertical axis
    colour_img = pli.fromarray(rgb_cube)
    colour_img.save('test.jpg')
    colour_img.show()
    return

if __name__ == '__main__':
    ref_image, images_to_align = select_images()

    identifications = alipy.ident.run(ref_image, images_to_align[1:], visu=False)
    outputshape = alipy.align.shape(ref_image)

    for id in identifications:
        if id.ok:
        # Variant 1, using only scipy and the simple affine transorm :
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

    aligned_images = sorted(glob.glob("alipy_out/*.fits"))

    rgb_list = read_aligned([ref_image]+aligned_images)

    create_colour(rgb_list)
