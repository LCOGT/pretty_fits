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
    m, imdata = detect_cosmics(data, readnoise=20., gain=1.4, sigclip=4., sigfrac=.2, objlim=5.)
    return imdata

def scale_data(data):
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
    # Recalculate the median because there is nothing less than the media now
    median = np.median(data)
    data-= median
    logging.warning('Max after median=%s' % data.max())
    max_val = np.percentile(data,99.5)
    logging.warning('99.5 =%s' % max_val)
    scaled = data*255./(max_val)
    scaled[scaled>255.]=255.
    logging.warning('Median of scaled=%s' % np.median(scaled))
    return scaled

if __name__ == '__main__':
    images_to_align = sorted(glob.glob("temp/*.fits"))
    ref_image = images_to_align[0]

    identifications = alipy.ident.run(ref_image, images_to_align[1:], visu=False)
    outputshape = alipy.align.shape(ref_image)
    # This is simply a tuple (width, height)... you could specify any other shape.

    for id in identifications:
        if id.ok:
        # Variant 1, using only scipy and the simple affine transorm :
            alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

    aligned_images = sorted(glob.glob("alipy_out/*.fits"))

    B_data, B_hdrs = fits.getdata(ref_image, header=True)
    V_data, V_hdrs = fits.getdata(aligned_images[0], header=True)
    R_data, R_hdrs = fits.getdata(aligned_images[1], header=True)
    # Scale the images
    rgb_list = [ scale_data(R_data), scale_data(V_data), scale_data(B_data)] # make sure these are in range 0-255
    rgb_cube = np.dstack(rgb_list).astype(np.uint8)[::-1, :, :]  # make cube, flip vertical axis
    color_img = pli.fromarray(rgb_cube)
    color_img.save('test.jpg')
    color_img.show()
