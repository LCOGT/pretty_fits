from astropy.io import fits
from datetime import datetime
from glob import glob
import os
import logging
import requests
from tempfile import mkstemp

from .clean import create_jupiter_image

logger = logging.getLogger(__name__)

def find_files(in_dir):
    path_match = "*.fits.fz"
    file_list = sorted(glob(os.path.join(in_dir, path_match)))
    return file_list

def write_temp_fits(data):
    tempcode, temp = mkstemp(suffix=".fits")
    hdu = fits.PrimaryHDU(data)
    hdul = fits.HDUList([hdu])
    hdul.writeto(temp, overwrite=True)
    return temp

def make_filename(filename, out_dir):
    hdu = fits.open(filename)
    date_obs = datetime.strptime(hdu[1].header['date-obs'][0:19],'%Y-%m-%dT%H:%M:%S')
    filename = "{}.jpg".format(date_obs.strftime('%Y%m%d%H%M%S'))
    outfile = os.path.join(out_dir, filename)
    return outfile

def clean_up_temps(temp_files):
    for tmp in temp_files:
        os.remove(tmp)
    return
