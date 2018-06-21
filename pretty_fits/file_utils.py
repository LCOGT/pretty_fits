from astropy.io import fits
from datetime import datetime
from glob import glob
import os
import logging
import requests

from pretty_fits.pretty_fits import create_jupiter_image

import settings

FILTER_ORDER = {
            'jupiter'   : {'up':1,'B':2,'zs':3},
            'mars'      : {},
            'saturn'    : {'zs':1}
                }

def lco_api_call(url, token):
    '''
    Get status of RequestID from the Valhalla API
    '''
    headers = {'Authorization': 'Token {}'.format(token)}
    try:
        r = requests.get(url, headers=headers, timeout=20.0)
    except requests.exceptions.Timeout:
        msg = "Observing portal API timed out"
        logging.error(msg)
        params['error_msg'] = msg
        return False, msg

    if r.status_code in [200,201]:
        logging.debug('Recieved data')
        return True, r.json()
    else:
        logging.error("Could not send request: {}".format(r.content))
        return False, r.content

def get_archive_data(out_dir, request_id):
    # Only look for data which has completed
    url = "{}{}".format(settings.PORTAL_REQUEST_API, request_id)
    state, r = lco_api_call(url, settings.VALHALLA_TOKEN)
    if not state:
        logging.debug('Failed')
        return
    for req in r['requests']:
        if req['state'] == 'COMPLETED':
            subreq_id = req['id']
            url = "{}?REQNUM={}&ordering=-id&RLEVEL=91".format(settings.ARCHIVE_FRAMES_API, subreq_id)
            success, r = lco_api_call(url, settings.ARCHIVE_TOKEN)
            if success:
                out_path = os.path.join(out_dir, request_id, str(subreq_id))
                dl_sort_data_files(r, out_path)
    return False

def download_file(out_file, url):
    logging.debug("Downloading: {}".format(out_file))
    dt = requests.get(url)
    f = open(out_file, 'wb')
    f.write(dt.content)
    f.close()
    return

def dl_sort_data_files(r, out_path):
    # Make subdirectory
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    for response in r['results']:
        out_file = os.path.join(out_path, response['filename'])
        download_file(out_file, response['url'])
    return

def sort_jupiter_data(out_path, filter_name):
    print("Sorting data in {}".format(out_path))
    files = {1:'',2:'', 3:''}
    if filter_name:
        path_match = "{}*.fits.fz".format(filter_name)
    else:
        path_match = "*.fits.fz"
    file_list = glob(os.path.join(out_path, path_match))
    for fits_file in file_list:
        hdr = fits.open(fits_file)
        new_filename = os.path.join(out_path,"{}-{}.fits.fz".format(hdr[1].header['FILTER'], hdr[1].header['SITEID']))
        os.rename(os.path.join(out_path, fits_file),  new_filename)
        filter_index = FILTER_ORDER['jupiter'][hdr[1].header['FILTER']]
        files[filter_index] = new_filename
    return files

def find_directories(out_dir, request_id, filter_name=None):
    print("Finding dirs: {} {}".format(out_dir, request_id))
    files_grouped = []
    for dirname, dirnames, filenames in os.walk(os.path.join(out_dir, request_id)):
        for subdirname in dirnames:
            filepath = os.path.join(dirname, subdirname)
            infiles = sort_jupiter_data(filepath, filter_name)
            if filter_name:
                filter_index = FILTER_ORDER['jupiter'][filter_name]
                input_files = [infiles[filter_index]]
            else:
                input_files = [infiles[1], infiles[2], infiles[3]]
            files_grouped.append(input_files)
    return files_grouped

def find_files(in_dir):
    path_match = "*.fits.fz"
    file_list = sorted(glob(os.path.join(out_path, path_match)))
    return file_list
