import pyfits as pf
import PIL.Image as pli
import numpy as np

import alipy
import glob

def reshape(data):
    return data.shape

def scale_data(data):
    data1 = data.reshape(data.shape[0]*data.shape[1])
    max_val = np.percentile(data1,99.95)
    scaled = data*256./max_val
    new_scaled = np.ma.masked_greater(scaled, 255.)
    new_scaled.fill_value=255.
    img_data = new_scaled.filled()
    small_data =np.resize(img_data, (2000,2000))
    return img_data

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

    B_data, B_hdrs = pf.getdata(ref_image, header=True)
    V_data, V_hdrs = pf.getdata(aligned_images[0], header=True)
    R_data, R_hdrs = pf.getdata(aligned_images[1], header=True)
    # do your scaling / matching / etc operations here
    print reshape(B_data), reshape(V_data), reshape(R_data)
    rgb_list = [ scale_data(R_data), scale_data(V_data), scale_data(B_data)] # make sure these are in range 0-255
    rgb_cube = np.dstack(rgb_list).astype(np.uint8)[::-1, :, :]  # make cube, flip vertical axis
    color_img = pli.fromarray(rgb_cube)
    color_img.save('test.jpg')
