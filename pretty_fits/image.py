from PIL import Image, ImageChops, ImageDraw, ImageFont
from astropy.io import fits



def add_label(filename, label_text, label_font):
    '''
        Writes label_text into Pillow Image image.
    :param image: Pillow Image Object - will be modified with the label_text
    :param label_text:
    :param label_font: The true type font on the system to use in the label. throws IOError if it can't be found.
    :return:
    '''
    image = Image.open(filename)
    (width, height) = image.size
    font_size = 20
    offset = 4
    font = ImageFont.truetype(label_font, font_size)
    (font_width, font_height) = font.getsize(label_text)
    while font_width > (int(width)-offset):
        font_size -= 1
        font = ImageFont.truetype(label_font, font_size)
        (font_width, font_height) = font.getsize(label_text)

    d = ImageDraw.Draw(image)
    d.text((offset, int(height) - offset - font_height), label_text, font=font, fill=255)
    image.save(filename)
    return

def crop_align(filename, offsetx, offsety):
    im = Image.open(filename)
    img = ImageChops.offset(im, offsetx, offsety)
    area = (200, 0, 1300, 1000)
    img.crop(area).save(filename)
    return

def align_files(file_group, out_dir, planet=False, stamp=False):
    offsets = calc_offsets(file_group)
    for idx, files in enumerate(file_group):
        hdu = fits.open(files[0])
        date_obs = datetime.strptime(hdu[1].header['date-obs'][0:19],'%Y-%m-%dT%H:%M:%S')
        filename = "{}.jpg".format(date_obs.strftime('%Y%m%d%H%M%S'))
        outfile = os.path.join(out_dir, filename)
        create_jupiter_image(files, outfile)
        (offsetx, offsety) = offsets[idx]
        if stamp:
            label_text = date_obs.strftime('%Y-%m-%d %H:%M:%S')
        crop_align(outfile, offsetx, offsety)
        if stamp:
            add_label(outfile, label_text, label_font='/Library/Fonts/Courier New.ttf')
    return
