import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from glob import glob


def white_edge_cutter(img):
    """
    cut white edge of image
    :param img: array, rgb value of a image
    :return: function to cut white edge
    """
    left = 0
    while all(img[:, left, :].sum(axis=1) == 3):
        left += 1
    right = img.shape[1] - 1
    while all(img[:, right, :].sum(axis=1) == 3):
        right -= 1
    top = 0
    while all(img[top, :, :].sum(axis=1) == 3):
        top += 1
    bottom = img.shape[0] - 1
    while all(img[bottom, :, :].sum(axis=1) == 3):
        bottom -= 1

    def cutter(img):
        return img[top-1:bottom, left-1:right, :]
    return cutter


def merge_img(imgs, ncol):
    """
    merge images
    :param ncol: number of images each row
    :param imgs: images. iterable
    :return: merged big image
    """
    BLANK_WIDTH = 0.05
    BLANK_HEIGHT = 0.01
    w, h = imgs[0].shape[1], imgs[0].shape[0]
    blankw = int(w * BLANK_WIDTH)
    blankh = int(h * BLANK_HEIGHT)
    rows = []
    for i in range(0, len(imgs), ncol):
        temp = imgs[i]
        for j in range(1, ncol):
            if i+j < len(imgs):
                temp = np.hstack([temp, np.ones([h, blankw, 3]), imgs[i+j]])
        rows.append(temp)
    merged_width = rows[0].shape[1]
    if rows[-1].shape[1] != merged_width:
        diff = merged_width - rows[-1].shape[1]
        conwid1 = diff//2
        conwid2 = diff - conwid1
        rows[-1] = np.hstack([np.ones([h, conwid1, 3]), rows[-1], np.ones([h, conwid2, 3])])
    result = rows[0]
    for i in range(1, len(rows)):
        result = np.vstack((result, np.ones([blankh, merged_width, 3]), rows[i]))
    return result


def quick_view(name, save=None, column=2):
    """
    quick view multiple plots
    :param name: png file name, support wild card, or list of filename
    :param save: save path for merged image
    :param column: columns of plots
    :return: merged image
    """
    if isinstance(name, list):
        files = name
    elif isinstance(name, str):
        files = glob(name)
    else:
        raise TypeError("name should be string or list")
    imgs = list(map(mpimg.imread, files))
    cutter = white_edge_cutter(imgs[0])
    imgs_cut = list(map(cutter, imgs))
    merged_img = merge_img(imgs_cut, column)
    plt.imshow(merged_img)
    if save:
        mpimg.imsave(save, merged_img)
    return merged_img


def get_cpt(image, cptout):
    """
    extract cptfile from a colorbar image file
    :param image: colorbar image file
    :param cptout: output cpt file path
    :return:
    """
    head = """# COLOR_MODEL = RGB
    """
    img = mpimg.imread(image)
    is_vertical = True
    colorbar = None
    for nrow in range(img.shape[0]):
        if len(set(np.sum(img[nrow, :, :], axis=1))) > 0.2 * img.shape[1]:
            is_vertical = False
            colorbar = img[nrow+5, :, :]
            break
    if is_vertical:
        for ncol in range(img.shape[1]):
            if len(set(np.sum(img[:, ncol, :], axis=1))) > 0.2 * img.shape[2]:
                colorbar = img[:, ncol+5, :]
                break
    if colorbar is None:
        raise ValueError('colorbar must account for most of the input image')

    show = []
    for _ in range(20):
        show.append(colorbar)
    show = np.array(show)
    plt.imshow(show)

    # TODO cut and write
    return colorbar


if __name__ == '__main__':
    colorbar = get_cpt('./testdata/colorbar.png', 'temp')
