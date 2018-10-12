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


def _view_row(row):
    show = []
    for _ in range(20):
        show.append(row)
    show = np.array(show)
    plt.imshow(show)


def get_cpt(image, cptout):
    """
    extract cptfile from a colorbar image file
    assumption: this image should have white background and black rectangle border
    :param image: colorbar image file
    :param cptout: output cpt file path
    :return:
    """
    head = """# COLOR_MODEL = RGB
    """
    border_color = np.array([0, 0, 0])  # black border

    def check_border(line):
        diff = abs(line[:, :3] - border_color)
        is_black = np.mean(diff, axis=1) < 0.5
        is_continuous = np.array([False for _ in range(len(is_black))])
        for j in range(len(line)-1):
            if np.mean(line[j] - line[j+1]) < 0.1:
                is_continuous[j] = True
        is_border = is_black & is_continuous
        max_len, beg, end = 0, 0, 0
        cur_beg = None
        for j in range(len(is_border)):
            if is_border[j]:
                if cur_beg is None:
                    cur_beg = j
                if (j == len(is_border) - 1) or (j < len(is_border) - 1 and (not is_border[j+1])):
                    if j - cur_beg > max_len:
                        beg, end = cur_beg, j
                        max_len = end - beg
                    cur_beg = None
        return [max_len, beg, end]

    img = mpimg.imread(image)
    result = []
    for nrow in range(img.shape[0]):
        result.append(check_border(img[nrow, :, :]))
    return result
    is_vertical = True
    colorbar = None
    for nrow in range(img.shape[0]):
        if len(set(np.sum(img[nrow, :, :], axis=1))) > 0.3 * img.shape[1]:
            print(nrow)
            is_vertical = False
            colorbar = img[nrow+5, :, :]
            break
    if is_vertical:
        for ncol in range(img.shape[1]):
            if len(set(np.sum(img[:, ncol, :], axis=1))) > 0.3 * img.shape[0]:
                colorbar = img[:, ncol+5, :]
                break
    if colorbar is None:
        raise ValueError('colorbar must account for most of the input image')
    prev, cur = None, None
    for i in range(1, len(colorbar)):  # find edge of colorbar
        if abs(sum(colorbar[i, :3]) - sum(colorbar[i-1, :3])) > 2.5:
            if prev is None:
                prev = i
            else:
                if cur is None:
                    cur = i
                else:
                    prev = cur
                    cur = i
                if cur - prev > 0.3 * len(colorbar):
                    print(prev, cur)
                    colorbar = colorbar[prev+1:cur-1]
                    break
    # TODO need more test and output cpt
    _view_row(colorbar)
    return colorbar


if __name__ == '__main__':
    colorbar = get_cpt('./testdata/colorbar.png', 'temp')
