"""
If using this code for your own research, 
please cite Mergny and Schmidt, 2024, MultIHeaTS
or contact author at cyril.mergny@universite-paris-saclay.fr
"""

from PIL import Image
import numpy as np


def load_image(filename):
    img = Image.open(filename)
    img.load()
    data = np.asarray(img, dtype="int32")
    return data


def convert2alb(pixel):
    """
    Convert pixels to albedo values.
    Min/Max albedo values are based on Rathbun2010.
    Min/Max pixel values are retrieved from the image
    by taking the mean of the 5% min/max.
    """
    amax = 0.8  # Rathbun2010
    amin = 0.4  # Rathbun2010
    pmax = 198.4  # Mean of top 5%
    pmin = 89.4  # Mean of bot 5% resolved.
    return (pixel - pmin) / (pmax - pmin) * (amax - amin) + amin


def downsample_img(omatrix, factor=10):
    """
    Downsample the original matrix by taking the mean of each block of pixels.
    Args:
    - omatrix (numpy.ndarray): Original matrix representing the image.
    - factor (int): Size of the block for downsampling. Default is 10.
    Returns:
    - dmatrix (numpy.ndarray): Downsampled matrix.
    """
    # Get the dimensions of the original matrix
    height, width = omatrix.shape
    new_height = height // factor
    new_width = width // factor
    dmatrix = np.zeros((new_height, new_width))

    for i in range(new_height):
        for j in range(new_width):
            dmatrix[i, j] = np.mean(
                omatrix[i * factor : (i + 1) * factor, j * factor : (j + 1) * factor]
            )

    return dmatrix


def retrieve_albedo(lat, long, albedos):
    """
    Convert position in lat/long into indexes in col/rows.
    Input: latitude, longitude in degrees
    Returns: ilat, ilong indexes
    """
    ny, nx = albedos.shape
    i_lat = int((90 - lat) / 180 * (ny - 1))
    i_long = int((360 - long) / 360 * (nx - 1))
    return albedos[i_lat, i_long]


def import_albedos():
    """
    1) Import photometric image of Europa
    2) Convert pixel values to albedo values
    3) Downsample albedo matrix
    Values not observed by Galileo are set to negative vals, badly coded.
    """
    img = load_image("europa_mosaic.jpg")
    img[img < 15] = -1000  # Values not resolved by Galileo
    albedos = convert2alb(img)
    albedos[albedos < 0] = -1000  # Values not resolved by Galileo
    return downsample_img(albedos, factor=50)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    lat, long = 0, 0
    dmatrix = import_albedos()
    alb = retrieve_albedo(lat, long, dmatrix)
    print(alb)

    amin, amax = 0.4, 0.8
    fig, ax = plt.subplots()
    im = ax.imshow(dmatrix, cmap="gray", vmin=amin, vmax=amax)
    fig.colorbar(im)
    plt.show()
