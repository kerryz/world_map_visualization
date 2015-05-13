# -*- coding: utf-8 -*-
from __future__ import division

import matplotlib.pyplot as plt
# from matplotlib.colors import hex2color

from projection_formulae import AzimuthalEquidistant, SansonFlamsteed


LANDSHAPE_FILEPATH = "../data/landshape.txt"

CENTER_COORDS = (45, 0.0)  # latitude, longitude
LAT_RANGE = (10, 80)  # latitude range


def main():
    # get coordinates
    landshape_coords = get_landshape_coords(LANDSHAPE_FILEPATH, lat_range=LAT_RANGE)

    # Azimuth equidistant projection
    azimuth = AzimuthalEquidistant(CENTER_COORDS)
    plot_map(azimuth, landshape_coords, "Azimuth Equidistant Projection")

    # Sanson-Flamsteed, a.k.a. sinusoidal projection
    sanson_flam = SansonFlamsteed(CENTER_COORDS[1])  # only requires the longitude, the meridian
    plot_map(sanson_flam, landshape_coords, "Sanson-Flamsteed Projection")

    # show plots
    plt.show()


def plot_map(projection, landshape_coords, title):
    """
    Parameters
    ----------
    projection : Projection
        object that should inheret the Projection abstract class
    landshape_coords : [(float, float)]
        [(latitude, longitude)]
    title : string
        the title of the map
    """
    projection_coords = [projection.project(coord) for coord in landshape_coords]
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.plot(
        [x for (x, y) in projection_coords],
        [y for (x, y) in projection_coords],
        'w.', markersize=0.04
    )
    # ax.set_axis_bgcolor((0.9686274509803922, 0.9450980392156862, 0.8745098039215686))
    ax.set_axis_bgcolor('black')


def plot_azimuthal(landshape_coords):
    azimuth = AzimuthalEquidistant(CENTER_COORDS)
    azimuth_coords = [azimuth.project(coord) for coord in landshape_coords]
    fig, ax = plt.subplots()
    ax.set_title("Azimuth Equidistant Projection")
    ax.plot(
        [x for (x, y) in azimuth_coords],
        [y for (x, y) in azimuth_coords],
        'w.', markersize=0.04
    )
    # ax.set_axis_bgcolor((0.9686274509803922, 0.9450980392156862, 0.8745098039215686))
    ax.set_axis_bgcolor('black')


def get_landshape_coords(landshape_filepath, lat_range=(-90, 90), long_range=(-180, 180)):
    """
    Returns
    -------
    landshape_coords : [(float, float)]
        [(latitude, longitude)]
    """
    landshape_coords = []
    with open(landshape_filepath) as file:
        for line in file:
            # landshape.txt stores coordinates in (longitude, latitude) order
            long, lat = line.strip().split(" ")
            lat, long = float(lat), float(long)
            if lat_range[0] <= lat <= lat_range[1] and long_range[0] <= long <= long_range[1]:
                landshape_coords.append((lat, long))
    return landshape_coords


if __name__ == "__main__":
    main()
