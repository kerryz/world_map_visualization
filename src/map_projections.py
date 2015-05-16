# -*- coding: utf-8 -*-
from __future__ import division

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Circle
# from matplotlib.colors import hex2color

from projection_formulae import (
    AzimuthalEquidistant, SansonFlamsteed, Mercator,
    get_region_projection_area
)


LANDSHAPE_FILEPATH = "../data/landshape.txt"
REGIONS_FILEPATH = "../data/regions.txt"
REGIONS_DIRPATH = "../data/regions/"
CITIES_FILEPATH = "../data/cities.txt"

# CITY_SIZE_FACTOR = 2*1e-7
CITY_SIZE_FACTOR = 2*1e-9
ECONOMIC_LEVELS = 7

CENTER_COORDS = (0.0, 0.0)  # latitude, longitude
LAT_RANGE = (10, 80)  # latitude range
ALL_LAT_RANGE = (-90, 90)
ALL_LONG_RANGE = (-180, 180)


def main():
    # get (latitude, longitude) coordinates
    landshape_coords = get_landshape_coords(LANDSHAPE_FILEPATH, lat_range=LAT_RANGE)

    # Azimuth equidistant projection
    azimuth = AzimuthalEquidistant(CENTER_COORDS)
    # plot_landshape(azimuth, landshape_coords, "Azimuth Equidistant Projection")

    # Sanson-Flamsteed, a.k.a. sinusoidal projection
    sanson_flam = SansonFlamsteed(CENTER_COORDS)  # only requires the longitude, the meridian
    # plot_landshape(sanson_flam, landshape_coords, "Sanson-Flamsteed Projection")

    # Mercator projection
    mercator = Mercator(CENTER_COORDS)
    # Mercator projects infinite y-values at the poles, remove these values
    # mercator_landshape_coords = [(lat, long) for (lat, long) in landshape_coords if -85 < lat < 85]
    # fig, ax = plot_landshape(mercator, mercator_landshape_coords, "Mercator projection")
    fig, ax = plot_landshape(mercator, landshape_coords, "Mercator projection")

    region_coords = get_landshape_coords("../data/regions/115/landshape.txt")
    region_parts = get_landparts("../data/regions/115/landparts.txt")
    # print get_region_projection_area(mercator, landshape_coords, region_coords, region_parts)
    plot_landparts(ax, mercator, region_coords, region_parts, "green")
    print get_region_projection_area((fig, ax), mercator, region_coords, region_parts)


    # fig, ax = plt.subplots()
    # plot_regions((fig, ax), sanson_flam, REGIONS_FILEPATH, REGIONS_DIRPATH, lat_range=LAT_RANGE)
    # plot_cities((fig, ax), sanson_flam, CITIES_FILEPATH, lat_range=LAT_RANGE)
    # ax.set_ylim([-2, 3])



    # show plots
    plt.show()


def plot_cities((fig, ax), projection, cities_filepath,
                lat_range=ALL_LONG_RANGE, long_range=ALL_LONG_RANGE):
    cities_info = []
    lat_min, lat_max = lat_range
    long_min, long_max = long_range
    with open(cities_filepath) as file:
        for line in file:
            city, size, long, lat = line.strip().split("|")
            size, long, lat = float(size), float(long), float(lat)
            if lat_min <= lat <= lat_max and long_min <= long <= long_max:
                cities_info.append((city, size, (lat, long)))
    color = "#0551BE"
    for _, size, coord in cities_info:
        x, y = projection.project(coord)
        ax.add_artist(Circle(
            xy=(x, y), radius=size*CITY_SIZE_FACTOR, alpha=0.5,
            facecolor=color, edgecolor='none'
        ))
        # ax.plot(x, y, 'o', color=color, markersize=size*CITY_SIZE_FACTOR, markeredgecolor=None)


def plot_regions((fig, ax), projection, regions_filepath, regions_dirpath,
                 lat_range=ALL_LAT_RANGE, long_range=ALL_LONG_RANGE):
    # get region economic prosperity levels
    econ_levels = []
    with open(regions_filepath) as file:
        for line in file:
            _, econ_level = line.strip().split("|")
            econ_levels.append(int(econ_level))

    # setup coloring
    ax.set_axis_bgcolor("#BFECFF")  # ocean blue
    # color map
    econ_cmap = colors.LinearSegmentedColormap.from_list(
        'econ_colors', [colors.hex2color("#DDFF28"), colors.hex2color("#D70A0A")])
    sm = plt.cm.ScalarMappable(cmap=econ_cmap, norm=plt.Normalize(vmin=1, vmax=7))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cbar = fig.colorbar(sm, ticks=[7, 6, 5, 4, 3, 2, 1])
    cbar.set_label("Economic level", rotation=270)
    cbar.ax.get_yaxis().labelpad = 15

    # optimization
    norm = 1/ECONOMIC_LEVELS
    # walk through regions directory and plot all the regions
    _, regions_dirs, _ = next(os.walk(regions_dirpath))
    regions_dirs = [d for d in regions_dirs if d[0] != "."]  # remove hidden subdirectories
    for reg_dir in regions_dirs:
        root, _, files = next(os.walk(os.path.join(regions_dirpath, reg_dir)))
        files = [f for f in files if f[0] != "."]  # remove hidden files
        for file in files:
            if file.startswith("landshape"):
                landshape_coords = get_landshape_coords(os.path.join(root, file))
            elif file.startswith("landparts"):
                landparts = get_landparts(os.path.join(root, file))

        if lat_range != ALL_LAT_RANGE or long_range != ALL_LONG_RANGE:
            landshape_coords, landparts = filter_landshape_landparts(
                landshape_coords, landparts, lat_range, long_range)
        # if everything was filtered out, continue to next region
        if len(landshape_coords) <= 2:
            continue

        color = econ_cmap(econ_levels[int(reg_dir)] * norm)
        color = colors.rgb2hex(color)
        plot_landparts(ax, projection, landshape_coords, landparts, color)


def plot_landparts(ax, projection, landshape_coords, landparts, color):
    projection_coords = [projection.project(coord) for coord in landshape_coords]

    """
    # TEST
    # Use below test to only draw landpart nr j to k
    j = 1
    k = 3
    projection_coords = projection_coords[landparts[j]:landparts[k]]
    landparts = landparts[j:k]
    print landparts
    norm = landparts[0]
    landparts = [x - norm for x in landparts]
    """

    xs = [x for (x, y) in projection_coords]
    ys = [y for (x, y) in projection_coords]

    n = len(landparts)
    # args = []

    for i in xrange(n):
        start_i = landparts[i]
        if i == n-1:
            # last element
            ax.fill(xs[start_i:], ys[start_i:], color=color, linewidth=0.04, edgecolor='white')
            # args += [xs[start_i:], ys[start_i:], color]

            # See start and end points of this landpart
            # ax.plot(xs[start_i], ys[start_i], 'go', markersize=10, alpha=0.8)
            # ax.plot(xs[-1], ys[-1], 'ro', markersize=10, alpha=0.8)
            break
        end_i = landparts[i+1]

        # See start and end points of this landpart
        # ax.plot(xs[start_i], ys[start_i], 'go', markersize=10, alpha=0.8)
        # ax.plot(xs[end_i-1], ys[end_i-1], 'ro', markersize=10, alpha=0.8)

        # args += [xs[start_i:end_i], ys[start_i:end_i], color]
        ax.fill(xs[start_i:end_i], ys[start_i:end_i],
                color=color, linewidth=0.04, edgecolor='white')

    # ax.fill(*args)


def plot_landshape(projection, landshape_coords, title):
    """
    Parameters
    ----------
    projection : Projection
        object that should inheret the Projection abstract class
    landshape_coords : [(float, float)]
        [(latitude, longitude)]
    title : string
        the title of the map

    Returns
    -------
    fig, ax
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
    return fig, ax


def get_landshape_coords(landshape_filepath, lat_range=ALL_LAT_RANGE, long_range=ALL_LONG_RANGE):
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


def get_landparts(landparts_filepath):
    landparts = []
    with open(landparts_filepath) as file:
        for line in file:
            landparts.append(int(line.strip()))
    return landparts


def filter_landshape_landparts(
        landshape, landparts,
        lat_range=ALL_LAT_RANGE, long_range=ALL_LONG_RANGE):
    """
    Filters the landshape to only include coordinates within `lat_range` and `long_range`.
    Updates (a copy of) `landparts` to correspond to filtered `landshape`.

    Parameters
    ----------
    landshape : [(float, float)]
        [(latitude, longitude)]
    landparts : numpy Jx1 array
        specifies which coordinates in `landshape` are connected
        Example:
            if landparts = [0, 40, 63]
            then coordinates landshape[0:40] are connected (is one land part)
            landshape[40:63] is a landpart and landshape[63:J] is a landpart
    lat_range : (float, float)
        lat_range = (latitude_min, latitude_max)
    long_range : (float, float)
        long_range = (longitude_min, longitude_max)

    Returns
    -------
    (new_landshape, landparts)
        updated, copied arrays
    """
    landparts = np.array(landparts)

    lat_min, lat_max = lat_range
    long_min, long_max = long_range

    new_landshape = []
    removed = 0
    last_j = 0
    for i, (lat, long) in enumerate(landshape):
        if not (lat_min <= lat <= lat_max and long_min <= long <= long_max):
            # in each iteration, landparts refers to the new filtered landshape
            # therefore use i-removed
            landparts, last_j = update_landparts(i - removed, landparts, last_j)
            removed += 1
        else:
            new_landshape.append((lat, long))
    # edge case: if all of the last land part is removed
    # the last element in `landparts` has to be removed
    if landparts[-1] >= len(new_landshape):
        landparts = landparts[:-1]
    return new_landshape, landparts


def update_landparts(i, landparts, last_j):
    """
    Parameters
    ----------
    i : int
        index of the coordinate in `landshape` to be removed
    landparts : [int]
        specifies which coordinates in `landshape` are connected
    last_j : int
        index of `landparts` from where we should start looking
        for landparts[j] < i < landparts[j+1]
        just for optimization, so we don't have to look through
        all of `landparts` every time
    """
    n = len(landparts)
    updated = False
    j = last_j - 1
    while not updated:
        j += 1
        if j == n - 1:
            # the coordinate to be removed is in the last landpart
            # `landparts` doesn't have to be changed: do nothing
            updated = True
        elif landparts[j] <= i < landparts[j+1]:
            # the coordinate to be removed is in the landpart that satisfies
            # landparts[j] <= i < landparts[j+1]
            # since we're removing a coordinate from this landpart,
            # all landparts that are listed later will have their index -= 1
            landparts[j+1:] -= 1
            if landparts[j+1] == landparts[j]:
                # after performing this update, there is a possibility that
                # landparts[j+1] == landparts[j]. In this case, we have to remove one of them
                landparts = np.delete(landparts, j+1)
            updated = True
    return landparts, j


if __name__ == "__main__":
    main()
