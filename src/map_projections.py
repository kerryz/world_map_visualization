# -*- coding: utf-8 -*-
from __future__ import division

import sys
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pow, sin, cos
from matplotlib import colors
from matplotlib.patches import Circle
# from matplotlib.colors import hex2color

from projection_formulae import (
    AzimuthalEquidistant, SansonFlamsteed, Mercator,
    get_region_projection_area
)

# can be changed if path to data directory is passed through command line
# see bottom of this page
DATA_DIR_PATH = os.path.join("..", "data")

LANDSHAPE_FILEPATH = os.path.join(DATA_DIR_PATH, "landshape.txt")
REGIONS_FILEPATH = os.path.join(DATA_DIR_PATH, "regions.txt")
CITIES_FILEPATH = os.path.join(DATA_DIR_PATH, "cities.txt")
REGIONS_DIRPATH = os.path.join(DATA_DIR_PATH, "regions")
# -----------------------------------

JAPAN_LANDSHAPE_FILEPATH = os.path.join(REGIONS_DIRPATH, "115", "landshape.txt")
JAPAN_LANDPARTS_FILEPATH = os.path.join(REGIONS_DIRPATH, "115", "landparts.txt")
TAIWAN_LANDSHAPE_FILEPATH = os.path.join(REGIONS_DIRPATH, "231", "landshape.txt")
TAIWAN_LANDPARTS_FILEPATH = os.path.join(REGIONS_DIRPATH, "231", "landparts.txt")

# CITY_SIZE_FACTOR = 2*1e-7
CITY_SIZE_FACTOR = 2*1e-9
ECONOMIC_LEVELS = 7

CENTER_COORDS = (0.0, 0.0)  # latitude, longitude
LAT_RANGE_2 = (10, 80)  # latitude range
ALL_LAT_RANGE = (-90, 90)
ALL_LONG_RANGE = (-180, 180)
# Mercator projects infinite y-values at the poles, don't use those values
LAT_RANGE_MERCATOR = (-80, 85)

BEIJING_COORDS = (39.9167, 116.3833)  # latitude, longitude
LA_COORDS = (34.0500, -118.2500)  # Los Angeles


def main():
    azimuth = AzimuthalEquidistant(CENTER_COORDS)
    sanson_flam = SansonFlamsteed(CENTER_COORDS)
    mercator = Mercator(CENTER_COORDS)
    print """
    Enter an option:

    1: Visualize the whole world
    2: Azimuthal equidistant, Sanson-Flamsteed (sinusoidal), and Mercator projections
    3: Geodesic path between Beijing and Los Angeles
    """

    command = raw_input("Enter option: ").strip()

    # Visualize whole world
    if command == "1":
        fig, ax = plot_regions(
            mercator, REGIONS_FILEPATH, REGIONS_DIRPATH,
            "Whole world visualization", lat_range=LAT_RANGE_MERCATOR)
        # plot_lat_long_grid(
        #     ax, mercator, LAT_RANGE_MERCATOR, ALL_LONG_RANGE, lat_grids=13, long_grids=23)
        plot_cities((fig, ax), mercator, CITIES_FILEPATH, lat_range=LAT_RANGE_MERCATOR)
        plt.show()
    # end option 1

    # Azimuthal equidistant, Sanson-Flamsteed (sinusoidal), and Mercator projections
    elif command == "2":
        japan_region_coords = get_landshape_coords(JAPAN_LANDSHAPE_FILEPATH)
        japan_region_parts = get_landparts(JAPAN_LANDPARTS_FILEPATH)
        taiwan_region_coords = get_landshape_coords(TAIWAN_LANDSHAPE_FILEPATH)
        taiwan_region_parts = get_landparts(TAIWAN_LANDPARTS_FILEPATH)

        # Azimuthal
        print
        print "Calculating Azimuthal equidistant projection..."
        fig, ax = plot_regions(
            azimuth, REGIONS_FILEPATH, REGIONS_DIRPATH,
            "Azimuthal equidistant projection", lat_range=LAT_RANGE_2)
        plot_lat_long_grid(ax, azimuth, LAT_RANGE_2, ALL_LONG_RANGE)
        plot_cities((fig, ax), azimuth, CITIES_FILEPATH, lat_range=LAT_RANGE_2)
        # ax.set_ylim([0.0, 2.5])
        japan_pixels_a = get_region_projection_area(
            (fig, ax), azimuth, japan_region_coords, japan_region_parts)
        taiwan_pixels_a = get_region_projection_area(
            (fig, ax), azimuth, taiwan_region_coords, taiwan_region_parts)

        # Sanson-Flamsteed (sinusoidal)
        print "Calculating sinusoidal projection..."
        fig, ax = plot_regions(
            sanson_flam, REGIONS_FILEPATH, REGIONS_DIRPATH,
            "Sinusoidal projection", lat_range=LAT_RANGE_2)
        plot_lat_long_grid(ax, sanson_flam, LAT_RANGE_2, ALL_LONG_RANGE, long_grids=9)
        ax.set_xlim([-3.5, 3.5])
        plot_cities((fig, ax), sanson_flam, CITIES_FILEPATH, lat_range=LAT_RANGE_2)
        japan_pixels_s = get_region_projection_area(
            (fig, ax), sanson_flam, japan_region_coords, japan_region_parts)
        taiwan_pixels_s = get_region_projection_area(
            (fig, ax), sanson_flam, taiwan_region_coords, taiwan_region_parts)

        # Mercator
        print "Calculating Mercator projection..."
        fig, ax = plot_regions(
            mercator, REGIONS_FILEPATH, REGIONS_DIRPATH,
            "Mercator projection", lat_range=LAT_RANGE_2)
        plot_lat_long_grid(ax, mercator, LAT_RANGE_2, ALL_LONG_RANGE)
        plot_cities((fig, ax), mercator, CITIES_FILEPATH, lat_range=LAT_RANGE_2)
        japan_pixels_m = get_region_projection_area(
            (fig, ax), mercator, japan_region_coords, japan_region_parts)
        taiwan_pixels_m = get_region_projection_area(
            (fig, ax), mercator, taiwan_region_coords, taiwan_region_parts)

        print
        print "----------------------------------"
        print "         Number of pixels"
        print "----------------------------------"
        print "                 |   Japan | Taiwan"
        print "Azimuthal        |     {0} | {1}".format(japan_pixels_a, taiwan_pixels_a)
        print "Sanson-Flamsteed |     {0} | {1}".format(japan_pixels_s, taiwan_pixels_s)
        print "Mercator         |     {0} | {1}".format(japan_pixels_m, taiwan_pixels_m)
        print "----------------------------------"

        plt.show()
    # end option 2

    # Geodesic path between Beijing and Los Angeles
    elif command == "3":
        landshape_coords = get_landshape_coords(LANDSHAPE_FILEPATH)
        # Azimuthal, center at (latitude=0, longitude=0)
        fig, ax = plot_landshape(azimuth, landshape_coords, "Azimuthal Equidistant Projection, Center at (0, 0)")
        plot_geodesic(ax, azimuth, LA_COORDS, BEIJING_COORDS)
        # Azimuthal, center at Beijing
        azimuth_bj = AzimuthalEquidistant(BEIJING_COORDS)
        fig, ax = plot_landshape(azimuth_bj, landshape_coords, "Azimuthal Equidistant Projection, Center at Beijing")
        plot_geodesic(ax, azimuth_bj, LA_COORDS, BEIJING_COORDS)
        # Sanson-Flamsteed, a.k.a. sinusoidal projection
        fig, ax = plot_landshape(sanson_flam, landshape_coords, "Sinusoidal projection")
        plot_geodesic(ax, sanson_flam, LA_COORDS, BEIJING_COORDS)
        # Mercator projection
        mercator_landshape_coords = [
            (lat, long) for (lat, long) in landshape_coords
            if LAT_RANGE_MERCATOR[0] < lat < LAT_RANGE_MERCATOR[1]
        ]
        fig, ax = plot_landshape(mercator, mercator_landshape_coords, "Mercator projection")
        plot_geodesic(ax, mercator, LA_COORDS, BEIJING_COORDS)

        plt.tight_layout()
        plt.show()

    else:
        print('Invalid Command.')


def plot_regions(projection, regions_filepath, regions_dirpath, title,
                 lat_range=ALL_LAT_RANGE, long_range=ALL_LONG_RANGE):
    """
    Returns
    -------
    fig, ax
    """
    fig, ax = plt.subplots()
    ax.set_title(title)
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

    fig.tight_layout()
    return fig, ax


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
    fig.tight_layout()
    return fig, ax


def plot_geodesic(ax, projection, start, end):
    """
    Assumes geodesic path will not go around the poles
    """

    coords = get_geodesic_line(start, end)

    projected_coords = [projection.project(coord) for coord in coords]
    xs = [x for (x, y) in projected_coords]
    ys = [y for (x, y) in projected_coords]

    # segment coordinates where they change signs (pos to neg, or neg to pos)
    # otherwise will connect lines from opposite sides of map
    longs = [long for (lat, long) in coords]
    longs_change = 0  # index where longs has changed sign
    for i, long in enumerate(longs[1:]):
        if long * longs[i] < 0:
            # only if the current and the previous values are of opposite signs
            # will the multiplication result in negative value
            longs_change = i + 1
            break

    if longs_change:
        ax.plot(xs[:longs_change], ys[:longs_change], 'b-')
        ax.plot(xs[longs_change:], ys[longs_change:], 'b-')
    else:
        ax.plot(xs, ys, 'b-')


def get_geodesic_line(start, end, nr_segments=200):
    (lat0, long0) = start
    (lat1, long1) = end

    lat0_rad, long0_rad = math.radians(lat0), math.radians(long0)
    lat1_rad, long1_rad = math.radians(lat1), math.radians(long1)

    frac_delta = (1.0/(nr_segments - 1))  # fractional increment

    # distance radius
    dist_rad = 2 * math.asin(
        sqrt(
            pow((sin((lat0_rad - lat1_rad) / 2)), 2)
            + cos(lat0_rad) * cos(lat1_rad)
            * pow((sin((long0_rad-long1_rad)/2)), 2)
        )
    )

    # start coordinates
    lats = [lat0]
    longs = [long0]

    # coordinates in between start and end
    f = frac_delta
    for _ in xrange(1, nr_segments - 1):
        # f is a fraction along the path from start to end
        A = sin((1-f)*dist_rad) / sin(dist_rad)
        B = sin(f*dist_rad) / sin(dist_rad)
        x = A*cos(lat0_rad)*cos(long0_rad) + B*cos(lat1_rad)*cos(long1_rad)
        y = A*cos(lat0_rad)*sin(long0_rad) + B*cos(lat1_rad)*sin(long1_rad)
        z = A*sin(lat0_rad) + B*sin(lat1_rad)
        lat_i = math.atan2(z, sqrt(pow(x, 2) + pow(y, 2)))
        long_i = math.atan2(y, x)
        lat_i_deg = math.degrees(lat_i)
        long_i_deg = math.degrees(long_i)
        lats.append(lat_i_deg)
        longs.append(long_i_deg)

        f += frac_delta

    # end coordinates
    lats.append(lat1)
    longs.append(long1)

    return zip(lats, longs)


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


def plot_lat_long_grid(
        ax, projection, lat_range, long_range,
        lat_grids=9, long_grids=9, color="#D1FFFF"):
    (lat_min, lat_max) = lat_range
    (long_min, long_max) = long_range

    lat_grid = np.linspace(lat_min, lat_max, num=lat_grids)
    long_grid = np.linspace(long_min, long_max, num=long_grids)

    lats = np.linspace(lat_min, lat_max, num=300)
    longs = np.linspace(long_min, long_max, num=300)

    for lat in lat_grid:
        projected_coords = [projection.project(coord) for coord in zip([lat]*300, longs)]
        xs = [x for (x, y) in projected_coords]
        ys = [y for (x, y) in projected_coords]
        ax.plot(xs, ys, '-', color=color, linewidth=0.5)

    for long in long_grid:
        projected_coords = [projection.project(coord) for coord in zip(lats, [long]*300)]
        xs = [x for (x, y) in projected_coords]
        ys = [y for (x, y) in projected_coords]
        ax.plot(xs, ys, '-', color=color, linewidth=0.5)


if __name__ == "__main__":

    if len(sys.argv) > 1:
        DATA_DIR_PATH = sys.argv[1]
        LANDSHAPE_FILEPATH = os.path.join(DATA_DIR_PATH, "landshape.txt")
        REGIONS_FILEPATH = os.path.join(DATA_DIR_PATH, "regions.txt")
        CITIES_FILEPATH = os.path.join(DATA_DIR_PATH, "cities.txt")
        REGIONS_DIRPATH = os.path.join(DATA_DIR_PATH, "regions")

    main()
