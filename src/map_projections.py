from __future__ import division

import matplotlib.pyplot as plt

from azimuthal_equidistant import AzimuthalEquidistant


LANDSHAPE_FILEPATH = "../data/landshape.txt"

CENTER_COORDS = (0.0, 0.0)


def main():
    # get coordinates
    landshape_coords = get_landshape_coords(LANDSHAPE_FILEPATH, lat_range=(-90, 90))

    # Azimuth equidistant projection
    azimuth = AzimuthalEquidistant(CENTER_COORDS)
    azimuth_coords = [azimuth.project(coord) for coord in landshape_coords]
    plt.plot(
        [x for (x, y) in azimuth_coords],
        [y for (x, y) in azimuth_coords],
        '.', markersize=1.0
    )
    plt.show()


def get_landshape_coords(landshape_filepath, lat_range=(-90, 90), long_range=(-90, 90)):
    landshape_coords = []
    with open(landshape_filepath) as file:
        lat_min, lat_max = 0, 0
        long_min, long_max = 0, 0
        for line in file:
            # landshape.txt stores coordinates in (longtitude, latitude) order
            long, lat = line.strip().split(" ")
            lat, long = float(lat), float(long)
            lat_min, lat_max = min(lat_min, lat), max(lat_max, lat)
            long_min, long_max = min(long_min, long), max(long_max, long)
            if lat_range[0] <= lat <= lat_range[1]:
                landshape_coords.append((lat, long))
    print "lat range", (lat_min, lat_max)
    print "long range", (long_min, long_max)
    return landshape_coords


if __name__ == "__main__":
    main()