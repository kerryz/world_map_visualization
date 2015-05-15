# -*- coding: utf-8 -*-
from __future__ import division
from math import pi, log, radians, cos, sin, tan, acos


CENTER_COORD = (45.0, 0.0)  # (latitude, longitude)


def main():
    azimuth = AzimuthalEquidistant(CENTER_COORD)
    print azimuth.project((45.0, 10.0))


class Projection(object):
    """
    Abstract projection class
    """
    def project(self, (lat, long)):
        """
        Projects the latitude and longitude to 2d (x, y) map coordinates.

        Returns
        -------
        (x, y) : (float, float)
            the projected 2d map coordinates
        """
        raise NotImplementedError()


class AzimuthalEquidistant(Projection):
    """
    Azimuthal equidistant projection.
    Formulas taken from
    http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
    """
    def __init__(self, (lat0, long0)):
        """
        Parameters
        ----------
        (lat0, long0) : (float, float)
            latitude and longitude (degrees) coordinates of the center of the projections,
            called phi_1 and lambda_0 in the Wolfram link
        """
        self.p1 = radians(lat0)  # phi_1
        self.l0 = radians(long0)  # lambda_0
        # optimization: storing these values because they are used
        # when calculating c, the angular distance from the center
        self.cos_p1 = cos(self.p1)
        self.sin_p1 = sin(self.p1)

    def project(self, (lat, long)):
        # if same coordinates as center point, then c = 0
        # and k = c / sin(c) will perform division by 0
        if (lat, long) == (self.p1, self.l0):
            return (0.0, 0.0)

        # convert to radians and rename variables to match Wolfram link
        p = radians(lat)  # phi
        l = radians(long)  # lambda
        # optimization, store repeated calculations
        sin_p = sin(p)
        cos_p = cos(p)
        cosp_cosll0 = cos_p * cos(l - self.l0)

        # calculate c, the angular distance from the center
        cos_c = self.sin_p1 * sin_p + self.cos_p1 * cosp_cosll0
        c = acos(cos_c)
        if c == 0:
            return (0.0, 0.0)
        # k'
        try:
            k = c / sin(c)
        except ZeroDivisionError:
            print lat, long

        # calculate Azimuthal Equidistant projection x, y coordinates
        x = k * cos_p * sin(l - self.l0)
        y = k * (self.cos_p1 * sin_p - self.sin_p1 * cosp_cosll0)

        return (x, y)


class SansonFlamsteed(Projection):
    """
    The Sanson–Flamsteed or the Mercator equal-area projection, or
    the sinusoidal projection is a pseudocylindrical equal-area map projection.

    date: 2015-05-13
    http://en.wikipedia.org/wiki/Sinusoidal_projection
    """
    def __init__(self, (lat0, long0)):
        """
        Parameters
        ----------
        (lat0, long0) : (float, float)
            latitude and longitude (degrees) coordinates of the center of the projections,
            called phi_1 and lambda_0 in the Wolfram link
        """
        # lambda_0, same notation as the Wikipedia page
        self.l0 = long0

    def project(self, (lat, long)):
        # convert to radians and rename variables to match Wikipedia page
        p = radians(lat)
        l = radians(long)

        x = (l - self.l0) * cos(p)
        y = p

        return (x, y)


class Mercator(Projection):
    """
    Mercator projection.
    Formulas taken from
    http://mathworld.wolfram.com/MercatorProjection.html
    """
    def __init__(self, (lat0, long0)):
        """
        Parameters
        ----------
        (lat0, long0) : (float, float)
            latitude and longitude (degrees) coordinates of the center of the projections,
            called phi_1 and lambda_0 in the Wolfram link
        """
        self.l0 = long0

    def project(self, (lat, long)):
        # TODO: be able to move central coord

        # y-values go to infinity at the poles, fix
        if lat <= -89:
            lat = -89
        if lat >= 89:
            lat = 89

        # convert to radians and rename variables to match Wolfram page
        p = radians(lat)
        l = radians(long)

        x = l - self.l0
        temp = tan(0.25*pi + 0.5*p)
        y = log(temp)
        return (x, y)


if __name__ == "__main__":
    main()
