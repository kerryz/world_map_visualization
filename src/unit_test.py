from __future__ import division

import unittest
import numpy as np

from map_projections import filter_landshape_landparts


class TestLandshapeLandparts(unittest.TestCase):

    def test_filter_beginning(self):
        landshape = [
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(10.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        correct_landparts = np.array([0, 2, 5])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_middle(self):
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(10.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        correct_landparts = np.array([0, 3, 5])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_end(self):
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0),
            (10.0, 90.0), (11.0, 90.0), (9.0, 90.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(10.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        correct_landparts = np.array([0, 3, 6])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_multiple_beginning(self):
        landshape = [
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(11.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        correct_landparts = np.array([0, 1, 4])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_multiple_middle(self):
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(11.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        correct_landparts = np.array([0, 3, 4])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_multiple_end(self):
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0),
            (11.0, 90.0), (9.0, 90.0), (10.0, 90.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(11.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        correct_landparts = np.array([0, 3, 6])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_multiple_everywhere(self):
        landshape = [
            (9.0, 0.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (9.0, 1.0), (9.0, 2.0),
            (11.0, 90.0), (9.0, 3.0), (9.0, 4.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(11.0, 80.0)
        )

        landshape.remove((9.0, 0.0))
        landshape.remove((9.0, 1.0))
        landshape.remove((9.0, 2.0))
        landshape.remove((9.0, 3.0))
        landshape.remove((9.0, 4.0))
        correct_landparts = np.array([0, 2, 3])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_whole_landpart_beginning(self):
        """
        When a whole landpart is removed, `landparts` should remove that element
        """
        landshape = [
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(20.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        landshape.remove((11.0, 90.0))
        correct_landparts = np.array([0, 3])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_whole_landpart_middle(self):
        """
        When a whole landpart is removed, `landparts` should remove that element
        """
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(20.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        landshape.remove((11.0, 90.0))
        correct_landparts = np.array([0, 3])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)

    def test_filter_whole_landpart_end(self):
        """
        When a whole landpart is removed, `landparts` should remove that element
        """
        landshape = [
            (30.0, 130.0), (31.0, 131.0), (32.0, 132.0),
            (40.0, 140.0), (41.0, 141.0), (42.0, 142.0),
            (9.0, 90.0), (10.0, 90.0), (11.0, 90.0)
        ]
        landparts = np.array([0, 3, 6])

        new_landshape, new_landparts = filter_landshape_landparts(
            landshape, landparts, lat_range=(20.0, 80.0)
        )

        landshape.remove((9.0, 90.0))
        landshape.remove((10.0, 90.0))
        landshape.remove((11.0, 90.0))
        correct_landparts = np.array([0, 3])

        np.testing.assert_array_equal(new_landshape, landshape)
        np.testing.assert_array_equal(new_landparts, correct_landparts)


if __name__ == "__main__":
    unittest.main()
