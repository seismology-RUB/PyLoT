import unittest
from pylot.core.pick.utils import get_quality_class


class TestQualityClassFromUncertainty(unittest.TestCase):
    """
    Test function that assigns a quality value [0...4] to a pick uncertainty.
    The pick uncertainty is compared to the error classes.
    A pick uncertainty that is below the first error class is assigned the best quality, quality 0.
    A pick uncertainty that is above the first error class but below the second is assigned quality 1 and so on.
    A pick uncertainty that is larger than the biggest error class is assigned quality 4.
    The upper border of a quality class is inclusive, the lower border exclusive. Meaning if a value is exactly on the
    border between two classes, it is assigned into the higher quality class (represented by the lower number).
    """

    def setUp(self):
        # entries hold upper/lower bound of error classes
        self.error_classes = [float(x) for x in range(1, 9, 2)]
        # [1.0, 3.0, 5.0, 7.0]

    def test_out_of_lower_bound(self):
        # Error out of lower bound of classes
        self.assertEqual(0, get_quality_class(0.5, self.error_classes))

    def test_out_of_upper_bound(self):
        # Error out of upper bound of error classes
        self.assertEqual(4, get_quality_class(14.7, self.error_classes))

    def test_on_lower_border(self):
        # Error exactly on lower bound
        self.assertEqual(0, get_quality_class(1., self.error_classes))

    def test_on_upper_border(self):
        # Error exactly on upper bound
        self.assertEqual(3, get_quality_class(7., self.error_classes))

    def test_on_middle_border_inclusive(self):
        # Error exactly between two classes, since lower bound is exclusive and upper bound is inclusive it should
        # fall into the class with better quality
        self.assertEqual(1, get_quality_class(3., self.error_classes))
        self.assertNotEqual(2, get_quality_class(3., self.error_classes))

    def test_in_class1(self):
        # Error exactly in class 1
        self.assertEqual(1, get_quality_class(1.5, self.error_classes))

    def test_in_class2(self):
        # Error exactly in class 2
        self.assertEqual(2, get_quality_class(3.5, self.error_classes))

    def test_in_class3(self):
        # Error exactly in class 3
        self.assertEqual(3, get_quality_class(5.6, self.error_classes))

if __name__ == '__main__':
    unittest.main()
