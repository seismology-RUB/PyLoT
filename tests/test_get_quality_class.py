import unittest
from pylot.core.pick.utils import getQualityFromUncertainty as get_quality_class


class TestQualityClassFromUncertainty(unittest.TestCase):
    """
    Test function that assigns a quality value [0...4] to a pick uncertainty.
    """

    def setUp(self):
        # entries hold upper/lower bound of error classes
        self.error_classes = [float(x) for x in range(1, 9, 2)]
        # [1.0, 3.0, 5.0, 7.0]

    def test_out_of_lower_bound(self):
        # Error out of lower bound of classes
        self.assertEqual(get_quality_class(0.5, self.error_classes), 0)

    def test_out_of_upper_bound(self):
        # Error out of upper bound of error classes
        self.assertEqual(get_quality_class(14.7, self.error_classes), 4)

    def test_on_lower_border(self):
        # Error exactly on lower bound
        self.assertEqual(get_quality_class(1., self.error_classes), 0)

    def test_on_upper_border(self):
        # Error exactly on upper bound
        self.assertEqual(get_quality_class(7., self.error_classes), 3)

    def test_on_middle_border_inclusive(self):
        # Error exactly between two classes, since lower bound is exclusive and upper bound is inclusive it should
        # fall into the class with better quality
        self.assertEqual(get_quality_class(3., self.error_classes), 1)
        self.assertNotEqual(get_quality_class(3., self.error_classes), 2)

    def test_in_class0(self):
        # Error exactly in class 1
        self.assertEqual(get_quality_class(1.5, self.error_classes), 1)

    def test_in_class2(self):
        # Error exactly in class 2
        self.assertEqual(get_quality_class(3.5, self.error_classes), 2)

    def test_in_class3(self):
        # Error exactly in class 3
        self.assertEqual(get_quality_class(5.6, self.error_classes), 3)

if __name__ == '__main__':
    unittest.main()
