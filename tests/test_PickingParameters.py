import unittest
from pylot.core.pick.autopick import PickingParameters

class TestPickingParameters(unittest.TestCase):

    def setUp(self):
        self.simple_dict = {'a': 3, 'b': 14}
        self.nested_dict = {'a': self.simple_dict, 'b': self.simple_dict}

    def assertParameterEquality(self, dic, instance):
        """Test wether all parameters given in dic are found in instance"""
        for key, value in dic.items():
            self.assertEqual(value, getattr(instance, key))

    def test_add_params_from_dict_simple(self):
        pickparam = PickingParameters()
        pickparam.add_params_from_dict(self.simple_dict)
        self.assertParameterEquality(self.simple_dict, pickparam)

    def test_add_params_from_dict_nested(self):
        pickparam = PickingParameters()
        pickparam.add_params_from_dict(self.nested_dict)
        self.assertParameterEquality(self.nested_dict, pickparam)

    def test_init(self):
        pickparam = PickingParameters(self.simple_dict)
        self.assertParameterEquality(self.simple_dict, pickparam)


if __name__ == '__main__':
    unittest.main()