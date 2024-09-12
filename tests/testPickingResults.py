import unittest

from pylot.core.pick.autopick import PickingResults


class TestPickingResults(unittest.TestCase):

    def setUp(self):
        self.pr = PickingResults()

    def test_non_existing_key_dot_access(self):
        """Accessing an attribute in the class that wasnt added to the dict should give a AttributeError"""
        with self.assertRaises(AttributeError):
            self.pr.doesntexist

    def test_non_existing_key_dict_access(self):
        """Accessing a missing attribute in a dictionary throws a KeyError"""
        with self.assertRaises(KeyError):
            self.pr['keydoesnotexist']

    def test_dot_member_creation(self):
        self.pr.x = 0
        self.assertEqual(self.pr.x, 0)
        self.pr.x += 42
        self.assertEqual(self.pr.x, 42)

    def test_dot_builtin_member(self):
        self.assertEqual(self.pr.weight, 4)
        self.pr.weight = 99
        self.assertEqual(self.pr.weight, 99)

    def test_key_access(self):
        self.pr['y'] = 11
        self.assertEqual(self.pr['y'], 11)

    def test_builtin_fields(self):
        self.assertEqual(self.pr['weight'], 4)

    def test_in(self):
        self.assertFalse('keydoesnotexist' in self.pr)
        self.pr['k'] = 0
        self.assertTrue('k' in self.pr)

    def test_keys_function(self):
        a = 99
        self.pr.newkey = a
        self.assertIn(a, self.pr.values())
        self.assertIn('newkey', self.pr.keys())

    def test_len_and_clear(self):
        self.pr.clear()
        self.assertEqual(len(self.pr), 0)
        self.pr.a = 6
        self.pr['b'] = 9
        self.assertEqual(len(self.pr), 2)

    def test_get_default(self):
        self.assertEqual(self.pr.get('keynotexisting', 42), 42)
        weight = self.pr.get('weight', -1)
        self.assertEqual(weight, 4)
        self.assertNotEqual(weight, -1)

    def test_dunder_attributes(self):
        """Storing Pythons special dunder method in a dictionary is valid and should not override the instances dunder
        methods"""
        prev_len = len(self.pr)
        try:
            self.pr['__len__'] = None
        except Exception:
            self.fail("test_dunder_attributes failed to add a dunder attribute to the dictionary keys")
        try:
            curr_len = len(self.pr)
        except Exception:
            self.fail("test_dunder_attributes overwrote an instance internal dunder method")
        self.assertEqual(prev_len + 1, curr_len)  # +1 for the added __len__ key/value-pair

        self.pr.__len__ = 42

        self.assertEqual(42, self.pr['__len__'])
        self.assertEqual(prev_len + 1, curr_len, msg="__len__ was overwritten")
