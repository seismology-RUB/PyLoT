import unittest
from pylot.core.pick.autopick import PickingResults

class TestPickingResults(unittest.TestCase):
    def setUp(self):
        self.pr = PickingResults()

    def test_dot_member_creation(self):
        self.pr.x = 0
        self.assertEqual(self.pr.x, 0)
        self.pr.x += 42
        self.assertEqual(self.pr.x, 42)

    def test_dot_builtin_member(self):
        self.assertEqual(self.pr.Pflag, 0)
        self.pr.Pflag = 99
        self.assertEqual(self.pr.Pflag, 99)

    def test_key_access(self):
        self.pr['y'] = 11
        self.assertEqual(self.pr['y'], 11)

    def test_builtin_fields(self):
        self.assertEqual(self.pr.Pflag, 0)

    def test_missing_attribute(self):
        # accessing a missing attribute in a dictionary throws a KeyError
        with self.assertRaises(KeyError):
            self.pr['keydoesnotexist']

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
        pflag = self.pr.get('Pflag', -1)
        self.assertEqual(pflag, 0)
        self.assertNotEqual(pflag, -1)