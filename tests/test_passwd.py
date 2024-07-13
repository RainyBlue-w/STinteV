# coding: utf-8

""" Tests for check. """

from ..utils import _check
import unittest


class TestCheck(unittest.TestCase):

    def test_length(self):
        rv = _check.password('1234567')
        self.assertTrue(repr(rv) == 'terrible')
        self.assertTrue(rv.message == 'Too short')
    
    def test_regular(self):
        rv = _check.password('qwertyuio')
        self.assertTrue(repr(rv) == 'simple')
        self.assertTrue(rv.message == 'Regular password')

    def test_by_step(self):
        rv = _check.password('abcdefghi')
        self.assertTrue(repr(rv) == 'simple')
        self.assertTrue(rv.message == 'Regular password')
        
    def test_common(self):
        rv = _check.password('password')
        self.assertTrue(repr(rv) == 'simple')
        self.assertTrue(rv.message == 'Common password')

    def test_medium(self):
        rv = _check.password('tdnwh001')
        self.assertTrue(repr(rv) == 'medium')
        self.assertTrue(rv.message == 'Risky password')

    def test_strong(self):
        rv = _check.password('tdnWwh001')
        self.assertTrue(repr(rv) == 'strong')
        self.assertTrue(rv.message == 'Strong password')

    def test_very_strong(self):
        rv = _check.password('tdnWwh001@')
        self.assertTrue(repr(rv) == 'strong')
        self.assertTrue(rv.message == 'Very strong password')

if __name__ == '__main__':
    unittest.main()
    rv = _check.password('tdnwh01')
    print(rv.strength)
    
    