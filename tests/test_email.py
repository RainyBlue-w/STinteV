from ..utils import _check
import unittest


class TestEmailCheck(unittest.TestCase):
    
    def test_email(self):
        rv = _check.email('1078497976@q.com')
        self.assertTrue()