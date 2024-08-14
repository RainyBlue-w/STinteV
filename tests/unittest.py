from stintev.utils import _check
import unittest
from stintev.forms import RequesetResetPwdForm, ResetPwdForm


class TestEmailCheck(unittest.TestCase):
    
    def test_email(self):
        rv = _check.email('1078497976@q.com')
        self.assertTrue()
        
class TestResetFormCheck(unittest.TestCase):
    
    def test_reset_form(self):
        form = RequesetResetPwdForm()
        form.email.data = '1078497976@qq.com'
        self.assertTrue(form.validate())
        
if __name__ == '__main__':
    unittest.main()