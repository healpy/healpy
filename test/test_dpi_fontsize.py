"""
Test DPI and font size functionality for visufunc and newvisufunc
"""

import unittest
import re


class TestFontSizeAndDPI(unittest.TestCase):
    """Test that font sizes are relative and DPI is supported"""

    def test_visufunc_no_hardcoded_fontsizes(self):
        """Test that visufunc.py has no hard-coded font sizes in text rendering"""
        with open('lib/healpy/visufunc.py', 'r') as f:
            content = f.read()
        
        # Check that text uses fontsize from dictionary, not hard-coded numbers
        # Pattern looks for fontsize=<number> not fontsize=fontsize[...]
        text_hardcoded = re.findall(
            r'\.text\([^)]*fontsize\s*=\s*\d+[^)]*\)', content, re.DOTALL
        )
        
        self.assertEqual(len(text_hardcoded), 0,
                        f"Found {len(text_hardcoded)} instances of hardcoded numeric fontsize in .text() calls")

    def test_visufunc_has_fontsize_parameter(self):
        """Test that visufunc functions have fontsize parameter"""
        with open('lib/healpy/visufunc.py', 'r') as f:
            content = f.read()
        
        # Check that main viewing functions have fontsize parameter
        functions = ['mollview', 'gnomview', 'cartview', 'orthview', 'azeqview']
        for func in functions:
            pattern = rf'def {func}\([^)]*fontsize\s*=\s*None'
            has_param = bool(re.search(pattern, content, re.DOTALL))
            self.assertTrue(has_param,
                           f"{func} function should have fontsize=None parameter")
        
        # Check that fontsize defaults are set
        self.assertIn("fontsize_defaults", content,
                     "Functions should set up fontsize_defaults dictionary")
        self.assertIn("'large'", content,
                     "Defaults should include 'large' relative size")
        self.assertIn("'medium'", content,
                     "Defaults should include 'medium' relative size")

    def test_projview_has_dpi_parameter(self):
        """Test that projview function has dpi parameter"""
        with open('lib/healpy/newvisufunc.py', 'r') as f:
            content = f.read()
        
        # Check function signature has dpi parameter
        has_dpi_param = bool(re.search(r'def projview\([^)]*dpi\s*=\s*None', content, re.DOTALL))
        self.assertTrue(has_dpi_param,
                       "projview function signature should include dpi=None parameter")

    def test_projview_dpi_documented(self):
        """Test that dpi parameter is documented in projview docstring"""
        with open('lib/healpy/newvisufunc.py', 'r') as f:
            content = f.read()
        
        # Check that dpi is documented
        self.assertTrue('dpi :' in content or 'dpi:' in content,
                       "dpi parameter should be documented in docstring")

    def test_projview_uses_dpi_parameter(self):
        """Test that projview uses dpi parameter when creating figures"""
        with open('lib/healpy/newvisufunc.py', 'r') as f:
            content = f.read()
        
        # Check that dpi is passed to figure creation
        dpi_usage = len(re.findall(r'dpi\s*=\s*dpi', content))
        
        # We expect at least 2 usages (one for main figure, one for subplot)
        self.assertGreaterEqual(dpi_usage, 2,
                               f"Expected dpi parameter to be used at least 2 times in figure creation, found {dpi_usage}")


if __name__ == '__main__':
    unittest.main()
