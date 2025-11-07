"""
Test DPI and font size functionality for visufunc and newvisufunc
"""

import unittest
import re


class TestFontSizeAndDPI(unittest.TestCase):
    """Test that font sizes are relative and DPI is supported"""

    def test_visufunc_no_hardcoded_fontsizes(self):
        """Test that visufunc.py has no hard-coded font sizes"""
        with open('lib/healpy/visufunc.py', 'r') as f:
            content = f.read()
        
        # Check for hard-coded fontsize=14 or fontsize=12
        hardcoded_14 = re.findall(r'fontsize\s*=\s*14\b', content)
        hardcoded_12 = re.findall(r'fontsize\s*=\s*12\b', content)
        
        self.assertEqual(len(hardcoded_14), 0,
                        f"Found {len(hardcoded_14)} instances of hardcoded fontsize=14")
        self.assertEqual(len(hardcoded_12), 0,
                        f"Found {len(hardcoded_12)} instances of hardcoded fontsize=12")

    def test_visufunc_uses_relative_fontsizes(self):
        """Test that visufunc.py uses relative font sizes"""
        with open('lib/healpy/visufunc.py', 'r') as f:
            content = f.read()
        
        # Check for relative font sizes
        large_count = len(re.findall(r'fontsize\s*=\s*["\']large["\']', content))
        medium_count = len(re.findall(r'fontsize\s*=\s*["\']medium["\']', content))
        
        # We expect at least 10 'large' and 1 'medium' based on the changes made
        self.assertGreaterEqual(large_count, 10,
                               f"Expected at least 10 'large' font sizes, found {large_count}")
        self.assertGreaterEqual(medium_count, 1,
                               f"Expected at least 1 'medium' font size, found {medium_count}")

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
