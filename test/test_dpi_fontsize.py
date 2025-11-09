"""
Test DPI and font size functionality for visufunc and newvisufunc
"""

import unittest
import inspect


class TestFontSizeAndDPI(unittest.TestCase):
    """Test that font sizes are relative and DPI is supported"""

    def test_visufunc_has_fontsize_parameter(self):
        """Test that visufunc functions have fontsize parameter"""
        import healpy as hp
        
        # Check that main viewing functions have fontsize parameter
        functions = {
            'mollview': hp.mollview,
            'gnomview': hp.gnomview,
            'cartview': hp.cartview,
            'orthview': hp.orthview,
            'azeqview': hp.azeqview
        }
        
        for func_name, func in functions.items():
            sig = inspect.signature(func)
            self.assertIn('fontsize', sig.parameters,
                         f"{func_name} function should have fontsize parameter")
            # Check default is None
            self.assertIsNone(sig.parameters['fontsize'].default,
                             f"{func_name} fontsize parameter should default to None")

    def test_graticule_has_fontsize_parameter(self):
        """Test that graticule function has fontsize parameter"""
        import healpy as hp
        
        sig = inspect.signature(hp.graticule)
        self.assertIn('fontsize', sig.parameters,
                     "graticule function should have fontsize parameter")

    def test_projview_has_dpi_parameter(self):
        """Test that projview function has dpi parameter"""
        import healpy as hp
        
        sig = inspect.signature(hp.projview)
        self.assertIn('dpi', sig.parameters,
                     "projview function should have dpi parameter")
        # Check default is None
        self.assertIsNone(sig.parameters['dpi'].default,
                         "projview dpi parameter should default to None")

    def test_projview_dpi_documented(self):
        """Test that dpi parameter is documented in projview docstring"""
        import healpy as hp
        
        docstring = hp.projview.__doc__
        self.assertIsNotNone(docstring, "projview should have a docstring")
        # Check that dpi is mentioned in the docstring
        self.assertIn('dpi', docstring.lower(),
                     "dpi parameter should be documented in docstring")


if __name__ == '__main__':
    unittest.main()
