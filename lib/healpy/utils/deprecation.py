# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Copied from `gammapy`

__all__ = [
    "HealpyDeprecationWarning",
    "deprecated",
    "deprecated_renamed_argument",
    "deprecated_attribute",
]


class HealpyDeprecationWarning(Warning):
    """
    The Healpy deprecation warning
    """


def deprecated(since, **kwargs):
    """
    Used to mark a function or class as deprecated.

    Reuses Astropy's deprecated decorator.
    Check arguments and usage in `~astropy.utils.decorator.deprecated`

    Parameters
    ----------
    since : str
        The release at which this API became deprecated.  This is required.
    """
    from astropy.utils import deprecated

    kwargs["warning_type"] = HealpyDeprecationWarning
    return deprecated(since, **kwargs)


def deprecated_renamed_argument(old_name, new_name, since, **kwargs):
    """Deprecate a _renamed_ or _removed_ function argument.

    Check arguments and usage in `~astropy.utils.decorator.deprecated_renamed_argument`
    """
    from astropy.utils import deprecated_renamed_argument

    kwargs["warning_type"] = HealpyDeprecationWarning
    return deprecated_renamed_argument(old_name, new_name, since, **kwargs)


def deprecated_attribute(name, since, **kwargs):
    """
    Used to mark a public attribute as deprecated.  This creates a
    property that will warn when the given attribute name is accessed.
    """
    from astropy.utils import deprecated_attribute

    kwargs["warning_type"] = HealpyDeprecationWarning
    return deprecated_attribute(name, since, **kwargs)
