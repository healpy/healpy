"""Various generic useful functions 
"""


def is_seq(o):
    """Check if the object is a sequence.

    Parameters
    ----------
    o : any object
      The object to check
    
    Returns
    -------
    is_seq : bool, scalar
      True if *o* is a sequence, False otherwise
    """
    return hasattr(o, "__len__")


def is_seq_of_seq(o):
    """Check if the object is a sequence of sequences. No check is done on
    the sizes of the sequences.

    Parameters
    ----------
    o : any object
      The object to check
    
    Returns
    -------
    is_seq_of_seq : bool
      True if *o* is a sequence of sequences, False otherwise.
    """
    if not is_seq(o):
        return False
    for s in o:
        if not is_seq(s):
            return False
    return True


def is_like2d(o):
    """Check if *o* is conformable to a 2d array.

    Parameters
    ----------
    o : any object
      The object to check
   
    Returns
    -------
    is_like2d : bool, scalar
      True if *o* is conformable to a 2d array, False otherwise.
    """
    if not is_seq(o):
        return False
    size = None
    for s in o:
        if not is_seq(s):
            return False
        if size is None:
            size = len(s)
        if len(s) != size:
            return False
    return True


def len_array_or_arrays(o):
    """Returns the length of a single array or list of arrays
    
    Parameters
    ----------
    o : either array or sequence of arrays

    Returns
    -------
    length : length of array
    """
    if is_seq_of_seq(o):
        length = len(o[0])
    else:
        length = len(o)
    return length
