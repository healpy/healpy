def test_is_seq():
    import numpy as np
    from healpy.cookbook import is_seq

    assert not is_seq(None)
    assert not is_seq(1)
    assert not is_seq(1.)
    assert not is_seq(np.array(1))
    assert is_seq((1, 2, 3))
    assert is_seq([1, 2, 3])
    assert is_seq(np.array([1, 2, 3]))
    assert is_seq(np.array([[1], [2], [3]]))
    assert is_seq(())
    assert is_seq([])
    assert is_seq(np.array([]))


def test_is_seq_of_seq():
    import numpy as np
    from healpy.cookbook import is_seq_of_seq

    assert not is_seq_of_seq(None)
    assert not is_seq_of_seq(1)
    assert not is_seq_of_seq(1.)
    assert not is_seq_of_seq(np.array(1))
    assert not is_seq_of_seq((1, 2, 3))
    assert not is_seq_of_seq([1, 2, 3])
    assert not is_seq_of_seq(np.array([1, 2, 3]))
    assert is_seq_of_seq(((1, 2, 3), (4, 5), (6,)))
    assert is_seq_of_seq([[1], [2, 3], [4, 5, 6]])
    assert is_seq_of_seq(np.array([[1], [2], [3]]))
    assert is_seq_of_seq(((1,), [2], np.array([3])))
    assert is_seq_of_seq(())
    assert is_seq_of_seq([])
    assert is_seq_of_seq(np.array([]))

    # allow None
    assert not is_seq_of_seq([[1], [2], None], False)
    assert is_seq_of_seq([[1], [2], None], True)


def test_len_array_or_arrays():
    import numpy as np
    from healpy.cookbook import len_array_or_arrays
    
    # Test with single array
    assert len_array_or_arrays(np.array([1, 2, 3])) == 3
    assert len_array_or_arrays([1, 2, 3, 4]) == 4
    
    # Test with list of arrays (all non-None)
    assert len_array_or_arrays([[1, 2, 3], [4, 5, 6]]) == 3
    assert len_array_or_arrays([np.array([1, 2, 3]), np.array([4, 5, 6])]) == 3
    
    # Test with None in list of arrays (regression test for bug)
    assert len_array_or_arrays([[1, 2, 3], [4, 5, 6], None]) == 3
    assert len_array_or_arrays([None, [1, 2, 3], [4, 5, 6]]) == 3
    assert len_array_or_arrays([[1, 2, 3], None, [4, 5, 6]]) == 3
    
    # Test with numpy arrays and None
    c_ee = np.linspace(0, 3e-6, 10000)
    c_ne = np.linspace(0, 1e-6, 10000)
    c_nn = np.linspace(0, 3e-5, 10000)
    c_ell = [c_nn, c_ne, c_ee, None]
    assert len_array_or_arrays(c_ell) == 10000
    
    # Test edge case: all None
    assert len_array_or_arrays([None, None, None]) == 3
