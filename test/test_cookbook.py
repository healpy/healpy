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
