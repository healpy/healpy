from functools import wraps

def is_singlemap(m):
    return not isinstance(m, (tuple, list))

def listify_first_argument(f):
    """If first argument is not list or tuple, makes a 1 element list"""
    @wraps(f)
    def wrapper(*args, **kwds):
        if not isinstance(args[0], (tuple, list)):
            args[0] = [args[0]] 
        return f(*args, **kwds)
    return wrapper
