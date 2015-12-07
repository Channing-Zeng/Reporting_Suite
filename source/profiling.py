import time
from functools import wraps
from source import info


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        info('Total time running %s: %s seconds' % (function.func_name, str(t1-t0)))
        return result

    return function_timer
