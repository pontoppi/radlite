##FILE:
##PURPOSE:


##Below Section: Import necessary functions
doutils = True
if doutils:
    import time


##DECORATOR: Timer
def func_timer(func):
    if not doutils:
        return func

    def wrapper(*args, **kwargs):
        timestart = time.time()
        funcres = func(*args, **kwargs)
        timeend = time.time()
        print("(Function {0} took {1:.2e}s.)".format(
                                        func.__name__, (timeend-timestart)))
        return funcres
    return wrapper
