import numpy as np
import pandas as pd
import scipy
import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint




def time_counter(args):
    def decorator(func):
        def inner():
            print(f"before do something {args}")
            result = func()
            return result
        return inner
    return decorator


@time_counter("Lee")
def do_something():
    print(f"doing something")

@time_counter("Lee")
def greet():
    print(f"Hello")


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    do_something()
    greet()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
