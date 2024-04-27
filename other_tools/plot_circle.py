import matplotlib.pyplot as plt
import numpy as np
import math

def area(r1, r2, d):
    s1 = r1**2 * math.acos((d**2 - r2**2 + r1**2)/(2 * d * r1))
    s2 = r2**2 * math.acos((d**2 - r1**2 + r2**2)/(2 * d * r2))
    s3 = math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d + r2 - r1) * (d + r1 + r2)) / 2
    return s1 + s2 - s3

def solve(r1, r2, targetArea, error):
    d1 = r1 + r2
    a1 = area(r1, r2, d1)
    d2 = abs(r1 - r2)
    a2 = area(r1, r2, d2)
    while True:
        d = (d1 + d2) / 2
        a = area(r1, r2, d)
        if abs(a - targetArea) <= error:
            return d
        if a > targetArea:
            d2 = d
            a2 = a
        else:
            d1 = d
            a1 = a

# input data
r1 = math.sqrt(1372/math.pi)
r2 = math.sqrt(1348/math.pi)
a = 1333
error = 0.00001

# find d for the given input 
d = solve(r1, r2, a, error)
print("The result is d=%.5f" % d)

# calculate area for the calculated value of d
a2 = area(r1, r2, d)

# print for verification
print("For r1=%.5f, r2=%.5f, d=%.5f the area is %.5f and the target area is %.5f"
      % (r1, r2, d, a2, a))