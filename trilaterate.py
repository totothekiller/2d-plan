import numpy as np

def trilaterate(points, distances):
    p1,p2,p3 = points
    r1,r2,r3 = distances
    def norm(v):
        return np.sqrt(np.sum(v**2))
    def dot(v1,v2):
        return np.dot(v1,v2)
    def cross(v1,v2):
        return np.cross(v1,v2)

    ex = (p2-p1) / norm(p2-p1)
    i = dot(ex, p3-p1)
    a = (p3-p1) - ex*i
    ey = a / norm(a)
    ez = cross(ex, ey)
    d = norm(p2-p1)
    j = dot(ey, p3-p1)
    x = (r1**2 - r2**2 + d**2) / (2*d)
    y = (r1**2 - r3**2 + i**2 + j**2) / (2*j) - (i/j) * x
    b = r1**2 - x**2 - y**2

    if (np.abs(b) < 0.0000000001):
        b = 0

    z = np.sqrt(b)
    if np.isnan(z):
        raise Exception('NaN met, cannot solve for z')

    a = p1 + ex*x + ey*y
    p4a = a + ez*z
    p4b = a - ez*z

    return p4a, p4b

def get_2d_intersections(x0, y0, r0, x1, y1, r1):
    d=np.sqrt((x1-x0)**2 + (y1-y0)**2)
    if d > r0 + r1 :
        return None
    if d < abs(r0-r1):
        return None
    if d == 0 and r0 == r1:
        return None
    else:
        a=(r0**2-r1**2+d**2)/(2*d)
        h=np.sqrt(r0**2-a**2)
        x2=x0+a*(x1-x0)/d   
        y2=y0+a*(y1-y0)/d   
        x3=x2+h*(y1-y0)/d     
        y3=y2-h*(x1-x0)/d 

        x4=x2-h*(y1-y0)/d
        y4=y2+h*(x1-x0)/d
        
        return (x3, y3, x4, y4)