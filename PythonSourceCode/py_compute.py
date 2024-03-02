from maya.OpenMaya import *




def bezier(t,p):
    return p[0]*(1-t)**3 + 3*p[1]*t*(1-t)**2 + 3*p[2]*t*t*(1-t) + p[3]*t**3

def bezier_t(v,p):
    min_t = 0.0
    max_t = 1.0
    while True:
        t = (min_t + max_t)/2
        tolerance = bezier(t,p)-v
        if tolerance > 0.0001:
            max_t= t
        elif tolerance < -0.0001:
            min_t = t
        else:
            return t

def createBezierPoints(num,p):
    px = p[0]
    py = p[1]
    j = 1/num
    points =[]
    for i in range(num+1):
        x = j*i

        t = bezier(x,px)
        y = bezier_t(t,py)
        points.append([x,y,0])
    return points

def createNormalBezierPoints(num,p):
    px = p[0]
    py = p[1]
    j = 1/num
    points = []
    for i in range(num+1):
        x = bezier(i*j,px)
        y = bezier(i*j,py)
        points.append([x,y,0])
    return points

def get_pointWeight(x,px,py):
    if x <= 0.0:
        return py[0]
    elif x >= 1.0:
        return py[3]
    t = bezier_t(x,px)
    return bezier(t,py)

def get_Weights(localPx,px,py):
    weights=[]
    for x in localPx:
        print(x)
        weight = get_pointWeight(x,px,py)
        weights.append(weight)
    return weights