import Solutions as sol

def avg_r(params):
    r_0 , q , mmax = params[0] , params[1] , params[2]
    r_1 = 0
    for m in range(1,mmax):
        r_1 += abs(r_m([r_0,q,m]))/m
    return r_1/np.pi

def avg_eccentricity(params):
    r_0 = params[0]
    return -avg_r(params)*np.pi/r_0
