import random
import time
from convexhull2 import computeHull
def JarvisMarch(points):
	'''
	Naive way of implementing computeHull
	Initialize p as leftmost point
	Do following while we do not come back to the first (or leftmost) point:
    - The next point q is the point such that the triplet (p, q, r) is counterclockwise for any other point r.
    - next[p] = q (Store q as next of p in the output convex hull).
    - p = q (Set p as q for next iteration).
	O(nh) with h being the number of points
	Taken from https://iq.opengenus.org/gift-wrap-jarvis-march-algorithm-convex-hull/
	'''
	if len(points) < 3:
		return points

	points.sort()
	p = 0

	hull = []

	hull.append(points[p])
	q = (p+1) % len(points)
	for i in range(len(points)):
		if (orientation(points[p], points[i], points[q]) == 2):
			q = i
	p = q

	while (p != 0):
		hull.append(points[p])
		q = (p+1) % len(points)
		for i in range(len(points)):
			if (orientation(points[p], points[i], points[q]) == 2):
				q = i
		p = q
	
	return hull

def orientation(p, q, r):
	val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
	ans = -1
	if (val == 0):
		ans = 0
	if (val > 0):
		ans = 1
	if (val < 0):
		ans = 2
	return ans
def funct_time(f, pts):
        start = time.time()
        f(pts)
        return time.time() - start
def generate_points(n):
    random.seed()
    xs = []
    ys = []
    pts = []
    for i in range(n):
        x = random.randint(-1000000, 1000000)
        y = random.randint(-1000000, 1000000)
        ys.append(y)
        while x in xs:
            x = random.randint(-1000000, 1000000)
        xs.append(x)
    for i in range(len(xs)):
        pts.append((xs[i],ys[i]))
    return pts
def test():
    # print(2)
    # pts = generate_points(2)
    # print("{} {}".format(funct_time(computeHull, pts),funct_time(JarvisMarch, pts)))

    # print(10)
    # pts = generate_points(10)
    # print("{} {}".format(funct_time(computeHull, pts),funct_time(JarvisMarch, pts)))

    # print(1000)
    # pts = generate_points(1000)
    # print("{} {}".format(funct_time(computeHull, pts),funct_time(JarvisMarch, pts)))

    # print(10000)
    # pts = generate_points(10000)
    # print("{} {}".format(funct_time(computeHull, pts),funct_time(JarvisMarch, pts)))

    print(100000)
    pts = generate_points(100000)
    print("{} {}".format(funct_time(computeHull, pts),funct_time(JarvisMarch, pts)))
test()