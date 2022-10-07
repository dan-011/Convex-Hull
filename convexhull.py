import math
import sys

EPSILON = sys.float_info.epsilon

'''
Given two points, p1 and p2,
an x coordinate, x,
and y coordinates y3 and y4,
compute and return the (x,y) coordinates
of the y intercept of the line segment p1->p2
with the line segment (x,y3)->(x,y4)
'''
def yint(p1, p2, x, y3, y4):
	x1, y1 = p1
	x2, y2 = p2
	x3 = x
	x4 = x
	px = ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / \
		 float((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4))
	py = ((x1*y2 - y1*x2)*(y3-y4) - (y1 - y2)*(x3*y4 - y3*x4)) / \
			float((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3-x4))
	return (px, py)

'''
Given three points a,b,c,
computes and returns the area defined by the triangle
a,b,c. 
Note that this area will be negative 
if a,b,c represents a clockwise sequence,
positive if it is counter-clockwise,
and zero if the points are collinear.
'''
def triangleArea(a, b, c):
	return (a[0]*b[1] - a[1]*b[0] + a[1]*c[0] \
                - a[0]*c[1] + b[0]*c[1] - c[0]*b[1]) / 2.0;

'''
Given three points a,b,c,
returns True if and only if 
a,b,c represents a clockwise sequence
(subject to floating-point precision)
'''
def cw(a, b, c):
	return triangleArea(a,b,c) < EPSILON;
'''
Given three points a,b,c,
returns True if and only if 
a,b,c represents a counter-clockwise sequence
(subject to floating-point precision)
'''
def ccw(a, b, c):
	return triangleArea(a,b,c) > EPSILON;

'''
Given three points a,b,c,
returns True if and only if 
a,b,c are collinear
(subject to floating-point precision)
'''
def collinear(a, b, c):
	return abs(triangleArea(a,b,c)) <= EPSILON

'''
Given a list of points,
sort those points in clockwise order
about their centroid.
Note: this function modifies its argument.
'''
def clockwiseSort(points):
	# get mean x coord, mean y coord
	xavg = sum(p[0] for p in points) / len(points)
	yavg = sum(p[1] for p in points) / len(points)
	angle = lambda p:  ((math.atan2(p[1] - yavg, p[0] - xavg) + 2*math.pi) % (2*math.pi))
	points.sort(key = angle)

'''
Replace the implementation of computeHull with a correct computation of the convex hull
using the divide-and-conquer algorithm
'''
def computeHull(points):
	'''
	This function sorts the set of points by x values and creates a hull using the _computehull helper function. The function then clockwise sorts the hull and returns
	it.
	'''
	points.sort()
	hull = _computeHull(points)
	clockwiseSort(hull)
	return hull
	
def _computeHull(points):
	'''
	_computehull is the helper function to the main function which breaks the set of points down to sizes of 2 or 3 points by dividing the set of points in half,
	creating lefthull and righthull with each recursion. The function then uses the lower_tangent and upper_tangent functions to find the points that create the lower
	and upper tangent lines of lefthull and righthull. The function then clockwise-sorts lefthull and righthull and loops through both sets of points. Starting from
	the upper tangent point of the lefthull moving along the lefthull to the lower tangent point, adding points to a new list pts along the way. Then, the function starts
	at the lower tangent point of the righthull and adds to pts all of the points up until the upper tangent point of the righthull. Finally, pts is sorted
	by x values and returned.
	'''
	if len(points) > 3:
		median = len(points)//2
		leftHull = _computeHull(points[:median])
		rightHull = _computeHull(points[median:])
	else:
		return points
	
	pts = []
	clockwiseSort(leftHull)
	clockwiseSort(rightHull)
	i_la,i_lb = lower_tangent(leftHull, rightHull)
	i_ua,i_ub = upper_tangent(leftHull, rightHull)
	for i in range(i_ua, i_la+1):
		pts.append(leftHull[i])
	if i_ub < i_lb:
		for i in range(i_lb, len(rightHull)):
			pts.append(rightHull[i])
	for i in range(i_ub+1):
		pts.append(rightHull[i])
	pts.sort()
	return pts
def lower_tangent(A, B):
	'''
	This function works by taking two sets, of points, A and B, and creating sets cws_A and cws_B that are the clockwise sorted sets of the points in A and B. The
	function starts with the rightmost element in A and the leftmost element in B and iterates through a loop which decrements the index of a and increments the index
	of b, tracing along the closkwise sorted sets until point a is in the lower tangent of cws_A, point b is in the lower tangent of cws_B, and a and b form the lower
	tangent of set cws_A and set cws_B. The function then returns the indices of the two points.
	'''
	a = A[-1]
	b = B[0]
	cws_A = []
	cws_B = []
	cws_A.extend(A)
	cws_B.extend(B)
	clockwiseSort(cws_A)
	clockwiseSort(cws_B)
	T = tangent(a,b)
	m = T[0]
	yin = T[1]
	i_a = cws_A.index(a)
	i_b = cws_B.index(b)
	before_a = (i_a - 1)%len(cws_A)
	before_b = (i_b - 1)%len(cws_B)
	after_a = (i_a + 1)%len(cws_A)
	after_b = (i_b + 1)%len(cws_B)
	
	while not(lower_tangent_good(before_a, after_a, yin, m, cws_A) and lower_tangent_good(before_b, after_b, yin, m, cws_B)):
		while not lower_tangent_good(before_a, after_a, yin, m, cws_A):
			i_a -= 1
			i_a = i_a%len(cws_A)
			a = cws_A[i_a]
			T = tangent(a,b)
			m = T[0]
			yin = T[1]
			before_a = (i_a - 1)%len(cws_A)
			after_a = (i_a + 1)%len(cws_A)
		while not lower_tangent_good(before_b, after_b, yin, m, cws_B):
			i_b += 1
			i_b = i_b%len(cws_B)
			b = cws_B[i_b]
			T = tangent(a,b)
			m = T[0]
			yin = T[1]
			before_b = (i_b - 1)%len(cws_B)
			after_b = (i_b + 1)%len(cws_B)
	return (i_a, i_b)

def lower_tangent_good(before, after, intercept, slope, pts):
	'''
	The input for this function is the index of the point before, the index of the point after, the intercept of the tangent line, the slope of the tangent line,
	and the set of clockwise sorted points. The purpose of this function is to check if the point is the lowermost point on the set of points. It does this by finding
	the y value of what the points before and after the point being checked are. The idea is that, if the point being checked is truly the uppermost point, then the y
	value of the points before and after it should be at or lie above the line created by the point and the point in the other set of points.
	'''
	if intercept is None: return True
	before_point = slope*pts[before][0] + intercept
	after_point = slope*pts[after][0] + intercept
	return before_point <= pts[before][1] and after_point <= pts[after][1]

def upper_tangent(A, B):
	'''
	This function works by taking two sets, of points, A and B, and creating sets cws_A and cws_B that are the clockwise sorted sets of the points in A and B. The
	function starts with the rightmost element in A and the leftmost element in B and iterates through a loop which increments the index of a and decrements the index
	of b, tracing along the closkwise sorted sets until point a is in the upper tangent of cws_A, point b is in the upper tangent of cws_B, and a and b form the upper
	tangent of set cws_A and set cws_B. The function then returns the indices of the two points.
	'''
	a = A[-1]
	b = B[0]
	cws_A = []
	cws_B = []
	cws_A.extend(A)
	cws_B.extend(B)
	clockwiseSort(cws_A)
	clockwiseSort(cws_B)
	T = tangent(a,b)
	m = T[0]
	yin = T[1]
	i_a = cws_A.index(a)
	i_b = cws_B.index(b)
	before_a = (i_a - 1)%len(cws_A)
	before_b = (i_b - 1)%len(cws_B)
	after_a = (i_a + 1)%len(cws_A)
	after_b = (i_b + 1)%len(cws_B)
	
	while not(upper_tangent_good(before_a, after_a, yin, m, cws_A) and upper_tangent_good(before_b, after_b, yin, m, cws_B)):
		while not upper_tangent_good(before_a, after_a, yin, m, cws_A):
			
			i_a += 1
			i_a = i_a%len(cws_A)
			a = cws_A[i_a]
			before_a = (i_a - 1)%len(cws_A)
			after_a = (i_a + 1)%len(cws_A)
			T = tangent(cws_A[i_a],cws_B[i_b])
			m = T[0]
			yin = T[1]
		while not upper_tangent_good(before_b, after_b, yin, m, cws_B):
			i_b -= 1
			i_b = i_b%len(cws_B)
			b = cws_B[i_b]
			before_b = (i_b - 1)%len(cws_B)
			after_b = (i_b + 1)%len(cws_B)
			T = tangent(a,b)
			m = T[0]
			yin = T[1]
	return (i_a, i_b)
def upper_tangent_good(before, after, intercept, slope, pts):
	'''
	The input for this function is the index of the point before, the index of the point after, the intercept of the tangent line, the slope of the tangent line,
	and the set of clockwise sorted points. The purpose of this function is to check if the point is the uppermost point on the set of points. It does this by finding the y
	value of what the points before and after the point being checked are. The idea is that, if the point being checked is truly the uppermost point, then the y value of the
	points before and after it should be at or lie below the line created by the point and the point in the other set of points.
	'''
	if intercept is None: return True
	before_point = slope*pts[before][0] + intercept
	after_point = slope*pts[after][0] + intercept
	return pts[before][1] <= before_point and pts[after][1] <= after_point

def tangent(p1, p2):
	'''
	Given two points, break down the points in to their x and y components. Use these components to evaluate the 
	slope and y intercept of the line the two points form. Included in this function is a try block that sets the y intercept to None if
	the two points are along the same y-axis point.
	'''
	x1 = p1[0]
	x2 = p2[0]
	y1 = p1[1]
	y2 = p2[1]

	m = (y2-y1)/(x2-x1)
	try:
		yin = (yint(p1, p2, 0, 1, 2))[1]
	except:
		yin = None
	return (m,yin)
