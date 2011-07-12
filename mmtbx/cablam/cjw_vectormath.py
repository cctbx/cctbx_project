#A set of basic vector math functions commonly used by Christopher Williams
#  programs, including but not limited to cablam.  Replace as convenient with
#  in-house functions
#Note: Assumes all vectors are 3-dimensional

import math

#Returns a vector connecting point p1 to point p2
def vectorize(p1, p2):
  v = [ p2[0]-p1[0] , p2[1]-p1[1] , p2[2]-p1[2] ]
  return v

#Returns the scalar length of a vector
def veclen(v):
  return math.sqrt( v[0]**2 + v[1]**2 + v[2]**2 )

#Dot product of two vectors
def dot(v1, v2):
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

#Cross product of two vectors
def cross(v1, v2):
  x = v1[1]*v2[2] - v1[2]*v2[1]
  y = v1[2]*v2[0] - v1[0]*v2[2]
  z = v1[0]*v2[1] - v1[1]*v2[0]
  return [x,y,z]

#Finds the line from a1 to a2, drops a perpendicular to it from b1, and returns
#  the point of intersection.
def perptersect(a1, a2, b1):
  #Find the slope of line A in each direction, A is in vector notation
  A = [a2[0]-a1[0], a2[1]-a1[1], a2[2]-a1[2]]
  #Solve the parametric equations (dot of perpendiculars=0). . .
  t = (A[0]*(b1[0]-a1[0]) + A[1]*(b1[1]-a1[1]) + A[2]*(b1[2]-a1[2])) / ((A[0]**2)+(A[1]**2)+(A[2]**2))
  # . . . and use the result to find the new point b2 on the line
  b2 = [a1[0]+A[0]*t, a1[1]+A[1]*t, a1[2]+A[2]*t]
  return b2

#Returns the perpendicular distance from point b1 to the a1-a2 line
def perpdist(a1, a2, b1):
  b2 = perptersect(a1, a2, b1)
  distance = veclen(vectorize(b1,b2))
  return distance
