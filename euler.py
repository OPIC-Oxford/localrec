#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/06/02 (SI)
# Modified: 2014/13/09 (JTH)

from math import *
import vector3


def euler_from_matrix(matrix):
    """converts a matrix to Eulers"""

    # this needs to be true for acos: -1.0 <= a < 1.0
    a = matrix.m[2][2]
    if a < -1:
        a = -1.000
    if a > 1:
        a = 1.000

    # tilt is always positive here, which means the final value will always be negative.
    # this results the other of the two equivalent options being written in the output file
    tilt = acos(a)

    if abs(tilt) < 0.00001:
        rot = radians(0.00)
        psi = atan2(-matrix.m[1][0], matrix.m[0][0])
    else:
        rot = atan2(matrix.m[2][1], matrix.m[2][0])
        psi = atan2(matrix.m[1][2], -matrix.m[0][2])

    return rot, tilt, psi


def euler_from_vector(vector):   
    """converts a view vector to Eulers that describe a rotation taking the reference vector [0,0,1] on the vector""" 

    vector.normalize()
    
    if abs(vector.x()) < 0.00001 and abs(vector.y()) < 0.00001:
        rot = radians(0.00)
        tilt = radians(0.00)
    else:
        rot = atan2(vector.y(), vector.x())
        tilt = acos(vector.z())

    psi = 0
    
    return rot, tilt, psi


def euler_print(rot, tilt, psi):
    """prints three angles in degrees"""

    print("%.6f\t%.6f\t%.6f"%(degrees(rot),degrees(tilt), degrees(psi)))



