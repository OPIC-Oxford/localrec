#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/06/02 (SI)
# Modified: 2014/09/12 (JTH)

import re
from math import *

class Vector3:
    """define Vector3 class and method to obtain individual vectors from lists with 3 values"""

    def __init__(self, v):
        if v == None:
            self.v = [0,0,1]
        else:
            self.v = v

    def set_vector(self, v):
        self.v = v

    def set_distance(self, d):
        self.distance = float(d)

    def print_vector(self):
        x, y, z = self.v
        print("[%.3f,%.3f,%.3f]"%(x, y, z)),

    def distance(self):
        d = self.distance
        return float(d)

    def length(self):
        x, y, z = self.v
        return sqrt(x**2 + y**2 + z**2)

    def normalize(self):
        try:
            magnitude = self.length()
            self.v[0] = self.v[0] / magnitude
            self.v[1] = self.v[1] / magnitude
            self.v[2] = self.v[2] / magnitude
        except ZeroDivisionError:
            self.v[0] = 0
            self.v[1] = 0
            self.v[2] = 0

    def x(self):
        return self.v[0]

    def y(self):
        return self.v[1]

    def z(self):
        return self.v[2]


def dot_product(v1, v2):
    """returns the dot product of two vectors"""

    x1, y1, z1 = v1.v
    x2, y2, z2 = v2.v

    return x1 * x2 + y1 * y2 + z1 * z2


def cross_product(v1, v2):
    """returns the cross product of two vectors"""

    x1, y1, z1 = v1.v
    x2, y2, z2 = v2.v

    x = y1 * z2 - y2 * z1
    y = z1 * x2 - z2 * x1
    z = x1 * y2 - x2 * y1

    return Vector3([x,y,z]) 


def vector_from_two_eulers(rot, tilt):
    """function that obtains a vector from the first two Euler angles"""

    x = sin(tilt)*cos(rot)
    y = sin(tilt)*sin(rot)
    z = cos(tilt)
    return Vector3([x,y,z]) 


def vectors_from_cmm(input_cmm, angpix):
    """function that obtains the input vector from a cmm file"""
    
    # coordinates in the CMM file need to be in Angstrom

    file_cmm = open(input_cmm, "r")
    vector_list = []
    counter=0

    for line in file_cmm.readlines():
        if 'marker id=' in line:
            line_values=line.split()
            for i in range(len(line_values)):
                if 'x=' in line_values[i]:
                    a = re.search('"(.*)"', line_values[i]).group(0)
                    x = float(a.translate(None, '""'))/angpix
                if 'y=' in line_values[i]:
                    b = re.search('"(.*)"', line_values[i]).group(0)
                    y = float(b.translate(None, '""'))/angpix
                if 'z=' in line_values[i]:
                    c = re.search('"(.*)"', line_values[i]).group(0)
                    z = float(c.translate(None, '""'))/angpix
       
            if counter != 0:
                vector = Vector3(None)
                x = x - x0
                y = y - y0
                z = z - z0
                vector.set_vector([x,y,z])
                vector_list.append(vector)
                counter = counter + 1
                continue
            else:
                x0 = x
                y0 = y
                z0 = z
                counter = counter + 1
                continue
        else:
             continue

    return vector_list
