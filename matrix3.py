#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/06/02 (SI)
# Modified: 2014/09/12 (JTH)

import os
from math import * 
from vector3 import *


class Matrix3:
    """define Matrix3 class and method to obtain individual matrices from lists with 9 values"""

    def __init__(self, m):
        if m == None:
            self.m = [[0,0,0],[0,0,0],[0,0,0]]
        else:
            self.m = [[0,0,0],[0,0,0],[0,0,0]]
            self.m[0][0] = m[0]
            self.m[0][1] = m[1]
            self.m[0][2] = m[2]
            self.m[1][0] = m[3]
            self.m[1][1] = m[4]
            self.m[1][2] = m[5]
            self.m[2][0] = m[6]
            self.m[2][1] = m[7]
            self.m[2][2] = m[8]

    def set_matrix(self, m):
        self.m[0][0] = m[0]
        self.m[0][1] = m[1]
        self.m[0][2] = m[2]
        self.m[1][0] = m[3]
        self.m[1][1] = m[4]
        self.m[1][2] = m[5]
        self.m[2][0] = m[6]
        self.m[2][1] = m[7]
        self.m[2][2] = m[8]

    def print_matrix(self):
        m = self.m
        print("%.6f\t"%(m[0][0])),
        print("%.6f\t"%(m[0][1])),
        print("%.6f\t"%(m[0][2]))
        print("%.6f\t"%(m[1][0])),
        print("%.6f\t"%(m[1][1])),
        print("%.6f\t"%(m[1][2]))
        print("%.6f\t"%(m[2][0])),
        print("%.6f\t"%(m[2][1])),
        print("%.6f\t"%(m[2][2]))


def matrix_from_euler(rot, tilt, psi):
    """create a rotation matrix from three Euler anges in ZYZ convention"""
    a = [0,0,0,0,0,0,0,0,0]
   
    a[0] =  cos(psi) * cos(tilt) * cos(rot) - sin(psi) * sin(rot)
    a[1] =  cos(psi) * cos(tilt) * sin(rot) + sin(psi) * cos(rot)
    a[2] = -cos(psi) * sin(tilt)
    a[3] = -sin(psi) * cos(tilt) * cos(rot) - cos(psi) * sin(rot)
    a[4] = -sin(psi) * cos(tilt) * sin(rot) + cos(psi) * cos(rot)
    a[5] =  sin(psi) * sin(tilt)
    a[6] =  sin(tilt)* cos(rot)
    a[7] =  sin(tilt)* sin(rot)
    a[8] =  cos(tilt)
    
    return Matrix3(a)	


def matrix_from_euler_zxz(rot, tilt, psi):
    """create a rotation matrix from three Euler anges in ZXZ convention"""
    a = [0,0,0,0,0,0,0,0,0]
   
    a[0] =  cos(psi) * cos(rot) - sin(psi) * cos(tilt) * sin(rot)
    a[1] = -cos(psi) * sin(rot) - sin(psi) * cos(tilt) * cos(rot)
    a[2] =  sin(psi) * sin(tilt)
    a[3] =  sin(psi) * cos(rot) + cos(psi) * cos(tilt) * sin(rot)
    a[4] = -sin(psi) * sin(rot) + cos(psi) * cos(tilt) * cos(rot)
    a[5] = -cos(psi) * sin(tilt)
    a[6] =  sin(tilt)* sin(rot)
    a[7] =  sin(tilt)* cos(rot)
    a[8] =  cos(tilt)
    
    return Matrix3(a)	


def matrix_multiply(m1, m2):
    a = [0,0,0,0,0,0,0,0,0]
    a[0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0]
    a[1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1]
    a[2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2]
    a[3] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0]
    a[4] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1]
    a[5] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2]
    a[6] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0]
    a[7] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1]
    a[8] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2]
    return Matrix3(a)


def matrix_transpose(m1):
    m2 = Matrix3(None)
    m2.m[0][0] = m1.m[0][0]
    m2.m[0][1] = m1.m[1][0]
    m2.m[0][2] = m1.m[2][0]
    m2.m[1][0] = m1.m[0][1]
    m2.m[1][1] = m1.m[1][1]
    m2.m[1][2] = m1.m[2][1]
    m2.m[2][0] = m1.m[0][2]
    m2.m[2][1] = m1.m[1][2]
    m2.m[2][2] = m1.m[2][2]
    return m2


def relion_create_symmetry_ops_file(symString, filename):
    """create a symmetry operator file by running relion_refine --print_symmetry_ops"""
    os.system("relion_refine --sym " + symString + " --print_symmetry_ops > " + " " + filename) 
    return 0 


def matrix_from_symmetry_ops_file(filename):
    """obtain the lists with 9 values for the set_matrix method by a symmetry operator file"""
   
    lines = open(filename, 'r').readlines()

    AllMatrices = []
    matrixLine = 0
    matrixFound = False
  
    for line in lines:
        values=line.split()
        if " R(" in line:
             matrixFound = True
             matrixLine = 0
             a = []
             continue
        if (matrixFound):
             a.append(float(values[0]))
             a.append(float(values[1]))
             a.append(float(values[2]))
             matrixLine = matrixLine + 1
             if (matrixLine == 3):
                 newMatrix = Matrix3(a)
                 AllMatrices.append(newMatrix)
                 matrixFound = False

    return AllMatrices

