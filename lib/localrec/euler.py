# **************************************************************************
# *
# * Authors:  Serban Ilca
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *
# * Oxford Particle Imaging Centre,
# * Division of Structural Biology, University of Oxford
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# **************************************************************************

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



