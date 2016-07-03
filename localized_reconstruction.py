#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created:  2014/06/02 (SI)
# Modified: 2014/10/13 (JTH)
# Modified: 2015/01/29 (SI)
# Modified: 2015/08/20 (JTH)
# Modified: 2015/09/04 (SI)
# Modified: 2015/09/22 (JTH)
# Modified: 2016/01/26 (JTH)

import os
import sys
import getopt
import re
import math
import random
import copy
import time
from itertools import izip

from star import *
from particle import * 
from matrix3 import *
from vector3 import *
from euler import *
from os.path import basename
from os.path import splitext
from distutils import spawn



def within_mindist(p1, p2, mindist):
    """ Returns True if two particles are closer to each other
    than the given distance in the projection. """

    x1 = p1.rlnCoordinateX
    y1 = p1.rlnCoordinateY
    x2 = p2.rlnCoordinateX
    y2 = p2.rlnCoordinateY
    distance_sqr = (x1 - x2)**2 + (y1 - y2)**2
    mindist_sqr = mindist**2

    return distance_sqr < mindist_sqr


def within_unique(p1, p2, unique):
    """ Returns True if two particles are closer to each other
    than the given angular distance. """

    v1 = vector_from_two_eulers(p1.rlnAnglePsi, p1.rlnAngleTilt)
    v2 = vector_from_two_eulers(p2.rlnAnglePsi, p2.rlnAngleTilt)

    dp = dot_product(v1, v2)/(v1.length()*v2.length())

    if dp < -1:
        dp = -1.000

    if dp > 1:
        dp = 1.000

    angle = math.acos(dp)

    return angle <= math.radians(unique)


def filter_unique(subparticles, subpart, unique):
    """ Return True if subpart is not close to any other subparticle
        by unique (angular distance). """
    for sp in subparticles:
        if (sp.rlnImageName[:6] != subpart.rlnImageName[:6] and
            within_unique(sp, subpart, unique)):
            return False

    return True


def filter_mindist(subparticles, subpart, mindist):
    """ Return True if subpart is not close to any other subparticle
    by mindist. """
    for sp in subparticles:
        if (sp.rlnImageName[:6] != subpart.rlnImageName[:6] and
            within_mindist(sp, subpart, mindist)):
            return False

    return True


def filter_side(subpart, side):
    return (abs(abs(subpart.rlnAngleTilt) - radians(90)) < side)


def filter_top(subpart, top):
    return (abs(abs(subpart.rlnAngleTilt) - radians(180)) < top)


def filter_subparticles(subparticles, filters):
    return [sp for sp in subparticles
            if all(f(subparticles, sp) for f in filters)]


def create_subparticles(particle, symmetry_matrices, subparticle_vector_list,
                        part_image_size, relax_symmetry, randomize, output,
                        unique, subparticles_total, align_subparticles,
                        subtract_masked_map, create_star, filters):
    """ Obtain all subparticles from a given particle and set
    the properties of each such subparticle. """

    part_index = particle.rlnImageName[0:6]
    part_prefix = splitext(particle.rlnImageName[7:])[0]
    part_filename = "%s_%s.mrc" % (part_prefix, part_index)
    part_stack = "%s_%s_%s.mrcs" % (part_prefix, output, part_index)

    # Euler angles that take particle to the orientation of the model
    rot  = -particle.rlnAnglePsi
    tilt = -particle.rlnAngleTilt
    psi  = -particle.rlnAngleRot
    matrix_particle = matrix_from_euler(rot, tilt, psi)

    subparticles = []
    subtracted = []
    subpart_id = 1
    subparticles_total += 1

    symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

    if randomize:
        # randomize the order of symmetry matrices, prevents preferred views
        random.shuffle(symmetry_matrix_ids)

    for subparticle_vector in subparticle_vector_list:
        rot, tilt, psi = euler_from_vector(subparticle_vector)
        matrix_from_subparticle_vector = matrix_from_euler(rot, tilt, psi)

        for symmetry_matrix_id in symmetry_matrix_ids:
            # symmetry_matrix_id can be later written out to find out
            # which symmetry matrix created this subparticle
            symmetry_matrix = symmetry_matrices[symmetry_matrix_id-1]

            subpart = copy.deepcopy(particle)

            m = matrix_multiply((matrix_multiply(matrix_from_subparticle_vector, symmetry_matrix)), matrix_particle)

            if align_subparticles:
                rotNew, tiltNew, psiNew = euler_from_matrix(m)
            else:
                m2 = matrix_multiply(symmetry_matrix, matrix_particle)
                rotNew, tiltNew, psiNew = euler_from_matrix(m2)

            # save Euler angles that take the model to the orientation of the subparticle
            subpart.setrlnAngleRot(-psiNew)
            subpart.setrlnAngleTilt(-tiltNew)
            subpart.setrlnAnglePsi(-rotNew)

            # subparticle origin
            x = -m.m[2][0] * subparticle_vector.distance + particle.rlnOriginX
            y = -m.m[2][1] * subparticle_vector.distance + particle.rlnOriginY
            z = -m.m[2][2] * subparticle_vector.distance

            # modify the subparticle defocus paramaters by its z location
            subpart.setrlnDefocusU(particle.rlnDefocusU + z)
            subpart.setrlnDefocusV(particle.rlnDefocusV + z)

            # save the subparticle coordinates (integer part) relative to the
            # user given image size and as a small shift in the origin (decimal part)
            x_d, x_i = math.modf(x)
            y_d, y_i = math.modf(y)
            subpart.setrlnCoordinateX(int(part_image_size/2) - x_i)
            subpart.setrlnCoordinateY(int(part_image_size/2) - y_i)
            subpart.setrlnOriginX(-x_d)
            subpart.setrlnOriginY(-y_d)

            overlaps = (unique >= 0 and
                        not filter_unique(subparticles, subpart, unique))

            if not overlaps:
                subpart.setrlnImageName("%06d@%s" % (subpart_id, part_stack))
                subpart.setrlnParticleName(str(subparticles_total))
                # Taken from out of this function (a bit of duplication)
                subpart.setrlnMicrographName(part_filename)
                subpart_filename = "%s_%s_subparticles.mrcs" % (part_prefix, part_index)
                subpart.setrlnImageName(subpart.rlnImageName[0:7] + subpart_filename)
                subparticles.append(subpart)
                subpart_id += 1
                subparticles_total += 1

            if relax_symmetry:
                # take just the first (random) one and finish
                break

    if subtract_masked_map:
        subtracted = clone_subtracted_subparticles(subparticles)

    # To preserve numbering, ALL sub-particles are written to STAR files before filtering
    if create_star:
        create_star(subparticles, "%s_%s.star" % (part_prefix, part_index))
        if subtract_masked_map:
            create_star(subtracted, "%s_subtracted_%s.star"
                        % (part_prefix, part_index))

    if filters:
        subparticles = filter_subparticles(subparticles, filters)
        if subtract_masked_map:
            subtracted = clone_subtracted_subparticles(subparticles)

    return subparticles, subtracted


def clone_subtracted_subparticles(subparticles):
    subparticles_subtracted = []

    for sp in subparticles:
        sp_new = copy.deepcopy(sp)
        sp_new.rlnImageName = sp_new.rlnImageName[:-24] + "subtracted" + sp_new.rlnImageName[-25:]
        subparticles_subtracted.append(sp_new)

    return subparticles_subtracted


def create_star(subparticles, star_filename):
    """ Create a Relion style STAR file for extracting (using relion_preprocess)
    all the subparticles for a given particle. """

    parameters = ["rlnMicrographName",
                  "rlnCoordinateX",
                  "rlnCoordinateY",
                  "rlnImageName"
                  ]
    write_star(subparticles, parameters, star_filename)


def run_command(command, output=""):
    if not output:
        print "+++ " + command
        sys.stdout.flush()
        os.system(command)
    else:
        os.system(command + " > " + output)


class ProgressBar():
    """ Implements a simple command line progress bar"""

    def __init__(self, width, percent, total):
        # setup toolbar
        self.width = width
        sys.stdout.write("%s[oo]" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width))
        self.count = 0 # total count
        self.c = 0 # progress count
        self.percent = percent
        self.timer = percent
        self.total = total

    def notify(self):
        if self.count == int(self.total * self.timer):
            sys.stdout.write("\b" * (self.c + 8))
            sys.stdout.write("." * self.c)
            sys.stdout.write("~~(,_,\">")
            sys.stdout.flush()
            self.timer += self.percent
            self.c += 1

        self.count += 1

