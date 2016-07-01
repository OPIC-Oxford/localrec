#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created:  2016/01/29 (SI)
# Modified: 2016/04/28 (JTH)

import os
import sys
import getopt
import re
import math
import random
import copy
import time

from star import *
from particle import * 
from matrix3 import *
from vector3 import *
from euler import *
from os.path import basename
from os.path import splitext
from distutils import spawn


def usage():
    print ""
    print "Usage: relion_reconstruct_relax.py [parameters] input.star"
    print "Obtains an asymmetrical reconstruction of a symmetrical object by relaxing its symmetry."
    print "This reconstruction can be used as a starting model for refinement without symmetry."
    print ""
    print "Parameters: "  
    print "--angpix 1.35                    Pixel size (A)."
    print "--relax_symmetry I1              Symmetry of the particle to be relaxed for the asymmetrical reconstruction."
    print "--output_map output.mrc          Name of the output map."
    print "--output_star output.star        Name of the output star file."    
    print "--maxres 5                       Maximum resolution."
    print "--j 8                            Number of threads."
    print ""


def run_command(command, output=""):
    if not output:
        print "+++ " + command
        os.system(command)
    else:
        os.system(command + " > " + output)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['angpix=',
                                                      'relax_symmetry=', 
                                                      'output_map=',
                                                      'output_star=',
                                                      'maxres=',
                                                      'j=',
                                                     ])
    except getopt.GetoptError as e:
        usage()
        print "Error: " + str(e)
        print " "
        sys.exit(2)

    if not (spawn.find_executable("relion_refine")):
        usage()
        print "Error: Relion not found."
        print "Make sure Relion programs are in $PATH."
        print " "
        sys.exit(2)

    relax_symmetry = "C1"
    angpix = -1
    maxres = -1
    nr_threads = 1
    output_star = ""
    output_map = ""

    for o, a in opts:
        if o == "--relax_symmetry":
            relax_symmetry = str(a)
        if o == "--angpix":
            angpix = float(a)
        if o == "--output_star":
            output_star = str(a)
        if o == "--output_map":
            output_map = str(a)
        if o == "--j":
            nr_threads = int(a)
        if o == "--maxres":
            maxres = str(a)

    if len(args) == 0:
        print " "
        print "Error: No input file given."
        usage()
        sys.exit(2)

    try:        
        open(args[0],"r") 
    except:
        print " "
        print "Error: Input file not found."
        usage()
        sys.exit(2)

    if angpix < 0:
        print " "
        print "Error: Pixel size not given."
        usage()
        sys.exit(2)

    if maxres < 0:
        print " "
        print "Warning: Maximum resolution not given, calculating the map to Nyquist."
        maxres = 2*angpix

    input_star_filename = str(args[0])
    matrices = []
    particles = []
    parameters = [] 
    
    # setup toolbar   
    toolbar_width = 70
    sys.stdout.write("%s[oo]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width))
    c = 0   
    timer = 0.01

    input_star = open(input_star_filename, "r")
    particles, parameters = read_star(input_star)
    input_star.close()

    relion_create_symmetry_ops_file(relax_symmetry, "relion_symops.tmp")
    symmetry_matrices = matrix_from_symmetry_ops_file("relion_symops.tmp")
    run_command("rm relion_symops.tmp", "")

    new_particles = []

    nparticle = 0

    for particle in particles:
        rot  = -particle.rlnAnglePsi
        tilt = -particle.rlnAngleTilt
        psi  = -particle.rlnAngleRot
        matrix_particle = matrix_from_euler(rot, tilt, psi)

        index = random.randint(0, len(symmetry_matrices)-1)

        symmetry_matrix = symmetry_matrices[index]

        new_particle = copy.deepcopy(particle)

        m = matrix_multiply(symmetry_matrix, matrix_particle)

        rotNew, tiltNew, psiNew = euler_from_matrix(m)

        new_particle.setrlnAngleRot(-psiNew)
        new_particle.setrlnAngleTilt(-tiltNew)
        new_particle.setrlnAnglePsi(-rotNew)

        new_particles.append(new_particle)

    write_star(new_particles, parameters, output_star)

    string = "relion_reconstruct --i " + str(output_star) + " --o " + str(output_map) + " --ctf --sym C1 --angpix " + str(angpix) + " --maxres " + str(maxres) + " --j " + str(nr_threads)
    print string
    os.system(string) 
    
    print "The output files have been written!"
    print " " 

    sys.exit(0)

if __name__ == "__main__":    
    main()

