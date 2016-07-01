#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2016/04/04 (JTH)

import os
import sys
import getopt
import re
import math
import random

from star import *

def usage():
    print ""
    print "Usage: relion_add_ctf_info.py [parameters]"
    print "Adds CTF info from a micrographs STAR file to particles."
    print ""
    print "Parameters: "
    print "--input particles.star          Input particles STAR file."
    print "--ctf micrographs_ctf.star      Input micrographs STAR file with CTF info."
    print "--output particles_ctf.star     Output particles STAR file."
    print ""


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['input=',
                                                      'ctf=',
                                                      'output=',
                                                     ])
    except getopt.GetoptError as e:
        usage()
        print "Error: " + str(e)
        print " "
        sys.exit(2)

    for o, a in opts:
        if o == "--input":
            input_filename = a
        if o == "--ctf":
            ctf_filename = a
        if o == "--output":
            output_filename = a

    try:
        input_star=open(input_filename, "r")
        ctf_star=open(ctf_filename, "r")
    except:
        usage()
        sys.exit(2)

    particles, particles_parameters = read_star(input_star)
    micrographs, micrographs_parameters = read_star(ctf_star)

    # merge the dictionary of tags from particles and micrographs
    parameters = particles_parameters.copy()
    parameters.update(micrographs_parameters)

    for particle in particles:
        for micrograph in micrographs:
            if (micrograph.rlnMicrographName == particle.rlnMicrographName):
                particle.setrlnVoltage(micrograph.rlnVoltage)
                particle.setrlnDefocusU(micrograph.rlnDefocusU)
                particle.setrlnDefocusV(micrograph.rlnDefocusV)
                particle.setrlnDefocusAngle(micrograph.rlnDefocusAngle)
                particle.setrlnSphericalAberration(micrograph.rlnSphericalAberration)
                particle.setrlnDetectorPixelSize(micrograph.rlnDetectorPixelSize)
                particle.setrlnCtfFigureOfMerit(micrograph.rlnCtfFigureOfMerit)
                particle.setrlnMagnification(micrograph.rlnMagnification)
                particle.setrlnAmplitudeContrast(micrograph.rlnAmplitudeContrast)
                particle.setrlnCtfImage(micrograph.rlnCtfImage)
                particle.setrlnPhaseShift(micrograph.rlnPhaseShift)
                continue

    write_star(particles, parameters, output_filename)

    input_star.close()

if __name__ == "__main__":    
    main()
