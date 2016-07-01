#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2015/08/13 (JTH)

import os
import sys
import getopt
import re
import math

from star import *


def usage():
    print ""
    print "Usage: relion_eulers_from_dynamo.py [parameters]"
    print "Adds Euler angles from a Dynamo table file into the STAR file."
    print ""
    print "Parameters: "
    print "--input particles.star          Input STAR file."
    print "--dynamotbl particles.tbl       Input Dynamo table file."
    print "--output particles_eulers.star  Output STAR file."
    print ""


def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], '', ['input=',
                                                              'dynamotbl=',
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
		if o == "--dynamotbl":
			tbl_filename = a
		if o == "--output":
			output_filename = a

	try:
		input_star=open(input_filename, "r")
		input_tbl=open(tbl_filename, "r")
	except:
		usage()
		sys.exit(2)

        particles, parameters = read_star(input_star)


        tbl_lines = input_tbl.readlines()

        i = 0
        for particle in particles:
            tbl_values = tbl_lines[i].split()
            rot, tilt, psi = eulers_from_dynamo()
            particle.setrlnAngleRot(math.radians(float(tbl_values[6])))
            particle.setrlnAngleTilt(math.radians(float(tbl_values[7])))
            particle.setrlnAnglePsi(math.radians(float(tbl_values[8])))

            i = i + 1

        new_parameters = {"rlnAngleRot": '', "rlnAngleTilt" : '', "rlnAnglePsi" : ''}
        parameters.update(new_parameters)
		
	write_star(particles, parameters, output_filename)

	input_star.close()
	input_tbl.close()

if __name__ == "__main__":    
    main()
