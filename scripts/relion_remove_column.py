#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2016/04/26 (JTH)

import os
import sys
import getopt
import re
import math
import random

from star import *


def usage():
    print ""
    print "Usage: relion_remove_column.py [parameters]"
    print "Removes a column from a STAR file."
    print ""
    print "Parameters: "
    print "--input particles.star          Input STAR file."
    print "--name rlnGroupName             Name of the column to be removed."
    print "--output particles_out.star     Output STAR file."
    print ""

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['input=',
                                                      'name=',
                                                      'output=',
                                                     ])
    except getopt.GetoptError as e:
        usage()
        print "Error: " + str(e)
        print " "
        sys.exit(2)

    micrograph = False

    for o, a in opts:
        if o == "--input":
            input_filename = a
        if o == "--name":
            name = a
        if o == "--output":
            output_filename = a
    try:
        input_star=open(input_filename, "r")
    except:
        usage()
        sys.exit(2)


    particles, parameters = read_star(input_star)

    if name in parameters:
        del parameters[name]

    write_star(particles, parameters, output_filename)

    input_star.close()

if __name__ == "__main__":    
    main()
