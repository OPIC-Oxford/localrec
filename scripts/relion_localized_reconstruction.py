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
from localized_reconstruction import *


def usage():
    print ""
    print "Usage: relion_create_subparticles.py [parameters] input.star"
    print "Calculates the coordinates and Euler angles for the subparticles defined by a vector and symmetry point group."
    print ""
    print "Parameters: "
    print "--split_stacks                   Split particle stacks (needs to be done once)."
    print "--masked_map file.mrc            Map with density that to be subtracted from particle images."
    print "--create_star                    Create new STAR files for extracting sub-particles."
    print "--extract_subparticles           Extract sub-particles from particle images."  
    print "--angpix 1.35                    Pixel size (A)."
    print "--sym I1                         Symmetry of the particle."
    print "--particle_size 512              Size of the particle box (pixels)."
    print "--subparticle_size 128           Size of the sub-particle box (pixels)."
    print "--randomize                      Randomize the order of the symmetry matrices. Useful for preventing preferred orientations (default: not)."
    print "--relax_symmetry                 Create one random subparticle for each particle (default: all symmetry related subparticles)."
    print "--vector 0,0,1                   Vector defining the location of the subparticle."
    print "--align_subparticles             Align subparticles to the standard orientation."
    print "--length 100                     Alternative length of the vector. Use to adjust the sub-particle center (default: length of the given vector; A)."
    print "--cmm vector.cmm                 A CMM file defining the location(s) of the subparticle(s) (use instead of --vector). Coordinates in Angstrom."
    print "--unique 1                       Keep only unique subparticles within angular distance (useful to remove overlapping sub-particles on symmetry axis)."
    print "--mindist 5                      Minimum distance between the subparticles in the image (all overlapping ones will be discarded; pixels)."
    print "--side 25                        Keep only particles within specified angular distance from side views (all others will be discarded; degrees)."
    print "--top 25                         Keep only particles within specified angular distance from top views (all others will be discarded; degrees)."
    print "--output subparticles            Output root for results."
    print "--j 8                            Number of threads."
    print "--np 4                           Number of MPI procs."
    print ""


def error(*msgs):
    usage()
    print "Error: " + '\n'.join(msgs)
    print " "
    sys.exit(2)



def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['sym=',
                                                      'randomize',
                                                      'relax_symmetry',
                                                      'vector=',
                                                      'align_subparticles',
                                                      'cmm=',
                                                      'length=',
                                                      'output=',
                                                      'create_star',
                                                      'unique=',
                                                      'mindist=',
                                                      'top=',
                                                      'side=',
                                                      'particle_size=',
                                                      'subparticle_size=',
                                                      'masked_map=',
                                                      'split_stacks',
                                                      'extract_subparticles',
                                                      'angpix=',
                                                      'output=',
                                                      'j=',
                                                      'np=',
                                                      ])
    except getopt.GetoptError as e:
        error(e)

    if not (spawn.find_executable("bimg")):
        error("Error: Bsoft not found.",
              "Make sure Bsoft programs are in $PATH.")

    if not (spawn.find_executable("relion_refine")):
        error("Error: Relion not found.",
              "Make sure Relion programs are in $PATH.")

    symString = "C1"
    output = ""
    subparticle_distances = []
    subparticle_vector_list = []
    input_cmm = ""
    part_image_size = -1
    relax_symmetry = False
    align_subparticles = False
    randomize = False
    star = False
    unique = -1
    mindist = -1
    minanglediff = -1
    side = -1
    top = -1
    masked_map = ""
    subtract_masked_map = False
    split_stacks = False
    extract_subparticles = False
    angpix = -1
    output_dir = ""
    nr_threads = 1
    nr_procs = 1

    vectors_string = ''
    distances_string = ''

    for o, a in opts:
        if o == "--randomize":
            randomize = True
        if o == "--relax_symmetry":
            relax_symmetry = True
        if o == "--sym":
            symString = str(a)
        if o == "--vector":
            vectors_string = a
        if o == "--align_subparticles":
            align_subparticles = True
        if o == "--length":
            distances_string = a
        if o == "--top":
            top = radians(float(a))
        if o == "--side":
            side = radians(float(a))
        if o == "--mindist":
            mindist = float(a)
        if o == "--unique":
            unique = float(a)
        if o == "--cmm":
            input_cmm = a
        if o == "--create_star":
            star = True
        if o == "--particle_size":
            part_image_size = int(a)
        if o == "--subparticle_size":
            subpart_image_size = int(a)
        if o == "--masked_map":
            masked_map = str(a)
            subtract_masked_map = True
        if o == "--split_stacks":
            split_stacks = True
        if o == "--extract_subparticles":
            extract_subparticles = True
        if o == "--angpix":
            angpix = float(a)
        if o == "--output":
            output = str(a)
        if o == "--j":
            nr_threads = int(a)
        if o == "--np":
            nr_procs = int(a)

    if len(args) == 0:
        error("Error: No input file given.")

    input_star_filename = str(args[0])

    if not os.path.exists(input_star_filename):
        error("Error: Input file '%s' not found." % input_star_filename)

    if part_image_size < 0:
        error("Error: Particle image size not given.")

    if subpart_image_size < 0:
        error("Error: Sub-particle image size not given.")

    if angpix < 0:
        error("Error: Pixel size not given.")

    matrices = []
    particles = []
    parameters = []
    subparticles = []
    all_subparticles = []
    all_subparticles_subtracted = []

    if input_cmm:
        subparticle_vector_list = vectors_from_cmm(input_cmm, angpix)
    else:
        subparticle_vector_list = vectors_from_string(vectors_string)

    if distances_string:
        # Change distances from A to pixel units
        subparticle_distances = [float(x) / angpix for x in distances_string.split(',')]

        if len(subparticle_distances) != len(subparticle_vector_list):
            error("Error: The number of distances doesn't match the number of vectors!")

        for vector, distance in izip(subparticle_vector_list, subparticle_distances):
            if distance > 0:
                vector.set_distance(distance)
            else:
                vector.compute_distance()
    else:
        for vector in subparticle_vector_list:
            vector.compute_distance()


    print "Creating subparticles using:"

    for subparticle_vector in subparticle_vector_list:
        print "Vector: ",
        subparticle_vector.print_vector()
        print ""
        print "Length: %.2f pixels" % subparticle_vector.distance
    print ""

    path = output + "/"
    output_subt = output + "_subtracted"

    run_command("mkdir -p " + output, "/dev/null")

    if split_stacks:
        print "Creating and splitting the particle stack..."
        print " Creating a normal stack from which nothing is subtracted..."
        run_command("beditimg -create " + str(part_image_size) + "," + str(part_image_size) + "," + str(part_image_size) + " -fill 0 dummy_mask.mrc", "")
        run_command("relion_project --i dummy_mask.mrc --o " + path + "particles" + " --ang " + input_star_filename + " --subtract_exp --angpix " + str(angpix), "")
        run_command("rm -f dummy_mask.mrc", "")
        run_command("bsplit -digits 6 -first 1 " + path + "particles" + ".mrcs:mrc " + path + "particles" + ".mrc", "")
        print "Finished splitting the particle stack!"
        print " "

    if subtract_masked_map:
        print "Creating and splitting the particle stack..."
        print " Creating a stack from which the projections of the masked particle have been subtracted..."
        run_command("relion_project --i " + masked_map + " --o " + path + "particles_subtracted" + " --ang " + input_star_filename + " --subtract_exp --ctf --angpix " + str(angpix))
        print "Finished subtracting the model projections!"
        print " "

        print "Splitting the subtracted particle stack..."
        run_command("bsplit -digits 6 -first 1 " + path + "particles_subtracted" + ".mrcs:mrc " + path + "particles_subtracted" + ".mrc", "")
        print "Finished splitting the subtracted particle stack!"
        print " "

    print "Creating subparticles..."

    # setup toolbar   
    toolbar_width = 70
    sys.stdout.write("%s[oo]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width))
    c = 0
    timer = 0.01

    input_star = open(path + "particles.star", "r")
    particles, parameters = read_star(input_star)
    input_star.close()

    relion_create_symmetry_ops_file(symString, "relion_symops.tmp")
    symmetry_matrices = matrix_from_symmetry_ops_file("relion_symops.tmp")
    run_command("rm relion_symops.tmp", "")

    nparticle = 0
    for particle in particles:
        subparticles = create_subparticles(particle, symmetry_matrices, subparticle_vector_list, part_image_size, relax_symmetry, randomize, "subparticles", unique, len(all_subparticles), align_subparticles)

        # to preserve numbering, ALL sub-particles are written to STAR files before filtering
        index = particle.rlnImageName[0:6]
        star_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + ".star"

        for subparticle in subparticles:
            particle_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + ".mrc"
            subparticle.setrlnMicrographName(particle_filename)
            subparticle_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + "_subparticles.mrcs"
            subparticle.setrlnImageName(subparticle.rlnImageName[0:7] + subparticle_filename)

        if subtract_masked_map:
            subparticles_subtracted = []
            star_filename_subtracted = splitext(particle.rlnImageName[7:])[0] + "_subtracted_" + index + ".star"
            for subparticle in subparticles:
                subparticle_new = copy.deepcopy(subparticle)
                subparticle_new.rlnImageName = subparticle_new.rlnImageName[:-24] + "subtracted" + subparticle_new.rlnImageName[-25:]
                subparticles_subtracted.append(subparticle_new)

        if star:
            create_star(subparticles, star_filename)
            if subtract_masked_map:
                create_star(subparticles_subtracted, star_filename_subtracted)

        if side > 0:
            subparticles = filter_subparticles_side(subparticles, side)
        if top > 0:
            subparticles = filter_subparticles_top(subparticles, top)
        if mindist > 0:
            subparticles = filter_subparticles_mindist(subparticles, mindist)

        if (side > 0 or top > 0 or mindist > 0) and subtract_masked_map:
            subparticles_subtracted = []
            for subparticle in subparticles:
                subparticle_new = copy.deepcopy(subparticle)
                subparticle_new.rlnImageName = subparticle_new.rlnImageName[:-24] + "subtracted" + subparticle_new.rlnImageName[-25:]
                subparticles_subtracted.append(subparticle_new)

        all_subparticles.extend(subparticles)
        if subtract_masked_map:
            all_subparticles_subtracted.extend(subparticles_subtracted)

        if nparticle == int(len(particles) * timer):
            sys.stdout.write("\b" * (c + 8))
            sys.stdout.write("." * c)
            sys.stdout.write("~~(,_,\">")
            sys.stdout.flush()
            timer = timer + 0.01
            c = c + 1

        nparticle = nparticle + 1
    sys.stdout.write("\n")

    print "Finished creating the subparticles!"
    print " "

    if extract_subparticles:
        print "Extracting the subparticles..."
        if nr_procs == 1:
            cmd = 'relion_preprocess '
        else:
            cmd = 'mpirun -np %s relion_preprocess_mpi ' % nr_procs

        suffix = '_substracted' if subtract_masked_map else ''

        def run_extract(suffix=''):
            args = '--extract --o subparticles --extract_size %s --coord_files "%sparticles%s_??????.star"' % (subpart_image_size, path, suffix)
            run_command(cmd + args)
            run_command('mv subparticles.star %s%s_preprocess.star' % (output, suffix))

        run_extract() # Run extraction without substracted density

        if subtract_masked_map:
            run_extract('_subtracted')

        run_command("mv Particles/%s/* %s/" % (output, output))
        run_command("rmdir Particles/" + output)
        print "Finished extracting the sub-particles!\n"


    write_star(all_subparticles, parameters, output + ".star")

    if subtract_masked_map:
        write_star(all_subparticles_subtracted, parameters, output_subt + ".star")

    print "The output files have been written!"
    print "   Parameters for subparticles: *** " + output + ".star **"

    if subtract_masked_map:
        print "   Parameters for subparticles after subtractions: *** " + output_subt + ".star ***"

    print " "

    sys.exit(0)

if __name__ == "__main__":    
    main()

