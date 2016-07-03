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
import argparse



class LocalizedReconstruction():
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Calculates the coordinates and Euler angles for the "
                        "subparticles defined by a vector and symmetry point group.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('input_star', help="Input star filename with particles.")
        add('--split_stacks', action='store_true',
            help="Split particle stacks (needs to be done once).")
        add('--masked_map',
            help="Map with density that to be subtracted from particle images.")
        add('--create_star', action='store_true',
            help="Create new STAR files for extracting sub-particles.")
        add('--extract_subparticles', action='store_true',
            help="Extract sub-particles from particle images.")
        addr('--angpix', type=float, help="Pixel size (A).", required=True)
        add('--sym', help="Symmetry of the particle.")
        addr('--particle_size', type=int, required=True,
            help="Size of the particle box (pixels).")
        addr('--subparticle_size', type=int, required=True,
            help="Size of the sub-particle box (pixels).")
        add('--randomize', action='store_true',
            help="Randomize the order of the symmetry matrices. \n"
                 "Useful for preventing preferred orientations (default: not).")
        add('--relax_symmetry', action='store_true',
            help="Create one random subparticle for each particle "
                 "(default: all symmetry related subparticles).")
        add('--vector', help="Vector defining the location of the subparticle.")
        add('--align_subparticles', action='store_true',
            help="Align subparticles to the standard orientation.")
        add('--length',
            help="Alternative length of the vector. Use to adjust the "
                 "sub-particle center (default: length of the given "
                 "vector; A).")
        add('--cmm',
            help="A CMM file defining the location(s) of the subparticle(s) "
                 "(use instead of --vector). Coordinates in Angstrom.")
        add('--unique', type=float, default=1,
            help="Keep only unique subparticles within angular distance "
                 "(useful to remove overlapping sub-particles on symmetry axis).")
        add('--mindist', type=float, default=5,
            help="Minimum distance between the subparticles in the image "
                 "(all overlapping ones will be discarded; pixels).")
        add('--side', type=float, default=25,
            help="Keep only particles within specified angular distance from side "
                 "views (all others will be discarded; degrees).")
        add('--top', type=float, default=25,
            help="Keep only particles within specified angular distance from top "
                 "views (all others will be discarded; degrees).")
        add('--output', default='subparticles',
            help="Output root for results.")
        add('--j', type=int, default=8, help="Number of threads.")
        add('--np', type=int, default=4, help="Number of MPI procs.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print "Error: " + '\n'.join(msgs)
        print " "
        sys.exit(2)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        if not (spawn.find_executable("bimg")):
            self.error("Error: Bsoft not found.",
                  "Make sure Bsoft programs are in $PATH.")

        if not (spawn.find_executable("relion_refine")):
            self.error("Error: Relion not found.",
                  "Make sure Relion programs are in $PATH.")

        if len(sys.argv) == 1:
            self.error("Error: No input file given.")

        input_star_filename = args.input_star

        if not os.path.exists(input_star_filename):
            self.error("Error: Input file '%s' not found." % input_star_filename)

        angpix = args.angpix
        part_image_size = args.particle_size
        subpart_image_size = args.subparticle_size
        output = args.output
        subtract_masked_map = args.masked_map is not None
        side = args.side
        top = args.top
        unique = args.unique
        mindist = args.mindist

        all_subparticles = []
        all_subparticles_subtracted = []

        if args.cmm:
            subparticle_vector_list = vectors_from_cmm(args.cmm, angpix)
        else:
            subparticle_vector_list = vectors_from_string(args.vector)

        if args.length:
            # Change distances from A to pixel units
            subparticle_distances = [float(x) / angpix for x in args.length.split(',')]

            if len(subparticle_distances) != len(subparticle_vector_list):
                self.error("Error: The number of distances doesn't match the number of vectors!")

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

        def create_stack(suffix, maskedFile, extraArgs=''):
            print "Creating and splitting the particle stack..."
            if suffix:
                print " Creating a stack from which the projections of the masked particle have been subtracted..."
            else:
                print " Creating a normal stack from which nothing is subtracted..."

            outputParticles = "%sparticles%s" % (path, suffix)
            run_command("relion_project --i %s --o %s --ang %s --subtract_exp --angpix %s %s"
                        % (maskedFile, outputParticles, input_star_filename, angpix, extraArgs))

            run_command("bsplit -digits 6 -first 1 %s.mrcs:mrc %s.mrc"
                        % (outputParticles, outputParticles))

            print "Finished splitting the particle stack!"
            print " "

        if args.split_stacks:
            maskedFile = 'dummy_mask.mrc'
            run_command("beditimg -create %s,%s,%s -fill 0 %s"
                        % (part_image_size, part_image_size, part_image_size, maskedFile))
            create_stack('', maskedFile)
            run_command("rm -f %s" % maskedFile)

        if subtract_masked_map:
            create_stack('_subtracted', args.masked_map, '--ctf')


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

        # Generate symmetry matrices with Relion convention
        symmetry_matrices = matrix_from_symmetry(args.sym)

        nparticle = 0

        # Define some conditions to filter subparticles
        filters = []

        if side > 0:
            filters.append(lambda x, y: filter_side(y, side))

        if top > 0:
            filters.append(lambda x, y: filter_top(y, top))

        if mindist > 0:
            filters.append(lambda x, y: filter_mindist(x, y, mindist))

        for particle in particles:
            subparticles, subtracted = create_subparticles(particle,
                                               symmetry_matrices,
                                               subparticle_vector_list,
                                               part_image_size,
                                               args.relax_symmetry,
                                               args.randomize,
                                               "subparticles", unique,
                                               len(all_subparticles),
                                               args.align_subparticles,
                                               subtract_masked_map,
                                               create_star, filters)

            all_subparticles.extend(subparticles)
            all_subparticles_subtracted.extend(subtracted)

            if nparticle == int(len(particles) * timer):
                sys.stdout.write("\b" * (c + 8))
                sys.stdout.write("." * c)
                sys.stdout.write("~~(,_,\">")
                sys.stdout.flush()
                timer += 0.01
                c += 1

            nparticle += 1

        sys.stdout.write("\n")

        print "Finished creating the subparticles!"
        print " "

        if args.extract_subparticles:
            print "Extracting the subparticles..."
            if args.np == 1:
                cmd = 'relion_preprocess '
            else:
                cmd = 'mpirun -np %s relion_preprocess_mpi ' % args.np

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


if __name__ == "__main__":    
    LocalizedReconstruction().main()

