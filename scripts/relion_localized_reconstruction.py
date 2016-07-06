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

from distutils import spawn
import argparse

from localized_reconstruction import *
from metadata import MetaData




class LocalizedReconstruction():

    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Calculates the coordinates and Euler angles for the "
                        "subparticles defined by a vector and symmetry point "
                        "group.")
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
        add('--unique', type=float, default=-1,
            help="Keep only unique subparticles within angular distance "
                 "(useful to remove overlapping sub-particles on symmetry axis).")
        add('--mindist', type=float, default=-1,
            help="Minimum distance between the subparticles in the image "
                 "(all overlapping ones will be discarded; pixels).")
        add('--side', type=float, default=-1,
            help="Keep only particles within specified angular distance from "
                 "side views (all others will be discarded; degrees).")
        add('--top', type=float, default=-1,
            help="Keep only particles within specified angular distance from "
                 "top views (all others will be discarded; degrees).")
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

    def validate(self, args):
        if args.extract_subparticles and not (spawn.find_executable("bimg")):
            self.error("Error: Bsoft not found.",
                       "Make sure Bsoft programs are in $PATH.")

        if not (spawn.find_executable("relion_refine")):
            self.error("Error: Relion not found.",
                       "Make sure Relion programs are in $PATH.")

        if len(sys.argv) == 1:
            self.error("Error: No input file given.")

        if not os.path.exists(args.input_star):
            self.error("Error: Input file '%s' not found."
                       % args.input_star)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        # Validate input arguments and required software (Relion and Bsoft)
        self.validate(args)

        particle_size = args.particle_size
        subpart_image_size = args.subparticle_size
        output = args.output
        subtract_masked_map = args.masked_map is not None

        # Load subparticle vectors either from Chimera CMM file or from
        # command line (command and semi-colon separated)
        # Distances can also be specified to modify vector lengths
        subparticle_vector_list = load_vectors(args.cmm, args.vector,
                                               args.length, args.angpix)

        run_command("mkdir -p " + output, "/dev/null")

        if args.extract_subparticles:
            create_initial_stacks(args.input_star, particle_size, args.angpix,
                                  args.split_stacks, args.masked_map, output)
            particles_star = output + "/particles.star"
        else:
            particles_star = args.input_star

        print "Creating subparticles..."


        md = MetaData(particles_star)

        # Initialize progress bar
        progressbar = ProgressBar(width=70, percent=0.01, total=len(md))

        # Generate symmetry matrices with Relion convention
        symmetry_matrices = matrix_from_symmetry(args.sym)

        # Define some conditions to filter subparticles
        filters = load_filters(args.side, args.top, args.mindist)

        # Compute all subparticles (included subtracted if masked_map given)
        mdOut = MetaData()
        mdOutSub = MetaData()

        for particle in md:
            subparticles, subtracted = create_subparticles(particle,
                                                   symmetry_matrices,
                                                   subparticle_vector_list,
                                                   particle_size,
                                                   args.relax_symmetry,
                                                   args.randomize,
                                                   "subparticles", args.unique,
                                                   len(mdOut),
                                                   args.align_subparticles,
                                                   subtract_masked_map,
                                                   args.create_star, filters)

            mdOut.addData(subparticles)
            mdOutSub.addData(subtracted)

            progressbar.notify()

        print "\nFinished creating the subparticles!\n"

        if args.extract_subparticles:
            extract_subparticles(subpart_image_size, args.np,
                                 args.masked_map, output)

        write_output_starfiles(md.getLabels(), mdOut, mdOutSub, output)


if __name__ == "__main__":    
    LocalizedReconstruction().main()

