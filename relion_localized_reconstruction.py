#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:  Serban Ilca
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
# *           J.M. de la Rosa Trevin
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


import os
import sys
import textwrap

from distutils import spawn
import argparse

from localrec import *
from pyrelion import MetaData


class LocalizedReconstruction():

    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
                    Localized reconstruction of subparticles. Normally the script is run in four steps:
                        1. Prepare particles.
                        2. Create subparticles.
                        3. Extract subparticles.
                        4. Reconstruct subparticles.
            '''))

        # Parameter groups
        general = self.parser.add_argument_group('General parameters')
        steps = self.parser.add_argument_group('Steps', 'Several steps can be combined in one run.')
        prepareParticles = self.parser.add_argument_group('Prepare particles')
        createSubparticles = self.parser.add_argument_group('Create subparticles')
        extractSubparticles = self.parser.add_argument_group('Extract subparticles')
        reconstructSubparticles = self.parser.add_argument_group('Reconstruct subparticles')

        # Shortcuts
        add = general.add_argument
        adds = steps.add_argument
        addpp = prepareParticles.add_argument
        addcs = createSubparticles.add_argument
        addes = extractSubparticles.add_argument
        addrs = reconstructSubparticles.add_argument

        # General parameters
        add('input_star', nargs='?', help="Input STAR filename with particles.")
        add('--output', default='subparticles',
              help="Output root for results.")

        # Parameters for "Steps" group
        adds('--prepare_particles', action='store_true',
              help="Prepare particles for extracting subparticles.")
        adds('--create_subparticles', action='store_true',
              help="Calculate the cooridnates and Euler angles for the subparticles.")
        adds('--extract_subparticles', action='store_true',
              help="Extract subparticles from particle images.")
        adds('--reconstruct_subparticles', action='store_true',
              help="Calculate a reconstruction of the subunit from the subparticles. "
                 "Subparticle coordinates are read from "
                 "[output]_subparticles.star and [output]_subparticles_subtracted.star")

        # Parameters for "Prepare particles" group
        addpp('--masked_map',
              help="Create another set of particles with partial signal subtraction using this map.")

        # Parameters for "Create subparticles" group
        addcs('--angpix', type=float, default=1, help="Pixel size (A). (default: 1)")
        addcs('--sym', default="C1", help="Symmetry of the particle. (default: C1)")
        addcs('--particle_size', type=int,
              help="Size of the particle box (pixels).")
        addcs('--randomize', action='store_true',
              help="Randomize the order of the symmetry matrices. "
                   "Useful for preventing preferred orientations.")
        addcs('--relax_symmetry', action='store_true',
              help="Create one random subparticle for each particle ")
        addcs('--vector', default="0,0,1", help="Vector defining the location of the subparticle. (default: 0,0,1)")
        addcs('--align_subparticles', action='store_true',
              help="Align subparticles to the standard orientation.")
        addcs('--length',
              help="Alternative length of the vector. Use to adjust the "
                 "subparticle center (A). (default: length of the given vector)")
        addcs('--cmm',
              help="A CMM file defining the location(s) of the subparticle(s) "
                 "(use instead of --vector). Coordinates should be in Angstrom.")
        addcs('--unique', type=float, default=-1,
              help="Keep only unique subparticles within angular distance "
                 "(useful to remove overlapping subparticles on symmetry axis).")
        addcs('--mindist', type=float, default=-1,
              help="Minimum distance between the subparticles in the image "
                 "(all overlapping ones will be discarded; pixels).")
        addcs('--side', type=float, default=-1,
              help="Keep only particles within specified angular distance from "
                 "side views (all others will be discarded; degrees).")
        addcs('--top', type=float, default=-1,
              help="Keep only particles within specified angular distance from "
                 "top views (all others will be discarded; degrees).")

        # Parameters for "Extract subparticles" group
        addes('--subparticle_size', type=int,
            help="Size of the subparticle box (pixels).")
        addes('--np', type=int, default=1, help="Number of MPI procs. (default: 1)")

        # Parameters for "Reconstruct subparticles" group
        addrs('--j', type=int, default=1, help="Number of threads.")
        addrs('--maxres', type=float, help="Maximum resolution of the reconstruction (A). (default: Nyquist)")
        addrs('--subsym', default="C1", help="Symmetry of the subparticle. (default: C1)")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print "\nError: " + '\n'.join(msgs)
        print " "
        sys.exit(2)

    def validate(self, args):
        # Check that required software is in PATH
        if not (spawn.find_executable("scipion")):
            self.error("Scipion not found.",
                       "Make sure Scipion is in $PATH.")

        if not (spawn.find_executable("relion_refine")):
            self.error("Relion not found.",
                       "Make sure Relion programs are in $PATH.")

        # Check that required parameters are given for each mode
        if args.prepare_particles or args.create_subparticles or args.extract_subparticles:
            if len(sys.argv) == 1:
                self.error("No input particles STAR file given.")
            if not os.path.exists(args.input_star):
                self.error("\nInput file '%s' not found." % args.input_star)

        if args.extract_subparticles:
            if not args.particle_size:
                self.error("Parameter --particle_size is required.")
            if not args.subparticle_size:
                self.error("Parameter --subparticle_size is required.")

        if args.reconstruct_subparticles:
            if not args.output:
                self.error("Parameter --output not specified. Cannot find STAR file with subparticles.")

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        # Validate input arguments and required software
        self.validate(args)

        run_command("mkdir -p " + args.output, "/dev/null")

        if args.prepare_particles:
            print "Preparing particles for extracting subparticles."
            create_initial_stacks(args.input_star, args.angpix, args.masked_map, args.output)
            print "\nFinished preparing the particles!\n"

        if args.create_subparticles:
            # Load subparticle vectors either from Chimera CMM file or from
            # command line (command and semi-colon separated)
            # Distances can also be specified to modify vector lengths
            subparticle_vector_list = load_vectors(args.cmm, args.vector,
                                                   args.length, args.angpix)

            particles_star = args.output + "/particles.star"

            if not os.path.exists(args.output + "/particles.star"):
                self.error("Input file '%s not found. "
                           "Run the script first with --prepare_particles option."
                           % particles_star)

            md = MetaData(particles_star)
            print "Creating subparticles..."

            # Initialize progress bar
            progressbar = ProgressBar(width=60, total=len(md))

            # Generate symmetry matrices with Relion convention
            symmetry_matrices = matrix_from_symmetry(args.sym)

            # Define some conditions to filter subparticles
            filters = load_filters(radians(args.side), radians(args.top), args.mindist)

            # Compute all subparticles (included subtracted if masked_map given)
            mdOut = MetaData()
            mdOutSub = MetaData()

            for particle in md:
                subparticles, subtracted = create_subparticles(particle,
                                                       symmetry_matrices,
                                                       subparticle_vector_list,
                                                       args.particle_size,
                                                       args.randomize,
                                                       args.output,
                                                       args.unique,
                                                       len(mdOut),
                                                       args.align_subparticles,
                                                       args.masked_map,
                                                       True,
                                                       filters)

                mdOut.addData(subparticles)
                mdOutSub.addData(subtracted)

                progressbar.notify()

            md.removeLabels('rlnOriginZ', 'rlnOriginalName')
            write_output_starfiles(md.getLabels(), mdOut, mdOutSub, args.output)

            print "\nFinished creating the subparticles!\n"

        if args.extract_subparticles:
            print "Extracting subparticles..."
            extract_subparticles(args.subparticle_size, args.np, args.masked_map, args.output, deleteParticles=True)
            print "\nFinished extracting the subparticles!\n"

        if args.reconstruct_subparticles:
            print "Reconstructing subparticles..."
            reconstruct_subparticles(args.j, args.output, args.maxres, args.subsym, args.angpix)
            print "\nFinished reconstructing the subparticles!\n"

        print "\nAll done!\n"

if __name__ == "__main__":
    LocalizedReconstruction().main()

