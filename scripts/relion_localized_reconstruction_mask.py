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

from distutils import spawn
import argparse

from localrec import *
from pyrelion import MetaData

class LocalizedReconstructionMask():

    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Creates a mask for localized reconstruction.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        addr('--angpix', type=float, help="Pixel size (A).", required=True)
        add('--sym', help="Apply symmetry on the mask.")
        addr('--particle_size', type=int, required=True,
            help="Size of the particle box (pixels).")
        addr('--radius', type=int, required=True,
            help="Radius of the mask (pixels).")
        add('--vector', help="Vector defining the location of the subparticle.")
        add('--length',
            help="Alternative length of the vector. Use to adjust the "
                 "subparticle center (default: length of the given "
                 "vector; A).")
        add('--cmm',
            help="A CMM file defining the location of the subparticle "
                 "(use instead of --vector). Coordinates in Angstrom."
		 "Only one vector is supported.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print "Error: " + '\n'.join(msgs)
        print " "
        sys.exit(2)

    def validate(self, args):
        if not (spawn.find_executable("bimg")):
            self.error("Error: Bsoft not found.",
                       "Make sure Bsoft programs are in $PATH.")

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        # Validate input arguments and required software (Bsoft)
        self.validate(args)

        # Load subparticle vectors either from Chimera CMM file or from
        # command line (command and semi-colon separated)
        # Distances can also be specified to modify vector lengths
        subparticle_vector_list = load_vectors(args.cmm, args.vector,
                                               args.length, args.angpix)

        mask_x = subparticle_vector_list[0].x()*subparticle_vector_list[0].distance() + args.particle_size / 2
        mask_y = subparticle_vector_list[0].y()*subparticle_vector_list[0].distance() + args.particle_size / 2
        mask_z = subparticle_vector_list[0].z()*subparticle_vector_list[0].distance() + args.particle_size / 2

        print "Creating a mask at the end point of the vector..."

        run_command("beditimg -create %s,%s,%s -sphere %s,%s,%s,%s -edge 3 -fill 1 mask.mrc"
                    % (args.particle_size, args.particle_size, args.particle_size,
                       mask_x, mask_y, mask_z, args.radius))

	if args.sym:
	        print "Symmetrizing the mask..."
		run_command("bsym -origin %s,%s,%s -sym %s mask.mrc mask_%s.mrc" %
			(int(args.particle_size/2), int(args.particle_size/2), int(args.particle_size/2),
			args.sym, args.sym))

        print "Inverting the mask(s)..."
	run_command("bar -multiply -1 -add 1 mask.mrc mask.mrc")

	if args.sym:
		run_command("bar -multiply -1 -add 1 mask_%s.mrc mask_%s.mrc" % (args.sym, args.sym))

	print "All done!"
	print " "
 	print "Use the following command to apply the mask on your particle:"
	print "bop -multiply 1,0 mask.mrc particle.mrc particle_masked.mrc"
	print "The masked particle can be used for partial signal subtraction for localized reconstruction." 
	print " "

if __name__ == "__main__":    
    LocalizedReconstructionMask().main()

