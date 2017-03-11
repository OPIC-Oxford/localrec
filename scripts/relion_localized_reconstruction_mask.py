#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:  Serban Ilca (serban@strubi.ox.ac.uk)
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

import localrec
from pyworkflow.em import runProgram
import pyworkflow.utils as pwutils


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

        add('--edge', type=int, default=3,
            help="Mask edge in pixels. Default: 3 pixels.")

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
        # self.validate(args)

        # Load subparticle vectors either from Chimera CMM file or from
        # command line (command and semi-colon separated)
        # Distances can also be specified to modify vector lengths
        subparticle_vector_list = localrec.load_vectors(args.cmm, args.vector,
                                                        args.length, args.angpix)

        size = args.particle_size
        half = size / 2
        v = subparticle_vector_list[0]
        d = v.distance()

        # Bsoft's way
        mask_x = v.x() * d + half
        mask_y = v.y() * d + half
        mask_z = v.z() * d + half

        print "Creating a mask at the end point of the vector..."

        runProgram("beditimg",
                   "-create %s,%s,%s -sphere %s,%s,%s,%s -edge %d -fill 1 mask.mrc"
                    % (size, size, size, mask_x, mask_y, mask_z,
                       args.radius, args.edge))

        if args.sym:
            print "Symmetrizing the mask..."
            runProgram("bsym", "-origin %s,%s,%s -sym %s mask.mrc mask_%s.mrc"
                       % (half, half, half, args.sym, args.sym))

        print "Inverting the mask(s)..."

        def invertMask(fn):
            runProgram("bar", "-multiply -1 -add 1 %s %s" % (fn, fn))

        invertMask("mask.mrc")

        if args.sym:
            invertMask("mask_%s.mrc" % args.sym)

        # Xmipp's way.
        # Let's create a phantom description file with the sphere in the
        # desired position and symmetrize it. We can set the background to 1
        # and subtract the sphere value (equivalent to mask invert)
        maskFn = 'xmask.mrc'
        v.scale(d)
        phantomFn = 'phantom.descr'
        phantomFile = open(phantomFn, 'w')
        phantomFile.write("%d %d %d 1 1\n" % (size, size, size))
        sym_matrices = localrec.matrix_from_symmetry(args.sym)

        for m in sym_matrices:
            vm = localrec.matrix_product(m, v)
            sx, sy, sz = vm.data()
            phantomFile.write("sph + -1 %d %d %d %d\n" % (sx, sy, sz, args.radius))

        phantomFile.close()

        runProgram('xmipp_phantom_create', '%s -o %s' % (phantomFn, maskFn))
        pwutils.cleanPath(phantomFn)

        runProgram('xmipp_transform_filter',
                   '%s --fourier real_gaussian %d ' % (maskFn, args.edge))


        print "All done!"
        print " "
        print "Use the following command to apply the mask on your particle:"
        print "bop -multiply 1,0 mask.mrc particle.mrc particle_masked.mrc"
        print "The masked particle can be used for partial signal subtraction for localized reconstruction."
        print " "


if __name__ == "__main__":    
    LocalizedReconstructionMask().main()

