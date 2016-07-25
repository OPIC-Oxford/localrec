#!/usr/bin/env python

# **************************************************************************
# *
# * Author:  Juha T. Huiskonen (juha@strubi.ox.ac.uk)
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

from pyrelion import MetaData
from localrec import *
import argparse


class CreateSymmetryRelatedParticles():
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Creates all symmetry copies for each particle as new lines in the input STAR file.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('input_star', help="Input STAR filename with particles.")
        add('--sym', default="C1", help="Symmetry of the particle.")
        addr('--output', required=True, help="Output STAR filename.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print "Error: " + '\n'.join(msgs)
        print " "
        sys.exit(2)

    def validate(self, args):
        if len(sys.argv) == 1:
            self.error("Error: No input file given.")

        if not os.path.exists(args.input_star):
            self.error("Error: Input file '%s' not found."
                       % args.input_star)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        print "Creating symmetry related particles..."

        md = MetaData(args.input_star)
        mdOut = MetaData()
        mdOut.addLabels(md.getLabels())
        
        new_particles = []

        symmetry_matrices = matrix_from_symmetry(args.sym)

        for particle in md:
            angles_to_radians(particle)
            new_particles.extend(create_symmetry_related_particles(particle,
                                                                   symmetry_matrices))
        mdOut.addData(new_particles)
        mdOut.write(args.output)

        print "All done!"


if __name__ == "__main__":

    CreateSymmetryRelatedParticles().main()
