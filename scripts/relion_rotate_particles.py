#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:  Serban Ilca (serban@strubi.ox.ac.uk)
# *           Juha T. Huiskonen (juha@strubi.ox.ac.uk)
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
            description="Rotates all particles in the input STAR file.")
        required = self.parser.add_argument_group('required arguments')
        add = self.parser.add_argument  # shortcut
        addr = required.add_argument

        add('input_star', help="Input STAR filename with particles.")
        add('--vector', default=None, help="Vector defining the additional rotation of the particles (x, y, z).")
        add('--angles', default=[], help="Euler angles defining the additional rotation of the particles (rot, tilt, psi).")
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

        print "Creating rotated particles..."

        md = MetaData(args.input_star)

        if args.vector and len(args.angles) != 0:
            print "Please only provide a vector or a triplet of Euler angles for the particle rotation."
            sys.exit(0)
        elif args.vector:
            vectors = load_vectors(None, args.vector, None, 1)
            rot_matrix = vectors[0].matrix()
        elif len(args.angles) != 0:
            angles = args.angles.split(',')
            if len(angles) != 3:
                print "Please provide exactly 3 Euler angles for the particle rotation."
                sys.exit(0)
            else:
                rot_rot = math.radians(float(angles[0]))
                rot_tilt = math.radians(float(angles[1]))
                rot_psi = math.radians(float(angles[2]))
                print rot_rot
                print rot_tilt
                print rot_psi
                rot_matrix = matrix_from_euler(rot_rot, rot_tilt, rot_psi)
        else:
            print "Please provide a vector or a triplet of Euler angles for the particle rotation."
            sys.exit(0)

        for particle in md:


            angles_to_radians(particle)
            new_particle = copy.deepcopy(particle)
            rot = particle.rlnAngleRot
            tilt = particle.rlnAngleTilt
            psi = particle.rlnAnglePsi

            matrix_particle = matrix_from_euler(rot, tilt, psi)

            m = matrix_multiply(matrix_particle, matrix_transpose(rot_matrix))
            rotNew, tiltNew, psiNew = euler_from_matrix(m)

            particle.rlnAngleRot = rotNew
            particle.rlnAngleTilt = tiltNew
            particle.rlnAnglePsi = psiNew

            print particle.rlnAngleRot
            print particle.rlnAngleTilt
            print particle.rlnAnglePsi
            sys.exit(0)

        md.write(args.output)

        print "All done!"
        print " "

if __name__ == "__main__":

    CreateSymmetryRelatedParticles().main()
