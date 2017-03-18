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
import argparse
from collections import OrderedDict


class RemoveOverlappingParticles():
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Removes overlapping particles within a group of particles read from the input STAR file."
                        "One particle is left and the rest are removed.")
        add = self.parser.add_argument
        add('input_star', help="Input STAR filename with particles.")
        add('--mindist', type=float, default=-1,
            help="Minimum distance allowed between particles (pixels).")
        add('--originalParticles', action='store_true',
            help="Use original particle names to define particle groups (default: micrograph names).")
        add('--output', help="Output STAR filename.")

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

    def get_micrographs(self, md):
        micrographs = []
        for particle in md:
            micrographs.append(particle.rlnMicrographName)
        return list(OrderedDict.fromkeys(micrographs))

    def get_originalParticleNames(self, md):
        originalParticleNames = []
        for particle in md:
            originalParticleNames.append(particle.rlnOriginalParticleName)
        return list(OrderedDict.fromkeys(originalParticleNames))

    def get_particles(self, md, micrograph):
        particles = []
        for particle in md:
            if (particle.rlnMicrographName == micrograph) or (particle.rlnOriginalParticleName == micrograph):
                particles.append(particle)
        return particles

    def overlaps(self, p1, p2, mindist):
        ''' returns True if overlap is found '''
        x1 = p1.rlnCoordinateX + p1.rlnOriginX
        y1 = p1.rlnCoordinateY + p1.rlnOriginY
        x2 = p2.rlnCoordinateX + p2.rlnOriginX
        y2 = p2.rlnCoordinateY + p2.rlnOriginY
        distance_sqr = (x1 - x2) ** 2 + (y1 - y2) ** 2
        mindist_sqr = mindist ** 2
        return distance_sqr < mindist_sqr

    def remove_overlapping_particles(self, particles, mindist):

        for i, particle in enumerate(particles):
            overlap_found = False
            for j, particle2 in enumerate(particles[i + 1:]):
                if mindist > 0:
                    if self.overlaps(particle, particle2, mindist):
                        particles.remove(particle2)

        return particles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        print "Remove overlapping particles..."

        md = MetaData(args.input_star)
        mdOut = MetaData()
        mdOut.addLabels(md.getLabels())

        new_particles = []

        # get list of tags 'micrographs' but can be anything that defines
        # a set of particles to compare
        if args.originalParticles:
            micrographs = self.get_originalParticleNames(md)
        else:
            micrographs = self.get_micrographs(md)

        for micrograph in micrographs:
            print "Processing micrograph %s..." % micrograph,
            particles = self.get_particles(md, micrograph)
            particles_unique = self.remove_overlapping_particles(particles, args.mindist)
            new_particles.extend(particles_unique)
            print "%s particle(s) kept." % len(particles_unique)

        mdOut.addData(new_particles)
        mdOut.write(args.output)

        print "All done!"


if __name__ == "__main__":
    RemoveOverlappingParticles().main()