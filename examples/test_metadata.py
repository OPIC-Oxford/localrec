#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:  J. M. de la Rosa Trevin (delarosatrevin@gmail.com)
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

from pyrelion import MetaData

"""
This script will show some of the basic usages of the MetaData class
to operate with Relion star files.
"""


def test_micrographs():
    # Create a new metadata and read micrograph star file

    md = MetaData("micrographs.star")

    # Just print its content
    md.printStar()

    # We can query the number of rows (micrographs in this case)  in the metadata

    print md.size() # should be 15 for the provided  micrographs.star file


    # We can checkout the columns (or labels) of the star file

    print md.getLabels()

    # We can remove the CTF related columns and write a new star file

    md.removeLabels('rlnDefocusU', 'rlnDefocusV', 'rlnDefocusAngle')


    # We can fill a column value with a constant
    # (if the column does not exists, it will be added)

    md.setLabels(rlnCtfFigureOfMerit=0, rlnDetectorPixelSize=7,
                 rlnRandomSubset=999)

    # Finally let's write our new micrographs star file
    md.write('micrograhs_noctf.star')


def test_particles():
    # We can also read particles star files

    md = MetaData("particles.star")

    # For example, we could remove particles meeting some criteria

    good = [particle for particle in md
            if float(particle.rlnParticleSelectZScore) > 0.5]

    # There are some Relion star labels that are read as string
    # so we need to care to convert them to the proper Python data type

    md.setData(good)

    print md.size() # should be 36 for the provided star file

    md.write('particles_good.star')


def test_ctf():
    # Copy all the values from the micrograph star file to particles

    input_filename = "particles.star"
    ctf_filename = "micrographs.star"

    inputMd = MetaData(input_filename)
    ctfMd = MetaData(ctf_filename)
    # Add labels from CTF star
    ctfLabels =  ctfMd.getLabels()
    inputMd.addLabels(ctfLabels)
    # Index micrographs by name in a dictionary to retrieve from particles
    micDict = {}

    for mic in ctfMd:
        micDict[mic.rlnMicrographName] = mic

    for particle in inputMd:
        mic = micDict[particle.rlnMicrographName]
        particle.copyValues(mic, *ctfLabels)

    inputMd.write('particles_ctf.star')



if __name__ == "__main__":

    test_micrographs()

    test_particles()

    test_ctf()