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

import math
import sys
import random
from itertools import izip

from matrix3 import *
from vector3 import *
from euler import *
from os.path import splitext

from pyrelion import MetaData


def within_mindist(p1, p2, mindist):
    """ Returns True if two particles are closer to each other
    than the given distance in the projection. """

    x1 = p1.rlnCoordinateX
    y1 = p1.rlnCoordinateY
    x2 = p2.rlnCoordinateX
    y2 = p2.rlnCoordinateY
    distance_sqr = (x1 - x2) ** 2 + (y1 - y2) ** 2
    mindist_sqr = mindist ** 2

    return distance_sqr < mindist_sqr


def within_unique(p1, p2, unique):
    """ Returns True if two particles are closer to each other
    than the given angular distance. """

    v1 = vector_from_two_eulers(p1.rlnAnglePsi, p1.rlnAngleTilt)
    v2 = vector_from_two_eulers(p2.rlnAnglePsi, p2.rlnAngleTilt)

    dp = dot_product(v1, v2) / (v1.length() * v2.length())

    if dp < -1:
        dp = -1.000

    if dp > 1:
        dp = 1.000

    angle = math.acos(dp)

    return angle <= math.radians(unique)


def filter_unique(subparticles, subpart, unique):
    """ Return True if subpart is not close to any other subparticle
        by unique (angular distance).
        For this function we assume that subpart is not contained
        inside."""
    for sp in subparticles:
        if within_unique(sp, subpart, unique):
            return False

    return True


def filter_mindist(subparticles, subpart, mindist):
    """ Return True if subpart is not close to any other subparticle
    by mindist. """
    for sp in subparticles:
        if (sp.rlnImageName[:6] != subpart.rlnImageName[:6] and
                within_mindist(sp, subpart, mindist)):
            return False

    return True


def filter_side(subpart, side):
    return (abs(abs(subpart.rlnAngleTilt) - radians(90)) < side)


def filter_top(subpart, top):
    return (abs(abs(subpart.rlnAngleTilt) - radians(180)) < top)


def filter_subparticles(subparticles, filters):
    return [sp for sp in subparticles
            if all(f(subparticles, sp) for f in filters)]


def angles_to_radians(particle):
    """ Convert the particle angles to radians. """
    particle.rlnAnglePsi = math.radians(particle.rlnAnglePsi)
    particle.rlnAngleTilt = math.radians(particle.rlnAngleTilt)
    particle.rlnAngleRot = math.radians(particle.rlnAngleRot)


def angles_to_degrees(particle):
    """ Convert the particle angles to radians. """
    particle.rlnAnglePsi = math.degrees(particle.rlnAnglePsi)
    particle.rlnAngleTilt = math.degrees(particle.rlnAngleTilt)
    particle.rlnAngleRot = math.degrees(particle.rlnAngleRot)


def create_subparticles(particle, symmetry_matrices, subparticle_vector_list,
                        part_image_size, randomize, output,
                        unique, subparticles_total, align_subparticles,
                        subtract_masked_map, do_create_star, filters):
    """ Obtain all subparticles from a given particle and set
    the properties of each such subparticle. """

    part_index = particle.rlnImageName[0:6]
    part_prefix = splitext(particle.rlnImageName[7:])[0]
    part_filename = "%s_%s.mrc" % (part_prefix, part_index)
    part_stack = "%s_%s_%s.mrcs" % (part_prefix, part_index, output)

    # We convert the particle angles to radian for further computations
    angles_to_radians(particle)

    # Euler angles that take particle to the orientation of the model
    rot = -particle.rlnAnglePsi
    tilt = -particle.rlnAngleTilt
    psi = -particle.rlnAngleRot
    matrix_particle = matrix_from_euler(rot, tilt, psi)

    subparticles = []
    subtracted = []
    subpart_id = 1
    subparticles_total += 1

    symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

    if randomize:
        # randomize the order of symmetry matrices, prevents preferred views
        random.shuffle(symmetry_matrix_ids)

    for subparticle_vector in subparticle_vector_list:
        matrix_from_subparticle_vector = subparticle_vector.matrix()

        for symmetry_matrix_id in symmetry_matrix_ids:
            # symmetry_matrix_id can be later written out to find out
            # which symmetry matrix created this subparticle
            symmetry_matrix = symmetry_matrices[symmetry_matrix_id - 1]

            subpart = particle.clone()

            m = matrix_multiply((matrix_multiply(matrix_from_subparticle_vector,
                                                 symmetry_matrix)), matrix_particle)

            if align_subparticles:
                rotNew, tiltNew, psiNew = euler_from_matrix(m)
            else:
                m2 = matrix_multiply(symmetry_matrix, matrix_particle)
                rotNew, tiltNew, psiNew = euler_from_matrix(m2)

            # save Euler angles that take the model to the orientation of the subparticle
            subpart.rlnAngleRot = -psiNew
            subpart.rlnAngleTilt = -tiltNew
            subpart.rlnAnglePsi = -rotNew

            # subparticle origin
            d = subparticle_vector.distance()
            x = -m.m[2][0] * d + particle.rlnOriginX
            y = -m.m[2][1] * d + particle.rlnOriginY
            z = -m.m[2][2] * d

            # modify the subparticle defocus paramaters by its z location
            if hasattr(particle, 'rlnDefocusU'):
                subpart.rlnDefocusU = particle.rlnDefocusU + z
                subpart.rlnDefocusV = particle.rlnDefocusV + z

            # save the subparticle coordinates (integer part) relative to the
            # user given image size and as a small shift in the origin (decimal part)
            x_d, x_i = math.modf(x)
            y_d, y_i = math.modf(y)
            subpart.rlnCoordinateX = int(part_image_size / 2) - x_i
            subpart.rlnCoordinateY = int(part_image_size / 2) - y_i
            subpart.rlnOriginX = -x_d
            subpart.rlnOriginY = -y_d

            overlaps = (unique >= 0 and
                        not filter_unique(subparticles, subpart, unique))

            if not overlaps:
                subpart.rlnImageName = "%06d@%s" % (subpart_id, part_stack)
                subpart.rlnParticleName = str(subparticles_total)
                subpart.rlnMicrographName = part_filename
                subparticles.append(subpart)
                subpart_id += 1
                subparticles_total += 1

    if subtract_masked_map:
        subtracted = clone_subtracted_subparticles(subparticles)

    # To preserve numbering, ALL sub-particles are written to STAR files before filtering
    if do_create_star:
        starfile = "%s_%s.star" % (part_prefix, part_index)
        create_star(subparticles, starfile)
        if subtract_masked_map:
            create_star(subtracted, add_suffix(starfile))

    if filters:
        subparticles = filter_subparticles(subparticles, filters)

        if subtract_masked_map:
            subtracted = clone_subtracted_subparticles(subparticles)

    return subparticles, subtracted


def create_symmetry_related_particles(particle, symmetry_matrices,
                                      keep_one=False):
    """ Return all related particles from the given symmetry matrices.
    If keep_one is True, randomly select only one of these equivalent
    particles.
    NOTE: Input particle should already contains angles in radians.
    """
    new_particles = []

    rot = -particle.rlnAnglePsi
    tilt = -particle.rlnAngleTilt
    psi = -particle.rlnAngleRot
    matrix_particle = matrix_from_euler(rot, tilt, psi)

    for symmetry_matrix in symmetry_matrices:
        m = matrix_multiply(symmetry_matrix, matrix_particle)
        rotNew, tiltNew, psiNew = euler_from_matrix(m)

        new_particle = particle.clone()
        new_particle.rlnAngleRot = -psiNew
        new_particle.rlnAngleTilt = -tiltNew
        new_particle.rlnAnglePsi = -rotNew
        angles_to_degrees(new_particle)
        new_particles.append(new_particle)

    if keep_one:
        new_particles = random.sample(new_particles, 1)

    return new_particles


def clone_subtracted_subparticles(subparticles):
    subparticles_subtracted = []

    for sp in subparticles:
        sp_new = sp.clone()
        sp_new.rlnImageName = add_suffix(sp.rlnImageName)
        sp_new.rlnMicrographName = add_suffix(sp.rlnMicrographName)
        subparticles_subtracted.append(sp_new)

    return subparticles_subtracted


def add_suffix(filename):
    return filename.replace('particles_', 'particles_subtracted_')


def create_star(subparticles, star_filename):
    """ Create a Relion style STAR file for extracting (using relion_preprocess)
    all the subparticles for a given particle. """

    md = MetaData()
    md.addLabels("rlnMicrographName", "rlnCoordinateX", "rlnCoordinateY",
                 "rlnImageName")
    md.addData(subparticles)
    md.write(star_filename)


def load_vectors(cmm_file, vectors_str, distances_str, angpix):
    """ Load subparticle vectors either from Chimera CMM file or from
    a vectors string. Distances can also be specified for each vector
    in the distances_str. """

    if cmm_file:
        subparticle_vector_list = vectors_from_cmm(cmm_file, angpix)
    else:
        subparticle_vector_list = vectors_from_string(vectors_str)

    if distances_str:
        # Change distances from A to pixel units
        subparticle_distances = [float(x) / angpix for x in
                                 distances_str.split(',')]

        if len(subparticle_distances) != len(subparticle_vector_list):
            raise Exception("Error: The number of distances does not match "
                            "the number of vectors!")

        for vector, distance in izip(subparticle_vector_list,
                                     subparticle_distances):
            if distance > 0:
                vector.set_distance(distance)
            else:
                vector.compute_distance()
    else:
        for vector in subparticle_vector_list:
            vector.compute_distance()

    print "Using vectors:"

    for subparticle_vector in subparticle_vector_list:
        print "Vector: ",
        subparticle_vector.normalize()
        subparticle_vector.compute_matrix()
        subparticle_vector.print_vector()
        print ""
        print "Length: %.2f pixels" % subparticle_vector.distance()
    print ""

    return subparticle_vector_list


def load_filters(side, top, mindist):
    """ Create some filters depending on the conditions imposed by the user.
    Each filter will return True if the subparticle will be kept in the
    subparticles list.
    """
    filters = []

    if side > 0:
        filters.append(lambda x, y: filter_side(y, side))

    if top > 0:
        filters.append(lambda x, y: filter_top(y, top))

    if mindist > 0:
        filters.append(lambda x, y: filter_mindist(x, y, mindist))

    return filters


def create_initial_stacks(input_star, particle_size, angpix,
                          split_stacks, masked_map, output):
    """ Create initial particle stacks (and star files) to be used
    for extraction of subparticles.
    If the subtracted_map is passed, another stack with subtracted
    particles will be created. """

    def create_stack(suffix, maskedFile, extra=''):
        print "Creating and splitting the particle stack..."
        if suffix:
            print(" Creating a stack from which the projections of the masked "
                  "particle have been subtracted...")
        else:
            print " Creating a normal stack from which nothing is subtracted..."

        outputParticles = "%s/particles%s" % (output, suffix)
        args = " --i %s --o %s --ang %s --subtract_exp --angpix %s " + extra
        run_command("relion_project" + args %
                    (maskedFile, outputParticles, input_star, angpix))

        run_command("bsplit -digits 6 -first 1 %s.mrcs:mrc %s.mrc"
                    % (outputParticles, outputParticles))

	run_command("rm -f %s.mrcs" % (outputParticles))

        print "Finished creating and splitting the particle stack!"
        print " "

    if split_stacks:
        maskedFile = 'dummy_mask.mrc'
        run_command("beditimg -create %s,%s,%s -fill 0 %s"
                    % (particle_size, particle_size, particle_size, maskedFile))
        create_stack('', maskedFile)
        run_command("rm -f %s" % maskedFile)

    if masked_map:
        create_stack('_subtracted', masked_map, '--ctf')


def extract_subparticles(subpart_size, np, masked_map, output):
    """ Extract subparticles images from each particle
    (Using 'relion_preprocess' as if the particle was a micrograph. """

    print "Extracting the subparticles..."

    if np == 1:
        cmd = 'relion_preprocess '
    else:
        cmd = 'mpirun -np %s relion_preprocess_mpi ' % np

    def run_extract(suffix=''):
        args = ('--extract --o subparticles --extract_size %s --coord_files '
                '"%s/particles%s_??????.star"') % (subpart_size, output, suffix)
        run_command(cmd + args)
        run_command("mv subparticles.star %s%s_preprocess.star"
                    % (output, suffix))
	run_command("rm -f %s/particles%s_??????.mrc" % (output, suffix))

    run_extract()  # Run extraction without subtracted density

    if masked_map:
        run_extract('_subtracted')

    run_command("mv Particles/%s/* %s/" % (output, output))
    run_command("rmdir Particles/" + output)

    print "Finished extracting the subparticles!\n"


def write_output_starfiles(labels, mdOut, mdOutSub, output):

    labels.extend(['rlnCoordinateX', 'rlnCoordinateY', 'rlnMicrographName'])
    print "Writing output STAR files."

    starfile1 = output + ".star"
    print "   Parameters for subparticles: \n      *** %s **" % starfile1
    # We convert back angles to degrees and write subparticles star file
    def _writeMd(md, starfile):
        for subpart in md:
            angles_to_degrees(subpart)
        md.addLabels(labels)
        md.write(starfile)

    _writeMd(mdOut, starfile1)

    if len(mdOutSub):
        starfile2 = starfile1.replace('.star', '_subtracted.star')
        print("   Parameters for subparticles after subtractions: \n"
              "      *** %s ***" % starfile2)
        _writeMd(mdOutSub, starfile2)

    print "The output files have been written!\n"


def reconstruct_subparticles(np, mdOutSub, output):

    print "Reconstructing subparticles."
    recfile = output + ".mrc"

    if np == 1:
        cmd = 'relion_reconstruct '
    else:
        cmd = 'mpirun -np %s relion_reconstruct_mpi' % np

    def run_reconstruct(suffix=''):
        args = ('--ctf --i %s%s.star --o subparticle%s_3d.mrc ') % (output, suffix, suffix)
        run_command(cmd + args)

    run_reconstruct()

    if len(mdOutSub):
        run_reconstruct('_subtracted')

    print "Subparticles have been reconsructed!\n"


def run_command(command, output=""):
    if not output:
        print "+++ " + command
        sys.stdout.flush()
        os.system(command)
    else:
        os.system(command + " > " + output)


class ProgressBar():
    """ Implements a simple command line progress bar"""

    def __init__(self, width, percent, total):
        # setup toolbar
        self.width = width
        sys.stdout.write("%s>->o" % (" " * width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (width))
        self.count = 0  # total count
        self.c = 0  # progress count
        self.percent = percent
        self.timer = percent
        self.total = total

    def notify(self):
        if self.count == int(self.total * self.timer):
            sys.stdout.write("\b" * (self.c + 8))
            sys.stdout.write("~" * self.c)
            sys.stdout.write(">))^)>")
            sys.stdout.flush()
            self.timer += self.percent
            self.c += 1

        self.count += 1
