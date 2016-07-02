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


def usage():
    print ""
    print "Usage: relion_create_subparticles.py [parameters] input.star"
    print "Calculates the coordinates and Euler angles for the subparticles defined by a vector and symmetry point group."
    print ""
    print "Parameters: "
    print "--split_stacks                   Split particle stacks (needs to be done once)."
    print "--masked_map file.mrc            Map with density that to be subtracted from particle images."
    print "--create_star                    Create new STAR files for extracting sub-particles."
    print "--extract_subparticles           Extract sub-particles from particle images."  
    print "--angpix 1.35                    Pixel size (A)."
    print "--sym I1                         Symmetry of the particle."
    print "--particle_size 512              Size of the particle box (pixels)."
    print "--subparticle_size 128           Size of the sub-particle box (pixels)."
    print "--randomize                      Randomize the order of the symmetry matrices. Useful for preventing preferred orientations (default: not)."
    print "--relax_symmetry                 Create one random subparticle for each particle (default: all symmetry related subparticles)."
    print "--vector 0,0,1                   Vector defining the location of the subparticle."
    print "--align_subparticles             Align subparticles to the standard orientation."
    print "--length 100                     Alternative length of the vector. Use to adjust the sub-particle center (default: length of the given vector; A)."
    print "--cmm vector.cmm                 A CMM file defining the location(s) of the subparticle(s) (use instead of --vector). Coordinates in Angstrom."
    print "--unique 1                       Keep only unique subparticles within angular distance (useful to remove overlapping sub-particles on symmetry axis)."
    print "--mindist 5                      Minimum distance between the subparticles in the image (all overlapping ones will be discarded; pixels)."
    print "--side 25                        Keep only particles within specified angular distance from side views (all others will be discarded; degrees)."
    print "--top 25                         Keep only particles within specified angular distance from top views (all others will be discarded; degrees)."
    print "--output subparticles            Output root for results."
    print "--j 8                            Number of threads."
    print "--np 4                           Number of MPI procs."
    print ""


def error(*msgs):
    usage()
    print "Error: " + '\n'.join(msgs)
    print " "
    sys.exit(2)


def within_mindist(p1, p2, mindist):
    """function that calculates whether two particles are closer to each other than the given distance in the projection"""

    x1 = p1.rlnCoordinateX
    y1 = p1.rlnCoordinateY
    x2 = p2.rlnCoordinateX
    y2 = p2.rlnCoordinateY
    distance_sqr = (x1 - x2)**2 + (y1 - y2)**2
    mindist_sqr = mindist**2

    if distance_sqr < mindist_sqr:
        return True
    else:
        return False


def within_unique(p1, p2, unique):
    """function that calculates whether two particles are closer to each other than the given angular distance"""

    v1 =  vector_from_two_eulers(p1.rlnAnglePsi, p1.rlnAngleTilt)
    v2 =  vector_from_two_eulers(p2.rlnAnglePsi, p2.rlnAngleTilt)

    dp =  dot_product(v1, v2)/(v1.length()*v2.length())

    if dp < -1:
        dp = -1.000
    if dp > 1:
        dp = 1.000

    angle = math.acos(dp)
    if angle <= math.radians(unique):
        return True
    else:
        return False


def filter_subparticles_mindist(subparticles, mindist):
    """function that filters all subparticles removing any overlaps"""

    new_subparticles = []

    n = 0
    while n < len(subparticles):
        subparticle=subparticles[n]
        overlaps = False

        m = 0
        while m < len(subparticles):
           subparticle2 = subparticles[m]
           if not subparticle.rlnImageName[:6] == subparticle2.rlnImageName[:6]:
               overlaps = within_mindist(subparticle, subparticle2, mindist)
           if overlaps:
               break

           m = m + 1              
        
        if not overlaps:
            new_subparticle = copy.copy(subparticle)
            new_subparticles.append(new_subparticle)

        n = n + 1

    return new_subparticles


def filter_subparticles_side(subparticles, side):
    new_subparticles = []

    n = 0
    while n < len(subparticles):
        subparticle=subparticles[n]

        if (abs(abs(subparticle.rlnAngleTilt) - radians(90)) < side):
            new_subparticle = copy.copy(subparticle)
            new_subparticles.append(new_subparticle)

        n = n + 1

    return new_subparticles


def filter_subparticles_top(subparticles, top):
    new_subparticles = []

    n = 0
    while n < len(subparticles):
        subparticle=subparticles[n]

        if (abs(abs(subparticle.rlnAngleTilt) - radians(180)) < top):
            new_subparticle = copy.copy(subparticle)
            new_subparticles.append(new_subparticle)

        n = n + 1

    return new_subparticles


def create_subparticles(particle, symmetry_matrices, subparticle_vector_list, part_image_size, relax_symmetry, randomize, output, unique, particles_total, align_subparticles):
    """function that obtains all subparticles from a given particle and sets the properties of each such subparticle"""

    subparticles = []

    # Euler angles that take particle to the orientation of the model
    rot  = -particle.rlnAnglePsi
    tilt = -particle.rlnAngleTilt
    psi  = -particle.rlnAngleRot
    matrix_particle = matrix_from_euler(rot, tilt, psi)

    particle_id = 1
    particles_total = particles_total + 1

    symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

    if randomize:
        # randomize the order of symmetry matrices, prevents preferred views
        random.shuffle(symmetry_matrix_ids)

    count_vectors = 0
    for subparticle_vector in subparticle_vector_list:
        rot, tilt, psi = euler_from_vector(subparticle_vector)
        matrix_from_subparticle_vector = matrix_from_euler(rot, tilt, psi)   

        for symmetry_matrix_id in symmetry_matrix_ids:
            overlaps = False

            # symmetry_matrix_id can be later written out to find out which symmetry matrix created this subparticle
            symmetry_matrix = symmetry_matrices[symmetry_matrix_id-1]

            subparticle = copy.deepcopy(particle)

            m = matrix_multiply((matrix_multiply(matrix_from_subparticle_vector, symmetry_matrix)), matrix_particle)

            if align_subparticles:
                rotNew, tiltNew, psiNew = euler_from_matrix(m)
            else:
                m2 = matrix_multiply(symmetry_matrix, matrix_particle)
                rotNew, tiltNew, psiNew = euler_from_matrix(m2)

            # save Euler angles that take the model to the orientation of the subparticle
            subparticle.setrlnAngleRot(-psiNew)
            subparticle.setrlnAngleTilt(-tiltNew)
            subparticle.setrlnAnglePsi(-rotNew)

            # subparticle origin
            x = -m.m[2][0] * subparticle_vector.distance + particle.rlnOriginX
            y = -m.m[2][1] * subparticle_vector.distance + particle.rlnOriginY
            z = -m.m[2][2] * subparticle_vector.distance             

            # modify the subparticle defocus paramaters by its z location
            subparticle.setrlnDefocusU(particle.rlnDefocusU + z)
            subparticle.setrlnDefocusV(particle.rlnDefocusV + z)
            
            # save the subparticle coordinates (integer part) relative to the user given image size and as a small shift in the origin (decimal part)

            subparticle_index = "{0:0>6}".format(str(particle_id))           
            particle_index = "{0:0>6}".format(str(int(particle.rlnImageName[0:6])))
            subparticle_image_name = splitext(particle.rlnImageName[7:])[0] + "_" + output + "_" + particle_index + ".mrcs"
            subparticle.setrlnImageName(subparticle_index + "@" + subparticle_image_name)
            subparticle.setrlnParticleName(str(particles_total))

            x_d, x_i = math.modf(x)
            y_d, y_i = math.modf(y)
            subparticle.setrlnCoordinateX(int(part_image_size/2) - x_i)
            subparticle.setrlnCoordinateY(int(part_image_size/2) - y_i)
            subparticle.setrlnOriginX(-x_d)
            subparticle.setrlnOriginY(-y_d)

            if unique >= 0:
                for subparticle2 in subparticles:
                    overlaps = within_unique(subparticle, subparticle2, unique)
                    if overlaps:
                        break

            if not overlaps:
                subparticles.append(subparticle)
                particle_id = particle_id + 1
                particles_total = particles_total + 1

            if relax_symmetry:
                # take just the first (random) one and finish
                break
        count_vectors = count_vectors + 1

    return subparticles


def create_star(subparticles, star_filename):
    '''function to create a Relion style STAR file for extracting (using relion_preprocess) all the subparticles for a given particle'''

    parameters = [
	"rlnMicrographName",
	"rlnCoordinateX",
	"rlnCoordinateY",
	"rlnImageName",
    ]  

    write_star(subparticles, parameters, star_filename)


def run_command(command, output=""):
    if not output:
        print "+++ " + command
        sys.stdout.flush()
        os.system(command)
    else:
        os.system(command + " > " + output)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['sym=',
                                                      'randomize',
                                                      'relax_symmetry',
                                                      'vector=',
                                                      'align_subparticles',
                                                      'cmm=',
                                                      'length=',
                                                      'output=',
                                                      'create_star',
                                                      'unique=',
                                                      'mindist=',
                                                      'top=',
                                                      'side=',
                                                      'particle_size=',
                                                      'subparticle_size=',
                                                      'masked_map=',
                                                      'split_stacks',
                                                      'extract_subparticles',
                                                      'angpix=',
                                                      'output=',
                                                      'j=',
                                                      'np=',
                                                      ])
    except getopt.GetoptError as e:
        error(e)

    if not (spawn.find_executable("bimg")):
        error("Error: Bsoft not found.",
              "Make sure Bsoft programs are in $PATH.")

    if not (spawn.find_executable("relion_refine")):
        error("Error: Relion not found.",
              "Make sure Relion programs are in $PATH.")

    symString = "C1"
    output = ""
    subparticle_distances = []
    subparticle_vector_list = []
    input_cmm = ""
    part_image_size = -1
    relax_symmetry = False
    align_subparticles = False
    randomize = False
    star = False
    unique = -1
    mindist = -1
    minanglediff = -1
    side = -1
    top = -1
    masked_map = ""
    subtract_masked_map = False
    split_stacks = False
    extract_subparticles = False
    angpix = -1
    output_dir = ""
    nr_threads = 1
    nr_procs = 1

    vectors_string = ''
    distances_string = ''

    for o, a in opts:
        if o == "--randomize":
            randomize = True
        if o == "--relax_symmetry":
            relax_symmetry = True
        if o == "--sym":
            symString = str(a)
        if o == "--vector":
            vectors_string = a
        if o == "--align_subparticles":
            align_subparticles = True
        if o == "--length":
            distances_string = a
        if o == "--top":
            top = radians(float(a))
        if o == "--side":
            side = radians(float(a))
        if o == "--mindist":
            mindist = float(a)
        if o == "--unique":
            unique = float(a)
        if o == "--cmm":
            input_cmm = a
        if o == "--create_star":
            star = True
        if o == "--particle_size":
            part_image_size = int(a)
        if o == "--subparticle_size":
            subpart_image_size = int(a)
        if o == "--masked_map":
            masked_map = str(a)
            subtract_masked_map = True
        if o == "--split_stacks":
            split_stacks = True
        if o == "--extract_subparticles":
            extract_subparticles = True
        if o == "--angpix":
            angpix = float(a)
        if o == "--output":
            output = str(a)
        if o == "--j":
            nr_threads = int(a)
        if o == "--np":
            nr_procs = int(a)

    if len(args) == 0:
        error("Error: No input file given.")

    input_star_filename = str(args[0])

    if not os.path.exists(input_star_filename):
        error("Error: Input file '%s' not found." % input_star_filename)

    if part_image_size < 0:
        error("Error: Particle image size not given.")

    if subpart_image_size < 0:
        error("Error: Sub-particle image size not given.")

    if angpix < 0:
        error("Error: Pixel size not given.")

    matrices = []
    particles = []
    parameters = []
    subparticles = []
    all_subparticles = []
    all_subparticles_subtracted = []

    if input_cmm:
        subparticle_vector_list = vectors_from_cmm(input_cmm, angpix)
    else:
        subparticle_vector_list = vectors_from_string(vectors_string)

    if distances_string:
        # Change distances from A to pixel units
        subparticle_distances = [float(x) / angpix for x in distances_string.split(',')]

        if len(subparticle_distances) != len(subparticle_vector_list):
            error("Error: The number of distances doesn't match the number of vectors!")

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

    if split_stacks:
        print "Creating and splitting the particle stack..."
        print " Creating a normal stack from which nothing is subtracted..."
        run_command("beditimg -create " + str(part_image_size) + "," + str(part_image_size) + "," + str(part_image_size) + " -fill 0 dummy_mask.mrc", "")
        run_command("relion_project --i dummy_mask.mrc --o " + path + "particles" + " --ang " + input_star_filename + " --subtract_exp --angpix " + str(angpix), "")
        run_command("rm -f dummy_mask.mrc", "")
        run_command("bsplit -digits 6 -first 1 " + path + "particles" + ".mrcs:mrc " + path + "particles" + ".mrc", "")
        print "Finished splitting the particle stack!"
        print " "

    if subtract_masked_map:
        print "Creating and splitting the particle stack..."
        print " Creating a stack from which the projections of the masked particle have been subtracted..."
        run_command("relion_project --i " + masked_map + " --o " + path + "particles_subtracted" + " --ang " + input_star_filename + " --subtract_exp --ctf --angpix " + str(angpix))
        print "Finished subtracting the model projections!"
        print " "

        print "Splitting the subtracted particle stack..."
        run_command("bsplit -digits 6 -first 1 " + path + "particles_subtracted" + ".mrcs:mrc " + path + "particles_subtracted" + ".mrc", "")
        print "Finished splitting the subtracted particle stack!"
        print " "

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

    relion_create_symmetry_ops_file(symString, "relion_symops.tmp")
    symmetry_matrices = matrix_from_symmetry_ops_file("relion_symops.tmp")
    run_command("rm relion_symops.tmp", "")

    nparticle = 0
    for particle in particles:
        subparticles = create_subparticles(particle, symmetry_matrices, subparticle_vector_list, part_image_size, relax_symmetry, randomize, "subparticles", unique, len(all_subparticles), align_subparticles)

        # to preserve numbering, ALL sub-particles are written to STAR files before filtering
        index = particle.rlnImageName[0:6]
        star_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + ".star"

        for subparticle in subparticles:
            particle_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + ".mrc"
            subparticle.setrlnMicrographName(particle_filename)
            subparticle_filename = splitext(particle.rlnImageName[7:])[0] + "_" + index + "_subparticles.mrcs"
            subparticle.setrlnImageName(subparticle.rlnImageName[0:7] + subparticle_filename)

        if subtract_masked_map:
            subparticles_subtracted = []
            star_filename_subtracted = splitext(particle.rlnImageName[7:])[0] + "_subtracted_" + index + ".star"
            for subparticle in subparticles:
                subparticle_new = copy.deepcopy(subparticle)
                subparticle_new.rlnImageName = subparticle_new.rlnImageName[:-24] + "subtracted" + subparticle_new.rlnImageName[-25:]
                subparticles_subtracted.append(subparticle_new)

        if star:
            create_star(subparticles, star_filename)
            if subtract_masked_map:
                create_star(subparticles_subtracted, star_filename_subtracted)

        if side > 0:
            subparticles = filter_subparticles_side(subparticles, side)
        if top > 0:
            subparticles = filter_subparticles_top(subparticles, top)
        if mindist > 0:
            subparticles = filter_subparticles_mindist(subparticles, mindist)

        if (side > 0 or top > 0 or mindist > 0) and subtract_masked_map:
            subparticles_subtracted = []
            for subparticle in subparticles:
                subparticle_new = copy.deepcopy(subparticle)
                subparticle_new.rlnImageName = subparticle_new.rlnImageName[:-24] + "subtracted" + subparticle_new.rlnImageName[-25:]
                subparticles_subtracted.append(subparticle_new)

        all_subparticles.extend(subparticles)
        if subtract_masked_map:
            all_subparticles_subtracted.extend(subparticles_subtracted)

        if nparticle == int(len(particles) * timer):
            sys.stdout.write("\b" * (c + 8))
            sys.stdout.write("." * c)
            sys.stdout.write("~~(,_,\">")
            sys.stdout.flush()
            timer = timer + 0.01
            c = c + 1

        nparticle = nparticle + 1
    sys.stdout.write("\n")

    print "Finished creating the subparticles!"
    print " "

    if extract_subparticles:
        print "Extracting the subparticles..."
        if nr_procs == 1:
            cmd = 'relion_preprocess '
        else:
            cmd = 'mpirun -np %s relion_preprocess_mpi ' % nr_procs

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

    sys.exit(0)

if __name__ == "__main__":    
    main()

