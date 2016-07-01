#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/09/23 (JTH)

import os
import sys
import getopt
import re
import copy

from star import *


def usage():
    print ""
    print "Usage: relion_create_movie_star.py [parameters] input.star [input.star]"
    print "Combines several movie.star files into one. Reads CTF info from the particle.star file"
    print ""
    print "Parameters: "
    print "--input particles.star         Input 'particles' STAR file with CTF info"
    print "--output particles_movie.star  Output 'movie' STAR file "
    print ""


def set_ctf_parameters_from_particles(movie_particles, particles):

	corresponding_particle = None
	for particle in particles:
		if movie_particles[0].rlnParticleName == particle.rlnImageName:
			corresponding_particle = particle
			break
			
	if corresponding_particle:
			
		for movie_particle in movie_particles:
			movie_particle.setrlnVoltage(corresponding_particle.rlnVoltage)
			movie_particle.setrlnDefocusU(corresponding_particle.rlnDefocusU)
			movie_particle.setrlnDefocusV(corresponding_particle.rlnDefocusV)
			movie_particle.setrlnDefocusAngle(corresponding_particle.rlnDefocusAngle)		
			movie_particle.setrlnSphericalAberration(corresponding_particle.rlnSphericalAberration)
			movie_particle.setrlnDetectorPixelSize(corresponding_particle.rlnDetectorPixelSize)
			movie_particle.setrlnMagnification(corresponding_particle.rlnMagnification)
			movie_particle.setrlnAmplitudeContrast(corresponding_particle.rlnAmplitudeContrast)
			movie_particle.setrlnCtfFigureOfMerit(corresponding_particle.rlnCtfFigureOfMerit)
	else:
		print "Error: Corresponding particle not found for movie particle " + movie_particles[0].rlnParticleName


def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], '', ['input=',
                                                      'output=',
                                                     ])
	except getopt.GetoptError as e:
		usage()
		print "Error: " + str(e)
		print " "
		sys.exit(2)

	for o, a in opts:
		if o == "--input":
			input_filename = a

		if o == "--output":
			output_filename = a

	try:
		output_star=open(output_filename, "w+")
		input_star=open(input_filename, "r")
	except:
		usage()
		sys.exit(2)

	if len(args) == 0:
		usage()
		sys.exit(2)
     
      
	movie_particles=[]
	particles, parameters = read_star(input_star)
     
	i = 0 
	while i < len(args):
		input_movie_star=open(args[i],"r")
		movie_particles_one_micrograph, parameters = read_star(input_movie_star)
		set_ctf_parameters_from_particles(movie_particles_one_micrograph, particles)
		
		movie_particles.extend(movie_particles_one_micrograph)
		
		i = i + 1
	 
	
	parameters = [
	"rlnMicrographName",
	"rlnCoordinateX",
	"rlnCoordinateY",
	"rlnImageName",
	"rlnDefocusU",
	"rlnDefocusV",
	"rlnDefocusAngle", 
	"rlnVoltage",
	"rlnSphericalAberration", 
	"rlnAmplitudeContrast",
	"rlnMagnification",
	"rlnDetectorPixelSize", 
	"rlnCtfFigureOfMerit",
	"rlnParticleName"]
		
	write_star(movie_particles, parameters, output_filename)

	input_star.close()
	output_star.close()


if __name__ == "__main__":    
    main()
