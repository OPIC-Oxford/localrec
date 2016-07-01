#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/06/02 (SI)
# Modified: 2014/10/20 (JTH)
# Modified: 2015/08/13 (JTH)

import particle, math

def read_star(input_star):
    """Function that reads all particles and their parameters from a star file"""
    star_lines = input_star.readlines()
    AllParticles = []
    parameters = {}
    index = 1

    for star_line in star_lines:
        star_values = star_line.split()

        if 1 <= len(star_values) <= 2:

            if star_values[0] == "_rlnVoltage":
                try: 
                    parameters['rlnVoltage'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnVoltage'] = index
                    index = index + 1   
            
            if star_values[0] == "_rlnDefocusU":
                try: 
                    parameters['rlnDefocusU'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnDefocusU'] = index
                    index = index + 1   
            
            if star_values[0] == "_rlnDefocusV":
                try: 
                    parameters['rlnDefocusV'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnDefocusV'] = index
                    index = index + 1   
            
            if star_values[0] == "_rlnDefocusAngle":
                try: 
                    parameters['rlnDefocusAngle'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnDefocusAngle'] = index
                    index = index + 1   

            if star_values[0] == "_rlnSphericalAberration":
                try: 
                    parameters['rlnSphericalAberration'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnSphericalAberration'] = index
                    index = index + 1   
            if star_values[0] == "_rlnDetectorPixelSize":
                try: 
                    parameters['rlnDetectorPixelSize'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnDetectorPixelSize'] = index
                    index = index + 1   
            if star_values[0] == "_rlnCtfFigureOfMerit":
                try: 
                    parameters['rlnCtfFigureOfMerit'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnCtfFigureOfMerit'] = index
                    index = index + 1   
            if star_values[0] == "_rlnMagnification":
                try: 
                    parameters['rlnMagnification'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnMagnification'] = index
                    index = index + 1   
            if star_values[0] == "_rlnAmplitudeContrast":
                try: 
                    parameters['rlnAmplitudeContrast'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnAmplitudeContrast'] = index
                    index = index + 1   

            if star_values[0] == "_rlnImageName":
                try: 
                    parameters['rlnImageName'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnImageName'] = index
                    index = index + 1   

            if star_values[0] == "_rlnCtfImage":
                try: 
                    parameters['rlnCtfImage'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnCtfImage'] = index
                    index = index + 1   

            if star_values[0] == "_rlnCoordinateX":
                try: 
                    parameters['rlnCoordinateX'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnCoordinateX'] = index
                    index = index + 1   

            if star_values[0] == "_rlnCoordinateY":
                try: 
                    parameters['rlnCoordinateY'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnCoordinateY'] = index
                    index = index + 1 

            if star_values[0] == "_rlnCoordinateZ":
                try: 
                    parameters['rlnCoordinateZ'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnCoordinateZ'] = index
                    index = index + 1 

            if star_values[0] == "_rlnNormCorrection":
                try: 
                    parameters['rlnNormCorrection'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnNormCorrection'] = index
                    index = index + 1   

            if star_values[0] == "_rlnMicrographName":
                try: 
                    parameters['rlnMicrographName'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnMicrographName'] = index
                    index = index + 1   

            if star_values[0] == "_rlnGroupName":
                try: 
                    parameters['rlnGroupName'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnGroupName'] = index
                    index = index + 1   

            if star_values[0] == "_rlnGroupNumber":
                try: 
                    parameters['rlnGroupNumber'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnGroupNumber'] = index
                    index = index + 1   

            if star_values[0] == "_rlnOriginX":
                try: 
                    parameters['rlnOriginX'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnOriginX'] = index
                    index = index + 1   

            if star_values[0] == "_rlnOriginY":
                try: 
                    parameters['rlnOriginY'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnOriginY'] = index
                    index = index + 1  

            if star_values[0] == "_rlnAngleRot":
                try: 
                    parameters['rlnAngleRot'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnAngleRot'] = index
                    index = index + 1   

            if star_values[0] == "_rlnAngleTilt":
                try: 
                    parameters['rlnAngleTilt'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnAngleTilt'] = index
                    index = index + 1   

            if star_values[0] == "_rlnAnglePsi":
                try: 
                    parameters['rlnAnglePsi'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnAnglePsi'] = index
                    index = index + 1   

            if star_values[0] == "_rlnClassNumber":
                try: 
                    parameters['rlnClassNumber'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnClassNumber'] = index
                    index = index + 1   

            if star_values[0] == "_rlnLogLikeliContribution":
                try: 
                    parameters['rlnLogLikeliContribution'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnLogLikeliContribution'] = index
                    index = index + 1   

            if star_values[0] == "_rlnRandomSubset":
                try: 
                    parameters['rlnRandomSubset'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnRandomSubset'] = index
                    index = index + 1   

            if star_values[0] == "_rlnParticleName":
                try: 
                    parameters['rlnParticleName'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnParticleName'] = index
                    index = index + 1   

            if star_values[0] == "_rlnOriginalParticleName":
                try: 
                    parameters['rlnOriginalParticleName'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnOriginalParticleName'] = index
                    index = index + 1   

            if star_values[0] == "_rlnNrOfSignificantSamples":
                try: 
                    parameters['rlnNrOfSignificantSamples'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnNrOfSignificantSamples'] = index
                    index = index + 1   

            if star_values[0] == "_rlnNrOfFrames":
                try: 
                    parameters['rlnNrOfFrames'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnNrOfFrames'] = index
                    index = index + 1   

            if star_values[0] == "_rlnMaxValueProbDistribution":
                try: 
                    parameters['rlnMaxValueProbDistribution'] = int(star_values[1][1:])
                    index = index + 1
                except:
                    parameters['rlnMaxValueProbDistribution'] = index
                    index = index + 1

        if len(star_values) > 2: 
            newParticle = particle.Particle()
            if 'rlnVoltage' in parameters: newParticle.setrlnVoltage(float(star_values[parameters['rlnVoltage']-1]))
            if 'rlnDefocusU' in parameters: newParticle.setrlnDefocusU(float(star_values[parameters['rlnDefocusU']-1]))
            if 'rlnDefocusV' in parameters: newParticle.setrlnDefocusV(float(star_values[parameters['rlnDefocusV']-1]))
            if 'rlnDefocusAngle' in parameters: newParticle.setrlnDefocusAngle(float(star_values[parameters['rlnDefocusAngle']-1]))
            if 'rlnSphericalAberration' in parameters: newParticle.setrlnSphericalAberration(float(star_values[parameters['rlnSphericalAberration']-1]))
            if 'rlnDetectorPixelSize' in parameters: newParticle.setrlnDetectorPixelSize(float(star_values[parameters['rlnDetectorPixelSize']-1]))
            if 'rlnCtfFigureOfMerit' in parameters: newParticle.setrlnCtfFigureOfMerit(float(star_values[parameters['rlnCtfFigureOfMerit']-1]))
            if 'rlnMagnification' in parameters: newParticle.setrlnMagnification(float(star_values[parameters['rlnMagnification']-1]))
            if 'rlnAmplitudeContrast' in parameters: newParticle.setrlnAmplitudeContrast(float(star_values[parameters['rlnAmplitudeContrast']-1]))
            if 'rlnImageName' in parameters: newParticle.setrlnImageName(star_values[parameters['rlnImageName']-1])
            if 'rlnCtfImage' in parameters: newParticle.setrlnCtfImage(star_values[parameters['rlnCtfImage']-1])
            if 'rlnCoordinateX' in parameters: newParticle.setrlnCoordinateX(float(star_values[parameters['rlnCoordinateX']-1]))
            if 'rlnCoordinateY' in parameters: newParticle.setrlnCoordinateY(float(star_values[parameters['rlnCoordinateY']-1]))
            if 'rlnCoordinateZ' in parameters: newParticle.setrlnCoordinateZ(float(star_values[parameters['rlnCoordinateZ']-1]))
            if 'rlnNormCorrection' in parameters: newParticle.setrlnNormCorrection(float(star_values[parameters['rlnNormCorrection']-1]))
            if 'rlnMicrographName' in parameters: newParticle.setrlnMicrographName(star_values[parameters['rlnMicrographName']-1])
            if 'rlnGroupName' in parameters: newParticle.setrlnGroupName(star_values[parameters['rlnGroupName']-1])
            if 'rlnGroupNumber' in parameters: newParticle.setrlnGroupNumber(star_values[parameters['rlnGroupNumber']-1])
            if 'rlnOriginX' in parameters: newParticle.setrlnOriginX(float(star_values[parameters['rlnOriginX']-1]))
            if 'rlnOriginY' in parameters: newParticle.setrlnOriginY(float(star_values[parameters['rlnOriginY']-1]))
            if 'rlnAngleRot' in parameters: newParticle.setrlnAngleRot(math.radians(float(star_values[parameters['rlnAngleRot']-1])))
            if 'rlnAngleTilt' in parameters: newParticle.setrlnAngleTilt(math.radians(float(star_values[parameters['rlnAngleTilt']-1])))
            if 'rlnAnglePsi' in parameters: newParticle.setrlnAnglePsi(math.radians(float(star_values[parameters['rlnAnglePsi']-1])))
            if 'rlnClassNumber' in parameters: newParticle.setrlnClassNumber(star_values[parameters['rlnClassNumber']-1])
            if 'rlnLogLikeliContribution' in parameters: newParticle.setrlnLogLikeliContribution(float(star_values[parameters['rlnLogLikeliContribution']-1]))
            if 'rlnRandomSubset' in parameters: newParticle.setrlnRandomSubset(star_values[parameters['rlnRandomSubset']-1])
            if 'rlnParticleName' in parameters: newParticle.setrlnParticleName(star_values[parameters['rlnParticleName']-1])
            if 'rlnOriginalParticleName' in parameters: newParticle.setrlnOriginalParticleName(star_values[parameters['rlnOriginalParticleName']-1])
            if 'rlnNrOfSignificantSamples' in parameters: newParticle.setrlnNrOfSignificantSamples(float(star_values[parameters['rlnNrOfSignificantSamples']-1]))
            if 'rlnNrOfFrames' in parameters: newParticle.setrlnNrOfFrames(star_values[parameters['rlnNrOfFrames']-1])
            if 'rlnMaxValueProbDistribution' in parameters: newParticle.setrlnMaxValueProbDistribution(float(star_values[parameters['rlnMaxValueProbDistribution']-1]))
            AllParticles.append(newParticle)

    return AllParticles, parameters


def write_star(particles, parameters, filename):
    """Function that writes all particles and their parameters to a star file"""

    counter = 1

    new_file = open(filename, 'w')
    new_file.write("\ndata_\n\nloop_\n")

    if 'rlnVoltage' in parameters: 
        new_file.write("_rlnVoltage #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnDefocusU' in parameters: 
        new_file.write("_rlnDefocusU #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnDefocusV' in parameters: 
        new_file.write("_rlnDefocusV #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnDefocusAngle' in parameters: 
        new_file.write("_rlnDefocusAngle #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnSphericalAberration' in parameters: 
        new_file.write("_rlnSphericalAberration #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnDetectorPixelSize' in parameters: 
        new_file.write("_rlnDetectorPixelSize #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnCtfFigureOfMerit' in parameters: 
        new_file.write("_rlnCtfFigureOfMerit #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnMagnification' in parameters: 
        new_file.write("_rlnMagnification #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnAmplitudeContrast' in parameters: 
        new_file.write("_rlnAmplitudeContrast #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnImageName' in parameters: 
        new_file.write("_rlnImageName #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnCtfImage' in parameters: 
        new_file.write("_rlnCtfImage #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnCoordinateX' in parameters: 
        new_file.write("_rlnCoordinateX #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnCoordinateY' in parameters: 
        new_file.write("_rlnCoordinateY #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnCoordinateZ' in parameters: 
        new_file.write("_rlnCoordinateZ #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnNormCorrection' in parameters: 
        new_file.write("_rlnNormCorrection #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnMicrographName' in parameters: 
        new_file.write("_rlnMicrographName #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnGroupName' in parameters: 
        new_file.write("_rlnGroupName #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnGroupNumber' in parameters: 
        new_file.write("_rlnGroupNumber #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnOriginX' in parameters: 
        new_file.write("_rlnOriginX #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnOriginY' in parameters: 
        new_file.write("_rlnOriginY #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnAngleRot' in parameters: 
        new_file.write("_rlnAngleRot #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnAngleTilt' in parameters: 
        new_file.write("_rlnAngleTilt #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnAnglePsi' in parameters: 
        new_file.write("_rlnAnglePsi #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnClassNumber' in parameters: 
        new_file.write("_rlnClassNumber #" + str(counter) + "\n")
        counter = counter + 1
    
    if 'rlnLogLikeliContribution' in parameters: 
        new_file.write("_rlnLogLikeliContribution #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnRandomSubset' in parameters: 
        new_file.write("_rlnRandomSubset #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnParticleName' in parameters: 
        new_file.write("_rlnParticleName #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnOriginalParticleName' in parameters: 
        new_file.write("_rlnOriginalParticleName #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnNrOfSignificantSamples' in parameters: 
        new_file.write("_rlnNrOfSignificantSamples #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnNrOfFrames' in parameters: 
        new_file.write("_rlnNrOfFrames #" + str(counter) + "\n")
        counter = counter + 1

    if 'rlnMaxValueProbDistribution' in parameters: 
        new_file.write("_rlnMaxValueProbDistribution #" + str(counter) + "\n")
        counter = counter + 1

    for i in range(len(particles)):
        if 'rlnVoltage' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnVoltage) + " \t")
        if 'rlnDefocusU' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnDefocusU) + " \t")
        if 'rlnDefocusV' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnDefocusV) + " \t")
        if 'rlnDefocusAngle' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnDefocusAngle) + " \t")
        if 'rlnSphericalAberration' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnSphericalAberration) + " \t")
        if 'rlnDetectorPixelSize' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnDetectorPixelSize) + " \t")
        if 'rlnCtfFigureOfMerit' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnCtfFigureOfMerit) + " \t")
        if 'rlnMagnification' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnMagnification) + " \t")
        if 'rlnAmplitudeContrast' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnAmplitudeContrast) + " \t")
        if 'rlnImageName' in parameters:
            new_file.write(str(particles[i].rlnImageName) + " \t")
        if 'rlnCtfImage' in parameters:
            new_file.write(str(particles[i].rlnCtfImage) + " \t")
        if 'rlnCoordinateX' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnCoordinateX) + " \t")
        if 'rlnCoordinateY' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnCoordinateY) + " \t")
        if 'rlnCoordinateZ' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnCoordinateZ) + " \t")
        if 'rlnNormCorrection' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnNormCorrection) + " \t")
        if 'rlnMicrographName' in parameters:
            new_file.write(str(particles[i].rlnMicrographName) + " \t")
        if 'rlnGroupName' in parameters:
            new_file.write(str(particles[i].rlnGroupName) + " \t")
        if 'rlnGroupNumber' in parameters:
            new_file.write("%d" %int(particles[i].rlnGroupNumber) + " \t")
        if 'rlnOriginX' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnOriginX) + " \t")
        if 'rlnOriginY' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnOriginY) + " \t")
        if 'rlnAngleRot' in parameters:
            new_file.write("%.6f" %(math.degrees(float(particles[i].rlnAngleRot))) + " \t")
        if 'rlnAngleTilt' in parameters:
            new_file.write("%.6f" %(math.degrees(float(particles[i].rlnAngleTilt))) + " \t")
        if 'rlnAnglePsi' in parameters:
            new_file.write("%.6f" %(math.degrees(float(particles[i].rlnAnglePsi))) + " \t")
        if 'rlnClassNumber' in parameters:
            new_file.write("%d" %int(particles[i].rlnClassNumber) + " \t")
        if 'rlnLogLikeliContribution' in parameters:
            new_file.write("%.7e" %float(particles[i].rlnLogLikeliContribution) + " \t")
        if 'rlnRandomSubset' in parameters:
            new_file.write("%d" %int(particles[i].rlnRandomSubset) + " \t")
        if 'rlnParticleName' in parameters:
            new_file.write(str(particles[i].rlnParticleName) + " \t")
        if 'rlnOriginalParticleName' in parameters:
            new_file.write(str(particles[i].rlnOriginalParticleName) + " \t")
        if 'rlnNrOfSignificantSamples' in parameters:
            new_file.write("%d" %int(particles[i].rlnNrOfSignificantSamples) + " \t")
        if 'rlnNrOfFrames' in parameters:
            new_file.write("%d" %int(particles[i].rlnNrOfFrames) + " \t")
        if 'rlnMaxValueProbDistribution' in parameters:
            new_file.write("%.6f" %float(particles[i].rlnMaxValueProbDistribution) + " \t")
        new_file.write("\n")
    new_file.write("\n")
    new_file.close()
