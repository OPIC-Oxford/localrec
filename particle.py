#!/usr/bin/env python

# Serban Ilca & Juha T. Huiskonen
# Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford

# Created: 2014/06/02 (SI)
# Modified: 2014/09/12 (JTH)
# Modified: 2015/08/13 (JTH)

import math

class Particle:
    """Class defining Particle and methods needed to set its necessary parameters"""

    def __init__(self):
         self.rlnVoltage = 0
         self.rlnDefocusU = 0
         self.rlnDefocusV = 0
         self.rlnDefocusAngle = 0
         self.rlnSphericalAberration = 0
         self.rlnDetectorPixelSize = 0
         self.rlnCtfFigureOfMerit = 0
         self.rlnMagnification = 0
         self.rlnAmplitudeContrast = 0
         self.rlnImageName = ''
         self.rlnCtfImage = ''
         self.rlnCoordinateX = 0
         self.rlnCoordinateY = 0
         self.rlnCoordinateZ = 0
         self.rlnNormCorrection = 0
         self.rlnMicrographName = ''
         self.rlnGroupName = ''
         self.rlnGroupNumber = 0
         self.rlnOriginX = 0
         self.rlnOriginY = 0
         self.rlnOriginZ = 0
         self.rlnAngleRot = 0
         self.rlnAngleTilt = 0
         self.rlnAnglePsi = 0
         self.rlnClassNumber = 0
         self.rlnLogLikeliContribution = 0
         self.rlnRandomSubset = 0
         self.rlnParticleName = 0
         self.rlnOriginalParticleName = 0
         self.rlnNrOfSignificantSamples = 0 
         self.rlnNrOfFrames = 0
         self.rlnMaxValueProbDistribution = 0

    def setrlnVoltage(self, x):
        self.rlnVoltage = x

    def setrlnDefocusU(self, x):
        self.rlnDefocusU = x

    def setrlnDefocusV(self, x):
        self.rlnDefocusV = x

    def setrlnDefocusAngle(self, x):
        self.rlnDefocusAngle = x

    def setrlnSphericalAberration(self, x):
        self.rlnSphericalAberration = x

    def setrlnDetectorPixelSize(self, x):
        self.rlnDetectorPixelSize = x

    def setrlnCtfFigureOfMerit(self, x):
        self.rlnCtfFigureOfMerit = x

    def setrlnMagnification(self, x):
        self.rlnMagnification = x

    def setrlnAmplitudeContrast(self, x):
        self.rlnAmplitudeContrast = x

    def setrlnImageName(self, x):
        self.rlnImageName = x

    def setrlnCtfImage(self, x):
        self.rlnCtfImage = x

    def setrlnCoordinateX(self, x):
        self.rlnCoordinateX = x

    def setrlnCoordinateY(self, x):
        self.rlnCoordinateY = x

    def setrlnCoordinateZ(self, x):
        self.rlnCoordinateZ = x

    def setrlnNormCorrection(self, x):
        self.rlnNormCorrection = x

    def setrlnMicrographName(self, x):
        self.rlnMicrographName = x

    def setrlnGroupName(self, x):
        self.rlnGroupName = x

    def setrlnGroupNumber(self, x):
        self.rlnGroupNumber = x

    def setrlnOriginX(self, x):
        self.rlnOriginX = x
  
    def setrlnOriginY(self, x):
        self.rlnOriginY = x

    def setrlnOriginZ(self, x):
        self.rlnOriginZ = x

    def setrlnAngleRot(self, x):
        self.rlnAngleRot = x

    def setrlnAngleTilt(self, x):
        self.rlnAngleTilt = x

    def setrlnAnglePsi(self, x):
        self.rlnAnglePsi = x

    def setrlnClassNumber(self, x):
        self.rlnClassNumber = x

    def setrlnLogLikeliContribution(self, x):
        self.rlnLogLikeliContribution = x

    def setrlnRandomSubset(self, x):
        self.rlnRandomSubset = x

    def setrlnParticleName(self, x):
        self.rlnParticleName = x
 
    def setrlnOriginalParticleName(self, x):
        self.rlnOriginalParticleName = x

    def setrlnNrOfSignificantSamples(self, x):
        self.rlnNrOfSignificantSamples = x

    def setrlnNrOfFrames(self, x):
        self.rlnNrOfFrames = x

    def setrlnMaxValueProbDistribution(self, x):
        self.rlnMaxValueProbDistribution = x

    def print_particle(self):
        print("%.6f" %float(self.rlnVoltage) + "\t")
  

