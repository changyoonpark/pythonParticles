from helpers import *
import numpy as np

class ParticleInitData:
	def __init__ (self,viscosity,particleMass,interactionLength):
		self.viscosity = viscosity
		self.particleMass =  particleMass
		self.interactionLength = interactionLength
		self.posDat = []
		self.velDat = []

		for i in range(0,12):
			for j in range(0,12):
				self.posDat.append(Vec2(0.5+self.interactionLength*i*(1.1+0.02*np.random.rand()),
					  					0.5+self.interactionLength*j*(1.1+0.02*np.random.rand())))
				self.velDat.append(Vec2(0.5+self.interactionLength*i*(1.1+0.02*np.random.rand()),
					  					0.5+self.interactionLength*j*(1.1+0.02*np.random.rand())))
