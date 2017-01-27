from helpers import *
import numpy as np

class ParticleInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):
	
		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		il = self.systemConstants["interactionlen"]
		
		self.posDat = []
		self.velDat = []

		for i in range(0,12):
			for j in range(0,12):
				self.posDat.append(Vec2(0.5+il*i*(1.1+0.02*np.random.rand()),
					  					0.5+il*j*(1.1+0.02*np.random.rand())))
				self.velDat.append(Vec2(0.5+il*i*(1.1+0.02*np.random.rand()),
					  					0.5+il*j*(1.1+0.02*np.random.rand())))
