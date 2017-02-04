from helpers import *
import numpy as np

class ParticleInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):
	
		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		d = self.systemConstants["diameter"]
		
		self.posDat = []
		self.velDat = []

		for i in range(0,12):
			for j in range(0,12):
				self.posDat.append(Vec2(d+d*i,
					  					d+d*j))
				self.velDat.append(Vec2(d+d*i,
					  					d+d*j))
