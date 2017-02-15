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

		for i in range(0,40):
			for j in range(0,40):
				self.posDat.append(Vec2(10.0+d+d*i,
					  					0.25+d+d*j))
				self.velDat.append(Vec2(0,
					  					0))
