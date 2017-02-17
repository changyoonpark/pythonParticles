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
					  					1.0+d+d*j))
				self.velDat.append(Vec2(0,
					  					0))

class BoundaryInitData:

	def __init__ (self,systemConstants):

		self.posDat = []
		padding = 1
		d = systemConstants["diameter"]

		numX = int((systemConstants["domain"][0] - padding *  2) / systemConstants["diameter"]) + 1 + 2
		numY = int((systemConstants["domain"][1] - padding *  2) / systemConstants["diameter"])

		self.wallParticleDiameter = (systemConstants["domain"][0]-2) / numX

		for i in range(0,numX):
			self.posDat.append(Vec2(padding - d + i * d, padding    ))
			self.posDat.append(Vec2(padding - d + i * d, padding - d))

		for j in range(0,numY):
			self.posDat.append(Vec2(padding - d, padding + d + j * d ))
			self.posDat.append(Vec2(padding    , padding + d + j * d ))
			self.posDat.append(Vec2(systemConstants["domain"][0] - padding    , padding + d + j * d ))
			self.posDat.append(Vec2(systemConstants["domain"][0] - padding + d, padding + d + j * d ))

