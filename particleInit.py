from helpers import *
import numpy as np
from math import sin,cos


class ParticleDiskInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):

		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		d = self.systemConstants["diameter"]

		self.posDat = []
		self.velDat = []
		self.pressureDat = []

		nr = 15
		pi = 3.141592
		c = Vec2(10,10);

		self.posDat.append(c);
		self.velDat.append(Vec2(0,0));


		for i in range(1,nr):
			ntheta = int((2*pi*(i*d)) / (d) + 0.5)		
			for j in range(0, ntheta):
				dtheta = 2.0 * pi / ntheta;
				self.posDat.append(c + Vec2( (i*d)*cos(dtheta * j) , (i*d)*sin(dtheta * j) ))
				self.velDat.append(Vec2(0,0))

	 
class ParticleInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):
	
		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		d = self.systemConstants["diameter"]
		
		self.posDat = []
		self.velDat = []

		for i in range(0,15):
			for j in range(0,15):

				self.posDat.append(Vec2(2.05+d+d*i,
					  					1.05+d+d*j))

				self.velDat.append(Vec2( 0 ,
					  					 0 ))

class BoundaryInitData:

	def __init__ (self,systemConstants):

		self.posDat = []
		padding = 1
		d = systemConstants["diameter"]

		numX = int((systemConstants["domain"][0] - padding *  2) / systemConstants["diameter"]) + 1 + 2
		numY = int((systemConstants["domain"][1] - padding *  2) / systemConstants["diameter"])

		# numX = int((4 - padding *  2) / systemConstants["diameter"]) + 1 + 2
		# numY = int((4 - padding *  2) / systemConstants["diameter"])

		# self.wallParticleDiameter = (systemConstants["domain"][0]-2) / numX

		# for i in range(0,numX):
		# 	self.posDat.append(Vec2(padding - d + i * d, padding    ))

		# for i in range(0,numX):
		# 	self.posDat.append(Vec2(padding - d + i * d, padding - d))

		# for j in range(0,numY-10):
		# 	self.posDat.append(Vec2(padding - d, padding + d + j * d ))
		# 	self.posDat.append(Vec2(padding    , padding + d + j * d ))
		# 	self.posDat.append(Vec2(systemConstants["domain"][0] - padding    , padding + d + j * d ))
		# 	self.posDat.append(Vec2(systemConstants["domain"][0] - padding + d, padding + d + j * d ))
		# 	# self.posDat.append(Vec2(4.5 - padding    , padding + d + j * d ))
		# 	# self.posDat.append(Vec2(4.5 - padding + d, padding + d + j * d ))

