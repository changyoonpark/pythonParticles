from helpers import *
import numpy as np
from math import cos,sin

class ParticleDiskInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):

		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		d = self.systemConstants["diameter"]

		self.posDat = []
		self.velDat = []
		self.a_external = []

		nr = 2
		pi = 3.141592
		c = Vec2(4,4)
		gr = 1000

		# self.posDat.append(c)
		# self.velDat.append(Vec2(0,0))
		# self.a_external.append(Vec2(0,0))

		for i in range(1,nr):
			ntheta = int((2*pi*(i*d)) / (d) + 0.5)		
			for j in range(0, ntheta):
				dtheta = 2.0 * pi / ntheta;
				direction = -Vec2(cos(dtheta * j), sin(dtheta * j))
				self.posDat.append(c + Vec2( (i*d)*cos(dtheta * j) , (i*d)*sin(dtheta * j) ))
				self.velDat.append(Vec2(0,0))

				self.a_external.append( gr * (c - self.posDat[-1]).dir())
				# self.a_external.append(Vec2(0,-10))				# self.f_external.append(Vec2(0,0))

class ParticleInitData:
	def __init__ (self,
		          systemConstants,
		          particleVariables):
	
		self.systemConstants = systemConstants
		self.particleVariables = particleVariables

		padding = 1
		d = systemConstants["diameter"]
		maxx = systemConstants["domain"][0]
		maxy = systemConstants["domain"][1]
		
		(self.posDat,pmass) = generateBlock(Vec2(padding,padding+d),Vec2(maxx-padding,maxy-padding - 2 * d),d)
		self.velDat = []
		self.a_external = []

		for _ in self.posDat:
			self.velDat.append(Vec2(0,0))
			self.a_external.append(Vec2(0,0))

		particleVariables["mass"] = pmass

class BoundaryInitData:

	def __init__ (self,systemConstants):

		self.posDat = []
		padding = 1
		d = systemConstants["diameter"]
		maxx = systemConstants["domain"][0]
		maxy = systemConstants["domain"][1]


		self.posDat += generateWall(Vec2(padding - d,padding), Vec2(maxx-padding + d,padding),d,2)
		self.posDat += generateWall(Vec2(maxx-padding,padding + d),Vec2(maxx-padding ,maxy-padding - d),d,2)
		# self.posDat += generateWall(Vec2(maxx-padding-d,maxy-padding-2*d),Vec2(padding - d,maxy-padding-2*d),d,2)
		self.posDat += generateWall(Vec2(padding - d,maxy-padding-d -d ), Vec2(padding - d,padding),d,2)

		for v in self.posDat:
			print(v)



def generateWall(startCoord,endCoord,d,layers):

	pos = []

	num = int((endCoord-startCoord).length() / d + 0.5) + 1
	stencil = (endCoord-startCoord) / num
	sdir = stencil.dir()

	s = sqrt(3)/2
	offset = Vec2(0,0)
	Pi = 3.141592

	rotmat = [[cos(-Pi/2),-sin(-Pi/2)],[sin(-Pi/2),cos(-Pi/2)]]


	snormaldir = Vec2(rotmat[0][0] * sdir.dir().x + rotmat[0][1] * sdir.dir().y,
			    	  rotmat[1][0] * sdir.dir().x + rotmat[1][1] * sdir.dir().y)
	
	offset = 0.5 * sdir * stencil.length() + s * snormaldir * stencil.length()
	# doOffset = false
	for i in range(0,num):
		pos.append(startCoord + stencil * i)

	for i in range(0,num - 1):
		pos.append(startCoord +  offset + stencil * i)

	return pos

def generateBlock(leftBottom,rightTop,d):

	pos = []
	nx = int((rightTop - leftBottom).x / (d))
	ny = int((rightTop - leftBottom).y / (d))

	# for i in range(0,nx):
	# 	for j in range(0,ny,2):
	# 		pos.append(leftBottom + Vec2( d * i * sqrt(3) / 2, d * j ))

	maxx = 0
	maxy = 0

	for i in range(0,nx):
		for j in range(0,ny,2):
			pos.append(leftBottom +
					   Vec2(d*i,d*sqrt(3)/2*j))

			if maxx < pos[-1].x :
				maxx = pos[-1].x 

			if maxy < pos[-1].y :
				maxy = pos[-1].y 

	for i in range(0,nx-1):
		for j in range(0,ny,2):
			pos.append(leftBottom +
					   Vec2(d*i,d*sqrt(3)/2*j) + 
					   Vec2(0.5,sqrt(3)/2)*d )

			if maxx < pos[-1].x :
				maxx = pos[-1].x 

			if maxy < pos[-1].y :
				maxy = pos[-1].y 

	particleMass = 1.022 ** 2 * ( maxx - leftBottom.x ) * ( maxy - leftBottom.y ) / len(pos) * 1000
	

	return (pos,particleMass)
