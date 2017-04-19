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

		nr = 10
		pi = 3.141592
		c = Vec2(4,4)


		discrad = nr * d


		self.posDat.append(c)
		self.velDat.append(Vec2(0,0))
		self.a_external.append(Vec2(0,0))


		rotmat = [[cos(pi/2),-sin(pi/2)],[sin(pi/2),cos(pi/2)]]

		omega = 0
		gr = 10

		for i in range(1,nr):
			ntheta = int((2*pi*(i*d)) / (d) + 0.5)		

			for j in range(0, ntheta):
				dtheta = 2.0 * pi / ntheta;
				direction = Vec2(cos(dtheta * j), sin(dtheta * j))
				radius = (i*d + (((i)/25 + 0.06)**2) );
				self.posDat.append( c + Vec2( radius*cos(dtheta * j), radius*sin(dtheta * j) ) )
				self.velDat.append(omega * radius * Vec2(rotmat[0][0] * direction.dir().x + rotmat[0][1] * direction.dir().y,
			    	  				            rotmat[1][0] * direction.dir().x + rotmat[1][1] * direction.dir().y))	

				pos = self.posDat[-1]
				self.a_external.append( gr * (c - pos).dir())
				# self.a_external.append(Vec2(0,0))

		particleVariables["mass"] = 0.91 * (pi * (discrad ** 2)) / len(self.posDat) * 1000.0

		# leftBottom = Vec2(0,0)
		# rightTop   = Vec2(systemConstants["domain"][1],systemConstants["domain"][1])

		# center = Vec2(systemConstants["domain"][1] / 2.0 ,
		# 			  systemConstants["domain"][1] / 2.0 )

		# rad = (systemConstants["domain"][0]) / 3

		# (tempPos,pmass) = generateBlockStaggered(leftBottom,rightTop,d)
		# particleVariables["mass"] = pmass

		# gr = 1000

		# for pos in tempPos:
		# 	if (pos - center).length() < rad:
		# 		self.posDat.append(pos)
		# 		self.velDat.append(Vec2(0,0))
		# 		self.a_external.append( gr * (center - pos).dir())

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
		
		# (self.posDat,pmass) = generateBlockSquare(Vec2(padding,padding),Vec2(maxx-padding,maxy-padding),d)

		# pillar
		# (self.posDat,pmass) = generateBlockSquare(Vec2(padding,padding+d + 0.027 * d),Vec2(maxx-padding,maxy-padding - 2 * d + 0.027 * d),d)

		# dam
		(self.posDat,pmass) = generateBlockStaggered(Vec2(padding+0.*d,padding+d),Vec2(d+maxx/2,maxy/1.5),d)
		self.velDat = []
		self.a_external = []
		print("Particle Mass : {}".format(pmass))

		for _ in self.posDat:
			self.velDat.append(Vec2(0,0))
			self.a_external.append(Vec2(0,0))

		self.particleVariables["mass"] = pmass

class BoundaryInitData:

	def __init__ (self,systemConstants):

		self.posDat = []
		padding = 1
		d = systemConstants["diameter"]
		maxx = systemConstants["domain"][0]
		maxy = systemConstants["domain"][1]

		slack = 0.08
		self.posDat += generateWall(Vec2(padding - d - slack * d ,padding), Vec2(maxx-padding + d + slack * d,padding),d,2)
		self.posDat += generateWall(Vec2(maxx-padding + slack * d,padding + d),Vec2(maxx-padding + slack * d,maxy-padding - d),d,2)
		self.posDat += generateWall(Vec2(maxx-padding + slack * d,maxy-padding-d), Vec2(padding - 2 * d - slack * d,maxy-padding-d),d,2)
		self.posDat += generateWall(Vec2(padding - d - slack * d,maxy-padding-d -d ), Vec2(padding - d - slack * d,padding),d,2)

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

	for i in range(0,num):
		pos.append(startCoord + stencil * i)

	for i in range(0,num - 1):
		pos.append(startCoord +  offset + stencil * i)

	return pos

def generateBlockStaggered(leftBottom,rightTop,d):

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

	particleMass = 1.025 ** 2 * ( maxx - leftBottom.x ) * ( maxy - leftBottom.y ) / len(pos) * 1000
	

	return (pos,particleMass)



def generateBlockSquare(leftBottom,rightTop,d):

	pos = []
	nx = int((rightTop - leftBottom).x / (d))
	ny = int((rightTop - leftBottom).y / (d))

	# for i in range(0,nx):
	# 	for j in range(0,ny,2):
	# 		pos.append(leftBottom + Vec2( d * i * sqrt(3) / 2, d * j ))

	maxx = 0
	maxy = 0

	for i in range(0,nx):
		for j in range(0,ny):

			pos.append(leftBottom +
					   Vec2(d*i,d*j))

			if maxx < pos[-1].x :
				maxx = pos[-1].x 

			if maxy < pos[-1].y :
				maxy = pos[-1].y 

	particleMass = 0.99 * ( maxx - leftBottom.x + d) * ( maxy - leftBottom.y + d ) / len(pos) * 1000
	# particleMass = 1.0


	return (pos,particleMass)
