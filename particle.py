import cProfile
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from matplotlib.animation import FuncAnimation
from multiprocessing import Process, Queue
from timer import *

class ParticlePair:
	def __init__ (self,pi,pj,relpos,relvel,dist,reldir):
		self.particlei = pi
		self.particlej = pj
		self.relpos = relpos
		self.relvel = relvel
		self.dist = dist
		self.reldir = reldir

	def __str__ (self):
		st = ""
		st += "Particle i : {}, Particle j : {}\n".format(self.particlei.pID,self.particlej.pID)
		st += "Rel Vel : {}\n".format(self.relvel)		
		st += "Rel Pos : {}\n".format(self.relpos)
		st += "Rel Dir : {}\n".format(self.reldir)		
		return st

class Vec2:

	def __init__ (self,x,y):
		self.x = x
		self.y = y

	def length(self):
		return sqrt(self.x*self.x+self.y*self.y)

	def dir(self):
		return self / self.length()
		# return Vec2(self.x / self.length(), self.y / self.length())

	def __mul__ (self,other):
		return Vec2(self.x * other, self.y * other)

	def __rmul__ (self,other):
		return Vec2(self.x * other, self.y * other)

	def __truediv__(self,other):
		return Vec2(self.x / other, self.y / other)

	def __add__ (self,other):
		return Vec2(self.x + other.x, self.y + other.y)

	def __sub__ (self,other):
		return Vec2(self.x - other.x, self.y - other.y)
		
	def __neg__ (self):
		return Vec2(-self.x,-self.y)

	def __str__ (self):
		return "Vector : ({},{})".format(self.x, self.y)

class Kernel:
	def __init__ (self,smoothingLength):
		self.h = smoothingLength

	def W(pairDat,h):
		q = pairDat.dist / h
		if q >= 0 and q <=1 :
			return 0.34104630662549*((2-q)**3-4*(1-q))/(h**2)
		elif q > 1 and q <= 2 :
			return 0.34104630662549*(2-q)**3/(h**2)
		else :
			return 0

	def gradW(pairDat,h):
		q = pairDat.dist / h
		if q >= 0 and q <= 1 :
			return 0.34104630662549*(4-3*(q-2)**2)/(pairDat.dist*h**2)*pairDat.reldir
		elif q >= 1 and q <= 2 :
			return 0.34104630662549*( -3*(2-q)**2)/(pairDat.dist*h**2)*pairDat.reldir
		else :
			return Vec2(0,0)

class Particle:
	def __init__ (self,
				  pID,
				  pos,vel,acc,
		          viscosity,
		          mass,
		          radius):
		self.pID = pID
		self.pos = pos
		self.vel = Vec2(0,0)
		self.acc = Vec2(0,0)
		self.fext = Vec2(0,0)
		# self.viscosity = np.float64(viscosity)
		# self.mass = np.float64(mass)
		# self.radius = np.float64(radius)
		self.viscosity = viscosity
		self.mass = mass
		self.radius = radius
		self.hash = (-1,-1)
		self.neighborList = []
		# self.foo = 0

	def initForce(self):
		self.fext = Vec2(0,0)

class ParticleSystem:
	def __init__ (self):
		self.q = Queue()		
		self.particleSet = []
		self.placeParticles(0.2)
		self.dimLim = Vec2(20,20)
		self.interactionLength = 0.5
		self.hashGridSize = self.interactionLength * 1.2
		self.gridLim = (int(self.dimLim.x // self.hashGridSize) + 1,
			            int(self.dimLim.y // self.hashGridSize) + 1)
		self.dt = 0.02
		self.walls = [(0.5,0.5),(19.5,19.5)]
		self.scat = None
		self.animation = None
		self.hashData = []
		self.
		for i in range(0,self.gridLim[0]):
			self.hashData.append([])
			for j in range(0,self.gridLim[1]):
				self.hashData[-1].append([])
		self.gridList = [[],[],[],[]]
		self.gridListAll = []
		self.pairsData = dict()
		gridNum = self.gridLim[0] * self.gridLim[1] // 4
		foo = 0
		for gridX in range(0,self.gridLim[0]):
			for gridY in range(0,self.gridLim[1]):
				self.gridList[foo].append((gridX,gridY))
				self.gridListAll.append((gridX,gridY))
				if len(self.gridList[foo]) % gridNum == 0 :
					foo += 1
					if foo > 3:
						foo = 3

	def getPoints(self):
		pointData = []
		for particle in self.particleSet:
			pointData.append(np.array([particle.pos.x,particle.pos.y]))
		return pointData

	def setPosAndVel(self,posVector,velVector):
		l = len(posVector)
		for i in range(0,l):
			self.particleSet[i].pos = posVector[i]
			self.particleSet[i].vel = velVector[i]

	def resetPairList(self):
		self.pairsData.clear()

	def calculateForces(self,withPos,withVel):
		self.setPosAndVel(withPos,withVel)
		self.resetForceBuffer()
		self.resetPairList()
		self.particleInteractions()
		self.boundaryInteractions()
		self.forceSum()
		return [particle.acc for particle in self.particleSet]

	def mult(self,vecList,k):
		return [vec * k for vec in vecList]

	def add(self,vecList1,vecList2):
		l = len(vecList1)
		return [vecList1[i]+vecList2[i] for i in range(0,l)]

	def update(self,t):
		tic()
		self.hash()
		posBeforeUpdate = [particle.pos for particle in self.particleSet]
		velBeforeUpdate = [particle.vel for particle in self.particleSet]
		self.rk4(posBeforeUpdate,velBeforeUpdate)
		# self.rk2(posBeforeUpdate,velBeforeUpdate)
		pdat = self.getPoints()
		print("time : {}".format(t))
		toc()
		if t % 3 == 0:
			self.scat.set_offsets(pdat)
		# if t == 100:

	def rk2(self,posBeforeUpdate,velBeforeUpdate):
		k1 = (velBeforeUpdate,
			  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.dt))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.dt))
		k2 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		l = len(k1)
		for i in range(0,l):
			self.particleSet[i].pos = posBeforeUpdate[i] + \
						self.dt*(k1[0][i])
			self.particleSet[i].vel = velBeforeUpdate[i] + \
						self.dt*(k1[1][i])

	def rk4(self,posBeforeUpdate,velBeforeUpdate):
		k1 = (velBeforeUpdate,
			  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.dt))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.dt))
		k2 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k2[0],0.5*self.dt))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k2[1],0.5*self.dt))
		k3 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k3[0],self.dt))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k3[1],self.dt))
		k4 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))
		l = len(k1)
		for i in range(0,l):
			self.particleSet[i].pos = posBeforeUpdate[i] + \
						(1.0/6.0)*self.dt*(k1[0][i]+2.0*k2[0][i]+2.0*k3[0][i]+k4[0][i])
			self.particleSet[i].vel = velBeforeUpdate[i] + \
						(1.0/6.0)*self.dt*(k1[1][i]+2.0*k2[1][i]+2.0*k3[1][i]+k4[1][i])

	def run(self):
		fig = plt.figure(figsize=(7, 7))
		ax = fig.add_axes([0, 0, 1, 1], frameon=True)
		ax.set_xlim(0, self.dimLim.x)
		ax.set_ylim(0, self.dimLim.y)
		pointData = self.getPoints()
		# xdat = np.array([np.float64(p[0]) for p in pointData])
		# ydat = np.array([np.float64(p[1]) for p in pointData])
		xdat = np.array([p[0] for p in pointData])
		ydat = np.array([p[1] for p in pointData])
		self.scat = ax.scatter(xdat,ydat,
						s=619*self.interactionLength**2,color = 'black',edgecolor= (1,1,1,0.5))
		animation = FuncAnimation(fig, self.update,interval=1,repeat=False)
		plt.show()

	def resetForceBuffer(self):
		for particle in self.particleSet:
			particle.fext = Vec2(0.,0.)

	def forceSum(self):
		for particle in self.particleSet:
			particle.acc = particle.fext / particle.mass
			particle.acc = particle.acc - Vec2(0,10.0)

	def boundaryInteractions(self):

		for particle in self.particleSet:

			if particle.pos.y < self.walls[0][1] : 
				particle.fext.y = particle.fext.y + 300*pow(self.walls[0][1] - particle.pos.y,2)
				particle.fext.y = particle.fext.y - 0.2 * particle.vel.y 
			elif particle.pos.y > self.walls[1][1] : 
				particle.fext.y = particle.fext.y - 300*pow(particle.pos.y - self.walls[1][1],2)
				particle.fext.y = particle.fext.y - 0.2 * particle.vel.y
			if particle.pos.x < self.walls[0][0] :
				particle.fext.x = particle.fext.x + 300*pow(self.walls[0][0] - particle.pos.x,2)
				particle.fext.x = particle.fext.x - 0.2 * particle.vel.x
			elif particle.pos.x > self.walls[1][0] : 
				particle.fext.x = particle.fext.x - 300*pow(particle.pos.x - self.walls[1][0],2)
				particle.fext.x = particle.fext.x - 0.2 * particle.vel.x

	def particleInteractionsWorker(self,gridPairs):
		for gridPair in gridPairs:
			(gridX,gridY) = gridPair
			for i in [-1,0,1]:
				for j in [-1,0,1]:
					if gridX + i < 0 or gridY + j < 0:
						continue
					if gridX + i >= self.gridLim[0] or gridY + j >= self.gridLim[1]:
						continue
					self.collideCells((gridX,gridY),(gridX+i,gridY+j))
		self.forAllParticlePairs(self.stupidSpring)

	def particleInteractions(self):
		self.particleInteractionsWorker(self.gridListAll)	

	def stupidSpring(self,pairData):
		dist = pairData.dist
		relvel = pairData.relvel
		reldir = pairData.reldir
		particlei = pairData.particlei
		particlei.fext = particlei.fext - 300 * pow(dist - self.interactionLength,1) * reldir
		particlei.fext = particlei.fext - 0.5 * relvel

	def forAllParticlePairs(self,doThis):
		# print("Iterating thru pairsdata")
		for pair in self.pairsData:
			# print(self.pairsData[pair])
			doThis(self.pairsData[pair])

	def collideParticles(self,particle1,particle2):
		if particle1 is particle2 : 
			return False
		relpos = (particle1.pos - particle2.pos)
		relvel = (particle1.vel - particle2.vel)
		dist = relpos.length()		
		reldir = relpos.dir()
		if dist <= self.interactionLength:
			if self.pairsData.get((particle1,particle2)) is None :
				self.pairsData[(particle1,particle2)]  \
				    = ParticlePair(particle1,particle2,relpos,relvel,dist,reldir)
				self.pairsData[(particle2,particle1)]  \
				    = ParticlePair(particle2,particle1,-relpos,-relvel,dist,-reldir)
			return True
		return False

	def collideCells(self,cell1,cell2):
		for particle1 in self.hashData[cell1[0]][cell1[1]]:
			for particle2 in self.hashData[cell2[0]][cell2[1]]:
				self.collideParticles(particle1, particle2)
				# if self.collideParticles(particle1,particle2):
					# particle1.neighborList.append(particle2)

	def wipeHash(self):
		for i in range(0,self.gridLim[0]):
			for j in range(0,self.gridLim[1]):
				if len(self.hashData[i][j]) is not 0:
					self.hashData[i][j] = []

	# def _hashWorker(self,particle):
	# 	particle.hash = (int(particle.pos.x // self.hashGridSize),
	# 				     int(particle.pos.y // self.hashGridSize))
	# 	self.hashData[particle.hash[0]][particle.hash[1]].append(particle)

	def hash(self):
		self.wipeHash()
		# self.pool.map(self._hashWorker,self.particleSet)
		for particle in self.particleSet:
			particle.hash = (int(particle.pos.x // self.hashGridSize),
						     int(particle.pos.y // self.hashGridSize))
			self.hashData[particle.hash[0]][particle.hash[1]].append(particle)


	def placeParticles(self,radius):
		num = 0
		for i in range(0,18):
			for j in range(0,18):
				newParticle = Particle(num,
									   Vec2(0.5+0.5*i*(1.1+0.02*np.random.rand()),
									   	    0.5+0.5*j*(1.1+0.02*np.random.rand())),
									   Vec2(0,0),Vec2(0,0),
									   0.1,0.1,0.1)
				self.particleSet.append(newParticle)
				num += 1
		print("Placed {} Particles.".format(num))

psys = ParticleSystem()
psys.run()
# def foo():
	# for i in range(0,20):
		# psys.update(i)
# cProfile.run('foo()')
