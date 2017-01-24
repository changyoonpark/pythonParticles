import cProfile
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from multiprocessing import Process
from timer import *

class Vec2:

	def __init__ (self,x,y):
		self.x = np.float16(x)
		self.y = np.float16(y)

	def length(self):
		return np.linalg.norm((self.x,self.y),2)

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

class Particle:
	def __init__ (self,
				  pID,
				  pos,vel,acc,
		          viscosity,
		          mass,
		          radius):
		self.q = Queue()		
		self.pID = pID
		self.pos = pos
		self.vel = Vec2(0,0)
		self.acc = Vec2(0,0)
		self.fext = Vec2(0,0)
		self.viscosity = np.float16(viscosity)
		self.mass = np.float16(mass)
		self.radius = np.float16(radius)
		self.hash = (-1,-1)
		self.neighborList = []




		
	def update(self,t):
		tic()
		self.hash()
		posBeforeUpdate = [particle.pos for particle in self.particleSet]
		velBeforeUpdate = [particle.vel for particle in self.particleSet]
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
		self.rk4(k1,k2,k3,k4,posBeforeUpdate,velBeforeUpdate)
		pdat = self.getPoints()
		print("time : {}".format(t))
		toc()
		if t % 2 == 0:
			self.scat.set_offsets(pdat)
		# if t == 100:
			

	def rk4(self,k1,k2,k3,k4,posBeforeUpdate,velBeforeUpdate):
		l = len(k1)
		for i in range(0,l):
			self.particleSet[i].pos = posBeforeUpdate[i] + \
						(1.0/6.0)*self.dt*(k1[0][i]+2.0*k2[0][i]+2.0*k3[0][i]+k4[0][i])
			self.particleSet[i].vel = velBeforeUpdate[i] + \
						(1.0/6.0)*self.dt*(k1[1][i]+2.0*k2[1][i]+2.0*k3[1][i]+k4[1][i])

	def run(self):
		fig = plt.figure(figsize=(7, 7))
		ax = fig.add_axes([0, 0, 1, 1], frameon=False)
		ax.set_xlim(0, self.dimLim.x)
		ax.set_ylim(0, self.dimLim.y)
		pointData = self.getPoints()
		xdat = np.array([np.float16(p[0]) for p in pointData])
		ydat = np.array([np.float16(p[1]) for p in pointData])
		self.scat = ax.scatter(xdat,ydat,
						s=619*self.interactionLength**2)
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

	def _particleInteractionsWorker(self,*gridPairs):
		for gridPair in gridPairs:
			print(gridPair)
			(gridX,gridY) = gridPair
			for i in [-1,0,1]:
				for j in [-1,0,1]:
					if gridX + i < 0 or gridY + j < 0:
						continue
					if gridX + i >= self.gridLim[0] or gridY + j >= self.gridLim[1]:
						continue
					self.collideCells((gridX,gridY),(gridX+i,gridY+j))
		self.q.put([particle.acc for particle in self.particleSet])

					


	def particleInteractions(self):
		print(self.gridList)
		self.resetForceBuffer()
		thread1 = Process(target = self._particleInteractionsWorker,args=(self.gridList[0]))
		thread2 = Process(target = self._particleInteractionsWorker,args=(self.gridList[1]))
		thread3 = Process(target = self._particleInteractionsWorker,args=(self.gridList[2]))
		thread4 = Process(target = self._particleInteractionsWorker,args=(self.gridList[3]))
		thread1.start(),thread2.start(),thread3.start(),thread4.start()
		forces = [self.q.get(),self.q.get(),self.q.get(),self.q.get()]
		thread1.join(),thread2.join(),thread3.join(),thread4.join()
		for i in range(len(self.ParticleSet)):
			self.particleSet[i].acc = forces[0][i] + \
			 						  forces[1][i] + \
			 						  forces[2][i] + \
			 						  forces[3][i]


	def doesCollide(self,particle1,particle2,dist,relpos,reldir,relvel):
		particle1.fext = particle1.fext - 300 * pow(dist - self.interactionLength,1) * reldir
		particle1.fext = particle1.fext - 0.5 * relvel

	def collideParticles(self,particle1,particle2):
		if particle1 is particle2 : 
			return False
		relpos = (particle1.pos - particle2.pos)
		relvel = (particle1.vel - particle2.vel)
		dist = relpos.length()		
		reldir = relpos.dir()
		# self.foo += 1
		if dist <= self.interactionLength:
			self.doesCollide(particle1,particle2,dist,relpos,reldir,relvel)
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
		for i in range(0,10):
			for j in range(0,10):
				newParticle = Particle(num,
									   Vec2(0.5+1.0*i*(1.01+0.02*np.random.rand()),
									   	    0.5+1.0*j*(1.01+0.02*np.random.rand())),
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
