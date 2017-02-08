import cProfile
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from multiprocessing import Process, Queue
from timer import *
from kernels import *
from helpers import *

class Particle:
	def __init__ (self,
				  pID,
				  pos,vel,acc,
		          particleVariables):
		self.pID = pID
		self.pos = pos
		self.vel = Vec2(0,0)
		self.acc = Vec2(0,0)
		self.fext = Vec2(0,0)
		self.particleVariables = particleVariables
		self.hash = (-1,-1)
		self.neighborList = []

	def initForce(self):
		self.fext = Vec2(0,0)

	def __str__(self):
		s =  "---------------------------------------\n"
		s += "Particle : {}\n".format(self.pID)
		s += "position : {}\n".format(self.pos)
		s += "velocity : {}\n".format(self.vel)
		s += "acceler  : {}\n".format(self.acc)
		s += "fext     : {}\n".format(self.fext)
		s += "	>>> Particle Variable Dictionary <<<\n"
		for key in self.particleVariables:
			s += "{} : {}\n".format(key,self.particleVariables[key])
		return s

class ParticleSystem:
	def __init__ (self,interactionAlgo,		       
		               particleInitData
		               ):

		self.systemConstants = dict()
		self.particleSet = []
		self.dimLim = Vec2(None,None,tup = particleInitData.systemConstants["domain"])
		self.scat = None
		self.animation = None
		self.hashData = []
		self.interactionOperations = interactionAlgo.getAlgoProcedure()

		if particleInitData is None : 
			self.systemConstants["interactionlen"] = 0.5			
			il = self.systemConstants["interactionlen"]
			num = 0
			for i in range(1,11):
				for j in range(1,11):
					newParticle = Particle(num,
										   Vec2(self.systemConstants["interactionlen"]+il*i,
										   	    self.systemConstants["interactionlen"]+il*j),
										   Vec2(0,0),Vec2(0,0),
										   0.1,0.1)
					self.particleSet.append(newParticle)
					num += 1
			print("Placed {} Particles.".format(num))

		else : 
			print("Reading initial particle data")
			num = 0
			self.systemConstants = particleInitData.systemConstants			
			posDat = particleInitData.posDat
			velDat = particleInitData.velDat
			for idx in range(0,len(posDat)) : 
				newParticle = Particle(num,posDat[idx],velDat[idx],Vec2(0,0),
					particleInitData.particleVariables
					)
				self.particleSet.append(newParticle)
				num += 1


		self.hashGridSize = self.systemConstants["interactionlen"] * 1.5
		self.gridLim = (int(self.dimLim.x // self.hashGridSize) + 1,
			            int(self.dimLim.y // self.hashGridSize) + 1)
		gridNum = self.gridLim[0] * self.gridLim[1] // 4
		for i in range(0,self.gridLim[0]):
			self.hashData.append([])
			for j in range(0,self.gridLim[1]):
				self.hashData[-1].append([])
		self.gridList = [[],[],[],[]]
		self.gridListAll = []
		self.pairsData = dict()

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
		self.constructContacts()

		for operation in self.operationFuncs:
			operation(self.systemConstants,self.pairsData,particle)

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

		self.resetPairList()
		self.constructContacts()	

		self.solveTimeStep()
		# posBeforeUpdate = [particle.pos for particle in self.particleSet]
		# velBeforeUpdate = [particle.vel for particle in self.particleSet]
		# self.rk4(posBeforeUpdate,velBeforeUpdate)
		# self.rk2(posBeforeUpdate,velBeforeUpdate)
		pdat = self.getPoints()
		print("time : {}".format(t))
		toc()
		if t % 3 == 0:
			self.scat.set_offsets(pdat)
		# if t == 100:

	def solveTimeStep(self):
		for operation in self.interactionOperations:
			operation(self.systemConstants,self.pairsData,self.particleSet)


	def rk2(self,posBeforeUpdate,velBeforeUpdate):
		k1 = (velBeforeUpdate,
			  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.systemConstants["dt"]))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.systemConstants["dt"]))
		k2 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		l = len(k1)
		for i in range(0,l):
			self.particleSet[i].pos = posBeforeUpdate[i] + \
						self.systemConstants["dt"]*(k1[0][i])
			self.particleSet[i].vel = velBeforeUpdate[i] + \
						self.systemConstants["dt"]*(k1[1][i])

	def rk4(self,posBeforeUpdate,velBeforeUpdate):
		k1 = (velBeforeUpdate,
			  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.systemConstants["dt"]))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.systemConstants["dt"]))
		k2 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k2[0],0.5*self.systemConstants["dt"]))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k2[1],0.5*self.systemConstants["dt"]))
		k3 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))

		intermediatePos = self.add(posBeforeUpdate,self.mult(k3[0],self.systemConstants["dt"]))
		intermediateVel = self.add(velBeforeUpdate,self.mult(k3[1],self.systemConstants["dt"]))
		k4 = (intermediateVel,
			  self.calculateForces(intermediatePos,intermediateVel))
		l = len(k1)
		for i in range(0,l):
			self.particleSet[i].pos = posBeforeUpdate[i] + \
						(1.0/6.0)*self.systemConstants["dt"]*(k1[0][i]+2.0*k2[0][i]+2.0*k3[0][i]+k4[0][i])
			self.particleSet[i].vel = velBeforeUpdate[i] + \
						(1.0/6.0)*self.systemConstants["dt"]*(k1[1][i]+2.0*k2[1][i]+2.0*k3[1][i]+k4[1][i])

	def run(self):
		fig = plt.figure(figsize=(7, 7))
		ax = fig.add_axes([0, 0, 1, 1], frameon=True)
		ax.set_xlim(0, self.dimLim.x)
		ax.set_ylim(0, self.dimLim.y)
		pointData = self.getPoints()
		xdat = np.array([p[0] for p in pointData])
		ydat = np.array([p[1] for p in pointData])
		self.scat = ax.scatter(xdat,ydat,
						s=619*self.systemConstants["interactionlen"]**2,color = 'black',edgecolor= (1,1,1,0.5))
		animation = FuncAnimation(fig,self.update,frames = 1,interval=1,repeat=False)
		plt.show()

	def resetForceBuffer(self):
		for particle in self.particleSet:
			particle.fext = Vec2(0.,0.)

	def forceSum(self):
		for particle in self.particleSet:
			particle.acc = particle.fext / particle.particleVariables["mass"]

	def constructContacts(self):
		self.wipeNeighborList()
		for gridPair in self.gridListAll:
			(gridX,gridY) = gridPair
			for i in [-1,0,1]:
				for j in [-1,0,1]:
					if gridX + i < 0 or gridY + j < 0:
						continue
					if gridX + i >= self.gridLim[0] or gridY + j >= self.gridLim[1]:
						continue
					self.collideCells((gridX,gridY),(gridX+i,gridY+j))

		# for particle in self.particleSet:
			# print("neighbor list for particle {}".format(particle.pID))
			# print(particle.neighborList)

	def collideParticles(self,particle1,particle2):
		if particle1 is particle2 :
			# Enable if you want to not include the particle as a neighbor for itself
			# return False
			particle1.neighborList.append(particle1)
			self.pairsData[(particle1,particle1)]  \
			    = ParticlePair(particle1,particle1,Vec2(0,0),Vec2(0,0),0,Vec2(0,0))
			return True

		relpos = (particle1.pos - particle2.pos)
		relvel = (particle1.vel - particle2.vel)
		dist = relpos.length()		
		reldir = relpos.dir()
		print("reldir:")
		print(reldir)
		if dist <= self.systemConstants["interactionlen"]:
			if self.pairsData.get((particle1,particle2)) is None :

				particle1.neighborList.append(particle2)
				particle2.neighborList.append(particle1)
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

	def wipeHash(self):
		for i in range(0,self.gridLim[0]):
			for j in range(0,self.gridLim[1]):
				if len(self.hashData[i][j]) is not 0:
					self.hashData[i][j] = []

	def wipeNeighborList(self):
		for particle in self.particleSet:
			particle.neighborList = []

	def hash(self):
		self.wipeHash()
		for particle in self.particleSet:
			particle.hash = (int(particle.pos.x // self.hashGridSize),
						     int(particle.pos.y // self.hashGridSize))
			print(particle.pos)
			print("----------")
			self.hashData[particle.hash[0]][particle.hash[1]].append(particle)


	# def placeParticles(self,radius):

