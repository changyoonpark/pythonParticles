import cProfile
import sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib.animation import FuncAnimation
from matplotlib import gridspec
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
		self.vel = vel
		self.acc = Vec2(0,0)
		self.particleVariables = particleVariables
		self.hash = (-1,-1)
		self.neighborList = []

	def __str__(self):
		s =  "---------------------------------------\n"
		s += "Particle : {}\n".format(self.pID)		
		s += "position : {}\n".format(self.pos)
		s += "velocity : {}\n".format(self.vel)
		s += "accelera : {}\n".format(self.particleVariables["f_p"] /  self.particleVariables["mass"])
		s += "	>>> Particle Variable Dictionary <<<\n"
		for key in self.particleVariables:
			s += "{} : {}\n".format(key,self.particleVariables[key])
		return s


class ParticleSystem:
	def __init__ (self,interactionAlgo,		       
		               particleInitData,
		               boundaryInitData):


		self.systemConstants = dict()
		self.particleSet = []
		self.dimLim = Vec2(None,None,tup = particleInitData.systemConstants["domain"])
		self.scat = None
		self.nodeScat = None
		self.animation = None
		self.hashData = []
		self.tstep = 0
		self.interactionAlgo = interactionAlgo


		print("Initializing matplotlib")

		# self.fig, (self.ax,self.ax2) = plt.subplots(figsize=(12, 8),nrows = 1,ncols = 1,gridspec_kw = {'height_ratios' : [2,1]})		
		self.fig, (self.ax, self.ax2) = plt.subplots(figsize=(12, 12),nrows = 2,ncols = 1)		
		# self.fig = plt.figure(1,)
		# self.fig = plt.subplots(nrows = 2,ncols = 1)		
		# self.ax = self.fig.add_subplot(gs[0])
		# self.ax2 = self.fig.add_subplot(gs[0])

		self.ax.set_xlim(0, self.dimLim.x)
		self.ax.set_ylim(0, self.dimLim.y)
		self.ax.set_aspect('equal')

		self.ax2.set_xlim(0, self.dimLim.x)
		self.ax2.set_ylim(0, self.dimLim.y)
		self.ax2.set_aspect('equal')
		# self.ax2 = self.fig2.add_axes([0, 0, 1, 1], frameon=True)		

		print("Reading initial particle data")
		num = 0
		self.systemConstants = particleInitData.systemConstants			
		posDat = particleInitData.posDat
		velDat = particleInitData.velDat


		for idx in range(0,len(posDat)) :
			newParticle = Particle(num,posDat[idx],velDat[idx],Vec2(0,0),
								   dict(particleInitData.particleVariables)
								   )
			newParticle.particleVariables["isBoundary"] = False
			# newParticle.particleVariables["pressure"] = newParticle.pos.x
			newParticle.particleVariables["a_external"] = particleInitData.a_external[idx]
			# newParticle.particleVariables["a_external"] = Vec2(0,0)
			self.particleSet.append(newParticle)
			num += 1


		# print("Initializing Matrix")
		# self.systemConstants["aMatrix"] = np.zeros((len(self.particleSet),len(self.particleSet)),dtype = float)
		# self.systemConstants["pVector"] = np.zeros((len(self.particleSet),1),dtype = float)
		# self.systemConstants["bVector"] = np.zeros((len(self.particleSet),1),dtype = float)


		print("Initializing boundary particle data")
		bPosDat = boundaryInitData.posDat
		for idx in range(0,len(bPosDat)) :
			newBoundaryParticle = Particle(num,bPosDat[idx],Vec2(0,0),Vec2(0,0),
								   dict())
			newBoundaryParticle.particleVariables["isBoundary"] = True
			newBoundaryParticle.vel = Vec2(0,0)
			self.particleSet.append(newBoundaryParticle)
			num += 1


		print("Hashing Particle Data")
		self.hashGridSize = self.systemConstants["interactionlen"] * 1.5
		self.gridLim = (int(self.dimLim.x // self.hashGridSize) + 1,
			            int(self.dimLim.y // self.hashGridSize) + 1)
		gridNum = self.gridLim[0] * self.gridLim[1] // 4

		for i in range(0,self.gridLim[0]):
			self.hashData.append([])
			for j in range(0,self.gridLim[1]):
				self.hashData[-1].append([])

		self.gridListAll = []
		self.pairsData = dict()

		for gridX in range(0,self.gridLim[0]):
			for gridY in range(0,self.gridLim[1]):
				# self.gridList[foo].append((gridX,gridY))
				self.gridListAll.append((gridX,gridY))
				# if len(self.gridList[foo]) % gridNum == 0 :
					# foo += 1
					# if foo > 3:
						# foo = 3



	def getPoints(self):
		pointData = []
		for particle in self.particleSet:
			pointData.append(np.array([particle.pos.x,particle.pos.y]))
		return pointData

	def getPressure(self):
		pData = []
		pdat = []

		pmax = -1
		pmin = 9E20

		cmap = cm.get_cmap('rainbow')

		for particle in self.particleSet:
			if particle.particleVariables["isBoundary"] is False :
				if particle.particleVariables["pressure"] > pmax:
					pmax = particle.particleVariables["pressure"]
				if particle.particleVariables["pressure"] < pmin:
					pmin = particle.particleVariables["pressure"]

				pData.append(particle.particleVariables["pressure"])
			else:
				pData.append(None)


		for p in pData:
			if p is None : 
				pdat.append((0.5,0.5,0.5))
			else:
				if p == 0:
					pdat.append((0,0,0))
				else:
					normalize = ((p-pmin)/(pmax-pmin+0.001))
					pdat.append(cmap(normalize))

		return pdat

	def getColorGrad(self):
		pData = []
		pdat = []

		gradmax = -1
		gradmin = 9E20

		cmap = cm.get_cmap('rainbow')

		for particle in self.particleSet:
			if particle.particleVariables["isBoundary"] is False :
				if particle.particleVariables.get("colorGrad") is None:
					break


				# if particle.particleVariables["colorGrad"] > gradmax:
				# 	gradmax = particle.particleVariables["colorGrad"]
				# if particle.particleVariables["colorGrad"] < gradmin:
				# 	gradmin = particle.particleVariables["colorGrad"]

				pData.append(particle.particleVariables["colorGradAfterGapClosing"])

				# pData.append(particle.particleVariables["colorGrad"])
			else:
				pData.append(None)


		for grad in pData:
			if grad is None : 
				pdat.append((0.5,0.5,0.5))
			else:
				if grad == 0:
					pdat.append((0,0,0))
				else:
					normalize = ((grad-gradmin)/(gradmax-gradmin+0.001))
					pdat.append(cmap(normalize))

		return pdat


	def getDensity(self):
		densData = []
		ddat = []

		dmax = -1
		dmin = 9999

		for particle in self.particleSet:
			if particle.particleVariables["isBoundary"] is False :
				if particle.particleVariables["rho"] > dmax:
					dmax = particle.particleVariables["rho"]
				if particle.particleVariables["rho"] < dmin:
					dmin = particle.particleVariables["rho"]

				densData.append(particle.particleVariables["rho"])
			else:
				densData.append(None)


		for d in densData:
			if d is None : 
				ddat.append((1,0,0))
			else:
				ddat.append((1-(d-dmin)/(dmax-dmin+0.01),1-(d-dmin)/(dmax-dmin+0.01),1-(d-dmin)/(dmax-dmin+0.01)))

		return ddat

	def resetPairList(self):
		self.pairsData.clear()

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

		pdat = self.getPoints()
		# ddat = self.getPressure()
		ddat = self.getColorGrad()
		# ddat = self.getDensity()

		print("time : {}".format(t))
		toc()

		# if self.tstep % 20 == 0 :
			# plt.savefig("plots/"+str(self.tstep//20)+".png");
			# print("Exporting Plot")
		# figg = plt.figure()
		# axx = figg.add_subplot(111)
		# axx.spy(self.systemConstants["aMatrix"])
		# print(self.systemConstants["aMatrix"].max())
		# plt.show()

		self.tstep += 1
		self.scat.set_offsets(pdat)
		self.scat.set_color(ddat)

		if self.systemConstants.get("grid") is not None:
			cmap = cm.get_cmap('rainbow')
			nodePDat = np.array([[self.systemConstants["grid"].nodes[nodeKey].nodePos.x,self.systemConstants["grid"].nodes[nodeKey].nodePos.y] for nodeKey in self.systemConstants["grid"].nodes ])
			nodeColorDat = np.array([self.systemConstants["grid"].nodes[nodeKey].quantities["color"] for nodeKey in self.systemConstants["grid"].nodes])
			# nodeColorDat = np.array([self.systemConstants["grid"].nodes[nodeKey].quantities["nonMaxSuppressedColorGradIntensity"] for nodeKey in self.systemConstants["grid"].nodes])
			# nodeColorDat = np.array([self.systemConstants["grid"].nodes[nodeKey].quantities["colorGradIntensity"] for nodeKey in self.systemConstants["grid"].nodes])
			nodeColorDat *= 1./nodeColorDat.max()
			nodeColorDat = [cmap(nodeColor) for nodeColor in nodeColorDat]
			self.nodeScat.set_offsets(nodePDat)
			self.nodeScat.set_color(nodeColorDat)

			# for nodekey in self.systemConstants["grid"].nodes:
			# 	dx = 0.03 * self.systemConstants["grid"].nodes[nodekey].quantities["colorGrad"].x
			# 	dy = 0.03 * self.systemConstants["grid"].nodes[nodekey].quantities["colorGrad"].y
			# 	self.ax.arrow(self.systemConstants["grid"].nodes[nodekey].nodePos.x - dx,
			# 				  self.systemConstants["grid"].nodes[nodekey].nodePos.y - dy,
			# 				  2 * dx,
			# 				  2 * dy,
			# 				  width = 0.02,
			# 				  head_length = 0.02,
			# 				  head_width = 0.05,
			# 				  fc='k', ec='k'
			# 				  )





		# self.ax2.set_xlim(-1,int(len(self.systemConstants["densityDeviation"])))
		# maxdev = 0
		# for dat in self.systemConstants["densityDeviation"]:
			# if maxdev < dat[1]:
				# maxdev = dat[1]   
		# self.ax2.set_ylim(0,maxdev*1.1)

		# plt.title('{} Iterations, Average Density Deviation {}'.format(len(self.systemConstants["densityDeviation"]),self.systemConstants["densityDeviation"][-1][1]))
		# self.densityDeviation.set_offsets(self.systemConstants["densityDeviation"])
		plt.savefig('plots/rotdisc/rotdisc_{}.png'.format(self.tstep))

	def solveTimeStep(self):
		for operation in self.interactionAlgo.getAlgoProcedure(self.tstep):
			operation(self.systemConstants,self.pairsData,self.particleSet)

	def run(self):

		pointData = self.getPoints()
		# densData = self.getDensity()
		densData = self.getPressure()

		xdat = np.array([p[0] for p in pointData])
		ydat = np.array([p[1] for p in pointData])

		self.scat = self.ax.scatter(xdat,ydat,c=densData,
						s=100*self.systemConstants["interactionlen"]**2,edgecolor= (1,1,1,0.5))

		self.nodeScat = self.ax2.scatter([],[],c=[],
						s=160*self.systemConstants["interactionlen"]**2, alpha=1.0, lw=0, marker='s')
		# self.nodeScat = self.ax.scatter([],[],c=[],
		# 				s=30*self.systemConstants["interactionlen"]**2, alpha=0.9, marker='s', lw=0)

		# self.densityDeviation = self.ax2.scatter([0],[1],linewidths = 0.1)

		# self.scat = self.ax.scatter(xdat,ydat,color='Black',
		# 				s=300*self.systemConstants["interactionlen"]**2,edgecolor= (1,1,1,0.5))
		animation = FuncAnimation(self.fig,self.update,frames = 1,interval=999999,repeat=False)
		# animation = FuncAnimation(self.fig,self.update)
		plt.show()

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
		if particle1 is particle2:
			if  self.pairsData.get((particle1,particle2)) is None:
				# Enable if you want to not include the particle as a neighbor for itself
				# return False
				particle1.neighborList.append(particle1)
				self.pairsData[(particle1,particle1)]  \
				    = ParticlePair(particle1,particle1,Vec2(0,0),Vec2(0,0),0,Vec2(0,0))
				return True
			else:
				return False

		relpos = (particle1.pos - particle2.pos)
		relvel = (particle1.vel - particle2.vel)
		dist = relpos.length()		
		reldir = relpos.dir()
		# print("reldir:")
		# print(reldir)
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
			# print(particle.pos)
			# print("----------")
			self.hashData[particle.hash[0]][particle.hash[1]].append(particle)





#  NOT USED.

	# def rk2(self,posBeforeUpdate,velBeforeUpdate):
	# 	k1 = (velBeforeUpdate,
	# 		  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

	# 	intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.systemConstants["dt"]))
	# 	intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.systemConstants["dt"]))
	# 	k2 = (intermediateVel,
	# 		  self.calculateForces(intermediatePos,intermediateVel))

	# 	l = len(k1)
	# 	for i in range(0,l):
	# 		self.particleSet[i].pos = posBeforeUpdate[i] + \
	# 					self.systemConstants["dt"]*(k1[0][i])
	# 		self.particleSet[i].vel = velBeforeUpdate[i] + \
	# 					self.systemConstants["dt"]*(k1[1][i])

	# def rk4(self,posBeforeUpdate,velBeforeUpdate):
	# 	k1 = (velBeforeUpdate,
	# 		  self.calculateForces(posBeforeUpdate,velBeforeUpdate))

	# 	intermediatePos = self.add(posBeforeUpdate,self.mult(k1[0],0.5*self.systemConstants["dt"]))
	# 	intermediateVel = self.add(velBeforeUpdate,self.mult(k1[1],0.5*self.systemConstants["dt"]))
	# 	k2 = (intermediateVel,
	# 		  self.calculateForces(intermediatePos,intermediateVel))

	# 	intermediatePos = self.add(posBeforeUpdate,self.mult(k2[0],0.5*self.systemConstants["dt"]))
	# 	intermediateVel = self.add(velBeforeUpdate,self.mult(k2[1],0.5*self.systemConstants["dt"]))
	# 	k3 = (intermediateVel,
	# 		  self.calculateForces(intermediatePos,intermediateVel))

	# 	intermediatePos = self.add(posBeforeUpdate,self.mult(k3[0],self.systemConstants["dt"]))
	# 	intermediateVel = self.add(velBeforeUpdate,self.mult(k3[1],self.systemConstants["dt"]))
	# 	k4 = (intermediateVel,
	# 		  self.calculateForces(intermediatePos,intermediateVel))
	# 	l = len(k1)
	# 	for i in range(0,l):
	# 		self.particleSet[i].pos = posBeforeUpdate[i] + \
	# 					(1.0/6.0)*self.systemConstants["dt"]*(k1[0][i]+2.0*k2[0][i]+2.0*k3[0][i]+k4[0][i])
	# 		self.particleSet[i].vel = velBeforeUpdate[i] + \
	# 					(1.0/6.0)*self.systemConstants["dt"]*(k1[1][i]+2.0*k2[1][i]+2.0*k3[1][i]+k4[1][i])

