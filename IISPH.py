from interactionAlgorithm import *


# Inherits from InteractionAlgorithm
class IISPH_Algorithm(InteractionAlogorithm):

# particleSystem is the list of all the particles.
# systemConstants stores the constants used in the simulation.
# pairsData stores the pre-computed data regarding the neighbors.
	def __init__(self,systemConstants,pairsData,particleSystem,W,gW,lW,omega = 0.5):
		super().__init__(systemConstants,pairsData,particleSystem,W,gW,lW)		
		# These operations are done in order, for the particlesystem.
		self.omega = omega
		self._proceduresInOrder = [self.calculateDensity,
		 		 				   self.calculate_dii,
		 	 					   self.calculate_sum_pj_dij,
		  						   self.calculateAdvectionForce,
		  						   self.calculateAdvectionVel,
		  						   self.calculateAdvectionDensity,
		  						   self.calculate_aii,
		  						   self.pressureSolver,
		 						   ]


	def getAlgoProcedure(self):
		return self._proceduresInOrder

	def calculateDensity(self):
		for particle in self.particleSystem:
		#   calculate density
			particle.particleVariables["rho"] = 0
			for neighbor in particle.neighborList:
				pairData = self.pairsData[(particle,neighbor)]
				particle.particleVariables["rho"] += \
					particle.particleVariables["mass"] * self.W(pairData,systemConstants["interactionlen"])

	def calculate_dii(self):
		for particle in self.particleSystem:
		#   calculate dii
			particle.particleVariables["d_ii"] = Vec2(0,0)
			for neighbor in particle.neighborList:
				pairData = self.pairsData[(particle,neighbor)]
				particle.particleVariables["d_ii"] -= \
					neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
					self.gW(pairData,systemConstants["interactionlen"])

			particle.particleVariables["d_ii"] = particle.particleVariables["d_ii"] * self.systemConstants["dt"]**2

	def calcaulate_sum_pj_dij(self):
		for particle in self.particleSystem:
			particle.particleVariables["sum_pj_dij"] = Vec2(0,0)
			for neigbor in particle.neighborList:
				pairData = self.pairsData[(particle,neighbor)]
				particle.particleVariables["sum_pj_dij"] += self.calcaulate_dij(pairData) * neighbor.particleVariables["pressure"]

	def calculate_dij(pairData):
		return -self.systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
				(pairData.particlej.particleVariables["rho"]**2) * self.gW(pairData,self.systemConstants["interactionlen"])
		#   calculate the sum of dij's

	def calcaulateAdvectionForce(self):
		for particle in self.particleSystem:
		#	calculate advection force (viscosity and gravity is applied.)
			particle.particleVariables["f_adv"] = Vec2(0,0)
			for neighbor in particle.neighborList:
				pairData = self.pairsData[(particle,neighbor)]
				particle.particleVariables["f_adv"] += \
					self.systemConstants["viscosity"] * \
					self.lW(pairData,self.systemConstants["interactionlen"]) * \
					neighbor.particleVariables["mass"] / \
					neighbor.particleVariables["rho"] * \
					neighbor.relvel
 		#   consider gravity.
			particle.particleVariables["f_adv"] -= \
				particle.particleVariables["mass"] * self.systemConstants["gravity"]

	def calcaulteAdvectionVel(self):
	#   calculate itermediate velocity
		for particle in self.particleSystem:
			particle.particleVariables["v_adv"] = \
				particle.vel + \
				particle.particleVariables["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]	

	def calculateAdvectionDensity(self):
		for particle in self.particleSystem:
		#   calculate advection density
			particle.particleVariables["rho_adv"] = particle.particleVariables["rho"]
			for neighbor in particle.neighborList:
				pairData = self.pairsData[(particle,neighbor)]
				v_adv_ij = particle.particleVariables["v_adv"] - neighbor.particleVariables["v_adv"]
				particle.particleVariables["rho_adv"] += \
					self.systemConstants["dt"] * neighbor.particleVariables["mass"] * \
					v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))

		#   initialize pressure estimate.
			particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"] 

	def calculate_aii(self):
		for particle in self.particleSystem:	
		# 	calculate a_ii for relaxed jacobi
			particle.particleVariables["a_ii"] = 0
			for neighbor in particle.neighborList:
				pairData  = self.pairsData[(particle,neighbor)]
				pairData2 = self.pairsData[(neighbor,particle)] 
				particle.particleVariables["a_ii"] += \
					neighbor.particleVariables["mass"] * \
					(particle.particleVariables["d_ii"] - self.calculate_dij(pairData2)) \
					.dot(self.gW(pairData,self.systemConstants["interactionlen"]))

	def calcAverageDensity(self):
		rhoAve = 0
		for particle in self.particleSystem:
			rhoAve += particle.particleVariables["rho"]
		return rhoAve / len(self.particleSystem)

	def estimatePressure(self,particle):
		junk = 0 
		for neighbor in particle.neighborList:
			pairData  = self.pairsData[(particle,neighbor)]
			pairData2 = self.pairsData[(neighbor,particle)]
			junk += neighbor.particleVariables["mass"] * \
					(  particle.particleVariables["sum_pj_dij"] - neighbor.particleVariables["d_ii"] * neighbor.particleVariables["pressure"]\
					 - neighbor.particleVariables["sum_pj_dij"] + self.calculate_dij(pairData2) * particle.particleVariables["pressure"])\
					 .dot(self.gW(pairData,self.systemConstants["interactionlen"]))

		 return (1-self.omega) * particle.particleVariables["pressure"] + \
		   		omega * (1/particle.particleVariables["a_ii"]) * \
		   		(self.systemConstants["rho0"] - particle.particleVariables["rho_adv"] - junk)


	def pressureSolver(self):
	# 	iteration counter
		l = 0
		while ( abs((self.calcAverageDensity() - self.systemConstants["rho0"])/self.systemConstants["rho0"]) > 0.02 ):
			for particle in self.particleSystem:
				self.estimatePressure(particle)












