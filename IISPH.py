# from interactionAlgorithm import *
from helpers import *

# Inherits from InteractionAlgorithm
class IISPH_Algorithm :

# particleSystem is the list of all the particles.
# systemConstants stores the constants used in the simulation.
# pairsData stores the pre-computed data regarding the neighbors.
	def __init__(self,W,gW,lW,omega = 0.5):
		# super().__init__(W,gW,lW)		
		# These operations are done in order, for the particlesystem.
		self.omega = omega
		self.W = W
		self.gW = gW
		self.lW = lW
		self._proceduresInOrder =  [
									   self.calculateDensity,
			  						   self.calculateAdvectionForce,
			  						   self.calculateWall,
			  						   self.calculateAdvectionVel,
			  						   self.calculate_dii,
			  						   self.calculateAdvectionDensity,
			  						   self.calculate_aii,
			  						   self.pressureSolver,
			  						   self.calculatePressureForce,
			  						   self.integration,
		 						   ]

	def getAlgoProcedure(self): 
		return self._proceduresInOrder

	def calculateDensity(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#   calculate density
			particle.particleVariables["rho"] = 0
			for neighbor in particle.neighborList:
				pairData = pairsData[(particle,neighbor)]
				particle.particleVariables["rho"] += \
					neighbor.particleVariables["mass"] * self.W(pairData,systemConstants["interactionlen"])
			# print("Density of particle {} : {}".format(particle.pID,particle.particleVariables["rho"]))

	def calculate_dii(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#   calculate dii
			particle.particleVariables["d_ii"] = Vec2(0,0)
			for neighbor in particle.neighborList:
				pairData = pairsData[(particle,neighbor)]
				particle.particleVariables["d_ii"] -= \
					neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
					self.gW(pairData,systemConstants["interactionlen"])

			particle.particleVariables["d_ii"] = particle.particleVariables["d_ii"] * systemConstants["dt"]**2

	def calculate_sum_pj_dij(self,particle,systemConstants,pairsData):
		particle.particleVariables["sum_pj_dij"] = Vec2(0,0)
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]
			print("dij : {}".format(self.calculate_dij(pairData,systemConstants)))
			particle.particleVariables["sum_pj_dij"] += self.calculate_dij(pairData,systemConstants) * neighbor.particleVariables["pressure"]
		if particle.pID == 10:
			print("sum_pj_dij for particle {} : {}".format(particle.pID,particle.particleVariables["sum_pj_dij"]))	

	def calculate_dij(self,pairData,systemConstants):		
		return (-systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
			   (pairData.particlej.particleVariables["rho"]**2)) * self.gW(pairData,systemConstants["interactionlen"])

	def calculateAdvectionForce(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#	calculate advection force (viscosity and gravity is applied.)
			particle.particleVariables["f_adv"] = Vec2(0,0)
			for neighbor in particle.neighborList:
				pairData  = pairsData[(particle,neighbor)]
				pairData2 = pairsData[(neighbor,particle)]
				particle.particleVariables["f_adv"] += \
					systemConstants["viscosity"] * \
					self.lW(pairData,systemConstants["interactionlen"]) * \
					(neighbor.particleVariables["mass"] / \
					 neighbor.particleVariables["rho"]) * \
					 pairData2.relvel
 		#   consider gravity.
			particle.particleVariables["f_adv"] -= \
				particle.particleVariables["mass"] * systemConstants["gravity"] * Vec2(0,-1)

	def calculateWall(self,systemConstants,pairsData,particleSet):
		# print(systemConstants)
		walls = systemConstants["walls"]
		domain = systemConstants["domain"]
		for particle in particleSet:			
			if particle.pos.y   < walls[1] : 
				particle.particleVariables["f_adv"].y = particle.fext.y + 300 * particle.particleVariables["mass"] * pow(walls[0][1] - particle.pos.y,2)
				particle.particleVariables["f_adv"].y = particle.fext.y - particle.particleVariables["mass"] * 0.2 * particle.vel.y 
			elif particle.pos.y > domain[1] - walls[1] : 
				particle.particleVariables["f_adv"].y = particle.fext.y - 300 * particle.particleVariables["mass"] * pow(particle.pos.y - walls[1][1],2)
				particle.particleVariables["f_adv"].y = particle.fext.y - particle.particleVariables["mass"] * 0.2 * particle.vel.y
			if particle.pos.x   < walls[0] :
				particle.particleVariables["f_adv"].x = particle.fext.x + 300 * particle.particleVariables["mass"] * pow(walls[0][0] - particle.pos.x,2)
				particle.particleVariables["f_adv"].x = particle.fext.x - particle.particleVariables["mass"] * 0.2 * particle.vel.x
			elif particle.pos.x > domain[1] - walls[0] : 
				particle.particleVariables["f_adv"].x = particle.fext.x - 300 * particle.particleVariables["mass"] * pow(particle.pos.x - walls[1][0],2)
				particle.particleVariables["f_adv"].x = particle.fext.x - particle.particleVariables["mass"] * 0.2 * particle.vel.x


	def calculateAdvectionVel(self,systemConstants,pairsData,particleSet):
	#   calculate itermediate velocity
		for particle in particleSet:
			particle.particleVariables["v_adv"] = \
				particle.vel + \
				systemConstants["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]	

	def calculateAdvectionDensity(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#   calculate advection density
			particle.particleVariables["rho_adv"] = particle.particleVariables["rho"]
			for neighbor in particle.neighborList:
				pairData = pairsData[(particle,neighbor)]
				v_adv_ij = particle.particleVariables["v_adv"] - neighbor.particleVariables["v_adv"]
				particle.particleVariables["rho_adv"] += \
					systemConstants["dt"] * neighbor.particleVariables["mass"] * \
					v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))

	def calculate_aii(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:	
		# 	calculate a_ii for relaxed jacobi
			particle.particleVariables["a_ii"] = 0
			for neighbor in particle.neighborList:
				pairData  = pairsData[(particle,neighbor)]
				pairData2 = pairsData[(neighbor,particle)] 
				particle.particleVariables["a_ii"] += \
					neighbor.particleVariables["mass"] * \
					(particle.particleVariables["d_ii"] - self.calculate_dij(pairData2,systemConstants)) \
					.dot(self.gW(pairData,systemConstants["interactionlen"]))
			# print("a_ii for particle {} : {}".format(particle.pID,particle.particleVariables["a_ii"]))

	def calcAverageDensity(self,particleSet):
		rhoAve = 0
		for particle in particleSet:
			rhoAve += particle.particleVariables["rho"]
		return rhoAve / len(particleSet)

	def estimatePressure(self,particle,systemConstants,pairsData,particleSet):
		junk = 0 
		for neighbor in particle.neighborList:
			pairData  = pairsData[(particle,neighbor)]
			pairData2 = pairsData[(neighbor,particle)]
			junk += neighbor.particleVariables["mass"] * \
					(  particle.particleVariables["sum_pj_dij"] - neighbor.particleVariables["d_ii"] * neighbor.particleVariables["pressure"]\
					 - neighbor.particleVariables["sum_pj_dij"] + self.calculate_dij(pairData2,systemConstants) * particle.particleVariables["pressure"])\
					 .dot(self.gW(pairData,systemConstants["interactionlen"]))

		return (1-self.omega) * particle.particleVariables["pressure"] + \
		   		self.omega * (1 / particle.particleVariables["a_ii"] ) * \
		   		(systemConstants["rho0"] - particle.particleVariables["rho_adv"] - junk)


	def pressureSolver(self,systemConstants,pairsData,particleSet):
	# 	iteration counter
		l = 0
		#   initialize pressure estimate.
		for particle in particleSet:
			particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"] 

		densityDeviation = 1
		while ( densityDeviation > 0.02 and l < 100):

			print("Performing Jacobi iteration for pressure field. Iteration : {} / Average Density Deviation : {}".format(l,densityDeviation))

			for particle in particleSet:
				self.calculate_sum_pj_dij(particle,systemConstants,pairsData)
				# print("sumpjdij")
				# print(particle.particleVariables["sum_pj_dij"])

			for particle in particleSet:
				particle.particleVariables["pressure_est"] = self.estimatePressure(particle,systemConstants,pairsData,particleSet)
				# print(particle.particleVariables["pressure_est"])					


			for particle in particleSet:
				particle.particleVariables["pressure"] = particle.particleVariables["pressure_est"]
				particle.particleVariables["rho"] = systemConstants["rho0"] + particle.particleVariables["pressure"] / 10000
				# if particle.pID == 10:
					# print(particle)

			densityDeviation = abs((self.calcAverageDensity(particleSet) - systemConstants["rho0"])/systemConstants["rho0"])
			l += 1

	def calculatePressureForce(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
			particle.particleVariables["fp_dt**2/m"] = particle.particleVariables["d_ii"] * particle.particleVariables["pressure"] + \
												particle.particleVariables["sum_pj_dij"]

	def integration(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
			particle.vel =  particle.particleVariables["v_adv"] + particle.particleVariables["fp_dt**2/m"] / systemConstants["dt"]
			particle.pos += systemConstants["dt"] * particle.vel











