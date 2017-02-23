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
		self.currentTime = 0
		self._proceduresInOrder =  [
									   self.calculateDensity,
			  						   self.calculateAdvectionForce,
			  						   self.calculateAdvectionVel,
			  						   self.calculate_dii,
			  						   self.calculateAdvectionDensity,
			  						   self.calculate_aii,
			  						   self.pressureSolver,
			  						   self.integration,
		 						   ]


	def getAlgoProcedure(self,t):
		if t == 0 : 
			return [self.initializeBoundaryParticleVariables] + self._proceduresInOrder
		else: 
			return self._proceduresInOrder

	def initializeBoundaryParticleVariables(self,systemConstants, pairsData, particleSet):
		for particle in particleSet:

			if particle.particleVariables["isBoundary"] is True:
				particle.particleVariables["nDens"] = 0

				for neighbor in particle.neighborList:
					if neighbor.particleVariables["isBoundary"] is True:
						pairData = pairsData[(particle,neighbor)]
						particle.particleVariables["nDens"] += self.W(pairData,systemConstants["interactionlen"])
				particle.particleVariables["psi"] = systemConstants["rho0"] / particle.particleVariables["nDens"]
				# print(particle.particleVariables["psi"])



	def calculateDensity(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
			# Dont include boundary particles.

			if particle.particleVariables["isBoundary"] is False: 
				particle.particleVariables["rho"] = 0
				for neighbor in particle.neighborList:
					pairData = pairsData[(particle,neighbor)]
					if neighbor.particleVariables["isBoundary"] is False:
						particle.particleVariables["rho"] += \
							neighbor.particleVariables["mass"] * self.W(pairData,systemConstants["interactionlen"])	
					else:
						particle.particleVariables["rho"] += \
							neighbor.particleVariables["psi"] * self.W(pairData,systemConstants["interactionlen"])

				# particle.particleVariables["rho"] = max(particle.particleVariables["rho"], systemConstants["rho0"])

	def calculate_dii(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#   calculate dii
			if particle.particleVariables["isBoundary"] is False: 

				particle.particleVariables["d_ii"] = Vec2(0,0)

				for neighbor in particle.neighborList:
					pairData = pairsData[(particle,neighbor)]

					if neighbor.particleVariables["isBoundary"] is False:
						particle.particleVariables["d_ii"] -= \
							neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
							self.gW(pairData,systemConstants["interactionlen"])
					else:
						particle.particleVariables["d_ii"] -= \
							neighbor.particleVariables["psi"] / (particle.particleVariables["rho"]**2) * \
							self.gW(pairData,systemConstants["interactionlen"])	

				particle.particleVariables["d_ii"] = particle.particleVariables["d_ii"] * systemConstants["dt"]**2

	def calculate_sum_pj_dij(self,particle,systemConstants,pairsData):
		if particle.particleVariables["isBoundary"] is False:

			particle.particleVariables["sum_pj_dij"] = Vec2(0,0)
			for neighbor in particle.neighborList:
				if neighbor.particleVariables["isBoundary"] is False:
					pairData = pairsData[(particle,neighbor)]
					particle.particleVariables["sum_pj_dij"] += self.calculate_dij(pairData,systemConstants) * neighbor.particleVariables["pressure"]

	def calculate_dij(self,pairData,systemConstants):		
		return (-systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
			   (pairData.particlej.particleVariables["rho"]**2)) * self.gW(pairData,systemConstants["interactionlen"])


	def calculateAdvectionForce(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#	calculate advection force (viscosity and gravity is applied.)
			if particle.particleVariables["isBoundary"] is False:

				particle.particleVariables["f_adv"] = Vec2(0,0)
				for neighbor in particle.neighborList:
					pairData  = pairsData[(particle,neighbor)]
					pairData2 = pairsData[(neighbor,particle)]

					if neighbor.particleVariables["isBoundary"] is False :
						particle.particleVariables["f_adv"] += \
							systemConstants["viscosity"] * \
							self.lW(pairData,systemConstants["interactionlen"]) * \
							(neighbor.particleVariables["mass"] / \
							 neighbor.particleVariables["rho"]) * \
							 pairData2.relvel

	 		#   consider gravity.
				particle.particleVariables["f_adv"] += \
					particle.particleVariables["mass"] * systemConstants["gravity"] * Vec2(0,-1)

	def calculateAdvectionVel(self,systemConstants,pairsData,particleSet):
	#   calculate itermediate velocity
		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				particle.vel += systemConstants["dt"] * particle.particleVariables["f_adv"] / particle.particleVariables["mass"]

	def calculateAdvectionDensity(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
		#   calculate advection density
			if particle.particleVariables["isBoundary"] is False:
				particle.particleVariables["rho_adv"] = particle.particleVariables["rho"]

				for neighbor in particle.neighborList:

					pairData = pairsData[(particle,neighbor)]
					v_adv_ij = particle.vel - neighbor.vel

					if neighbor.particleVariables["isBoundary"] is False:
						particle.particleVariables["rho_adv"] += \
							systemConstants["dt"] * neighbor.particleVariables["mass"] * \
							v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))
					else:
						particle.particleVariables["rho_adv"] += \
							systemConstants["dt"] * neighbor.particleVariables["psi"] * \
							v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))


	def calculate_aii(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:	
			if particle.particleVariables["isBoundary"] is False:
			# 	calculate a_ii for relaxed jacobi
				particle.particleVariables["a_ii"] = 0
				for neighbor in particle.neighborList:
					if neighbor.particleVariables["isBoundary"] is False:
						pairData  = pairsData[(particle,neighbor)]
						pairData2 = pairsData[(neighbor,particle)] 
						particle.particleVariables["a_ii"] += \
						  	 neighbor.particleVariables["mass"] * \
							(particle.particleVariables["d_ii"] - self.calculate_dij(pairData2,systemConstants)) \
							.dot(self.gW(pairData,systemConstants["interactionlen"]))

				# particle.particleVariables["a_ii"] = max(particle.particleVariables["a_ii"],0)


	def calcAverageDensity(self,particleSet):
		rhoAve = 0
		numP = 0
		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				rhoAve += particle.particleVariables["rho_est"]
				numP += 1

		return rhoAve / numP

	def estimateDensity(self,particle,pairsData,systemConstants):
		# junk = 0
		if particle.particleVariables["isBoundary"] is False:
			if particle.particleVariables["pressure_est"] > 0:
				particle.particleVariables["rho_est"] = 0

				for neighbor in particle.neighborList:
					# pairData2 = pairsData[(neighbor,particle)]

					# if neighbor.particleVariables["isBoundary"] is False:
					# 	junk +=    neighbor.particleVariables["mass"] * \
					# 			(  particle.particleVariables["sum_pj_dij"] - neighbor.particleVariables["d_ii"] * neighbor.particleVariables["pressure"]\
					# 			 - neighbor.particleVariables["sum_pj_dij"] + self.calculate_dij(pairData2,systemConstants) * particle.particleVariables["pressure"])\
					# 			 .dot(self.gW(pairData,systemConstants["interactionlen"]))
					# else:
					# 	junk +=    neighbor.particleVariables["psi"] * \
					# 			  (particle.particleVariables["sum_pj_dij"]).dot(self.gW(pairData,systemConstants["interactionlen"]))

				# particle.particleVariables["rho_est"] = (particle.particleVariables["a_ii"] * particle.particleVariables["pressure"] + junk)\
														 # + particle.particleVariables["rho_adv"]

					pairData  = pairsData[(particle,neighbor)]
					if neighbor.particleVariables["isBoundary"] is False:
						particle.particleVariables["rho_est"] += neighbor.particleVariables["mass"] * (particle.particleVariables["f_p"] / particle.particleVariables["mass"] - neighbor.particleVariables["f_p"] / neighbor.particleVariables["mass"])\
						.dot(self.gW(pairData,systemConstants["interactionlen"]))
					else:
						particle.particleVariables["rho_est"] += neighbor.particleVariables["psi"] * (particle.particleVariables["f_p"] / particle.particleVariables["mass"])\
						.dot(self.gW(pairData,systemConstants["interactionlen"]))

				particle.particleVariables["rho_est"] = particle.particleVariables["rho_est"] * (systemConstants["dt"] ** 2)
				particle.particleVariables["rho_est"] += particle.particleVariables["rho_adv"]

			else:
				particle.particleVariables["rho_est"] = systemConstants["rho0"]


	def estimatePressure(self,particle,systemConstants,pairsData,particleSet):
		junk = 0
		if particle.particleVariables["isBoundary"] is False: 
			for neighbor in particle.neighborList:
				pairData  = pairsData[(particle,neighbor)]
				pairData2 = pairsData[(neighbor,particle)]

				if neighbor.particleVariables["isBoundary"] is False:
					junk +=    neighbor.particleVariables["mass"] * \
							(  particle.particleVariables["sum_pj_dij"] - neighbor.particleVariables["d_ii"] * neighbor.particleVariables["pressure"]\
							 - neighbor.particleVariables["sum_pj_dij"] + self.calculate_dij(pairData2,systemConstants) * particle.particleVariables["pressure"])\
							 .dot(self.gW(pairData,systemConstants["interactionlen"]))
				else:
					junk +=    neighbor.particleVariables["psi"] * \
							  (particle.particleVariables["sum_pj_dij"]).dot(self.gW(pairData,systemConstants["interactionlen"]))

			pres = (1-self.omega) * particle.particleVariables["pressure"] + \
		   	    	self.omega * (1 / particle.particleVariables["a_ii"] ) * \
		   	    	(systemConstants["rho0"] - particle.particleVariables["rho_adv"] - junk)


			# print(particle.particleVariables["a_ii"])
			if abs(particle.particleVariables["a_ii"]) > 1E-9 :
				pres =  max(0,pres)
			else:
				pres = 0
			# print(pres)

			particle.particleVariables["pressure_est"] = pres



	def pressureSolver(self,systemConstants,pairsData,particleSet):
	# 	iteration counter
		l = 0
		#   initialize pressure estimate.
		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"]

		densityDeviation = 1
		foo = 0
		while ( (densityDeviation > 0.001 and l < 100)):
			if (abs(foo - densityDeviation) < 1E-6):
				break

			for particle in particleSet:
				self.calculate_sum_pj_dij(particle,systemConstants,pairsData)

			for particle in particleSet:
				self.estimatePressure(particle,systemConstants,pairsData,particleSet)	

			for particle in particleSet:
				self.calculatePressureForce(particle,pairsData,systemConstants)

			for particle in particleSet:
				self.estimateDensity(particle,pairsData,systemConstants)

			for particle in particleSet:
				if particle.particleVariables["isBoundary"] is False:
					particle.particleVariables["pressure"] = particle.particleVariables["pressure_est"]

			foo = densityDeviation
			densityDeviation = abs((self.calcAverageDensity(particleSet) - systemConstants["rho0"])/systemConstants["rho0"])
			print("Performing Jacobi iteration for pressure field. Iteration : {} / Average Density Deviation : {}".format(l,densityDeviation))
			l += 1

		# for particle in particleSet:
			# if particle.particleVariables["isBoundary"] is False:
				# if particle.particleVariables["rho"] < systemConstants["rho0"] :
					# print(particle.particleVariables["pressure"])


					
	def calculatePressureForce(self,particle,pairsData,systemConstants):

		# return (-systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
		# 	   (pairData.particlej.particleVariables["rho"]**2)) * self.gW(pairData,systemConstants["interactionlen"])



		# if particle.particleVariables["isBoundary"] is False:


		# 	particle.particleVariables["f_p"] = Vec2(0,0)
		# 	for neighbor in particle.neighborList:

		# 		pairData = pairsData[(particle,neighbor)]

		# 		if neighbor.particleVariables["isBoundary"] is False:
		# 			particle.particleVariables["f_p"] += (self.calculate_dij(pairData,systemConstants) * neighbor.particleVariables["pressure"])

		# 	particle.particleVariables["f_p"] += particle.particleVariables["d_ii"] * particle.particleVariables["pressure"]
		# 	particle.particleVariables["f_p"] = particle.particleVariables["f_p"] * particle.particleVariables["mass"] / (systemConstants["dt"])**2
		# 	print("----")
		# 	print(particle.particleVariables["f_p"])

		if particle.particleVariables["isBoundary"] is False:

			particle.particleVariables["f_p"] = Vec2(0,0)
			for neighbor in particle.neighborList:

				pairData = pairsData[(particle,neighbor)]

				if neighbor.particleVariables["isBoundary"] is False:
					particle.particleVariables["f_p"] += \
					              -particle.particleVariables["mass"] * neighbor.particleVariables["mass"] * \
					              (particle.particleVariables["pressure"] / ( particle.particleVariables["rho"] ** 2 ) + \
					               neighbor.particleVariables["pressure"] / ( neighbor.particleVariables["rho"] ** 2 )) \
					              *self.gW(pairData,systemConstants["interactionlen"])

				else:					
					particle.particleVariables["f_p"] += \
					              -particle.particleVariables["mass"] * neighbor.particleVariables["psi"] * \
					              (particle.particleVariables["pressure"]/(particle.particleVariables["rho"]**2))*self.gW(pairData,systemConstants["interactionlen"])



	def integration(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:
			self.calculatePressureForce(particle,pairsData,systemConstants)

		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				particle.vel += particle.particleVariables["f_p"] / particle.particleVariables["mass"] * systemConstants["dt"]
				particle.pos += systemConstants["dt"] * particle.vel











