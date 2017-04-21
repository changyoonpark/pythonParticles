from helpers import *
from grid import *
import scipy as sp


class IISPH_Algorithm :

# particleSystem is the list of all the particles.
# systemConstants stores the constants used in the simulation.
# pairsData stores the pre-computed data regarding the neighbors.
	def __init__(self,W,gW,lW,omega = 0):
		# super().__init__(W,gW,lW)		
		# These operations are done in order, for the particlesystem.

		self.omega = omega
		self.W = W
		self.gW = gW
		self.lW = lW
		self.currentTime = 0
		self._proceduresInOrder =  [
									   self.initializeGrid,
									   self.rasterizeGrid,
									   self.calculateDensity,
			  						   self.calculateAdvectionVel,
			  						   self.calculateSourceTerm,
			  						   self.calculateDiagonal,
			  						   self.solveForPressure,
			  						   self.integration,

		 						   ]

	def initializeGrid(self,systemConstants, pairsData, particleSet):
		systemConstants["grid"] = Grid(particleSet, systemConstants["interactionlen"])


	def rasterizeGrid(self,systemConstants, pairsData, particleSet):
    	# Init nodes
		for key in systemConstants["grid"].nodes:
			node = systemConstants["grid"].nodes[key]
			node.quantities["color"] = 0
			node.quantities["pressure"] = 0
			node.quantities["density"] = 0
			node.quantities["colorGrad"] = Vec2(0,0)
    	# Iterate through the particles and rasterize onto nodes
		for particle in particleSet:
			hashVal = systemConstants["grid"].hashFunction(particle.pos)

			for i in range(-1,3):
				for j in range(-1,3):

					nodes = systemConstants["grid"].nodes
					node = nodes.get((hashVal[0]+i,hashVal[1]+j))
					if node is None :
						continue
					node.quantities["color"]     += node.W(particle.pos)
					# node.quantities["color"]     += 0.1
					node.quantities["colorGrad"] += node.gW(particle.pos)
					# node.quantities["pressure"]  += node.W(particle.pos) * particle.particleVariables["pressure"]
					# node.quantities["density"]   += node.W(particle.pos) * particle.particleVariables["rho"]

		for key in systemConstants["grid"].nodes:
			print(systemConstants["grid"].nodes[key].quantities["color"])

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



	def getAlgoProcedure(self,t):
		if t == 0 : 
			return [self.initializeBoundaryParticleVariables] + self._proceduresInOrder
		else: 
			return self._proceduresInOrder


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


				# print(particle.particleVariables["rho"])



	def calculateAdvectionAccel(self,systemConstants,pairsData,particleSet):
		for particle in particleSet:

			#	calculate advection force (viscosity and gravity is applied.)
			if particle.particleVariables["isBoundary"] is False:

				particle.particleVariables["a_adv"] = Vec2(0,0)

				for neighbor in particle.neighborList:
					pairData  = pairsData[(particle,neighbor)]
					pairData2 = pairsData[(neighbor,particle)]

					if neighbor.particleVariables["isBoundary"] is False and neighbor is not particle:						
						particle.particleVariables["a_adv"] += \
							( 8.0 * systemConstants["viscosity"] / particle.particleVariables["rho"] ) * \
							( neighbor.particleVariables["mass"] / neighbor.particleVariables["rho"] ) * \
							( pairData.relvel.dot(pairData.relpos) / pairData.relpos.dot(pairData.relpos) ) * \
							( self.gW(pairData,systemConstants["interactionlen"]) )

				particle.particleVariables["a_adv"] += systemConstants["gravity"] * Vec2(0,-1)

				particle.particleVariables["a_adv"] += particle.particleVariables["a_external"]


	def calculateAdvectionVel(self,systemConstants,pairsData,particleSet):

		self.calculateAdvectionAccel(systemConstants,pairsData,particleSet)

		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				particle.vel += systemConstants["dt"] * particle.particleVariables["a_adv"]




	def calculateSourceTerm(self,systemConstants, pairsData, particleSet): 
		for particle in particleSet:

			if particle.particleVariables["isBoundary"] is False:

				sum1 = 0
				sum2 = 0

				for neighbor in particle.neighborList:
					pairData = pairsData[(particle,neighbor)]
					relvel = particle.vel - neighbor.vel
					
					if neighbor.particleVariables["isBoundary"] is False:
						mff = neighbor.particleVariables["mass"]
						sum1 += mff * (relvel).dot(self.gW(pairData,systemConstants["interactionlen"]))

					if neighbor.particleVariables["isBoundary"] is True:
						mfb = neighbor.particleVariables["psi"]
						sum2 += mfb * (relvel).dot(self.gW(pairData,systemConstants["interactionlen"]))

				particle.particleVariables["source"] = systemConstants["rho0"] - particle.particleVariables["rho"]
				particle.particleVariables["source"] -= systemConstants["dt"] * sum1 
				particle.particleVariables["source"] -= systemConstants["dt"] * sum2




	def calculateDiagonal(self,systemConstants, pairsData, particleSet): 

		gamma = systemConstants["gamma"]	
		dt = systemConstants["dt"]

		for particle in particleSet:

			if particle.particleVariables["isBoundary"] is False:

				particle.particleVariables["A_ff"] = 0

				sum1 = Vec2(0,0)
				sum2 = Vec2(0,0)

				mf = particle.particleVariables["mass"]
				rhof0 = systemConstants["rho0"]


				# Compute the sums inside the sums.
				# ---------------------------------------------------------
				for neighbor in particle.neighborList:
					pairData = pairsData[(particle,neighbor)]
					gradW = self.gW(pairData,systemConstants["interactionlen"])

					if neighbor.particleVariables["isBoundary"] is False:
						mff   = neighbor.particleVariables["mass"]
						sum1 += mff / (rhof0**2) * gradW

					else:
						mfb   = neighbor.particleVariables["psi"]
						sum2 += mfb / (rhof0**2) * gradW 
				
				sum1 = sum1 * (-1)
				sum2 = sum2 * (- 2 * gamma)					

				# Compute the sums.
				# ---------------------------------------------------------
				for neighbor in particle.neighborList:

					pairData = pairsData[(particle,neighbor)]
					gradW = self.gW(pairData,systemConstants["interactionlen"])

					if neighbor.particleVariables["isBoundary"] is False:
						mff   = neighbor.particleVariables["mass"]
						particle.particleVariables["A_ff"] += mff * (sum1 + sum2).dot(gradW)
						particle.particleVariables["A_ff"] += mff * (mf / (rhof0 ** 2)) * (- gradW).dot(gradW)

					else:
						mfb  = neighbor.particleVariables["psi"]
						particle.particleVariables["A_ff"] += mfb * (sum1 + sum2).dot(gradW)

				particle.particleVariables["A_ff"] = (dt ** 2) * particle.particleVariables["A_ff"] 


	def calculatePressureAccel(self,particle,systemConstants, pairsData):

		sum1 = Vec2(0,0)
		sum2 = Vec2(0,0)
		gamma = systemConstants["gamma"]

		rhof0 = systemConstants["rho0"]
		pf = particle.particleVariables["pressure"]

		# print("pressure : {}, density_est : {}".format(pf,particle.particleVariables["rho"]))
		
		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]			
			gradW = self.gW(pairData,systemConstants["interactionlen"])

			if neighbor.particleVariables["isBoundary"] is False:
				mff = neighbor.particleVariables["mass"]
				pff = neighbor.particleVariables["pressure"]
				sum1 += mff * (pf/(rhof0**2)+pff/(rhof0**2)) * gradW

			else:
				mfb = neighbor.particleVariables["psi"]
				sum2 += 2 * mfb * (pf/(rhof0**2)) * gradW
		return (- sum1 - gamma * sum2)

	def calculateDivAccel(self,particle,systemConstants, pairsData):

		sum1 = 0 

		rhof0 = systemConstants["rho0"]
		pf    = particle.particleVariables["pressure"]
		apf   = particle.particleVariables["a_p"]
		dt    = systemConstants["dt"]

		for neighbor in particle.neighborList:
			pairData = pairsData[(particle,neighbor)]			
			gradW = self.gW(pairData,systemConstants["interactionlen"])

			if neighbor.particleVariables["isBoundary"] is False:

				mff  = neighbor.particleVariables["mass"]
				pff  = neighbor.particleVariables["pressure"]
				apff = neighbor.particleVariables["a_p"]
				sum1 += mff * (apf - apff).dot(gradW)

			else:

				mfb = neighbor.particleVariables["psi"]
				sum1 += mfb * apf.dot(gradW)

		return (dt ** 2) * (sum1)


	def solveForPressure(self,systemConstants, pairsData, particleSet): 
		for particle in particleSet:
			particle.particleVariables["pressure"] = 0

		l=0
		rhoError = 1

		while ( (rhoError > 0.1 and l < 100) or l < 3):

			numFluidParticles = 0
			for particle in particleSet:
				if particle.particleVariables["isBoundary"] is False:
					numFluidParticles += 1
					particle.particleVariables["a_p"] = self.calculatePressureAccel(particle,systemConstants,pairsData)

			rhoError = 0

			for particle in particleSet:

				if particle.particleVariables["isBoundary"] is False:

					pf  = particle.particleVariables["pressure"]
					sf  = particle.particleVariables["source"]
					aff = particle.particleVariables["A_ff"]
					particle.particleVariables["Ap_f"] = self.calculateDivAccel(particle,systemConstants, pairsData)

					if abs(aff) > 1E-9:
						apf = particle.particleVariables["Ap_f"]
						particle.particleVariables["pressure"] = max(pf + self.omega * (sf-apf)/aff,0)
						# print(sf,apf,aff)
						# print(pf + self.omega * (sf-apf)/aff)

					rhoError += (apf - sf)

			rhoError = rhoError / (systemConstants["rho0"] * numFluidParticles) * 100.0

			print("Performing Jacobi iteration for pressure field. Iteration : {} / Average Density Deviation : {}".format(l,rhoError))
			l += 1




	def integration(self,systemConstants,pairsData,particleSet):

		maxVel = 0
		epsilon = 0.2
		for particle in particleSet:
			if particle.particleVariables["isBoundary"] is False:
				particle.vel += particle.particleVariables["a_p"] * systemConstants["dt"]
				# Apply XSPH
				for neighbor in particle.neighborList:
					if neighbor.particleVariables["isBoundary"] is False:
						pairData = pairsData[(particle,neighbor)]
						rhoab = 0.5 * (particle.particleVariables["rho"] + neighbor.particleVariables["rho"])				
						particle.vel += epsilon * (neighbor.particleVariables["mass"] * (- pairData.relvel / (rhoab) ) * self.W(pairData,systemConstants["interactionlen"]))

				if particle.vel.length() > maxVel :
					maxVel = particle.vel.length()

				particle.pos += systemConstants["dt"] * particle.vel


		# print("max Vel : {}".format(maxVel))
		print("Timestep estimate : {}".format(systemConstants["diameter"]/(maxVel+0.001)))



