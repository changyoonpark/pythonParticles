from helpers import *
from grid import *
import scipy as sp
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse import linalg

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
		self.pixelSmoothing = 0.25
		self._proceduresInOrder =  [
									   self.initializeGrid,
									   self.rasterizeGrid,
									   self.detectBoundaries,
									   self.calculateDensity,
			  						   self.calculateAdvectionVel,

			  						   self.solvePPE,
			  						   # self.calculateSourceTerm,

			  						   # self.calculateDiagonal,
			  						   # self.solveForPressure,

			  						   self.integration,

		 						   ]

	def initializeGrid(self,systemConstants, pairsData, particleSet):
		systemConstants["grid"] = Grid(particleSet, self.pixelSmoothing * systemConstants["interactionlen"])


	def rasterizeGrid(self,systemConstants, pairsData, particleSet):
    	# Init nodes
		for key in systemConstants["grid"].nodes:
			node = systemConstants["grid"].nodes[key]
			node.quantities["color"] = 0
			node.quantities["pressure"] = 0
			node.quantities["density"] = 0
			node.quantities["colorGrad"] = Vec2(0,0)
			node.quantities["colorGradIntensity"] = 0
			node.quantities["colorGradX"] = 0
			node.quantities["colorGradY"] = 0
			node.quantities["edgeType"] = None

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





	def traverseParticle(self, particle):		

		if particle.particleVariables["traversed"] is True :
			return

		if particle.particleVariables["edgeType"] == None :
			return


		particle.particleVariables["colorGradAfterGapClosing"] = 1.0
		particle.particleVariables["isFreeSurface"] = True
		particle.particleVariables["traversed"] = True

		for neighbor in particle.neighborList:
			if neighbor.particleVariables["edgeType"] == "Strong" :
				continue
			self.traverseParticle(neighbor)

	def detectBoundaries(self, systemConstants, pairsData, particleSet):
		k1 = 0.35
		k2 = 0.7
		for key in systemConstants["grid"].nodes:
			node = systemConstants["grid"].nodes[key]
			node.quantities["colorGrad"] = systemConstants["grid"].sampleScalarGradientFromGrid(node.nodePos,"color")
			node.quantities["colorGradIntensity"] = node.quantities["colorGrad"].length()

		for key in systemConstants["grid"].nodes:
			node = systemConstants["grid"].nodes[key]
			gradDir = node.quantities["colorGrad"]

			evalPoint1 = node.nodePos + gradDir * self.pixelSmoothing * systemConstants["interactionlen"]
			r = systemConstants["grid"].sampleScalarFromGrid(evalPoint1,"colorGradIntensity")

			evalPoint2 = node.nodePos - gradDir * self.pixelSmoothing * systemConstants["interactionlen"]
			p = systemConstants["grid"].sampleScalarFromGrid(evalPoint2,"colorGradIntensity") 

			if node.quantities["colorGradIntensity"] < r or node.quantities["colorGradIntensity"] < p :
				node.quantities["nonMaxSuppressedColorGradIntensity"] = 0.0
			else:
				node.quantities["nonMaxSuppressedColorGradIntensity"] = node.quantities["colorGradIntensity"]

			# print(node.quantities["nonMaxSuppressedColorGradIntensity"])

			# else:
			# 	if node.quantities["colorGradIntensity"] > k1 :
			# 		node.quantities["nonMaxSuppressedColorGradIntensity"] = node.quantities["colorGradIntensity"]
			# 		if node.quantities["colorGradIntensity"] < k2 :
			# 			node.quantities["edgeType"] = "weak"
			# 			node.quantities["nonMaxSuppressedColorGradIntensity"] = 0.5
			# 		else:
			# 			node.quantities["edgeType"] = "strong"
			# 			node.quantities["nonMaxSuppressedColorGradIntensity"] = 1.0
			# 	else:
			# 		node.quantities["nonMaxSuppressedColorGradIntensity"] = 0.0

		for particle in particleSet:
			particle.particleVariables["traversed"] = False
			particle.particleVariables["isFreeSurface"] = False
			particle.particleVariables["edgeType"] = None

			hashVal = systemConstants["grid"].hashFunction(particle.pos)
			particle.particleVariables["colorGrad"] = systemConstants["grid"].sampleScalarFromGrid(particle.pos,"nonMaxSuppressedColorGradIntensity")
			# print(particle.particleVariables["colorGrad"])
			if particle.particleVariables["colorGrad"] > k1 :
				if particle.particleVariables["colorGrad"] < k2 :
					particle.particleVariables["colorGrad"] = 0.5
					particle.particleVariables["edgeType"] = "Weak"
				else:
					particle.particleVariables["colorGrad"] = 1.0
					particle.particleVariables["edgeType"] = "Strong"
			else:
				particle.particleVariables["colorGrad"] = 0.0

			particle.particleVariables["colorGradAfterGapClosing"] = particle.particleVariables["colorGrad"] 


		for particle in particleSet:
			if particle.particleVariables["edgeType"] == "Strong" : 
				self.traverseParticle(particle)
			



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


	# def calculateSourceTerm(self,systemConstants, pairsData, particleSet): 
	# 	for particle in particleSet:

	# 		if particle.particleVariables["isBoundary"] is False:

	# 			sum1 = 0
	# 			sum2 = 0

	# 			for neighbor in particle.neighborList:
	# 				pairData = pairsData[(particle,neighbor)]
	# 				relvel = particle.vel - neighbor.vel
					
	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					mff = neighbor.particleVariables["mass"]
	# 					sum1 += mff * (relvel).dot(self.gW(pairData,systemConstants["interactionlen"]))

	# 				if neighbor.particleVariables["isBoundary"] is True:
	# 					mfb = neighbor.particleVariables["psi"]
	# 					sum2 += mfb * (relvel).dot(self.gW(pairData,systemConstants["interactionlen"]))

	# 			particle.particleVariables["source"] = systemConstants["rho0"] - particle.particleVariables["rho"]
	# 			particle.particleVariables["source"] -= systemConstants["dt"] * sum1 
	# 			particle.particleVariables["source"] -= systemConstants["dt"] * sum2


	def calculateOmega(self, particle, systemConstants, pairsData):		
		sum1 = Vec2(0,0)
		sum2 = Vec2(0,0)
		for neighbor in particle.neighborList:
			pairIJ = pairsData[(particle,neighbor)]
			if neighbor.particleVariables["isBoundary"] is False:
				sum1 += neighbor.particleVariables["mass"] * self.gW(pairIJ, systemConstants["interactionlen"])
			else:
				sum2 += neighbor.particleVariables["psi"] * self.gW(pairIJ, systemConstants["interactionlen"])

		return (sum1,sum2)

	def solvePPE(self, systemConstants, pairsData, particleSet):

		systemConstants["bVector"] = np.zeros(len(particleSet),dtype = float)
		systemConstants["pVector"] = np.zeros(len(particleSet),dtype = float)

		row = []
		col = []
		entry = []

		# Pre calculate Omega
		for particle in particleSet:
			(particle.particleVariables["omegaF"], particle.particleVariables["omegaB"]) = self.calculateOmega(particle, systemConstants, pairsData)

		# Set the RHS, And the corresponding Entries For the A Matrix
		for particle in particleSet:
			i = particle.pID

			if particle.particleVariables["isBoundary"] is False: 

				if particle.particleVariables["isFreeSurface"] is False:

					sum1 = 0
					sum2 = 0

					for neighbor in particle.neighborList:
						pairIJ = pairsData[(particle,neighbor)]
						relvel = particle.vel - neighbor.vel
						# if neighbor.particleVariables["isFreeSurface"] is False:
						if neighbor.particleVariables["isBoundary"] is False:
							sum1 += neighbor.particleVariables["mass"] * relvel.dot(self.gW(pairIJ,systemConstants["interactionlen"]))
						else:
							sum2 += neighbor.particleVariables["psi"] * relvel.dot(self.gW(pairIJ,systemConstants["interactionlen"]))
						# else:

					# print(len(particleSet), len())
					systemConstants["bVector"][i] = (1./(systemConstants["dt"] ** 2)) * (systemConstants["rho0"] - particle.particleVariables["rho"] - systemConstants["dt"] * (sum1 + sum2))

					# DIAGONAL TERM
					diagonal = 0

					# Term 1.1
					diagonal += ( - particle.particleVariables["omegaF"].dot(particle.particleVariables["omegaF"])/(systemConstants["rho0"]**2))
					# Term 1.3
					diagonal += ( - 2 * systemConstants["gamma"] * particle.particleVariables["omegaF"].dot(particle.particleVariables["omegaB"])/(systemConstants["rho0"]**2))
					# Term 3.1
					diagonal += ( particle.particleVariables["omegaB"].dot(particle.particleVariables["omegaF"])/(systemConstants["rho0"]**2))
					# Term 3.3
					diagonal += ( 2 * systemConstants["gamma"] * particle.particleVariables["omegaB"].dot(particle.particleVariables["omegaB"])/(systemConstants["rho0"]**2))

					row.append(i)
					col.append(i)
					entry.append(diagonal)

					# NEIGHBOR PARTICLE TERMS
					for neighbor in particle.neighborList:

						pairIJ = pairsData[(particle,neighbor)]
						j = neighbor.pID 
						if neighbor.particleVariables["isBoundary"] is False:							
							if neighbor.particleVariables["isFreeSurface"] is False:

								value = 0
								# Term 1.2
								value += ( particle.particleVariables["omegaF"].dot( - neighbor.particleVariables["mass"] * self.gW(pairIJ,systemConstants["interactionlen"] )) / (systemConstants["rho0"] ** 2))
								# Term 3.2
								value += ( particle.particleVariables["omegaB"].dot(   neighbor.particleVariables["mass"] * self.gW(pairIJ,systemConstants["interactionlen"] )) / (systemConstants["rho0"] ** 2))
								# Term 2.1
								value += neighbor.particleVariables["mass"] * self.gW(pairIJ,systemConstants["interactionlen"]).dot(neighbor.particleVariables["omegaF"]) / (systemConstants["rho0"] ** 2)
								# Term 2.3
								value += 2 * systemConstants["gamma"] * neighbor.particleVariables["mass"] * self.gW(pairIJ,systemConstants["interactionlen"]).dot(neighbor.particleVariables["omegaB"]) / (systemConstants["rho0"] ** 2)

								row.append(i)
								col.append(j)
								entry.append(value)

						# NEIGHBOR OF NEIGHBOR TERMS
								for nneighbor in neighbor.neighborList :
									pairJK = pairsData[(neighbor,nneighbor)]
									k = nneighbor.pID
									if nneighbor.particleVariables["isBoundary"] is False:
										if nneighbor.particleVariables["isFreeSurface"] is False:
											# Term 2.3
											value = neighbor.particleVariables["mass"] * self.gW(pairIJ, systemConstants["interactionlen"]).dot( nneighbor.particleVariables["mass"] * self.gW(pairJK, systemConstants["interactionlen"]) ) / (systemConstants["rho0"] ** 2)
											row.append(i)
											col.append(k)
											entry.append(value)						
				else:
					row.append(i)
					col.append(i)
					entry.append(1)
					systemConstants["bVector"][i] = 0

			else:
				row.append(i)
				col.append(i)
				entry.append(1)
				systemConstants["bVector"][i] = 0


		systemConstants["AMatrix"] = csr_matrix((entry, (row, col)))

		# for i in range(0,len(particleSet)):
			# print(systemConstants["AMatrix"][i,i])

		fig, ax = plt.subplots(figsize=(4, 4),nrows = 1,ncols = 1)		
		ax.spy(systemConstants["AMatrix"], markersize = 1)
		ax.figure.show()

		# systemConstants["pVector"] = linalg.spsolve(systemConstants["AMatrix"], systemConstants["bVector"])
		systemConstants["pVector"] = linalg.bicg(systemConstants["AMatrix"], systemConstants["bVector"])[0]
		print(len(systemConstants["pVector"]), len(systemConstants["bVector"]))
		for particle in particleSet:
			particle.particleVariables["pressure"] = systemConstants["pVector"][particle.pID]

		for particle in particleSet:
			# particle.particleVariables["a_p"] = self.calculatePressureAccel(particle,systemConstants,pairsData)
			particle.particleVariables["a_p"] = Vec2(0,0)

		# 	# Enforce pressure dirichlet condition on particles with no neighbors, or particles on the free surface 
		# 	if particle.particleVariables["isFreeSurface"] or (len(particle.neighborList) == 1):
		# 		pID = particle.pID
		# 		print("Free surface : {}".format(pID))
		# 		systemConstants["AMatrix"][:,pID] = 0
		# 		systemConstants["AMatrix"][pID,:] = 0
		# 		systemConstants["AMatrix"][pID,pID] = 1
		# 		systemConstants["bVector"][pID] = 0


		# print("Condition Number : {}".format(np.linalg.cond(systemConstants["AMatrix"].todense())))
		# print("foo")
		# print("eigs : {}".format(linalg.eigs(systemConstants["AMatrix"],k=3)))
		



	# def constructMatrix(self,systemConstants, pairsData, particleSet):

	# 	rhs = []

	# 	row = []
	# 	col = []
	# 	data = []

	# 	for particle in particleSet:

	# 		if particle.particleVariables["isBoundary"] is False:

	# 			i = particle.pID

	# 			# Set the RHS.
	# 			systemConstants["bVector"][i] += systemConstants["rho0"] - particle.particleVariables["rho_adv"]
	# 			omega = Vec2(0,0)

	# 			for neighbor in particle.neighborList:					
	# 				if neighbor.particleVariables["isBoundary"] is False:				
	# 					pairIJ = pairsData[(particle,neighbor)]
	# 					omega += neighbor.particleVariables["mass"] * self.gW(pairIJ,systemConstants["interactionlen"])

	# 			# Diagonal
	# 			systemConstants["aMatrix"][i,i] += particle.particleVariables["d_ii"].dot(omega)

	# 			# "j" Columns
	# 			for neighbor in particle.neighborList:
	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					j = neighbor.pID
	# 					pairIJ = pairsData[(particle,neighbor)]

	# 					systemConstants["aMatrix"][i,j] -= (neighbor.particleVariables["d_ii"] * neighbor.particleVariables["mass"])\
	# 	            						.dot(self.gW(pairIJ,systemConstants["interactionlen"]))

	# 			# "k" Columns
	# 			for neighbor in particle.neighborList:
	# 				if neighbor.particleVariables["isBoundary"] is False:										
	# 					pairIJ = pairsData[(particle,neighbor)]

	# 					for nneighbor in neighbor.neighborList:					
	# 						if nneighbor.particleVariables["isBoundary"] is False:
	# 							k = nneighbor.pID
	# 							pairJK = pairsData[(neighbor,nneighbor)]

	# 							systemConstants["aMatrix"][i,k] -= (self.calculate_dij(pairJK,systemConstants) * neighbor.particleVariables["mass"])\
	# 			            						  .dot(self.gW(pairIJ,systemConstants["interactionlen"]))
				
	# 			# "l" Columns
	# 			for neighbor in particle.neighborList:
	# 				if neighbor.particleVariables["isBoundary"] is False:								
	# 					l = neighbor.pID
	# 					pairIL = pairsData[(particle,neighbor)]

	# 					systemConstants["aMatrix"][i,l] += self.calculate_dij(pairIL,systemConstants).dot(omega)







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



