
	# def calcVelocityLaplacian(self,systemConstants, pairsData, particleSet, plot = None):
	# 	for particle in particleSet:
	# 		particle.particleVariables["vel_div"] = 0;
	# 		particle.particleVariables["vel_lap"] = Vec2(0,0)

	# 		if particle.particleVariables["isBoundary"] is False: 
	# 			for neighbor in particle.neighborList:
	# 				if neighbor.particleVariables["isBoundary"] is False and particle is not neighbor:
	# 					pairData = pairsData[(particle,neighbor)]
	# 					particle.particleVariables["vel_div"] -= (neighbor.particleVariables["mass"] / neighbor.particleVariables["rho"]) * pairData.relvel.dot(self.gW(pairData,systemConstants["interactionlen"]))

	# 					particle.particleVariables["vel_lap"] += \
	# 						( 8.0 ) * \
	# 						( neighbor.particleVariables["mass"] / neighbor.particleVariables["rho"] ) * \
	# 						( pairData.relvel.dot(pairData.relpos) / pairData.relpos.dot(pairData.relpos) ) * \
	# 						( self.gW(pairData,systemConstants["interactionlen"]) )


	# def checkDii(self,systemConstants, pairsData, particleSet, plot = None):
	# 	for particle in particleSet:
	# 		print("IDX : {}, D_II : {}".format(particle.pID,particle.particleVariables["d_ii"]))

	# def calcPressuregrad(self,systemConstants, pairsData, particleSet, plot = None):


	# 	for particle in particleSet:
	# 		if particle.particleVariables["isBoundary"] is False:
	# 			particle.particleVariables["press_grad"] = Vec2(0,0)

	# 			for neighbor in particle.neighborList:
	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					# print(particle.particleVariables["pressure"])
	# 					pairData = pairsData[(particle,neighbor)]
	# 					particle.particleVariables["press_grad"] += \
	# 						(neighbor.particleVariables["mass"] * particle.particleVariables["rho"]) * \
	# 						(particle.particleVariables["pressure"] / (particle.particleVariables["rho"] ** 2) + \
	# 						 neighbor.particleVariables["pressure"] / (neighbor.particleVariables["rho"] ** 2) ) * \
	# 						 self.gW(pairData,systemConstants["interactionlen"])


	# def constructMatrix(self,systemConstants, pairsData, particleSet, plot = None):

	# 	for particle in particleSet:

	# 		if particle.particleVariables["isBoundary"] is False:

	# 			i = particle.pID

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












	# def calculate_dii(self,systemConstants,pairsData,particleSet, plot = None):
	# 	for particle in particleSet:
	# 	#   calculate dii
	# 		if particle.particleVariables["isBoundary"] is False: 

	# 			particle.particleVariables["d_ii"] = Vec2(0,0)

	# 			for neighbor in particle.neighborList:
	# 				pairData = pairsData[(particle,neighbor)]


	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					particle.particleVariables["d_ii"] -= \
	# 						neighbor.particleVariables["mass"] / (particle.particleVariables["rho"]**2) * \
	# 						self.gW(pairData,systemConstants["interactionlen"])
	# 				else:
	# 					particle.particleVariables["d_ii"] -= \
	# 						neighbor.particleVariables["psi"] / (particle.particleVariables["rho"]**2) * \
	# 						self.gW(pairData,systemConstants["interactionlen"])	

	# 			particle.particleVariables["d_ii"] = particle.particleVariables["d_ii"] * systemConstants["dt"]**2

	# def calculate_sum_pj_dij(self,particle,systemConstants,pairsData, plot = None):
	# 	if particle.particleVariables["isBoundary"] is False:

	# 		particle.particleVariables["sum_pj_dij"] = Vec2(0,0)
	# 		for neighbor in particle.neighborList:
	# 			if neighbor.particleVariables["isBoundary"] is False:
	# 				pairData = pairsData[(particle,neighbor)]
	# 				particle.particleVariables["sum_pj_dij"] += self.calculate_dij(pairData,systemConstants) * neighbor.particleVariables["pressure"]


	# def calculate_dij(self,pairData,systemConstants):	
	# 	# dpi = pairData.particlei.particleVariables["mass"] / (pairData.particlej.particleVariables["rho"]**2)

	# 	return (-systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
	# 		   (pairData.particlej.particleVariables["rho"]**2)) * self.gW(pairData,systemConstants["interactionlen"])


	# def calculateAdvectionDensity(self,systemConstants,pairsData,particleSet, plot = None):
	# 	for particle in particleSet:
	# 	#   calculate advection density
	# 		if particle.particleVariables["isBoundary"] is False:
	# 			particle.particleVariables["rho_adv"] = particle.particleVariables["rho"]

	# 			for neighbor in particle.neighborList:

	# 				pairData = pairsData[(particle,neighbor)]
	# 				v_adv_ij = particle.vel - neighbor.vel

	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					particle.particleVariables["rho_adv"] += \
	# 						systemConstants["dt"] * neighbor.particleVariables["mass"] * \
	# 						v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))


	# 				else:
	# 					particle.particleVariables["rho_adv"] += \
	# 						systemConstants["dt"] * neighbor.particleVariables["psi"] * \
	# 						v_adv_ij.dot(self.gW(pairData,systemConstants["interactionlen"]))


	# def calculate_aii(self,systemConstants,pairsData,particleSet, plot = None):
	# 	for particle in particleSet:	
	# 		if particle.particleVariables["isBoundary"] is False:
	# 		# 	calculate a_ii for relaxed jacobi
	# 			particle.particleVariables["a_ii"] = 0
	# 			for neighbor in particle.neighborList:
	# 				pairData  = pairsData[(particle,neighbor)]
	# 				pairData2 = pairsData[(neighbor,particle)] 

	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					particle.particleVariables["a_ii"] += \
	# 					  	 neighbor.particleVariables["mass"] * \
	# 						(particle.particleVariables["d_ii"] - self.calculate_dij(pairData2,systemConstants)) \
	# 						.dot(self.gW(pairData,systemConstants["interactionlen"]))
	# 				else:
	# 					particle.particleVariables["a_ii"] += \
	# 					  	 neighbor.particleVariables["psi"] * \
	# 						(particle.particleVariables["d_ii"] - self.calculate_dij(pairData2,systemConstants)) \
	# 						.dot(self.gW(pairData,systemConstants["interactionlen"]))


	# def calcAverageDensity(self,particleSet):
	# 	rhoAve = 0
	# 	numP = 0
	# 	for particle in particleSet:
	# 		if particle.particleVariables["isBoundary"] is False:
	# 			rhoAve += particle.particleVariables["rho_est"]
	# 			numP += 1

	# 	return rhoAve / numP



	# def estimateDensity(self,particle,pairsData,systemConstants):
	# 	# junk = 0
	# 	if particle.particleVariables["isBoundary"] is False:
	# 		if particle.particleVariables["pressure_est"] > 0:
	# 			particle.particleVariables["rho_est"] = 0

	# 			for neighbor in particle.neighborList:
	# 				# pairData2 = pairsData[(neighbor,particle)]

	# 				# if neighbor.particleVariables["isBoundary"] is False:
	# 				# 	junk +=    neighbor.particleVariables["mass"] * \
	# 				# 			(  particle.particleVariables["sum_pj_dij"] - neighbor.particleVariables["d_ii"] * neighbor.particleVariables["pressure"]\
	# 				# 			 - neighbor.particleVariables["sum_pj_dij"] + self.calculate_dij(pairData2,systemConstants) * particle.particleVariables["pressure"])\
	# 				# 			 .dot(self.gW(pairData,systemConstants["interactionlen"]))
	# 				# else:
	# 				# 	junk +=    neighbor.particleVariables["psi"] * \
	# 				# 			  (particle.particleVariables["sum_pj_dij"]).dot(self.gW(pairData,systemConstants["interactionlen"]))

	# 			# particle.particleVariables["rho_est"] = (particle.particleVariables["a_ii"] * particle.particleVariables["pressure"] + junk)\
	# 													 # + particle.particleVariables["rho_adv"]

	# 				pairData  = pairsData[(particle,neighbor)]
	# 				if neighbor.particleVariables["isBoundary"] is False:
	# 					particle.particleVariables["rho_est"] += neighbor.particleVariables["mass"] * (particle.particleVariables["f_p"] / particle.particleVariables["mass"] - neighbor.particleVariables["f_p"] / neighbor.particleVariables["mass"])\
	# 					.dot(self.gW(pairData,systemConstants["interactionlen"]))
	# 				else:
	# 					particle.particleVariables["rho_est"] += neighbor.particleVariables["psi"] * (particle.particleVariables["f_p"] / particle.particleVariables["mass"])\
	# 					.dot(self.gW(pairData,systemConstants["interactionlen"]))

	# 			particle.particleVariables["rho_est"] = particle.particleVariables["rho_est"] * (systemConstants["dt"] ** 2)

	# 			particle.particleVariables["rho_est"] += particle.particleVariables["rho_adv"]

	# 		else:
	# 			particle.particleVariables["rho_est"] = systemConstants["rho0"]


	# def estimatePressure(self,particle,systemConstants,pairsData,particleSet):
	# 	junk = 0
	# 	if particle.particleVariables["isBoundary"] is False: 
	# 		for neighbor in particle.neighborList:

	# 			pairData  = pairsData[(particle,neighbor)]
	# 			pairData2 = pairsData[(neighbor,particle)]

	# 			if neighbor.particleVariables["isBoundary"] is False:

	# 				dji 	   = self.calculate_dij(pairData2,systemConstants)
	# 				sum_pj_dij = particle.particleVariables["sum_pj_dij"]
	# 				sum_pk_djk = neighbor.particleVariables["sum_pj_dij"]
	# 				pj_djj     = neighbor.particleVariables["pressure"] * neighbor.particleVariables["d_ii"] 
	# 				pi_dji 	   = particle.particleVariables["pressure"] * dji					
	# 				gradKernel = self.gW(pairData,systemConstants["interactionlen"])
	# 				mj         = neighbor.particleVariables["mass"]

	# 				junk += mj * (sum_pj_dij - pj_djj - sum_pk_djk + pi_dji).dot(gradKernel)

	# 			else:
	# 				psi_j 	   = neighbor.particleVariables["psi"]
	# 				sum_pj_dij = particle.particleVariables["sum_pj_dij"]
	# 				gradKernel = self.gW(pairData,systemConstants["interactionlen"])

	# 				junk += psi_j * sum_pj_dij.dot(gradKernel)




	# 		source = systemConstants["rho0"] - particle.particleVariables["rho_adv"];


	# 		# print("source/ junk")
	# 		# print(source,junk)		   										        	
	# 		# print(particle.particleVariables["a_ii"])			
	# 		# pres = 0

	# 		if abs(particle.particleVariables["a_ii"]) > 1E-9 :
	# 			pres = (1-self.omega) * particle.particleVariables["pressure"] + \
	# 		   	    	self.omega * (1 / particle.particleVariables["a_ii"] ) * \
	# 		   										        	(source - junk)
	# 			pres =  max(0,pres)
	# 		else:
	# 			pres = 0

	# 		particle.particleVariables["pressure_est"] = pres



	# def pressureSolver(self,systemConstants,pairsData,particleSet, plot):
	# # 	iteration counter
	# 	l = 0
	# 	#   initialize pressure estimate.
	# 	for particle in particleSet:
	# 		if particle.particleVariables["isBoundary"] is False:
	# 			particle.particleVariables["pressure"] = 0.5 * particle.particleVariables["pressure"]

	# 	densityDeviation = 1
	# 	foo = 0		

	# 	systemConstants["densityDeviation"] = []
	# 	while ( (densityDeviation > 0.005 and l < 100) or l < 3):
	# 		# if (abs(foo - densityDeviation) < 1E-6):
	# 			# break

	# 		for particle in particleSet:
	# 			self.calculate_sum_pj_dij(particle,systemConstants,pairsData)

	# 		for particle in particleSet:
	# 			self.estimatePressure(particle,systemConstants,pairsData,particleSet)	

	# 		for particle in particleSet:
	# 			self.calculatePressureForce(particle,pairsData,systemConstants)

	# 		for particle in particleSet:
	# 			self.estimateDensity(particle,pairsData,systemConstants)

	# 		for particle in particleSet:
	# 			if particle.particleVariables["isBoundary"] is False:
	# 				particle.particleVariables["pressure"] = particle.particleVariables["pressure_est"]


	# 		foo = densityDeviation
	# 		densityDeviation = abs((self.calcAverageDensity(particleSet) - systemConstants["rho0"])/systemConstants["rho0"])
	# 		print("Performing Jacobi iteration for pressure field. Iteration : {} / Average Density Deviation : {}".format(l,densityDeviation * 100))
	# 		systemConstants["densityDeviation"].append([l,densityDeviation * 100])
	# 		l += 1

	# 		# print(plot)			
	# 		# plot.set_offsets(systemConstants["densityDeviation"])
	# 		# plot.draw()
	# 	# for particle in particleSet:
	# 	# 	if particle.particleVariables["isBoundary"] is False and particle.pID == 24:
	# 	# 		print(particle)
	# 	# 		print(" Neighboring boundary particls : {}".format(len(particle.neighborList)))
	# 	# 		for neighbor in particle.neighborList:
	# 	# 			if neighbor.particleVariables["isBoundary"] is True :
	# 	# 				print("boundary particle {}, psi : {}.".format(neighbor.pID,neighbor.particleVariables["psi"]))

	# def tensileInstabCorrection(self, pairData, systemConstants):
	# 	epsilon = 0.2

	# 	if (pairData.particlej.particleVariables["isBoundary"]) :
	# 		pa = pairData.particlei.particleVariables["pressure"];
	# 		pb = pairData.particlei.particleVariables["pressure"];

	# 		rhoa = pairData.particlei.particleVariables["rho"];
	# 		rhob = pairData.particlei.particleVariables["rho"];
	# 	else:			
	# 		pa = pairData.particlei.particleVariables["pressure"];
	# 		pb = pairData.particlej.particleVariables["pressure"];

	# 		rhoa = pairData.particlei.particleVariables["rho"];
	# 		rhob = pairData.particlej.particleVariables["rho"];

	# 	ra = rb = 0

	# 	if pa < 0 :
	# 		ra = epsilon * abs(pa / (rhoa **2 ) )
	# 	else :
	# 		ra = 0.01 * abs(pa / (rhoa **2 ) )

	# 	if pb < 0:			
	# 		rb = epsilon * abs(pb / (rhob **2 ) )
	# 	else :
	# 		rb = 0.01 * abs(pb / (rhob **2 ) )

	# 	r = ra + rb
	# 	n = 4
	# 	W_deltaP = self.W( ParticlePair(None,None,None,None,systemConstants["diameter"],None), systemConstants["interactionlen"] )
	# 	fabn = pow(self.W(pairData,systemConstants["interactionlen"]) / W_deltaP,n)

	# 	return r*fabn


					
	# def calculatePressureForce(self,particle,pairsData,systemConstants):

	# 	# return (-systemConstants["dt"]**2 * pairData.particlej.particleVariables["mass"] / \
	# 	# 	   (pairData.particlej.particleVariables["rho"]**2)) * self.gW(pairData,systemConstants["interactionlen"])



	# 	# if particle.particleVariables["isBoundary"] is False:


	# 	# 	particle.particleVariables["f_p"] = Vec2(0,0)
	# 	# 	for neighbor in particle.neighborList:

	# 	# 		pairData = pairsData[(particle,neighbor)]

	# 	# 		if neighbor.particleVariables["isBoundary"] is False:
	# 	# 			particle.particleVariables["f_p"] += (self.calculate_dij(pairData,systemConstants) * neighbor.particleVariables["pressure"])

	# 	# 	particle.particleVariables["f_p"] += particle.particleVariables["d_ii"] * particle.particleVariables["pressure"]
	# 	# 	particle.particleVariables["f_p"] = particle.particleVariables["f_p"] * particle.particleVariables["mass"] / (systemConstants["dt"])**2

	# 	if particle.particleVariables["isBoundary"] is False:

	# 		particle.particleVariables["f_p"] = Vec2(0,0)
	# 		for neighbor in particle.neighborList:

	# 			pairData = pairsData[(particle,neighbor)]


	# 			# tiCorrection = self.tensileInstabCorrection(pairData,systemConstants)
	# 			tiCorrection = 0

	# 			if neighbor.particleVariables["isBoundary"] is False:
	# 				particle.particleVariables["f_p"] += \
	# 				              -particle.particleVariables["mass"] * neighbor.particleVariables["mass"] * \
	# 				              (particle.particleVariables["pressure"] / ( particle.particleVariables["rho"] ** 2 ) + \
	# 				               neighbor.particleVariables["pressure"] / ( neighbor.particleVariables["rho"] ** 2 ) + tiCorrection ) \
	# 				              *self.gW(pairData,systemConstants["interactionlen"])

	# 			else:					
	# 				particle.particleVariables["f_p"] += \
	# 				              -particle.particleVariables["mass"] * neighbor.particleVariables["psi"] * \
	# 				              (particle.particleVariables["pressure"]/(particle.particleVariables["rho"]**2) + tiCorrection ) * self.gW(pairData,systemConstants["interactionlen"])

	# 		# print("After pressure esimtate : ")
	# 		# print(particle)




