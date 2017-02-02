from kernels import W_4, gW_4, lW_v

class InteractionAlogorithm : 
	def __init__ (systemConstants,pairsData,particleSystem,W,gW,lW):

		self.systemConstants = systemConstants
		self.pairsData = pairsData
		self.particleSystem = particleSystem
		# Kernel functions
		self.W =   W
		self.gW = gW
		self.lW = lW

