class Kernel:
	def __init__ (self,smoothingLength):
		self.h = smoothingLength

	def W(pairDat,h):
		q = pairDat.dist / h
		if q >= 0 and q <=1 :
			return 0.34104630662549*((2-q)**3-4*(1-q))/(h**2)
		elif q > 1 and q <= 2 :
			return 0.34104630662549*(2-q)**3/(h**2)
		else :
			return 0

	def gradW(pairDat,h):
		q = pairDat.dist / h
		if q >= 0 and q <= 1 :
			return 0.34104630662549*(4-3*(q-2)**2)/(pairDat.dist*h**2)*pairDat.reldir
		elif q >= 1 and q <= 2 :
			return 0.34104630662549*( -3*(2-q)**2)/(pairDat.dist*h**2)*pairDat.reldir
		else :
			return Vec2(0,0)
