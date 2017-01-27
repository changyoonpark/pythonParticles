def W_4(pairDat,h):
	q = pairDat.dist / (2*h)
	if q >= 0 and q <=0.5 :
		
	elif q > 0.5 and q <= 1 :

	else :
		return 0

def gradW_4(pairDat,h):
	q = pairDat.dist / (2*h)
	if q >= 0 and q <= 0.5 :
		return * pairDat.reldir
	elif q >= 0.5 and q <= 1 :
		return * pairDat.reldir
	else :
		return Vec2(0,0)
