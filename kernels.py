def W_4(pairDat,h):
	q = pairDat.dist / (h)
	if q >= 0 and q <=0.5 :
		return (1.8189136353359467 / (h**2)) * \
		       (1 - 6 * q**2 + 6 * q**3)
	elif q > 0.5 and q <= 1 :
		return (1.8189136353359467 / (h**2)) * \
		       (2 * (1 - q)**3)
	else :
		return 0

def gW_4(pairDat,h):
	q = pairDat.dist / (h)
	
	if q >= 0 and q <= 0.5 :
		return (1.8189136353359467 / (h**3)) * \
			   (- 12 * q + 18 * q**2) * pairDat.reldir
	elif q >= 0.5 and q <= 1 :
		return (1.8189136353359467 / (h**3)) * \
		       (-6 * (1 - q)**2) * pairDat.reldir
	else :
		return Vec2(0,0)

def lW_v(pairDat,h):
	q = pairDat.dist / (h)
	if q >= 0 and q <= 1:
		return (12.732395447351628 / (h**4)) * (1 - q)
	else : 
		return 0
