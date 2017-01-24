class number:
	def __init__(self,n):
		self.value = n

def addone(number):
	number.value += 1

def forthisset(theset,dothis):
	for number in theset:
		dothis(number)

theset = [number(1),number(2),number(3)]

forthisset(theset,addone)

for el in theset:
	print(el.value)

