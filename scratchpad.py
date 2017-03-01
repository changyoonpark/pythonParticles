def foo(a,b):
	print("dfafd")

def foobar(a,b,c = 1):
	print(c)

bar = [foo,foobar]

for func in bar:
	func(1,2,3)
