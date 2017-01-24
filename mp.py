from multiprocessing import Process
def foo(*f):
	print(f)

t1 = Process(target = foo,args = ([1,2,3]))
t1.start()
t1.join()
print("End")