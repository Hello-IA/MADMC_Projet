import matplotlib.pyplot as plt
import numpy as np
import pdb



# import matplotlib.pyplot as plt
# import numpy as np

# def readFile(filename,w,v):
# 	f = open(filename, "r")
# 	i=0
# 	for line in f:
# 		if line[0]=="i":
# 			data = line.split()
# 			w[i]=int(data[1])
# 			v[i,0]=int(data[2])
# 			v[i,1]=int(data[3])
# 			i=i+1
# 		else:
# 			if line[0]=="W":
# 				data = line.split()	
# 				W=int(data[1])
# 	f.close()
# 	return W


# def readPoints(filename,p):
# 	f = open(filename, "r")
# 	nbPND = 0
# 	for line in f:
# 		nbPND += 1

# 	YN = np.zeros((nbPND,p))
# 	f = open(filename, "r")
# 	i=0
# 	for line in f:
# 		data = line.split()
# 		for j in range(p):
# 			YN[i][j]=int(data[j])
# 		i=i+1
# 	f.close()
# 	return YN

# w = np.zeros(100, dtype = int)
# v = np.zeros((100,2), dtype = int)

# filename = "data\\100_items\\2KP100-TA-0.dat"

# W = readFile(filename, w, v)
# print(W)

# file_y = "data\\100_items\\2KP100-TA-0.eff"

# Yn = readPoints(file_y, 2)
# print(Yn)



def readFile(filename,w,v):
	f = open(filename, "r")
	i=0
	for line in f:
		if line[0]=="i":
			data = line.split()
			w[i]=int(data[1])
			v[i,0]=int(data[2])
			v[i,1]=int(data[3])
			i=i+1 
		
		else:
			if line[0]=="W":
				data = line.split()	
				W=int(data[1])
		
    	
	f.close()
	return W, w,v



# def readPoints(filename,p):
# 	f = open(filename, "r")
# 	nbPND = 0
# 	for line in f:
# 		nbPND += 1
  
# 	YN = np.zeros((nbPND,p)) 
# 	f = open(filename, "r")
# 	i=0
# 	for line in f:
# 		data = line.split()
# 		for j in range(p):
# 			YN[i][j]=int(data[j])
# 		i=i+1
# 	f.close()
# 	return YN
def readPoints(filename,p):
	f = open(filename, "r")
	nbPND = 0
	for line in f:
			nbPND += 1

	YN = np.zeros((nbPND,p))
	f = open(filename, "r")
	i=0
	for line in f:
		data = line.split()
		for j in range(p):
			YN[i][j]=int(data[j])
		i=i+1
	f.close()
	return YN



"""

w = np.zeros(100,dtype=int)
v = np.zeros((100,2),dtype=int)
filename = "Data/100_items/2KP100-TA-0.dat"
capacity, weights, values = readFile(filename,w,v) # items 
p = 2
file_eff =  "data\\100_items\\2KP100-TA-0.eff"
YN = readPoints(file_eff,p)
print(YN)

# print("First 10 weights (w):", w[:10].reshape(-1,1))
# print("First 10 profits (v):", v[:10])
f = open(filename, "r")
lines = f.readlines()
f.close()

# Find the first 10 lines with item data (starting with 'i')
count = 0
for line in lines:
    if line[0] == "i" and count < 10:
        print(line.strip())  # Print the raw line
        count += 1 
        
"""