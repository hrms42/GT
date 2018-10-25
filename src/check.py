import sys

f = open(sys.argv[1])

s = f.readline()
s = f.readline()
s = s.strip().split()
numv = int(s[0])
print(numv)

numf = int(s[1])
print(numf)

for i in range(numv):
	f.readline()

for i in range(numf):
	s = f.readline()
	s = s.strip().split()
	nv = int(s[0])
	if nv==4:
		print("Quad found")


