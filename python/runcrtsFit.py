import sys
import os

for a in os.listdir(os.environ['CRTSDATADIR']): 
	if a.endswith('.txt'):
		a = a[:-4]
		os.system('python crtsFit.py -pMax 1 -maxEvals 10 -nsteps 100 -n %s -z 1' %(a)) #-z [search NED]