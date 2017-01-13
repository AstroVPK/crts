import os as os
import sys as sys

def splitDoc(name, path = None, **kwargs):
	
	if path is None:
		try:
			path = os.environ['CRTSDATADIR']
		except KeyError:
			raise KeyError('Environment variable "CRTSDATADIR" not set! Please set "CRTSDATADIR" to point where all CRTS data should live first...')
	extension = kwargs.get('extension', '.txt')
	fullPath = os.path.join(path, '+Data', name + extension)
	with open(fullPath, 'rb') as fileOpen:
		allLines = fileOpen.readlines()
	allLines = [line.rstrip('\n') for line in allLines]
	inputID = '0'
	for i in xrange(1, len(allLines)-1):
		splitLine = allLines[i].split(',')
		if splitLine[0] != inputID:
			inputID = splitLine[0]
			newFile = os.path.join(os.environ['CRTSDATADIR'], inputID + '.txt')
			file = open(newFile, 'a')
			file.write("MasterID,Mag,Magerr,RA,Dec,MJD,Blend\n")
			for j in range (1, len(splitLine)-1):
				file.write(splitLine[j] + ',')
			file.write(splitLine[-1])
			file.write('\n')
		else:
			for j in range (1, len(splitLine)-1):
				file.write(splitLine[j] + ',')
			file.write(splitLine[-1])
			file.write('\n')