def splitDoc(name, path = None, **kwargs):
	import os as os
	import sys as sys
	
	if path is None:
		try:
			path = os.environ['CRTSDATADIR']
		except KeyError:
			raise KeyError('Environment variable "CRTSDATADIR" not set! Please set "CRTSDATADIR" to point where all CRTS data should live first...')
	extension = kwargs.get('extension', '.txt')
	fullPath = os.path.join(path, name + extension)
	with open(fullPath, 'rb') as fileOpen:
		allLines = fileOpen.readlines()
	allLines = [line.rstrip('\n') for line in allLines]
	ID = '0'
	for i in xrange(1, len(allLines)-1):
		splitLine = allLines[i].split(',')
		if splitLine[1] != ID:
			inputID = splitLine[0]
			ID = splitLine[1]
			newFile = os.path.join(os.environ['CRTSDATADIR'], inputID + '.txt')
			file = open(newFile, 'a')
			file.write(allLines[i])
			file.write('\n')
		file.write(allLines[i])
		file.write('\n')