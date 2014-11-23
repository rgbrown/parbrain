from struct import unpack
import os
import sys

def printList(l, rowSize=10):
    length = len(l)
    numRows = length/rowSize

    for i in range(numRows):
        print "%3d: " %i,
        for j in range(rowSize):
            print "%7.4f " %l[i*rowSize+j],
        print ""
    print "%3d: " %numRows,
    for j in range(length-numRows*rowSize):
        print "%7.4f " %l[numRows*rowSize+j],
    print ""


def readTissueBlock(dirName, n_blocks, eqs_per_block, m_local, n_local, m_global, n_global):
	values = []
	fileName = '%s/tissueBlocks.dat' %(dirName)
	bf = open(fileName,'rb')
	while(bf.tell() < os.fstat(bf.fileno()).st_size):
		value = bf.read(8)
		lx=unpack('d',value)
		values.append(lx[0])

	bf.close()

	xCoords=values[:n_blocks]
	yCoords=values[n_blocks:2*n_blocks]
	print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< COORDNIATEs >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

	for i in range(n_blocks):
		print "Block #%.2d  (%7.4f, %7.4f)" %(i, xCoords[i], yCoords[i])
	"""
	n_processors : 2
	n_blocks : 8 (per processor)
	eqs_per_block : 34
	m_local : 2
	n_local : 4
	m_global : 2
	n_global : 1

	"""
	blockSize = eqs_per_block*n_local*m_local*m_global*n_global
	j=0
	for i in range(2*n_blocks,len(values)):
		if j==0: 
			print "======= TimeStamp = %f ========" %values[i]
		else:
			print "%7.4f "%values[i],
		if j<blockSize:
			j+=1
		else:
			j=0
			print ""
		

#	print values[2*n_blocks:]
	print "Len=%d" %len(values[2*n_blocks:])



def readInfo(dirName):
	fileName = '%s/info.dat' %(dirName)

	af = open(fileName, 'r')
	lines=af.readlines()
	header=[x.rstrip('\n') for x in lines[0].split(' ') if x!='' and x!='\n']
	values=[int(x.rstrip('\n')) for x in lines[1].split(' ') if x!='' and x!='\n']
	print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   I N F O   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	n_processors=values[0]
	n_blocks=values[1]*values[0] 
	eqs_per_block=values[2]
	m_local = values[3]
	n_local = values[4]
	m_global=values[5]
	n_global=values[6]
	for i in range(len(header)):
		print "%s : %d" %(header[i],values[i])


	readTissueBlock(dirName, n_blocks,eqs_per_block, m_local,n_local, m_global, n_global)


	

def readFlow(dirName,numLevels):
	values = []	
	fileName = '%s/flow.dat' %(dirName)

	bf = open(fileName,'rb')
	#	af.write("Lon,Lat,Amplitude\n")
	while(bf.tell() < os.fstat(bf.fileno()).st_size):
		value =bf.read(8) #read 8 bytes of data (skip spaces)
		lx=unpack('d',value) #
		values.append(lx)

	bf.close()

	print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   F L O W   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	lenEachTimeStampList = 2**numLevels
	i=0
	while i*lenEachTimeStampList < len(values):
		timeStamp = values[i*lenEachTimeStampList]
		oneTimeStampList = values[i*lenEachTimeStampList+1:(i+1)*lenEachTimeStampList]
		i+=1
		print "======= TimeStamp = %f ========" %timeStamp
		printList(oneTimeStampList)

def readPressure(dirName,numLevels):
	values = []	
	fileName = '%s/pressure.dat' %(dirName)

	bf = open(fileName,'rb')
	while(bf.tell() < os.fstat(bf.fileno()).st_size):
		value =bf.read(8) #read 8 bytes of data (skip spaces)
		lx=unpack('d',value) #
		values.append(lx)

	bf.close()

	print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   PRESSURE  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	lenEachTimeStampList = 2**numLevels/2
	i=0
	while i*lenEachTimeStampList < len(values):
		timeStamp = values[i*lenEachTimeStampList]
		oneTimeStampList = values[i*lenEachTimeStampList+1:(i+1)*lenEachTimeStampList]
		i+=1
		print "======= TimeStamp = %f ========" %timeStamp
		printList(oneTimeStampList)


if __name__ == '__main__':

	if len(sys.argv) == 3:
		dirName=sys.argv[1]
		numLevels = int(sys.argv[2])
	else:
		print "Usage: %s dirName numLevels" %sys.argv[0]
		sys.exit()

	readInfo(dirName)
	readFlow(dirName,numLevels)
	readPressure(dirName,numLevels)


