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


if __name__ == '__main__':

	if len(sys.argv) == 3:
		binaryNamePrefix=sys.argv[1]
		numLevels = int(sys.argv[2])
	else:
		print "Usage: %s binaryFileName numLevels" %sys.argv[0]
		sys.exit()
	values = []	
	binaryName = '%s' %(binaryNamePrefix)
	asciiName = '%s.txt' %(binaryName)
	bf = open(binaryName,'rb')
	while(bf.tell() < os.fstat(bf.fileno()).st_size):
		value =bf.read(8) #read 8 bytes of data (skip spaces)
		lx=unpack('d',value) #
		values.append(lx)

	bf.close()

	lenEachTimeStampList = 2**numLevels/2
	i=0
	while i*lenEachTimeStampList < len(values):
		timeStamp = values[i*lenEachTimeStampList]
		oneTimeStampList = values[i*lenEachTimeStampList+1:(i+1)*lenEachTimeStampList]
		i+=1
		print "======= TimeStamp = %f ========" %timeStamp
		printList(oneTimeStampList)

