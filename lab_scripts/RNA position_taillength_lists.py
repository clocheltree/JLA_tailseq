#Generates datafiles for position, taillength and composition from Taildata files


path = raw_input("tRNA position/taillength: Enter the directory to be analyzed (e.g. Data/):")

from glob import glob
import csv

limit = -200 #sets a limit to the accepted truncation of RNAs to include in the analysis

taildata_filenames = glob(path + '*_Taildata.csv')

for datafile in taildata_filenames:

	fcsv = open(datafile)

	Masterlist = []

	for line in csv.reader(fcsv):
		Masterlist.append(line)

	fcsv.close()


	#Sums up total number of mappable reads from the Masterlist
	totalreads = 0
	totalfulllengthreads = 0
	totaltailedreads = 0
	
	for line in Masterlist[1:]:
		if line[2] != '[]' and int(line[3]) >= limit:
			totalreads += int(line[1])
			if int(line[3]) > -1:
				totalfulllengthreads += int(line[1])
			if int(line[4]) > 0:
				totaltailedreads += int(line[1])



	#Generates the Positionlist and Taillists
	Endpositionlist = [['Position','%Reads','%Tailed (of total)','%Tail>=2','%Tail>=3','%Tail>=4','%Tail>=5',
		'%Tail>=6','%Tail>=7','%Tail>=8','%Tail>=9','%Tail>=10','Cumulative 3-end','Cumulative tailposition', 				'%Tailed (for each position)','Cumulative total length']]

	l = -200
	while l < 51:
		Endpositionlist.append([l, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		l += 1

	Taillengthlist = [['Tail length','%Encoded tails','%Unencoded tails']]
	Tailcomposition = [['Tail Position','A','T','C','G']]

	l = 0
	while l < 26:
		Taillengthlist.append([l, 0, 0])
		Tailcomposition.append([l+1, 0, 0, 0, 0])
		l += 1


	for line in Masterlist[1:]:
		if line[2] != '[]' and int(line[3]) >= limit: #Requires an identified gene

			Endpositionlist[int(line[3])+201][1] += 100*float(line[1])/totalreads

			x = int(line[3])
			while x < 51: #Generates the cumulative 3'end list
				Endpositionlist[x + 201][12] += 100*float(line[1])/totalreads
				x += 1

			x = int(line[3]) + int(line[4])
			while x < 51: #Generates the cumulative total length list
				Endpositionlist[x + 201][15] += 100*float(line[1])/totalreads
				x += 1

			if int(line[4]) < 26:
				Taillengthlist[int(line[4])+1][2] += 100*float(line[1])/totalreads

			if int(line[4]) > 0:
				x = 0
				while x < int(line[4]) and x < 10: #Calculates tail lengths
					Endpositionlist[int(line[3])+201][x+2] += 100*float(line[1])/totalreads
					x += 1

				x = int(line[3])
				while x < 51: #Generates the cumulative tailposition list
					Endpositionlist[x + 201][13] += 100*float(line[1])/totaltailedreads
					x += 1
	
				x = 0
				for n in line[5]: #Calculates average tail composition
					x += 1
					if x < 27:
						if n == 'A':
							Tailcomposition[x][1] += 100*float(line[1])/totalreads
						if n == 'T':
							Tailcomposition[x][2] += 100*float(line[1])/totalreads
						if n == 'C':
							Tailcomposition[x][3] += 100*float(line[1])/totalreads
						if n == 'G':
							Tailcomposition[x][4] += 100*float(line[1])/totalreads
							
			if int(line[3]) > -1:
				if int(line[3]) < 26:
					Taillengthlist[int(line[3])+1][1] += 100*float(line[1])/totalfulllengthreads
				else:
					Taillengthlist[26][1] += 100*float(line[1])/totalfulllengthreads

	
	l = 0
	for line in Endpositionlist[1:]:
		l += 1
		if float(line[1]) >= 0.1:
			Endpositionlist[l][14] += 100*float(line[2])/line[1]

			



	savefile = datafile.replace("Taildata", "Positiondata")
	f = open(savefile, 'w')
	csv.writer(f).writerows(Endpositionlist)
	f.close()

	savefile = datafile.replace("Taildata", "Taillengthdata")
	f = open(savefile, 'w')
	csv.writer(f).writerows(Taillengthlist)
	f.close()

	savefile = datafile.replace("Taildata","Tailcomposition")
	f = open(savefile, 'w')
	csv.writer(f).writerows(Tailcomposition)
	f.close()

	print datafile.split('/')[-1] + ": n= " + str(totalreads)

