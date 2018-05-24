#Compiles data from files analyzed by the tailseq analyzer


path = raw_input("Compiled Data: Enter the directory to be compiled (e.g. Data/): ")



from glob import glob
import csv


compiled = []


#Compiles the cumulative 3'end position data
datafilenames = glob(path + '*_Positiondata.csv')
nd = -1

for datafile in datafilenames:

	fcsv = open(datafile)

	nd += 1
	l = 0

	for line in csv.reader(fcsv):

		if nd < 1:
			compiled.append([line[0]])

		if l > 0:
			compiled[l].append(line[12])

		l += 1


	compiled[0].append(datafile.replace("_Positiondata.csv","_3'end"))

	fcsv.close()


#Compiles the cumulative total length data
datafilenames = glob(path + '*_Positiondata.csv')
nd = -1

for datafile in datafilenames:

	fcsv = open(datafile)

	nd += 1
	l = 0

	for line in csv.reader(fcsv):

		if nd < 1:
			compiled[l].append('')
			compiled[l].append(line[0])

		if l > 0:
			compiled[l].append(line[15])

		l += 1


	compiled[0].append(datafile.replace("_Positiondata.csv","_TotalLength"))

	fcsv.close()



#Compiles the cumulative unencoded tail length
datafilenames = glob(path + '*_Tailcomposition.csv')
nd = -1

for datafile in datafilenames:

	fcsv = open(datafile)

	nd += 1
	l = 0

	for line in csv.reader(fcsv):

		if nd < 1:
			compiled[l].append('')
			compiled[l].append(line[0])

		if l > 0:
			compiled[l].append(float(line[1])+float(line[2])+float(line[3])+float(line[4]))

		l += 1


	compiled[0].append(datafile.replace("_Tailcomposition.csv",""))

	fcsv.close()


savefile = path + path.replace("/","") + '_compileddata.csv'
f = open(savefile, 'w')
csv.writer(f).writerows(compiled)
f.close()




	


