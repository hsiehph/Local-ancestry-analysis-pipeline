
from operator import itemgetter
import itertools
import pandas as pd
import sys, os


if __name__ == "__main__":
	inputMSPfile = sys.argv[1]
	cutoff_segLength = float(".".join(sys.argv[2].split("_")))

	with open(inputMSPfile) as fin:
		header1 = fin.readline()
		header2 = fin.readline()
		
	dat = pd.read_csv(inputMSPfile, delimiter="\t", skiprows=[0])

	new_df1 = None

	for k, g in itertools.groupby(dat.iterrows(), key=lambda row: row[1][6]):
		list_seg = [ list(elem[1]) for elem in list(g) ]
		seg_genLength = list_seg[len(list_seg)-1][4] - list_seg[0][3]
		if new_df1 is None:
			new_df1 = list_seg
		else:
			if seg_genLength < cutoff_segLength:
				ancCode_last = new_df1[len(new_df1)-1][6]
				for l in list_seg:
					l[6] = ancCode_last

				new_df1.extend(list_seg)
			else:
				new_df1.extend(list_seg)


	new_df2 = None
	for k, g in itertools.groupby(dat.iterrows(), key=lambda row: row[1][7]):
		list_seg = [ list(elem[1]) for elem in list(g) ]
		seg_genLength = list_seg[len(list_seg)-1][4] - list_seg[0][3]
		if new_df2 is None:
			new_df2 = list_seg
		else:
			if seg_genLength < cutoff_segLength:
				ancCode_last = new_df2[len(new_df2)-1][7]
				for l in list_seg:
					l[7] = ancCode_last

				new_df2.extend(list_seg)
			else:
				new_df2.extend(list_seg)

	new_df = []
	for i in range(len(new_df1)):
		l_tmp = new_df1[i][:7]
		l_tmp.append(new_df2[i][7])
		new_df.append(l_tmp)


	with sys.stdout as fout:
		fout.write(header1)
		fout.write(header2)
		for elem in new_df:
			fout.write("\t".join([str(x) for x in elem]) + "\n")

	


