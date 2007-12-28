import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

class QualityControl:
	"""
	2007-12-19
		an abstract class for more specific classes to inherit
		the functionality is to do quality control between different types of SNP data
	"""
	def __init__(self, **keywords):
		if keywords.has_key('debug'):
			self.debug = int(keywords['debug'])
	
	def calculate_row_NA_rate(self, strain_acc_list, category_list, data_matrix):
		"""
		2007-12-14
		"""
		sys.stderr.write("Calculating row NA_rate ...")
		ecotypeid_duplicate_NA_rate_ls = []
		for i in range(len(strain_acc_list)):
			ecotypeid = int(strain_acc_list[i])
			duplicate = int(category_list[i])
			no_of_NAs = 0
			row = data_matrix[i]
			no_of_SNPs = len(row)
			for call in row:
				if call==0:
					no_of_NAs += 1
			NA_rate = no_of_NAs/float(no_of_SNPs)
			ecotypeid_duplicate_NA_rate_ls.append([ecotypeid, duplicate, NA_rate])
		sys.stderr.write("Done.\n")
		return ecotypeid_duplicate_NA_rate_ls
	
	def calculate_col_NA_rate(self, col_id_ls, data_matrix):
		"""
		2007-12-18
			data_matrix is a list embedded with list
		"""
		sys.stderr.write("Calculating col NA_rate ...")
		col_id_NA_rate_ls = []
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			no_of_NAs = 0
			no_of_rows = len(data_matrix)
			for j in range(no_of_rows):
				call = data_matrix[j][i]
				if call==0:
					no_of_NAs += 1
			NA_rate = no_of_NAs/float(no_of_rows)
			col_id_NA_rate_ls.append([col_id, NA_rate])
		sys.stderr.write("Done.\n")
		return col_id_NA_rate_ls
	
	def get_col_matching_dstruc(self):
		pass
	
	def get_row_matching_dstruc(self):
		pass
	
	def cmp_row_wise(self, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
		2007-12-18
			strain wise
		2007-12-20
			make it more generic
		2007-12-21 even more generic
		"""
		sys.stderr.write("Comparing row-wise for mismatches ...\n")
		row_id2NA_mismatch_rate = {}
		for row_id1, row_id2 in row_id12row_id2.iteritems():
			row_index1 = row_id2row_index1[row_id1]
			row_index2 = row_id2row_index2[row_id2]
			no_of_mismatches = 0
			no_of_non_NA_pairs = 0
			no_of_NAs = 0
			no_of_totals = 0
			for col_id1, col_index1 in col_id2col_index1.iteritems():
				if col_id1 in col_id12col_id2:
					col_id2 = col_id12col_id2[col_id1]
					col_index2 = col_id2col_index2[col_id2]
					no_of_totals += 1
					if data_matrix1[row_index1][col_index1] == 0:
						no_of_NAs += 1
					if data_matrix1[row_index1][col_index1] != 0 and data_matrix2[row_index2][col_index2]!=0:
						no_of_non_NA_pairs += 1
						if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
							no_of_mismatches += 1
			if no_of_totals >0 and no_of_non_NA_pairs>0:
				NA_rate = no_of_NAs/float(no_of_totals)
				mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
				row_id2NA_mismatch_rate[row_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs]
			else:
				sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
		sys.stderr.write("Done.\n")
		return row_id2NA_mismatch_rate
	
	def _cal_pairwise_dist(self, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
		2007-12-21
		"""
		sys.stderr.write("Calculating pairwise distance ...")
		row_id2pairwise_dist = {}
		for row_id1, row_index1 in row_id2row_index1.iteritems():
			pairwise_dist = []
			for row_id2, row_index2 in row_id2row_index2.iteritems():
				no_of_mismatches = 0
				no_of_non_NA_pairs = 0
				for col_id1, col_id2 in col_id12col_id2.iteritems():
					col_index1 = col_id2col_index1[col_id1]
					col_index2 = col_id2col_index2[col_id2]
					if data_matrix1[row_index1][col_index1]!=0 and data_matrix2[row_index2][col_index2]!=0:
						no_of_non_NA_pairs += 1
						if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
							no_of_mismatches += 1
				if no_of_non_NA_pairs>0:
					mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
					pairwise_dist.append([mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs])
				else:
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
			pairwise_dist.sort()
			row_id2pairwise_dist[row_id1] = pairwise_dist
		sys.stderr.write("Done.\n")
		return row_id2pairwise_dist
	
	def trim_row_id2pairwise_dist(self, row_id2pairwise_dist, min_no_of_non_NA_pairs=10):
		"""
		2007-12-26
			used to throw away unreliable pairwise comparisons
		"""
		new_row_id2pairwise_dist = {}
		for row_id, pairwise_dist_ls in row_id2pairwise_dist.iteritems():
			new_pairwise_dist_ls = []
			for row in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = row
				if no_of_non_NA_pairs>=min_no_of_non_NA_pairs:
					new_pairwise_dist_ls.append(row)
			new_row_id2pairwise_dist[row_id] = new_pairwise_dist_ls
		return new_row_id2pairwise_dist

	def cmp_col_wise(self, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
		2007-12-18
			SNP wise
		2007-12-20
			make it more generic
		"""
		sys.stderr.write("Comparing col-wise for mismatches ...\n")
		col_id2NA_mismatch_rate = {}
		no_of_rows1 = len(data_matrix1)
		for col_id1, col_id2 in col_id12col_id2.iteritems():
			col_index1 = col_id2col_index1[col_id1]
			col_index2 = col_id2col_index2[col_id2]
			no_of_mismatches = 0
			no_of_non_NA_pairs = 0
			no_of_NAs = 0
			no_of_totals = 0
			for row_id1, row_index1 in row_id2row_index1.iteritems():
				if row_id1 in row_id12row_id2:
					row_id2 = row_id12row_id2[row_id1]
					row_index2 = row_id2row_index2[row_id2]
					no_of_totals += 1
					if data_matrix1[row_index1][col_index1] == 0:
						no_of_NAs += 1
					if data_matrix1[row_index1][col_index1] != 0 and data_matrix2[row_index2][col_index2] != 0:
						no_of_non_NA_pairs += 1
						if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
							no_of_mismatches += 1
			if no_of_totals >0 and no_of_non_NA_pairs>0:
				NA_rate = no_of_NAs/float(no_of_totals)
				mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
				col_id2NA_mismatch_rate[col_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs]
			else:
				sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(col_id1, col_id2))
		sys.stderr.write("Done.\n")
		return col_id2NA_mismatch_rate
	
	def load_dstruc(self):
		"""
		2007-12-21
			let children classes to fill in details
			need to load:
				self.data_matrix1, self.data_matrix2
				self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2
				self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2
		"""
		pass

	def get_row_id2info(self, row_id_ls, curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype'):
		"""
		2007-12-16
		2007-12-20
			make it more generic
		"""
		row_id2info = {}
		for row_id in row_id_ls:
			ecotypeid, duplicate = row_id
			curs.execute("SELECT c.comment, e.nativename, e.stockparent FROM %s c, %s e where e.id=c.ecotypeid and e.id=%s and c.duplicate=%s"%(calls_250k_duplicate_comment_table, ecotype_table, ecotypeid, duplicate))
			rows = curs.fetchall()
			if rows:
				comment, nativename, stockparent = rows[0]
				directory = os.path.split(comment)[0]	#take the 1st
				directory = os.path.split(directory)[-1]	#take the last
				row_id2info[row_id] = '%s,%s'%(nativename,directory)
				row_id2info[row_id] = row_id2info[row_id].decode('utf-8', 'ignore')
				#ecotypeid_duplicate2info[key_pair] = ecotypeid_duplicate2info[key_pair].decode('latin10')
			else:
				row_id2info[row_id] = '%s'%row_id
		return row_id2info
	
	def on_click_row(self, event):
		"""
		2007-12-14
		2007-12-16
			need all_ls and ecotypeid_duplicate2info to be global variable
		2007-12-18
			replace all_ls with ecotypeid_duplicate2NA_mismatch_rate
		2007-12-20
			make it more generic
			replace ecotypeid_duplicate2info with row_id2info
			replace ecotypeid_duplicate2NA_mismatch_rate with row_id2NA_mismatch_rate
		"""
		# get the x and y coords, flip y from top to bottom
		import pylab
		x, y = event.x, event.y
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for key, value in self.row_id2NA_mismatch_rate.iteritems():
					NA_rate, mismatch_rate = value[:2]
					if abs(NA_rate-event.xdata)<0.005 and abs(mismatch_rate-event.ydata)<0.005:
						pylab.text(event.xdata, event.ydata, self.row_id2info[key], size=8)
						print "row id: %s, NA_mismatch data: %s, info: %s"%(key, value, self.row_id2info[key])
	
	def on_click_col(self, event):
		"""
		2007-12-18
			need ins.snpid2NA_mismatch_rate
		2007-12-20
			make it more generic
			replace snpid2NA_mismatch_rate with col_id2NA_mismatch_rate
		"""
		# get the x and y coords, flip y from top to bottom
		import pylab
		x, y = event.x, event.y
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for key, value in self.col_id2NA_mismatch_rate.iteritems():
					NA_rate, mismatch_rate = value[:2]
					if abs(NA_rate-event.xdata)<0.005 and abs(mismatch_rate-event.ydata)<0.005:
						pylab.text(event.xdata, event.ydata, key, size=8)
						print "col id: %s, NA_mismatch data: %s"%(key, value)
	
	def plot_NA_mismatch_rate(self, NA_mismatch_rate_ls, on_click_func, title=''):
		"""
		2007-12-14
		"""
		NA_rate_ls = []
		mismatch_rate_ls = []
		for row in NA_mismatch_rate_ls:
			NA_rate, mismatch_rate = row[:2]
			NA_rate_ls.append(NA_rate)
			mismatch_rate_ls.append(mismatch_rate)
		import pylab
		pylab.clf()
		pylab.plot(NA_rate_ls, mismatch_rate_ls, '.')
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min(NA_rate_ls), min(mismatch_rate_ls))-0.1
		diagonal_end = max(max(NA_rate_ls), max(mismatch_rate_ls))+0.1
		pylab.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])
		if title:
			pylab.title(title)
		pylab.xlabel('NA rate')
		pylab.ylabel('mismatch rate')
		pylab.show()
		pylab.connect('button_press_event', on_click_func)
	
	def plot_row_NA_mismatch_rate(self, title=''):
		"""
		2007-12-20
		"""
		self.row_id2NA_mismatch_rate = self.cmp_row_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 )
		self.row_id2info = self.get_row_id2info(self.row_id2NA_mismatch_rate.keys(), self.curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype')
		self.plot_NA_mismatch_rate(self.row_id2NA_mismatch_rate.values(), self.on_click_row, title=title)
	
	def plot_col_NA_mismatch_rate(self, title=''):
		"""
		2007-12-20
		"""
		self.col_id2NA_mismatch_rate = self.cmp_col_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		self.plot_NA_mismatch_rate(self.col_id2NA_mismatch_rate.values(), self.on_click_col, title=title)
	
	def cal_row_id2pairwise_dist(self):
		self.row_id2pairwise_dist = self._cal_pairwise_dist(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
