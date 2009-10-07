import os, sys

def process_function_arguments(keywords, argument_default_dict, error_doc='', class_to_have_attr=None, howto_deal_with_required_none=1, default_value_in_list=0):
	"""
	2008-04-28
		add argument default_value_in_list to differentiate argument_default_dict whose value is a list which potentially
		could include option instructions (short_option, has_argument, description_for_option).
	2008-04-10
		add argument: howto_deal_with_required_none. if =1, prompt user to enter. else, terminate the program.
	2008-04-09
		if required argument is not given, prompt the user to get it.
	2008-04-02
		add class_to_have_attr to assign the argument and values
	2008-02-28
		example of argument_default_dict:
		
			argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('experimenter',0,):None,\
								('mapping_file',1,):None,\
								('input_dir',1, ):None,\
								('output_dir',1, ):'/Network/Data/250k/db/raw_data/',\
								('array_data_table',0, ):'array_data',\
								('probes_table',0, ):'probes',\
								('strain_info_table',0, ):'strain_info',\
								('array_info_table',1, ):'array_info',\
								('commit',0, int):0,\
								('debug',0, int):0,\
								('report',0, int):0}

		process arguments from keywords based on argument_default_dict
		
		keywords is a dictionary of given arguments. 'argument_name':'argument_value'
		
		argument_default_dict is a dictionary of default arguments.
			the key is a tuple, ('argument_name', is_argument_required, argument_type) and argument_type is optional.
			the value is could just be the default_value or a list = [default_value, 'short_option', 'long_option', has_argument, description_for_option]
				'short_option', has_argument, description_for_option are all optional.
	"""
	import sys
	ad = {}
	argument_key_ls = argument_default_dict.keys()
	argument_key_ls.sort()
	argument_key_ls.reverse()	#to keep 'user' appearing in front of 'password'.
	for argument_key in argument_key_ls:
		if default_value_in_list:
			default_value = argument_default_dict[argument_key][0]
		else:
			default_value = argument_default_dict[argument_key]
		argument, is_argument_required = argument_key[0:2]
		if len(argument_key)>2:
			argument_type = argument_key[2]
		else:
			argument_type = None
		if keywords.has_key(argument):
			if keywords[argument]!='' and keywords[argument]!=None:	#only when keywords has this argument and it's not nothing. change default_value
				if argument_type!=None:	#cast to the desired type
					default_value = argument_type(keywords[argument])
				else:
					default_value = keywords[argument]
		if is_argument_required==1 and (default_value=='' or default_value==None):
			if howto_deal_with_required_none==1:
				if argument=='passwd' or argument=='Passwd' or argument=='password':
					import getpass
					default_value = getpass.getpass("%s: "%argument)
				else:
					default_value = raw_input("%s: "%argument)
			else:
				if error_doc:
					sys.stderr.write(error_doc)
				sys.stderr.write('Error: %s is required but %s given.\n'%(argument, default_value))
				sys.exit(2)
		if argument_type!=None:	#cast to the desired type
			default_value = argument_type(default_value)
		if class_to_have_attr:
			setattr(class_to_have_attr, argument, default_value)
		ad[argument] = default_value
	return ad

def process_options(argv_list, option_default_dict, error_doc=''):
	"""
	2008-04-25
		add option 'help' into option_default_dict
		program exits if '-h' or '--help' is set.
		throw away 'try ... except ...' around getopt
	2008-04-20
		wraps the option handling (getopt), the usual block underneath "if __name__ == '__main__':"
	example of option_default_dict:
	
		option_default_dict = {('z', 'hostname', 1, 'hostname of the db server', 1, ): 'papaya.usc.edu',\
							('d', 'dbname', 1, '', 1, ): 'stock_250k',\
							('u', 'user', 1, '', 1, ):None,\
							('p', 'passwd', 1, '', 1, ):None,\
							('t', 'call_info_table', 1, '', 1, ): 'call_info',\
							('i', 'cmp_data_filename', 1, 'the data file to be compared with.', 1, ): None,\
							('q', 'call_QC_table', 1, '', 1, ): 'call_QC',\
							('m', 'QC_method_id', 1, 'id in table QC_method', 1, int): None,\
							('c', 'commit', 0, 'commit db transaction', 0, int):0,\
							('b', 'debug', 0, 'toggle debug mode', 0, int):0,\
							('r', 'report', 0, 'toggle report, more verbose stdout/stderr.', 0, int):0}
	
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	import sys
	if len(argv_list) == 1:
		print error_doc
		sys.exit(2)
	option_default_dict[('h', 'help', 0, 'Display this documentation', 0, int)] =0
	
	#prepare long_options_list and short_options_str for getopt
	long_options_list = []
	short_options_list = []
	short_option2long_option = {}
	long_option2has_argument = {}
	for option_key in option_default_dict:
		short_option, long_option, has_argument = option_key[0:3]
		if has_argument:
			long_options_list.append('%s='%long_option)
		else:
			long_options_list.append('%s'%long_option)
		if has_argument:
			short_options_list.append('%s:'%short_option)
		else:
			short_options_list.append(short_option)
		short_option2long_option[short_option] = long_option
		long_option2has_argument[long_option] = has_argument
	short_options_str  = ''.join(short_options_list)
	
	#handle options
	import getopt, traceback
	opts, args = getopt.getopt(argv_list[1:], short_options_str, long_options_list)
	opts_dict = {}	#a dictionary 
	for opt, arg in opts:
		if opt[1]=='-':	#it's long option
			long_option = opt[2:]
		else:	#it's short option
			short_option = opt[1:]
			long_option = short_option2long_option[short_option]
		if long_option=='help':
			print error_doc
			sys.exit(2)
		if long_option2has_argument[long_option]:
			opts_dict[long_option] = arg
		else:	#toggle the bit for options which don't have arguments
			opts_dict[long_option] = 1
	"""
	except:
		traceback.print_exc()
		print sys.exc_info()
		print error_doc
		sys.exit(2)
	"""
	return opts_dict

class ProcessOptions(object):
	"""
	2008-04-28
		a class to replace process_options(), process_function_arguments(), generate_program_doc(), turn_option_default_dict2argument_default_dict()
	"""
	getopt_ready = 0
	def __init__(self, argv_list=None, option_default_dict=None, error_doc='', program_name=None):
		self.argv_list = argv_list
		self.option_default_dict = option_default_dict
		self.error_doc = error_doc
		if program_name!=None and program_name!='':
			self.program_name = program_name
		elif (program_name==None or program_name=='') and type(argv_list)==list and len(argv_list)>0:
			self.program_name = argv_list[0]
		else:
			self.program_name = ''
	
	def prepare_for_getopt(self, option_default_dict):
		"""
		2008-09-09
			in the case that short_option is not provided,
				get the has_argument and description_for_option if they are provided.
		2008-05-02
			correct a bug in deal with the rest with no short_option specified
		2008-04-28
			2nd-generation of process_options()
			option_default_dict becomes an enlargement of and compatible to argument_default_dict.
			
			option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server',],\
								('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
								('user', 1, ): [None, 'u', 1, ],\
								('passwd', 1, ): [None, 'p', ],\
								('commit', 0, int): [0, 'c', 0, 'commit db transaction', ],\
								('debug', 0, int): [0, 'b', 0, 'toggle debug mode', ],\
								('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.', ]}
		
			option_default_dict is a dictionary for option handling, including argument_default_dict info
				the key is a tuple, ('argument_name', is_argument_required, argument_type) and argument_type is optional.
				
				the value could just be the default_value or a list = [default_value, 'short_option', has_argument, description_for_option]
				'short_option', has_argument, description_for_option are all orderly optional. You can just supply 'short_option' or 'short_option', has_argument.
		"""
		import sys
		option_default_dict[('help', 0, int)] = [0, 'h', 0, 'Display this documentation', ]
		std_option_default_dict = option_default_dict.copy()	#for generate_program_doc()
		
		#prepare long_options_list and short_options_str for getopt
		long_options_list = []
		short_options_list = []
		short_option2long_option = {}
		long_option2has_argument = {}
		options_with_no_short_option = []
		for option_key, option_value in option_default_dict.iteritems():
			long_option = option_key[0]
			if type(option_value)==list and len(option_value)>=4:
				default_value, short_option, has_argument, description_for_option = option_value[:4]
			elif type(option_value)==list and len(option_value)>=3:
				default_value, short_option, has_argument = option_value[:3]
				description_for_option = ''
			elif type(option_value)==list and len(option_value)==2:
				default_value, short_option = option_value[:2]
				has_argument = 1	#not specified, assuming this option has argument.
				description_for_option = ''
			else:	#04/28/08 deal with these later. use short options that are not used by others
				options_with_no_short_option.append(option_key)
				continue
			
			if short_option==None or short_option=='':
				options_with_no_short_option.append(option_key)
				continue
			#prepare short options
			if short_option in short_option2long_option:
				sys.stderr.write("Error: short option %s already used by %s.\n"%(short_option, short_option2long_option[short_option]))
				sys.exit(3)
			short_option2long_option[short_option] = long_option
			if has_argument:
				short_options_list.append('%s:'%short_option)
			else:
				short_options_list.append(short_option)
			
			#prepare long options
			if long_option in long_option2has_argument:
				sys.stderr.write("Error: long option %s already used.\n"%(long_option))
				sys.exit(3)
			long_option2has_argument[long_option] = has_argument
			if has_argument:
				long_options_list.append('%s='%long_option)
			else:
				long_options_list.append('%s'%long_option)
			
			std_option_default_dict[option_key] = [default_value, short_option, has_argument, description_for_option]
		#deal with the rest with no short_option specified
		#always assume has_argument=1
		chr_ordinal_ls = range(97, 97+26) + range(65, 65+26)	#potential candidates
		eng_letter_ls = map(chr, chr_ordinal_ls)
		eng_letters = ''.join(eng_letter_ls)
		for option_key in options_with_no_short_option:
			long_option = option_key[0]
			option_value = option_default_dict[option_key]
			#2008-09-09 get the has_argument and description_for_option if they are provided.
			if type(option_value)==list and len(option_value)>=4:
				has_argument = option_value[2]
				description_for_option = option_value[3]
			elif type(option_value)==list and len(option_value)>=3:
				has_argument = option_value[2]
				description_for_option = ''
			else:
				has_argument = 1
				description_for_option = ''
			short_option = None
			for letter in long_option + eng_letters:	#start with the long_option, then try all english letters
				if letter not in eng_letter_ls:	#this character has to be english letters. long_option might contain non-letters. like '_'.
					continue
				if letter not in short_option2long_option:
					short_option2long_option[letter] = long_option
					short_options_list.append('%s:'%letter)
					short_option = letter
					break
			if long_option in long_option2has_argument:
				sys.stderr.write("Error: long option %s already used.\n"%(long_option))
				sys.exit(3)
			long_option2has_argument[long_option] = 1
			long_options_list.append('%s='%long_option)
			std_option_default_dict[option_key] = [option_value[0], short_option, has_argument, description_for_option]
			
		self.short_options_str  = ''.join(short_options_list)
		self.short_option2long_option = short_option2long_option
		self.long_option2has_argument = long_option2has_argument
		self.long_options_list = long_options_list
		self.std_option_default_dict = std_option_default_dict
		self.getopt_ready = 1
	
	#@property	not in 2.3
	def program_doc(self):
		"""
		2008-05-02
			fix a bug in toggle need_star
		2008-04-28
			sort the options in short_option ascending order. structure of std_option_default_dict changed.
		2008-04-28
			2nd-generation of generate_program_doc()
			option_default_dict becomes an enlargement of and compatible to argument_default_dict.
		"""
		if self.getopt_ready==0:
			self.prepare_for_getopt(self.option_default_dict)
		
		usage_str = 'Usage: %s [OPTIONS]'%self.program_name
		argument_list_str_ls = ['Argument List:']
		
		#sort the options in short_option ascending order
		option_key_ls_to_be_sorted = []
		for option_key, option_value in self.std_option_default_dict.iteritems():
			if option_key[0] == 'help':	#help would be appended in the end
				continue
			option_key_ls_to_be_sorted.append((option_value[1], option_key))	#option_value[1] = short_option
		option_key_ls_to_be_sorted.sort()
		
		for short_option, option_key in option_key_ls_to_be_sorted:
			
			option_value = self.std_option_default_dict[option_key]
			default_value, short_option, has_argument, description_for_option = option_value
			long_option, is_option_required = option_key[0:2]
			need_star = 0	#required and default_value is empty
			this_argument_ls = ['\t']
			if is_option_required and (default_value==None or default_value==''):
				usage_str += ' -%s %s'%(short_option, long_option.upper())
				need_star = 1
			
			if has_argument:
				this_argument_ls.append('-%s ...,\t--%s=...'%(short_option, long_option))
			else:
				this_argument_ls.append('-%s,\t--%s'%(short_option, long_option))
			
			if need_star:
				this_argument_ls.append('*')
			if description_for_option:
				if description_for_option[-1]=='.':
					this_argument_ls.append('\t%s '%description_for_option)
				else:
					this_argument_ls.append('\t%s. '%description_for_option)
			else:
				this_argument_ls.append('\t')
			if has_argument and default_value!=None and default_value!='':
				this_argument_ls.append('"%s"(default)'%default_value)
			argument_list_str_ls.append(''.join(this_argument_ls))
		this_argument_ls = ['\t', '-h,\t--help', '\tDisplay this documentation']
		argument_list_str_ls.append(''.join(this_argument_ls))
		argument_list_str_ls.append("\nFor required(*) options, if no argument is given, you'll be prompted.")
		program_doc = '\n%s\n\n'%usage_str
		program_doc += '\n'.join(argument_list_str_ls)
		program_doc += '\n'
		return program_doc
	
	program_doc = property(program_doc)
	
	#@property
	def long_option2value(self):
		#handle options
		if self.getopt_ready==0:
			self.prepare_for_getopt(self.option_default_dict)
		
		error_doc = self.program_doc + self.error_doc
		
		import sys
		if len(self.argv_list) == 1:
			print error_doc
			sys.exit(2)
		
		import getopt, traceback
		opts, args = getopt.getopt(self.argv_list[1:], self.short_options_str, self.long_options_list)
		opts_dict = {}	#a dictionary 
		for opt, arg in opts:
			if opt[1]=='-':	#it's long option
				long_option = opt[2:]
			else:	#it's short option
				short_option = opt[1:]
				long_option = self.short_option2long_option[short_option]
			if long_option=='help':
				print error_doc
				sys.exit(2)
			if self.long_option2has_argument[long_option]:
				opts_dict[long_option] = arg
			else:	#toggle the bit for options which don't have arguments
				opts_dict[long_option] = 1
		"""
		except:
			traceback.print_exc()
			print sys.exc_info()
			print error_doc
			sys.exit(2)
		"""
		return opts_dict
	long_option2value = property(long_option2value)

	def process_function_arguments(cls, keywords, argument_default_dict, error_doc='', class_to_have_attr=None, howto_deal_with_required_none=1):
		"""
		2009-10-07
			"is not ''" works for character numpy.array, while "!=''" doesn't in if-condition.
		2008-10-25 
			if default_value is None, and no value is given by a user, no casting into type specified.
		2008-07-09
			not just ='passwd', if the option name (variable 'argument' here) contains 'passwd', use getpass().
		2008-05-06
			return ad as well, just like the standalone function
		2008-04-28
			largely copied from the standalone process_function_arguments()
			argument_default_dict = option_default_dict in prepare_for_getopt()
			The only difference is that value of this argument_default_dict could be non-list, which is plain default_value. To be compatible with previous argument_default_dict.
		"""
		import sys
		argument_key_ls = argument_default_dict.keys()
		argument_key_ls.sort()
		argument_key_ls.reverse()	#to keep 'user' appearing in front of 'password'.
		ad = {}
		for argument_key in argument_key_ls:
			if type(argument_default_dict[argument_key])!=list:
				default_value = argument_default_dict[argument_key]
			else:
				default_value = argument_default_dict[argument_key][0]
			argument, is_argument_required = argument_key[0:2]
			if len(argument_key)>2:
				argument_type = argument_key[2]
			else:
				argument_type = None
			if keywords.has_key(argument):
				if keywords[argument] is not '' and keywords[argument] is not None:	# 2009-10-07 "is not ''" works for character numpy.array, while "!=''" doesn't in if-condition.
					#only when keywords has this argument and it's not nothing. change default_value
					if argument_type!=None:	#cast to the desired type
						default_value = argument_type(keywords[argument])
					else:
						default_value = keywords[argument]
			if is_argument_required==1 and (default_value=='' or default_value==None):
				if howto_deal_with_required_none==1:
					if argument.find('passwd')!=-1 or argument.find('Passwd')!=-1 or argument.find('password')!=-1:
						import getpass
						default_value = getpass.getpass("%s: "%argument)
					else:
						default_value = raw_input("%s: "%argument)
				else:
					if error_doc:
						sys.stderr.write(error_doc)
					sys.stderr.write('Error: %s is required but %s given.\n'%(argument, default_value))
					sys.exit(2)
			if argument_type!=None and default_value is not None:	#cast to the desired type 2008-10-25 default_value is not None
				default_value = argument_type(default_value)
			if class_to_have_attr:
				setattr(class_to_have_attr, argument, default_value)
			ad[argument] = default_value
		return ad
	process_function_arguments = classmethod(process_function_arguments)
	
def turn_option_default_dict2argument_default_dict(option_default_dict):
	"""
	2008-04-20
		option_default_dict contains the info of argument_default_dict.
		to avoid repetitive code
	"""
	argument_default_dict = {}
	for option_key, default_value in option_default_dict.iteritems():
		short_option, long_option, has_argument, description_for_option, is_option_required = option_key[0:5]
		if len(option_key)==6:
			argument_type = option_key[5]
			argument_key = (long_option, is_option_required, argument_type)
		else:
			argument_key = (long_option, is_option_required, )
		argument_default_dict[argument_key] = default_value
	return argument_default_dict

def generate_program_doc(program_name, option_default_dict):
	"""
	2008-04-25
		cosmetic change. if description_for_option already has . in the end. No '.' appended.
	2008-04-25
		add -h, --help to program_doc
	2008-04-24
		generates better doc
	2008-04-21
		automatically generate documentation like "Usage: ... \n Argument List: ... " from option_default_dict
	
	example of option_default_dict:
	
		option_default_dict = {('z', 'hostname', 1, 'hostname of the db server', 1, ): 'papaya.usc.edu',\
							('d', 'dbname', 1, '', 1, ): 'stock_250k',\
							('u', 'user', 1, '', 1, ):None,\
							('p', 'passwd', 1, '', 1, ):None,\
							('t', 'call_info_table', 1, '', 1, ): 'call_info',\
							('i', 'cmp_data_filename', 1, 'the data file to be compared with.', 1, ): None,\
							('q', 'call_QC_table', 1, '', 1, ): 'call_QC',\
							('m', 'QC_method_id', 1, 'id in table QC_method', 1, int): None,\
							('c', 'commit', 0, 'commit db transaction', 0, int):0,\
							('b', 'debug', 0, 'toggle debug mode', 0, int):0,\
							('r', 'report', 0, 'toggle report, more verbose stdout/stderr.', 0, int):0}
	
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	usage_str = 'Usage: %s [OPTIONS]'%program_name
	argument_list_str_ls = ['Argument List:']
	option_key_ls = option_default_dict.keys()
	option_key_ls.sort()
	for option_key in option_key_ls:
		default_value = option_default_dict[option_key]
		short_option, long_option, has_argument, description_for_option, is_option_required = option_key[0:5]
		need_star = 0
		this_argument_ls = ['\t']
		if is_option_required and default_value==None or default_value=='':
			usage_str += ' -%s %s'%(short_option, long_option.upper())
			need_star = 1
		if has_argument:
			this_argument_ls.append('-%s ...,\t--%s=...'%(short_option, long_option))
		else:
			this_argument_ls.append('-%s,\t--%s'%(short_option, long_option))
		if need_star:
			this_argument_ls.append('*')
		if description_for_option:
			if description_for_option[-1]=='.':
				this_argument_ls.append('\t%s '%description_for_option)
			else:
				this_argument_ls.append('\t%s. '%description_for_option)
		else:
			this_argument_ls.append('\t')
		if has_argument and default_value!=None and default_value!='':
			this_argument_ls.append('"%s"(default)'%default_value)
		argument_list_str_ls.append(''.join(this_argument_ls))
	this_argument_ls = ['\t', '-h,\t--help', '\tDisplay this documentation']
	argument_list_str_ls.append(''.join(this_argument_ls))
	argument_list_str_ls.append("\nFor required(*) options, if no argument is given, you'll be prompted.")
	program_doc = '\n%s\n\n'%usage_str
	program_doc += '\n'.join(argument_list_str_ls)
	program_doc += '\n'
	return program_doc
