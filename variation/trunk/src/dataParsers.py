""" 
This library offers functions to parse different types of SNPs data from multiple formats into a snpsData objects.

Bjarni Vilhjalmsson, bvilhjal@usc.edu
"""
import time 
import MySQLdb

from snpsdata import *

# this should be fixed
homedir = "/Users/bjarni/"
resultsdir = ""




accessionName2Id = {'Ms-0': 92, 'Mt-0': 79, 'HR-10': 34, 'Knox-18': 4, 'Spr1-6': 20, 'Knox-10': 3, 'Spr1-2': 19, 'Nok-3': 80, 'Wa-1': 81, 'RRS-10': 2, 'C24': 57, 'Wt-5': 74, 'Ra-0': 69, 'Gu-0': 54, 'Mz-0': 73, 'Tsu-1': 78, 'Fei-0': 82, 'Bur-0': 93, 'Omo2-3': 22, 'Pu2-23': 30, 'Rmx-A180': 6, 'Kondara': 88, 'Tamm-2': 41, 'CS22491': 58, 'Zdr-6': 26, 'Ren-11': 48, 'Zdr-1': 25, 'Ren-1': 47, 'Ler-1': 55, 'Fab-2': 13, 'Yo-0': 61, 'Wei-0': 59, 'Got-7': 45, 'Tamm-27': 42, 'Kas-2': 75, 'Ts-5': 85, 'Ts-1': 84, 'Pu2-7': 29, 'Mr-0': 77, 'Ei-2': 53, 'Mrk-0': 72, 'Lz-0': 52, 'Bil-7': 16, 'Bil-5': 15, 'Sq-8': 38, 'Fab-4': 14, 'Sq-1': 37, 'Omo2-1': 21, 'Var2-1': 17, 'Var2-6': 18, 'Shahdara': 89, 'Uod-7': 50, 'Uod-1': 49, 'Lov-5': 12, 'Lov-1': 11, 'Gy-0': 68, 'Col-0': 62, 'Kin-0': 91, 'NFA-8': 35, 'Nd-1': 56, 'Got-22': 46, 'Br-0': 65, 'HR-5': 33, 'Ull2-3': 24, 'Ull2-5': 23, 'Est-1': 66, 'CIBC-17': 40, 'Ct-1': 76, 'Cvi-0': 51, 'Oy-0': 95, 'LL-0': 87, 'Bor-4': 28, 'Bor-1': 27, 'Pna-10': 8, 'Pna-17': 7, 'Ga-0': 71, 'Bay-0': 70, 'Eden-2': 10, 'Eden-1': 9, 'Pro-0': 86, 'Kz-1': 43, 'RRS-7': 1, 'Kz-9': 44, 'Edi-0': 94, 'An-1': 63, 'CIBC-5': 39, 'Ws-0': 60, 'Ws-2': 96, 'Van-0': 64, 'Rmx-A02': 5, 'Se-0': 83, 'Lp2-2': 31, 'Lp2-6': 32, 'NFA-10': 36, 'Ag-0': 67, 'Sorbo': 90}

accessionId2Name = {1: 'RRS-7', 2: 'RRS-10', 3: 'Knox-10', 4: 'Knox-18', 5: 'Rmx-A02', 6: 'Rmx-A180', 7: 'Pna-17', 8: 'Pna-10', 9: 'Eden-1', 10: 'Eden-2', 11: 'Lov-1', 12: 'Lov-5', 13: 'Fab-2', 14: 'Fab-4', 15: 'Bil-5', 16: 'Bil-7', 17: 'Var2-1', 18: 'Var2-6', 19: 'Spr1-2', 20: 'Spr1-6', 21: 'Omo2-1', 22: 'Omo2-3', 23: 'Ull2-5', 24: 'Ull2-3', 25: 'Zdr-1', 26: 'Zdr-6', 27: 'Bor-1', 28: 'Bor-4', 29: 'Pu2-7', 30: 'Pu2-23', 31: 'Lp2-2', 32: 'Lp2-6', 33: 'HR-5', 34: 'HR-10', 35: 'NFA-8', 36: 'NFA-10', 37: 'Sq-1', 38: 'Sq-8', 39: 'CIBC-5', 40: 'CIBC-17', 41: 'Tamm-2', 42: 'Tamm-27', 43: 'Kz-1', 44: 'Kz-9', 45: 'Got-7', 46: 'Got-22', 47: 'Ren-1', 48: 'Ren-11', 49: 'Uod-1', 50: 'Uod-7', 51: 'Cvi-0', 52: 'Lz-0', 53: 'Ei-2', 54: 'Gu-0', 55: 'Ler-1', 56: 'Nd-1', 57: 'C24', 58: 'CS22491', 59: 'Wei-0', 60: 'Ws-0', 61: 'Yo-0', 62: 'Col-0', 63: 'An-1', 64: 'Van-0', 65: 'Br-0', 66: 'Est-1', 67: 'Ag-0', 68: 'Gy-0', 69: 'Ra-0', 70: 'Bay-0', 71: 'Ga-0', 72: 'Mrk-0', 73: 'Mz-0', 74: 'Wt-5', 75: 'Kas-2', 76: 'Ct-1', 77: 'Mr-0', 78: 'Tsu-1', 79: 'Mt-0', 80: 'Nok-3', 81: 'Wa-1', 82: 'Fei-0', 83: 'Se-0', 84: 'Ts-1', 85: 'Ts-5', 86: 'Pro-0', 87: 'LL-0', 88: 'Kondara', 89: 'Shahdara', 90: 'Sorbo', 91: 'Kin-0', 92: 'Ms-0', 93: 'Bur-0', 94: 'Edi-0', 95: 'Oy-0', 96: 'Ws-2'}


def get2010DataFromDb(host="papaya.usc.edu",chromosomes=[1,2,3,4,5], db = "at", dataVersion="3", user = None,passwd = None):
    """
    Retrieve 2010 data from DB.  Returns a list of RawSnpsData objects. 

    """
    rt = time.time()
    decoder = RawDecoder()  #Other unused informative letters are ['R','Y','S','M','K','W']:
    
    print "Connecting to db."
    if not user:
        import sys
        print "Username:"
        user = sys.stdin.readline().rstrip()
    if not passwd:
        import getpass
        passwd = getpass.getpass()
    try:
        conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = db)
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)
    cursor = conn.cursor ()
        #Get distinct accessions and their id.
        #Generate an internal dictionary using their id.
    numRows = int(cursor.execute("select distinct g.accession, acc.name from at.genotype g, at.allele al, at.accession acc, at.locus l, at.alignment an where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version>="+dataVersion+" and l.offset=0 and g.accession=acc.id and acc.id<97 order by acc.name"))

    dict = {}
    accessions = []
    i = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        dict[int(row[0])]=i
        accessions.append(int(row[0]))
        i = i+1


    print "Fetching 2010 data:"
    snpsds=[]
    for chromosome in chromosomes:
        print "    Chromosome",chromosome
        numRows = int(cursor.execute("select distinct l.position, g.accession, acc.name, al.base from at.genotype g, at.allele al, at.accession acc, at.locus l, at.alignment an where g.allele=al.id and l.id=al.locus and l.alignment=an.id  and an.version=3 and l.offset=0 and g.accession=acc.id and acc.id<97 and l.chromosome="+str(chromosome)+" order by l.position, acc.name"))
        print "    ",numRows,"rows retrieved."
        positions = []
        snps = []
        if numRows > 0:
            row = cursor.fetchone()
            newPosition = int(row[0])
            while(1):
                if not row:
                    break;
                positions.append(newPosition)
                oldPosition = newPosition
                snp = ['NA']*96  #Initialize to missing data.
                while(oldPosition==newPosition):
                    snp[dict[int(row[1])]]=decoder[row[3]]
                    row = cursor.fetchone()
                    if not row:
                        break;
                    newPosition = int(row[0])
                snps.append(snp)

        snpsd = RawSnpsData(snps,positions,accessions=accessions)
        snpsds.append(snpsd)
                
    cursor.close ()
    conn.close ()
    dif = int(time.time() - rt)
    print "It took "+str(dif/60)+" min. and "+str(dif%60)+" sec. to fetch data."

    return snpsds



def getPerlgenDataFromDb(host="papaya.usc.edu", db = "chip", chromosomes=[1,2,3,4,5], user = None,passwd = None):
    """
    Retrieve Perlgen data from DB.  Returns a list of RawSnpsData objects. 

    """

    rt = time.time()
    perlgenTo2010 = {'bay-0':'Bay-0','bor-4':'Bor-4','br-0':'Br-0','bur-0':'Bur-0',
                     'c24':'C24','col-0':'Col-0','cvi-0':'Cvi-0','est-1':'Est-1',
                     'fei-0':'Fei-0','got-7':'Got-7','ler-1':'Ler-1','lov-5':'Lov-5',
                     'nfa-8':'NFA-8','rrs-10':'RRS-10','rrs-7':'RRS-7','sha':'Shahdara',
                     'tamm-2':'Tamm-2','ts-1':'Ts-1','tsu-1':'Tsu-1','van-0':'Van-0'}
    
    perlgenAccessionToId = {}  #Get accession id's (or ecotype id?)  from DB
    

    decoder = {'N':'NA','A':'A','C':'C','G':'G','T':'T'} 
    print "Connecting to db."
    if not user:
        import sys
        print "Username:"
        user = sys.stdin.readline().rstrip()
    if not passwd:
        import getpass
        passwd = getpass.getpass()
    try:
        conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = db)
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)

    cursor = conn.cursor()
    numRows = int(cursor.execute("select distinct d.ecotype from snp_combined_may_9_06_no_van d order by d.ecotype"))
    dict = {}
    accessions = []
    i = 0
    while(1):
        row = cursor.fetchone()
        if not row:
            break;
        dict[row[0]]=i
        accessions.append(accessionName2Id[perlgenTo2010[row[0]]])
        i = i+1

    print "Fetching Perlgen data:"
    snpsds=[]
    for chromosome in chromosomes:
        print "    Chromosome",chromosome
        numRows = int(cursor.execute("select distinct d.position, d.ecotype, d.mbml98 from snp_combined_may_9_06_no_van d where d.chromosome="+str(chromosome)+" and d.mbml98 is not null order by d.position, d.ecotype;"))
        print "    ",numRows,"rows retrieved."
        positions = []
        snps = []
        if numRows > 0:
            row = cursor.fetchone()
            newPosition = int(row[0])
            while(1):
                if not row:
                    break;
                positions.append(newPosition)
                oldPosition = newPosition
                snp = ['NA']*20  #Initialize to missing data.
                while(oldPosition==newPosition):
                    snp[dict[row[1]]]=decoder[row[2]]
                    row = cursor.fetchone()
                    if not row:
                        break;
                    newPosition = int(row[0])
                snps.append(snp)
        snpsd = RawSnpsData(snps,positions,accessions=accessions)
        snpsds.append(snpsd)

    cursor.close ()
    conn.close ()
    dif = int(time.time() - rt)
    print "It took "+str(dif/60)+" min. and "+str(dif%60)+" sec. to fetch data."
    return snpsds

def parse2010Data(datafile=None):
    """
    Returns 2010 Data.
    
    Loads it from a file.
    """

    #Imputed data
    datadir = homedir+"Projects/data/2010-Zhao_et_al/"
    if not datafile:
        datafile = datadir+"SNPs.csv"
    
    #Initialization
    accessions = []
    positions = [] #list[chr][position_index]
    genotypes = [] #list[chr][position_index][acces]
    for i in range(0,5):
        positions.append([])
    
    
    chromosomes = []
    f = open(datafile, 'r')
    lines = f.readlines()
    f.close()
    for chr in (lines[0].rstrip()).split(",")[1:]:
        chromosomes.append(int(chr)-1)
    
    line = (lines[1].rstrip()).split(",")[1:]
    for i in range(0,len(chromosomes)):
        positions[chromosomes[i]].append(int(line[i]))

    for j in range(0,5):
        l = []
        for k in range(0,len(positions[j])):
            l.append([])
        genotypes.append(l)
            

    accessions = []
    for i in range(2,len(lines)):
        line = (lines[i].rstrip()).split(",")
        accessions.append(line[0])
        line = line[1:]
        for j in range(0,5):
            for k in range(0,len(positions[j])):
                genotypes[j][k].append(line[i])

    print accessions

    import random
    #Converting genotype and filtering.
    countAll = 0
    countGood = 0 
    countStupid = 0
    decoder = {'.':'NA'}
    newgenotypes = [[],[],[],[],[]]
    newpositions = [[],[],[],[],[]]
    for i in range(0,len(genotypes)):
        for j in range(0,len(positions[i])):
            countAll = countAll+1
            k = 0
            ntl = [] #list of observed nucleotides.
            for nt in ['0','1','2','3']:
                if nt in genotypes[i][j]:
                    decoder[nt]=k
                    ntl.append(nt)
                    k = k+1
            if k==2:
                countGood = countGood + 1
                l = []
                for nt in genotypes[i][j]:
                    l.append(decoder[nt])
                newgenotypes[i].append(l)
                newpositions[i].append(positions[i][j])
            else:
                if k==1:
                    countStupid = countStupid+1
    print countAll," SNPs in all"
    print countGood," SNPs used"
    print countStupid," Stupid SNPs thrown away"

    for i in range(0,5):
        print newpositions[i][-1]
    del positions
    del genotypes
    positions = newpositions
    genotypes = newgenotypes


    chromasomes = []
    for i in range(0,5):
        chromasomes.append(SnpsData(genotypes[i],positions[i],accessions=accessions))
            
    #print positions[4]
    return(chromasomes)

def parse250DataRaw(imputed = True):
    """
    Returns 250K Data in a list of RawSnpsData objects (one for each chromosome).

    Set imputed to False, if un-imputed data is preferred.

    """

    accessions250To2010 = { "Ag0":"Ag-0", "An1":"An-1", "Bay0B":"Bay-0", "Bay0":"Bay-0", "Bil5":"Bil-5", "Bil7":"Bil-7", "Bor1":"Bor-1", "Bor4":"Bor-4", "Bor4B":"Bor-4", "Br0":"Br-0", "BroA":"Br-0", "BurOB":"Bur-0", "Bur0":"Bur-0", "BUR0A":"Bur-0", "C24A":"C24", "C24B":"C24", "CIBC17":"CIBC-17", "Cibc17":"CIBC-17", "CIBC5":"CIBC-5", "Cibc5":"CIBC-5", "CS2249":"CS22491", "Col0":"Col-0", "LY1":"Col-0", "Ct1":"Ct-1", "Cvi0":"Cvi-0", "CviOA":"Cvi-0", "Eden1":"Eden-1", "Eden2":"Eden-2", "Edi0":"Edi-0", "El2":"Ei-2", "Ei2":"Ei-2", "Est1B":"Est-1", "Est1":"Est-1", "EST1A":"Est-1", "Fab2":"Fab-2", "Fab4":"Fab-4", "Fei0":"Fei-0", "FeiOB":"Fei-0", "Ga0":"Ga-0", "Got22":"Got-22", "Got7":"Got-7", "Got7A":"Got-7", "Gu0":"Gu-0", "Gy0":"Gy-0", "HR10":"HR-10", "HR5":"HR-5", "Kas2":"Kas-2", "Kas1":"Kas-2", "Kin0":"Kin-0", "Knox10":"Knox-10", "Kno10":"Knox-10", "Knox18":"Knox-18", "Kno18":"Knox-18", "Kondara":"Kondara", "Kz1":"Kz-1", "Kz9":"Kz-9", "Kz_9":"Kz-9", "LL0":"LL-0", "Ler1A":"Ler-1", "Ler1":"Ler-1","LY3":"Ler-1", "Lov1":"Lov-1", "Lov_1":"Lov-1", "Lov5":"Lov-5", "Lov5B":"Lov-5", "LOV5A":"Lov-5", "Lp22":"Lp2-2", "LP22":"Lp2-2", "Lp26":"Lp2-6", "LP26":"Lp2-6", "Lz0":"Lz-0", "Mr0":"Mr-0", "Mrk0":"Mrk-0", "Ms0":"Ms-0", "Ms_0":"Ms-0", "Mt0":"Mt-0", "Mz0":"Mz-0", "NFA10":"NFA-10", "Nfa10":"NFA-10", "NFA8":"NFA-8", "Nafa8B":"NFA-8", "Nd1":"Nd-1", "Nok3":"Nok-3", "Nok_3":"Nok-3", "Omo21":"Omo2-1", "Omo23":"Omo2-3", "Omo2_3":"Omo2-3", "Oy0":"Oy-0", "Pna10":"Pna-10", "Pna_10":"Pna-10", "Pna17":"Pna-17", "Pna_17":"Pna-17", "Pro0":"Pro-0", "Pu223":"Pu2-23", "Pu27":"Pu2-7", "Pu2_7":"Pu2-7", "RRS10":"RRS-10", "Rrs10B":"RRS-10", "RRS7":"RRS-7", "Rrs7B":"RRS-7", "Ra0":"Ra-0", "Ren1":"Ren-1", "Ren11":"Ren-11", "Ren_11":"Ren-11", "RmxA02":"Rmx-A02", "RmxA180":"Rmx-A180", "Se0":"Se-0", "Se0_2":"Se-0", "ShaA":"Shahdara", "Sorbo":"Sorbo", "Spr12":"Spr1-2", "Spr1_2":"Spr1-2", "Spr16":"Spr1-6", "Spr1_6":"Spr1-6", "Sq1":"Sq-1", "Sq_1":"Sq-1", "Sq8":"Sq-8", "Tamm2":"Tamm-2", "Tamm2B":"Tamm-2", "Tamm27":"Tamm-27", "Tamm_27":"Tamm-27", "Ts1":"Ts-1", "Ts1A":"Ts-1", "Ts5":"Ts-5", "Ts_5":"Ts-5", "Tsu1":"Tsu-1", "Tsu1B":"Tsu-1", "Ts1B":"Tsu-1", "Ull23":"Ull2-3", "Ull25":"Ull2-5", "Uod1":"Uod-1", "Uod_1":"Uod-1", "Uod7":"Uod-7", "Van0":"Van-0", "VanOA":"Van-0", "Var21":"Var2-1", "Var26":"Var2-6", "Var2_6":"Var2-6", "Wa1":"Wa-1","Wa_1":"Wa-1", "Wei0":"Wei-0", "Ws0":"Ws-0", "Ws2":"Ws-2", "Ws_2":"Ws-2", "Wt5":"Wt-5", "Yo0":"Yo-0", "Zdr1":"Zdr-1", "Zdr_1":"Zdr-1", "Zdr6":"Zdr-6"}


    if imputed:
        datadir = homedir+"Projects/data/NPUTE_results/"
        datafile = datadir+"250K_snp2acc_NPUTE.csv"
        datafilecol = datadir+"250K_snp2acc_NPUTE.colnames.csv"
        datafilerow = datadir+"250K_snp2acc_NPUTE.rownames.csv"
    else:
        datadir = homedir+"Projects/data/020808dump/"
        datafile = datadir+"250K_snp2acc_GOOD.csv"
        datafilecol = datadir+"250K_snp2acc_GOOD.colnames.csv"
        datafilerow = datadir+"250K_snp2acc_GOOD.rownames.csv"
        
    positions = [] #list[chr][position_index]
    genotypes = [] #list[chr][position_index][acces]
    accessions = []
    
    phenotype = None
    filteredPhenotypes = [] #list[acces]
    individuals = [] #list[acces][chr][position_index]  Filtered for specific phenotypes.
    
    #Reading column data
    f = open(datafilecol, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if imputed:
            acc = line.rstrip()
        else:
            acc = (line.rstrip().split("."))[1]            
        accessions.append(accessionName2Id[accessions250To2010[acc]])

    #Reading row data
    f = open(datafilerow, 'r')
    lines = f.readlines()
    f.close()
    positions = [[],[],[],[],[]] #1 list per chromasome
    for line in lines:
        if imputed:
            l = line.rstrip().split(".")
            chr = int(l[0])
            pos = int(l[1])
        else:
            pos = int(line.rstrip())
            chr = pos/100000000
            pos = pos%100000000
        positions[chr-1].append(pos)

    #Reading genotype data
    rawgenotypes = [[],[],[],[],[]] #1 list per chromasome
    f = open(datafile, 'r')
    for i in range(0,5):
        for p in positions[i]:
            line = f.readline()
            l = (line.rstrip()).split(",")
            if l == []:
                raise Exception("Data problem") 
            rawgenotypes[i].append(l)
    f.close()
    print "raw genotype read"

    import random
    #Converting genotype and filtering.
    decoder = RawDecoder()
    genotypes = [[],[],[],[],[]]
    newpositions = [[],[],[],[],[]]
    for i in range(0,len(rawgenotypes)):
        for j in range(0,len(positions[i])):
            l = []
            for nt in rawgenotypes[i][j]:
                l.append(decoder[nt])
            genotypes[i].append(l)
            newpositions[i].append(positions[i][j])
        print "Chromosome",i+1,":",len(newpositions[i]),"SNPs."

    """
    for i in range(0,5):
        print newpositions[i][-1]
    """
    del positions
    del rawgenotypes
    positions = newpositions
    
    chromasomes = []
    for i in range(0,5):
        chromasomes.append(RawSnpsData(genotypes[i],positions[i],accessions=accessions))
            
    return(chromasomes)




def parse250KDataFiles(imputed = True):
    """
    Returns 250K Data as a list of SnpsData objects (not RawSnpsData).
    
    Set imputed to False, if un-imputed data is preferred.

    """
    
    if imputed:
        datadir = homedir+"Projects/data/NPUTE_results/"
        datafile = datadir+"250K_snp2acc_NPUTE.csv"
        datafilecol = datadir+"250K_snp2acc_NPUTE.colnames.csv"
        datafilerow = datadir+"250K_snp2acc_NPUTE.rownames.csv"
    else:
        datadir = homedir+"Projects/data/020808dump/"
        datafile = datadir+"250K_snp2acc_GOOD.csv"
        datafilecol = datadir+"250K_snp2acc_GOOD.colnames.csv"
        datafilerow = datadir+"250K_snp2acc_GOOD.rownames.csv"
        
    positions = [] #list[chr][position_index]
    genotypes = [] #list[chr][position_index][acces]
    accessions = []
    
    phenotype = None
    filteredPhenotypes = [] #list[acces]
    individuals = [] #list[acces][chr][position_index]  Filtered for specific phenotypes.
    
    #Reading column data
    f = open(datafilecol, 'r')
    lines = f.readlines()
    f.close()
    for line in lines:
        if imputed:
            acc = line.rstrip()
        else:
            acc = (line.rstrip().split("."))[1]            
        accessions.append(acc)

    #Reading row data
    f = open(datafilerow, 'r')
    lines = f.readlines()
    f.close()
    positions = [[],[],[],[],[]] #1 list per chromasome
    for line in lines:
        if imputed:
            l = line.rstrip().split(".")
            chr = int(l[0])
            pos = int(l[1])
        else:
            pos = int(line.rstrip())
            chr = pos/100000000
            pos = pos%100000000
        positions[chr-1].append(pos)

    #Reading genotype data
    rawgenotypes = [[],[],[],[],[]] #1 list per chromasome
    f = open(datafile, 'r')
    for i in range(0,5):
        for p in positions[i]:
            line = f.readline()
            l = (line.rstrip()).split(",")
            if l == []:
                raise Exception("Data problem") 
            rawgenotypes[i].append(l)
    f.close()
    print "raw genotype read"

    import random
    #Converting genotype and filtering.
    countAll = 0
    countGood = 0 
    countStupid = 0
    decoder = {'NA':'NA'}
    genotypes = [[],[],[],[],[]]
    newpositions = [[],[],[],[],[]]
    for i in range(0,len(rawgenotypes)):
        for j in range(0,len(positions[i])):
            countAll = countAll+1
            k = 0
            ntl = [] #list of observed nucleotides.
            for nt in ['A','C','G','T']:
                if nt in rawgenotypes[i][j]:
                    decoder[nt]=k
                    ntl.append(nt)
                    k = k+1
            if k==2:
                countGood = countGood + 1
                l = []
                for nt in rawgenotypes[i][j]:
                    l.append(decoder[nt])
                genotypes[i].append(l)
                newpositions[i].append(positions[i][j])
            else:
                if k==1:
                    countStupid = countStupid+1
    print countAll," SNPs in all"
    print countGood," SNPs used"
    print countStupid," Stupid SNPs thrown away"

    for i in range(0,5):
        print newpositions[i][-1]
    del positions
    del rawgenotypes
    positions = newpositions
    
    chromasomes = []
    for i in range(0,5):
        chromasomes.append(SnpsData(genotypes[i],positions[i],accessions=accessions))
            
    return(chromasomes)



def parseRawDataFile(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()


"""
def parseDataFile(filename, missingDataLimit = 0, filterProb=1):
    snps = []
    positions = []
    f = open(filename, 'r')
    lines =	f.readlines()
    f.close()
    accessions = (lines[0].rstrip()).split()
    accessions = accessions[1:len(accessions)-1]
    num = len(accessions)	
    decoder	= {'N':-1}
    i =1
    while i	< len(lines) and lines[i][0].isdigit():
        if random.random()<=filterProb:
            line = (lines[i].rstrip()).split()
            snp = line[1:num+1]
            if (snp.count('N') <= float(num)*missingDataLimit) and int('A' in snp)+ int('C' in snp) + int('G' in snp) + int('T' in snp)==2:
                snp = []
                k = 0
                
                for lt in ['A','C','G','T']:
                    if lt in line[1:num+1]:
                        decoder[lt]=k
                        k = k+1
                for j in range(1,num+1):
                    snp.append(decoder[line[j]])
                snps.append(snp)
                positions.append(int(line[0]))
		
        i = i+1
    print "File loaded, filtered, and converted."		
    snpsd = SnpsData(snps,positions)
    return snpsd
"""

def parseMSFile(filename, baseScale=1000000):
    f = open(filename, 'r')
    lines =	f.readlines()
    f.close()
    i=0
    data = []
    invBaseScale = 1/float(baseScale)
    while i	< len(lines):
        line = lines[i]
        if line.startswith("//"):
            num =0
            positions = []
            snps = []
            i = i+1			
            while i < len(lines) and not lines[i].startswith("//"):
                line = lines[i]	
                if line.startswith("segsites:"):
                    num = int(line[9:].rstrip())
                if line.startswith("positions:"):
                    l1 = line[10:].rstrip().split()
                    l2 = [0.0]*len(l1)
                    for j in range(0,len(l1)):
                        l2[j]=float(l1[j])
                        snps.append([])	    #Initializing the snps.
                    positions = l2
                if line[0].isdigit():
                    line = line.rstrip()
                    snps[0].append(int(line[0]))
                    for j in range(1,len(positions)):
                        snps[j].append(int(line[j]))
                                    
                i = i+1
            newSnps = []
            newPositions = []
            if len(positions)>0:
                newSnps.append(snps[0])
                newPositions.append(positions[0])
            k = 0 
            for j in range(1,len(positions)):
                if (positions[j]-positions[j-1])>=invBaseScale:
                    newSnps.append(snps[j])
                    newPositions.append(positions[j])
                else:
                    if baseScale>10000 and k<=j:  #The precision of ms is 1/10000 (sometimes we have many more snps)
                        k = j+1
                        while(k < len(positions) and (positions[k]-positions[j-1])<invBaseScale):
                            k = k+1
                        last = 0.0
                        for h in range(0,k-j): # Order statistic 
                            last = random.betavariate(1,k-j-h)*(1-last)+last
                            positions[j+h] = positions[j+h]+last/10000
                        if (positions[j]-positions[j-1])>=invBaseScale:
                            newSnps.append(snps[j])
                            newPositions.append(positions[j])
			
            data.append(SnpsData(newSnps,newPositions,baseScale))
        else:
            i = i+1
    return data


def parseMSFileFilter(filename, baseScale=1000000, fixPos=True, filterProb=1.0):
    f = open(filename, 'r')
    lines =	f.readlines()
    f.close()
    i=0
    data = []
    invBaseScale = 1/float(baseScale)
    while i	< len(lines):
        line = lines[i]
        if line.startswith("//"):
            num =0
            positions = []
            snps = []
            i = i+1			
            while i	< len(lines) and not lines[i].startswith("//"):
                line = lines[i]	
                if line.startswith("segsites:"):
                    num = int(line[9:].rstrip())
                if line.startswith("positions:"):
                    l1 = line[10:].rstrip().split()
                    l2 = [0.0]*len(l1)
                    for j in range(0,len(l1)):
                        l2[j]=float(l1[j])
                        snps.append([])	    #Initializing the snps.
                    positions = l2
                if line[0].isdigit():
                    line = line.rstrip()
                    snps[0].append(int(line[0]))
                    for j in range(1,len(positions)):
                        snps[j].append(int(line[j]))
                                    
                i = i+1
            newSnps = []
            newPositions = []
            l = 0
            pos = -1
            for j in range(0,len(positions)):
                if random.random()<=filterProb:
                    npos = positions[j]
                    if (npos-pos)>=invBaseScale:
                        newSnps.append(snps[j])
                        newPositions.append(npos)
                    else:
                        if fixPos and baseScale>10000 and l<=j :  #The precision of ms
                            l = j+1
                            while(l<len(positions) and (positions[l]-pos)<invBaseScale):
                                l = l+1
                            last = 0.0
                            for h in range(0,l-j): # Order statistic 
                                last = random.betavariate(1,l-j-h)*(1-last)+last
                                positions[j+h] = positions[j+h]+last/10000
                                                
                            if (npos-pos)>=invBaseScale:
                                newSnps.append(snps[j])
                                newPositions.append(npos)
                    pos = positions[j]
            data.append(SnpsData(newSnps,newPositions,baseScale))
        else:
            i = i+1
    return data



def parseMSData(filename, baseScale=1000000,sampleNum=None, fixPos=True):
    """ 
    Parses randomly chosen ms outputs from previously calculated files.  Depends on the datafiles!
    """
    f = open(filename, 'r')
    lines =	f.readlines()
    totSampleNum = int(lines[0].split()[2])
    f.close()
    i=0
    klist = []
    if sampleNum:
        k = 0
        while k < sampleNum:
            val = int(random.random()*totSampleNum+1)
            if not val in klist:
                klist.append(val)
                k = k+1
    else:
        klist = range(0,totSampleNum)
    k = -1
    data = []
    invBaseScale = 1/float(baseScale)
    while i < len(lines):
        line = lines[i]
        if line.startswith("//"):
            k = k+1
            i = i+1	
            if k in klist:
                num =0
                positions = []
                snps = []						
                while i	< len(lines) and not lines[i].startswith("//"):
                    line = lines[i]	
                    if line.startswith("segsites:"):
                        num = int(line[9:].rstrip())
                    if line.startswith("positions:"):
                        l1 = line[10:].rstrip().split()
                        l2 = [0.0]*len(l1)
                        for j in range(0,len(l1)):
                            l2[j]=float(l1[j])
                            snps.append([])	    #Initializing the snps.
                        positions = l2
                    if line[0].isdigit():
                        line = line.rstrip()
                        snps[0].append(int(line[0]))
                        for j in xrange(1,len(positions)):
                            snps[j].append(int(line[j]))
                                            
                    i = i+1
                newSnps = []
                newPositions = []
                if len(positions)>0:
                    newSnps.append(snps[0])
                    newPositions.append(positions[0])
                    pos = positions[0]
                l = 0
                for j in range(1,len(positions)):
                    npos = positions[j]
                    if (npos-pos)>=invBaseScale:
                        newSnps.append(snps[j])
                        newPositions.append(npos)
                    else:
                        if fixPos and baseScale>10000 and l<=j :  #The precision of ms
                            l = j+1
                            while(l<len(positions) and (positions[l]-pos)<invBaseScale):
                                l = l+1
                            last = 0.0
                            for h in range(0,l-j): # Order statistic 
                                last = random.betavariate(1,l-j-h)*(1-last)+last
                                positions[j+h] = positions[j+h]+last/10000
                                                    
                            if (npos-pos)>=invBaseScale:
                                newSnps.append(snps[j])
                                newPositions.append(npos)
                    pos = positions[j]
                data.append(SnpsData(newSnps,newPositions,baseScale))
        else:
            i = i+1
    return data
    

def parseMSDataFilter(filename, baseScale=1000000,sampleNum=None, fixPos=True, filterProb = 1.0):
    """ 
    Parses randomly chosen ms outputs from previously calculated files.  Depends on the datafiles!
    """
	#print filename
    f = open(filename, 'r')
    lines =	f.readlines()
    totSampleNum = int(lines[0].split()[2])
    f.close()
    i=0
    klist = []
    if sampleNum:
        k = 0
        while k < sampleNum:
            val = int(random.random()*totSampleNum+1)
            if not val in klist:
                klist.append(val)
                k = k+1
    else:
        klist = range(0,totSampleNum)
    k = -1
    data = []
    invBaseScale = 1/float(baseScale)
    while i	< len(lines):
        line = lines[i]
        if line.startswith("//"):  #New data
            k = k+1
            i = i+1	
            if k in klist:
                num =0
                positions = []
                snps = []						
                while i	< len(lines) and not lines[i].startswith("//"):
                    line = lines[i]	
                    if line.startswith("segsites:"):
                        num = int(line[9:].rstrip())
                    if line.startswith("positions:"):
                        l1 = line[10:].rstrip().split()
                        l2 = [0.0]*len(l1)
                        for j in range(0,len(l1)):
                            l2[j]=float(l1[j])
                            snps.append([])	    #Initializing the snps.
                        positions = l2
                    if line[0].isdigit():
                        if len(positions)==len(line.rstrip()):  #Hack added to fix ms generated errors in files!
                            line = line.rstrip()
						        #snps[0].append(int(line[0]))
                            for j in xrange(0,len(line)):
                                snps[j].append(int(line[j]))
                        else:
                            print "line nr.",i
                            print "positions:",positions
                            print "line:",line.rstrip()
                            print "old line:",line
                            raise Exception()
                    i = i+1
                newSnps = []
                newPositions = []
                l = 0
                pos = -1
                debug1 = positions[:]
                debug2 = snps[:]
                r = 0				
                for j in range(0,len(positions)):
                    if random.random()<=filterProb:
                        r = r + 1
                        npos = positions[j]
                        if (npos-pos)>=invBaseScale:
                            newSnps.append(snps[j])
                            newPositions.append(npos)
                        else:
                            if fixPos and baseScale>10000 and l<=j :  #The precision of ms
                                l = j+1
                                while(l<len(positions) and (positions[l]-pos)<invBaseScale):
                                    l = l+1
                                last = 0.0
                                for h in range(0,l-j): # Order statistic 
                                    last = random.betavariate(1,l-j-h)*(1-last)+last
                                    positions[j+h] = positions[j+h]+last/10000
                
                                if (npos-pos)>=invBaseScale:
                                    newSnps.append(snps[j])
                                    newPositions.append(npos)
                        pos = positions[j]
                if r ==0 :
                    print "No luck"
                    print len(positions)
                if len(newSnps)==0 or len(newSnps[0])==0:
                    print "data is empty!!!"
                    print "It succeded ",r,"number of times."
                    print debug1
                    print newPositions
                    print newSnps
                    print debug2
                snpsd = SnpsData(newSnps,newPositions,baseScale)
                if len(snpsd.snps)==0 or len(snpsd.snps[0])==0:
                    print snpsd.positions
                    print snpsd.snps
                data.append(snpsd)
        else:
            i = i+1
    if data[0]==[]:
        print 
    return data
    





#--------------------------------------------------------------------------------#

def _testDBParser_():
    
    snpsDatas250 = parse250DataRaw()


    #snpsDatasPerlgen = get2010DataFromDb()
    snpsDatasPerlgen = getPerlgenDataFromDb()
    print len(snpsDatasPerlgen)
    for snpsd in snpsDatasPerlgen:
        print len(snpsd.positions)
    
    
    snpsDatas250[0].compareWith(snpsDatasPerlgen[0])

    merged250nPerlgen = []
    for i in range(0,len(snpsDatas250)):
        merged250nPerlgen.append(snpsDatasPerlgen[i].mergeData(snpsDatas250[i]))
    
    writeRawSnpsDatasToFile("Perlgen_250.out",merged250nPerlgen)

    snpsDatas2010 = get2010DataFromDb()
    print len(snpsDatas2010)
    for snpsd in snpsDatas2010:
        print len(snpsd.positions)


    snpsDatas2010[0].compareWith(snpsDatasPerlgen[0])

    merged250n2010nPerlgen = []
    for i in range(0,len(merged250nPerlgen)):
        merged250n2010nPerlgen.append(snpsDatas2010[i].mergeData(merged250nPerlgen[i]))

    writeRawSnpsDatasToFile("2010_Perlgen_250.out",merged250n2010nPerlgen)

   
    """  
    snpErrorRate = []
    for chr in [0,1,2,3,4]:
        snpErrorRate +=snpsds2[chr].compareWith(snpsds1[chr])
    print "Average Genome Wide Snp Error:",sum(snpErrorRate)/float(len(snpErrorRate))
    """

if __name__ == "__main__":
    pass
    #_testDBParser_()
