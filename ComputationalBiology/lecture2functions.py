# RETURNS REVERSE COMPLEMENT SEQUENCE
# IF/ELIF/ELSE STATEMENTS

# to use functions, need to import file as XXX
# EX: import lecture2functions as f2

def compbase(base):
    '''Returns complement of base'''
    if base == "A":
        compbase = "T"
    elif base == "C":
        compbase = "G"
    elif base == "G":
        compbase = "C"
    elif base == "T":
        compbase = "A"
    else:
        compbase = base
        # if there's an N in the sequence, for example, keep the N
    return compbase


def compseq(seq):
    '''Returns complement of seq'''
    comp = []  # empty list
    for base in seq:
        nuc = compbase(base)
        comp.append(nuc)
    compstr = "".join(comp)
    return compstr


def rev(seq):
    '''Returns the reverse of seq'''
    return seq[::-1]


def revcomp(seq):
    '''Returns the reverse complement of seq'''
    rv = rev(seq)
    rvcp = compseq(rv)
    return rvcp


# COUNT THE NUMBER OF START CODONS ("ATG") ON FORWARD STRAND

def countSTART(seq):
    '''Count number start codons on forward strand in seq'''
    cnt = 0
    for i in range(len(seq)):
        if seq[i:(i + 3)] == "ATG":
            cnt = cnt + 1
    return cnt


def countSTART2(seq):
    '''Count number start codons on forward strand in seq'''
    return seq.count("ATG")  # Python has many shortcuts


# COUNT THE NUMBER OF STOP CODONS ("TAG","TAA", OR "TGA") ON FORWARD STRAND

def countSTOP(seq):
    '''Count number stop codons on forward strand in seq'''
    stopcod = ["TAG", "TAA", "TGA"]
    cnt = 0
    for i in range(len(seq)):
        if seq[i:(i + 3)] in stopcod:  # nice shortcut
            cnt = cnt + 1
    return cnt


# IMPORTANT FUNCTION FOR YOUR FIRST HW ASSIGNMENT

# to use file, first need to name variable
# EX: ace2 = f2.loadFASTA("ACE2.fasta")

def loadFASTA(filename):
    '''Outputs a sequence string from the FASTA file named filename'''
    infile = open(filename)  # opens the file
    seqlist = []
    temp = infile.readline()  # reads a single line
    # we don't want the first line
    # of the FASTA file ">..."
    for line in infile:  # reads one line at a time
        # when it gets to the end of the file it stops
        temp = line.replace("\n", "")  # removes \n
        temp = temp.replace("\r", "")  # removes \r too
        seqlist.append(temp)
    infile.close()  # closes the file
    seq = "".join(seqlist)  # combines the list into a string
    return seq


# there are many different ways to do the same thing in Python
# the advantage of using "with" below is that you don't have to remember
# to close the file when you are done with it

def loadFASTA2(filename):
    '''Outputs a sequence string from the FASTA file named filename'''
    seqlist = []
    with open(filename, "r") as infile:
        temp = infile.readline()  # reads a single line
        # we don't want the first line
        # of the FASTA file ">..."
        for line in infile:  # reads one line at a time
            # when it gets to the end of the file it stops
            temp = line.replace("\n", "")  # removes \n
            temp = temp.replace("\r", "")  # removes \r too
            seqlist.append(temp)
    seq = "".join(seqlist)  # combines the list into a string
    return seq


# DICTIONARIES

def dicDNA2aa():
    '''Returns DNA to amino acid dictionary'''
    aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
          '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
          'G']
    # aa is list of amino acids
    codons = [['TTT', 'TTC'],
              ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
              ['ATT', 'ATC', 'ATA'],
              ['ATG'],
              ['GTT', 'GTC', 'GTA', 'GTG'],
              ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
              ['CCT', 'CCC', 'CCA', 'CCG'],
              ['ACT', 'ACC', 'ACA', 'ACG'],
              ['GCT', 'GCC', 'GCA', 'GCG'],
              ['TAT', 'TAC'],
              ['TAA', 'TAG', 'TGA'],
              ['CAT', 'CAC'],
              ['CAA', 'CAG'],
              ['AAT', 'AAC'],
              ['AAA', 'AAG'],
              ['GAT', 'GAC'],
              ['GAA', 'GAG'],
              ['TGT', 'TGC'],
              ['TGG'],
              ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
              ['GGT', 'GGC', 'GGA', 'GGG']]
    # codons is a list of lists: each element is a list of codons
    # the first element has 2 codons which code for the first element in aa
    # the second element has 6 codons which code for the second element in aa
    # etc.
    # it is important that aa and codons are in the same order

    # create aaDIC dictionary, codons are keys, amino acids are values
    aaDIC = {}  # empty dictionary
    for i in range(len(codons)):
        for j in range(len(codons[i])):
            aaDIC[codons[i][j]] = aa[i]
    return aaDIC


def dna2aa(seq, aaDIC):
    '''Returns amino acid sequence coded by seq (assume in frame)'''
    pro = []
    for i in range(0, len(seq), 3):  # increment by 3
        cod = seq[i:(i + 3)]
        acid = aaDIC.get(cod, "?")
        # if the codon is not in the dictionary returns ?

        pro += acid  # does NOT create a new list
    prostring = "".join(pro)  # makes a string
    return prostring


# DON'T PEEK
# CLASSWORK SOLUTIONS BELOW
# COMPARE TWO SEQUENCES
# Assume they are aligned and of the same length
# Calculate the number nucleotide differences, number syn and nonsyn changes

#
# My attempts
#

#bob = "ATGTTTTATCCAGGTGATTAA"
#tom = "ATGTTCTACCCAGGTAATTAA"
#bobPro = f2.dna2aa(bob, aaDIC)
#tomPro = f2.dna2aa(tom, aaDIC)

def compareNuc(seq1, seq2):
    '''Compares 2 sequences and returns # of nucleotide differences'''
    nucDiff = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            nucDiff += 1
    return nucDiff


def compareAA(seq1, seq2, aaDIC):
    '''Compares 2 sequences and calculates # syn & nonsyn AA'''
    syn = 0
    nonsyn = 0
    seq1pro = dna2aa(seq1, aaDIC)
    seq2pro = dna2aa(seq2, aaDIC)
    for i in range(len(seq1pro)):
        if seq1pro[1] == seq2pro[1]:
            syn += 1
        else:
            nonsyn += 1
    print("# of synonymous AA: " + syn)
    print("# of non-synonymous AA: " + nonsyn)






#
# ANSWERS BELOW
#

def countdiff(seq1, seq2):
    '''Returns number of differences in seq1, seq2 (assume aligned, same length)'''
    cnt = 0
    for i in range(len(seq1)):
        if (seq1[i] != seq2[i]):
            cnt = cnt + 1
    return cnt


def countNONSYN(seq1, seq2, aaDIC):
    '''Returns number amino acid differences coded by seq1, seq2 (assume aligned, same length, in frame)'''
    pro1 = dna2aa(seq1, aaDIC)
    pro2 = dna2aa(seq2, aaDIC)
    return countdiff(pro1, pro2)
