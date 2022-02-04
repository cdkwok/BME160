#!/usr/bin/env python3
# Name: Christie Kwok (cdkwok)

class ProteinParam :
    '''
    Calculates the total amino acid count, molecular weight, Molar extinction coefficient, 
    Mass Extinction Coefficient, pI, and amino acid composition when given an input string.
    
    Author: Christie Kwok
    Date: April 20, 2020
    Input is a string of characters.
    Initalized: input string reprenting the protein sequence, and a dictionary of amino acid
                counts of the intialized input string.
    Class Attributes: mwH2O, aaNterm, aaCterm 
                    Dictionaries aa2mw,  aa2abs280, aa2chargePos, aa2chargeNeg
    Instance Attributes: protein, aaCountDictionary
    Methods: __init__(), aaCount(), aaComposition(), molecularWeight(), _charge_(), pI(), 
            molarExtinction(), massExtinction().
    '''
# Constant Defintion:
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        Initializes and capitalizes the object(input string) and a dictionary of 
        amino acid counts from the input. 
        '''
        
        self.protein = protein.upper() #instance attribute of capitalized input string
        
        
        # instance attribute of a  dictionary of each amino acid and 
        #the counts of the amount of time it shows up in the input
        self.aaCountDictionary = {
            'A' : self.protein.count('A'), 'C' : self.protein.count('C'), 'D' : self.protein.count('D'), 
            'E' : self.protein.count('E'), 'F' : self.protein.count('F'), 'G' : self.protein.count('G'), 
            'H' : self.protein.count('H'), 'I' : self.protein.count('I'), 'L' : self.protein.count('L'), 
            'K' : self.protein.count('K'), 'M' : self.protein.count('M'), 'N' : self.protein.count('N'), 
            'P' : self.protein.count('P'), 'Q' : self.protein.count('Q'), 'R' : self.protein.count('R'), 
            'S' : self.protein.count('S'), 'T' : self.protein.count('T'), 'V' : self.protein.count('V'), 
            'Y' : self.protein.count('Y'), 'W' : self.protein.count('W')
        }
        # I tried removing Dictionary from the name and everywhere else it shows up in my code,
        # and it compltely breaks my code and I don't understand why
        

    def aaCount (self):
        '''The input is the initalized object.
           The output is the total amino acid count.
        '''
        return sum(self.aaCountDictionary.values()) #sums the values of aaCountDictionary

    def pI (self):
        '''The input is the net_charge from _charge_ 
           The output is the pH at which the protein is a zwitterion.
        
        Logic: I had wanted to find the pH at which pI = 0 for the protein. So I call charge to find that pH.
        The pH covers a range of 0 to 14, with everything in between. So I used a range of 0 to 1400 divide by 100
        to cover for the in between pHs.
        
        note: received help from lecture(David Bernick) and from LSS(Mohammad Abdulqader)
        '''
        maxPossibleCharge = 10000000 #arbitrary value for comparison
        maxPossiblepH = 10000000 #set to initial value
        
        #this for loop checks and if statement continuously calculates the charge of the 
        #protein at each pH in the range
        for pH in range(0, 1401):
            floatpH = pH / 100 
                #divides pH by 100 to get a float number for all possible values, and sets it to floatpH
            proteinCharge = abs(self._charge_(floatpH)) 
                #calls _charge_ and sets it to proteinCharge
            if proteinCharge < maxPossibleCharge:
                maxPossibleCharge = proteinCharge
                maxPossiblepH = floatpH
        return(maxPossiblepH)

    def aaComposition (self) :
        '''Input is is the object, returns the counts of the dictionary(aaCountDictionary)'''
        return self.aaCountDictionary

    def _charge_ (self,pH):
        '''Input is from the floatpH from pI.
           Output is the net charge of the protein.
        
        Logic:
        First, the individual charges for the N-terminus and C-terminus were calculated.
        Then I used for loops to do the summation of the postitively charged and negatively 
        charged amino acids. To get the net charge, the postive charge and negative charge were
        subtracted. This method takes the pH given by pI(), and the values from the dictionaries
        as part of its calculations.
        
        Got help from LSS(Mohammad Abdulqader)'''
        
        posCharge = ((10 ** self.aaNterm)/(10 ** (pH) + 10 ** self.aaNterm))
            #this calculates the postivie charge of the N-terminus 
            
        negCharge =((10 ** pH)/(10 ** (pH) + 10 ** self.aaCterm))
            #this calculates the negative charge of the C-terminus 
            
        for aa in self.aa2chargePos:
            posCharge += self.aaCountDictionary[aa] * ((10 ** self.aa2chargePos[aa])/(10 ** (pH) + 10 ** self.aa2chargePos[aa]))
            #calculates the total positive charge of the protein
            #iterates through aa2chargePos and pulls the count of the positively charged aa
            # from aaCountDictionary
            
        for aa in self.aa2chargeNeg:
            negCharge += self.aaCountDictionary[aa] * ((10 ** pH)/(10 ** (pH) + 10 ** self.aa2chargeNeg[aa]))
            #calculates the total negative charge of the protein
            #iterates through aa2chargeNeg and pulls the count of the negatively charged aa
            # from aaCountDictionary
            
        netCharge = posCharge - negCharge #total charge
        
        return(netCharge)

    def molarExtinction (self):
        '''Input is the counts of Trp, Tyr, and Cys from aaCountDictionary, 
           Output is the molar extinction value
           
           Logic:
           I found the individual molar extinction values for Tyr, Trp, and Cys by
           getting their individual counts and multiplying it their coefficient values.
           Then the molar extinction values were summed to get a total molar extinction value.
           '''
        
        tyrME = self.aaCountDictionary.get('Y') * self.aa2abs280.get('Y')
            #calculates the molar extinction of Tyr by pulling the values from the respective dictionaries
            
        trpME = self.aaCountDictionary.get('W') * self.aa2abs280.get('W')
            #calculates the molar extinction of Trp by pulling the values from the respective dictionaries
            
        cysME = self.aaCountDictionary.get('C') * self.aa2abs280.get('C')
            #calculates the molar extinction of Cys by pulling the values from the respective dictionaries
            
        totalME = tyrME + trpME + cysME 
            #sums them together
        
        return(totalME)
        

    def massExtinction (self):
        '''Input is the molecular weight of the protein 
           Output is the mass extinction'''
        
        myMW =  self.molecularWeight()
            #sets the molecular weight of the protein to myMW
            
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''input is the amino acid counts from aaCountDictionary 
           Output is the molecular weight of the protein
           
           Logic:
            I used a for loop to for the summation of the amino acids. This is specific to each
            amino acid(its counts and its molecular weight). The weight of water was subtracted 
            to represent the loss of of water mass during peptide bond formation. A single water
            was added after the loop, as adding that water during the loop caused water to be added
            for every amino acid(even the ones that had zero counts).'''
        
        aaMW = 0
        MW = 0
        for aa in self.aaCountDictionary:
            aaMW += ((self.aaCountDictionary[aa]) * (self.aa2mw[aa] - self.mwH2O))
            #calculates the total molecular weight of the amino acids

        if aaMW != 0:
	#this checks for single non-amino acid characters   
            MW = aaMW + self.mwH2O   
            #adds the weight of water to the total aa molecular weight
            
        return(MW)


class NucParams:
    '''Calcuates amino acid composition, codon composition, nucleotide composition, and nucleotide count of the input.
        
        Author: Christie Kwok
        Date: April 27, 2020
        Input: a string of characters
        Initialized: the input, dictionaries: amino acid composition, codon composition, 
                nucleotide composition, and nucleotide count
        Class Attributes: Dictionaries: rnaCodonTable, Set: validCharSet
        Instance Attributes: Dictionaries: aaComp, nucComp, codonComp, codonList
        Methods: __init__(), addSequence(), aaComposition(), nucComposition(), codonComposition(), nucCount()
        '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    validCharSet = {'A', 'G', 'C', 'T', 'U', 'N'} #set of valid nucleotides; 
        #I had to make this a class attribute because it kept giving me errors when initialized

    def __init__ (self, inString=''):
        ''' Initializes dictionaries, and the input.'''
        
        self.sequence = ''
            #initaliezes sequence as nothing (first)
            
        self.sequence = self.addSequence(inString)
            #initialzes self.sequence to be cleaned input from addSequence()
        
        self.aaComp = {aa : 0 for aa in NucParams.rnaCodonTable.values()}
            #initialzes a dictionary of amino acid counts
        
        self.nucComp =  {dnaNuc : 0 for dnaNuc in 'ATCGUN'}
            #initialzes  a dictionary of nucleotide counts
     
        self.codonComp = {codon : 0 for codon in NucParams.rnaCodonTable.keys()}
            #initialzes a  dictionary of codon counts
        
        self.codonList = []
            #initialzes a initalizes a list for later use
            
        
    def addSequence (self, inSeq):
        '''Input: inSeq
            Output: cleaned version of inSeq
            
            Logic: this takes the input, removes all invalid characters.
                It is then initalized in __init__ and appends to self.sequence to be utilized by other methods.
                I removed .upper() and .replace(), as FastAreader handles this. (compared to my last version).
                
            Note: I learned from LSS that it was ideal that we add to the dictionaries
            in addSeqeunce(), but I had coded everything before LSS.
        '''
        
        
        for char in inSeq:
            if char not in NucParams.validCharSet:
                inSeq.replace(char, "")
                #replaces non-valid characters with nothing
        
        #I had a problem where self.sequence would return NoneType, so this is to mitigate that
        #this forces self.sequence to equal inSeq, and continue to concatonate input
        if self.sequence == None:
            self.sequence = inSeq
        else:
            self.sequence = self.sequence + inSeq
        
        
        
        
    def aaComposition(self): 
        '''Input: a list of codons (in the input)
            Output: a dictionary of amino acids (in the input)
            
            Logic: checks if each index(codon) is in the rnaCodonTable. If true, then increments the count of that codon.
        '''
        
        for codon in self.codonList:
            #iterates through codonList
            aminoAcid = NucParams.rnaCodonTable.get(codon) 
                #gets the amino acid value using the codon key
            self.aaComp[aminoAcid] += 1
                #increments aacomp
                
        return self.aaComp
    
    def nucComposition(self): 
        '''Input: the cleaned up sequence
            Output: a dictionary of nucleotides (in the input)
            
            Logic: checks if the nucleotide is valid, then increments the counts of the nucleotude.
        '''
        
        for nuc in self.sequence:
            #iterates through the input
            if nuc in NucParams.validCharSet:
                #checks if each base in the input is valid
                self.nucComp[nuc] += 1
                #increments nucComp
    
        return self.nucComp
    
    def codonComposition(self):
        '''Input: cleaned up sequence(input)
            Output: a dictionary of codons (in the input)
            
            Logic: replaces all T's with U's, to represent RNA. This is then converted to a list, 
                with each index containing three characters of the string. If the index contains an invalid
                character(s), that index is ignored. If the index contains  valid character(s), than the
                codon count is incremented, and that index is appended to codonList.
        '''
        
        rnaString = self.sequence.replace("T", "U")
        
        
        rnaList = [rnaString[i:i+3] for i in range(0, len(rnaString), 3)]
            #this slices the string every 3 characters to make a codon and puts it in a list
            #cite: https://pythonexamples.org/python-split-string-into-specific-length-chunks/
                # mostly for syntax
       
        for codon in rnaList:
            if codon in NucParams.rnaCodonTable.keys():
                #checks if the codon is in the rnaCodonTable keys
                self.codonComp[codon] += 1
                    #increments codon
                self.codonList.append(codon)
                    #adds it to codonList
            
        return self.codonComp
    
    def nucCount(self):
        '''Input: dictionary of nucleotides
            Output: the total nucleotide count'''
        
        return sum(self.nucComp.values())
    
    
    



import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


class ORFfinder:
    '''Calcuates the frame, start position of the frame, end position of the frame, and length of the input sequence.
        
        Author: Christie Kwok
        Date: May 11, 2020
        Input: a string of characters, in the form of a DNA or RNA sequence. The RNA sequence will have 'U' replaced
            with 'T'. (Input: ACTGU)
        Initialized: the input with replacement of U to T, start and end codons,  short(est) and longest genes.
        Class Attributes: none
        Instance Attributes: sequence, startCodons, stopCodons, reverseSequence, header, minGene, longestGene 
        Methods: __init__(), orf(), reverse(), calculateORFs(),
        '''
    def __init__(self, sequence, startCodons, stopCodons , minGene, longestGene):
        ''' Initalize components'''
        self.sequence = sequence.replace('U', 'T') 
            # to insure we have a DNA sequence
        self.startCodons = startCodons
        self.stopCodons = stopCodons
        self.reverseSequence = self.reverse()
        self.minGene = minGene
        self.longestGene = longestGene
    
    def orf(self):
        ''' Find all the orfs of the sequence. 
            Input: sequence of characters from input (self.sequence)
            Output: all possible open reading frames
            
            Logic: First I set up a list of ORFs, and all the possible frames. Then I iterate through the possible frames
            and their nucleotides(by codon length) of the sequence. The sequence is sepearted into codons. If the codon
            is a start codon, then its position is added to the startSet. It then checks for a stop codon. When a stop codon
            is found, an calculateORF is called and a new ORF is made. It is added to orfList. The same happens for the
            reverse strand; reverse() is called to find the reverse complement of the original sequence.
        '''
        
        orfList= list()
            # a list that will contain all possible open reading frames
        codonLength = 3
        startSet = {0} 
            # handles dangling start
        frames = [i for i in range(1, 4)] 
            # frames 1, 2, and 3

        # for forward strand 
        for frame in frames:
            for nuc in range(frame - 1, len(self.sequence), codonLength):
                #for a nucleotide in the input sequence and in the correct frame
                codon = self.sequence[nuc: nuc + codonLength]
                
                    #seperates the sequence into codons

                if codon in self.startCodons:
                    startSet.add(nuc) 
                    
                        #add position of start codon
                elif codon in self.stopCodons:
                    
                    stopPosition = (nuc + 2)
                    newORF = self.calculateORFs(frame, startSet, stopPosition)
                        #add position of stop codon and make a new ORF
                    orfList += newORF
                        #add new ORF to the list
                    startSet = set()
            orfList += self.calculateORFs(frame, startSet, len(self.sequence)-1)
            startSet = {0}

        # for reverse complment
        #the same as above but for reverse sequence
        for frame in frames:
            for nuc in range(frame - 1, len(self.sequence), codonLength):
                codon = self.reverseSequence[nuc : nuc + codonLength]
                
                if codon in self.startCodons:
                    startSet.add(nuc)
                    
                elif codon in self.stopCodons:
                    stopPosition = (nuc + 2)
                    revORF = self.calculateORFs(frame * -1, startSet, stopPosition)
                    orfList += revORF
                    startSet = set()
            orfList += self.calculateORFs(frame * -1, startSet, len(self.sequence)-1)
            startSet = {0}

        return orfList
    
    def reverse(self):
        ''' Reverse compliment self.sequence.
            Input: self.sequence
            Output: reverse complement of self.sequence
            
            Logic: make the orignal sequence lowercase for ease of replacement. 
                (visually easier to determine what to replace)
        '''

        lowerSeq = self.sequence.lower()
        revSeq = lowerSeq.replace('a','T') #replace a with T
        revSeq = revSeq.replace('t', 'A') #replace t with A
        revSeq = revSeq.replace('c', 'G') #replace c with G
        revSeq = revSeq.replace('g', 'C') #replace g with C
        revSeq = revSeq[::-1]
        
        return revSeq

    def calculateORFs(self, frame, startSet, stopPosition):
        ''' caluclate the frame, start, stop and length
            
            Input: frame, startSet, stopPosition from orf()
            Output: the (+/-) frame number, start codon position, stop codon position, and gene length
            
            Logic: checks to see if there are start codons. Calculates the length of the gene(adding one because 
            it was off by one). Returns the ORF in the wanted output and appends that to a list.'''
        
        
        listORF = []

        
        if self.longestGene and len(startSet) > 0: 
            startSet = {min(startSet)}

        for startPosition in startSet:
            #position of the start codon(s)
            length = stopPosition - startPosition + 1
                #calculates the length, adding 1 because it was off by one - this fixes that
            if length >= self.minGene:
                #makes sure the gene length is greater than minGene
                if frame > 0:
                    orf = [frame, startPosition + 1, stopPosition + 1, abs(length)]
                    listORF.append(orf)
                if frame < 0:
                    orf = [frame, len(self.sequence) - stopPosition,
                           len(self.sequence)- startPosition, abs(length)]
                        #to get the correct start and stop position for reverse complement
                        #I substracted the position from the length of the sequence and swapped them
                    listORF.append(orf)
                
        return listORF
