#!/usr/bin/env python3
# Name: Christie Kwok
# Group Members: Priyanka Khaware (pkhaware)
# got help from: LSS(Tutor: Mohammad Abdulqader)
'''
Receives a FastQ (extension .fastq) file that is translated from one format to another. The possible input formats are
as follows: PHRED64, PHRED33, Solexa PHRED64, and PHRED64 with a string of B's at the end. The possible output formats
are PHRED64 and PHRED33.

Example:  From the sequence in illumina1.3.fastq file.
    Input format: -P64
    Output format: -P33
    
    Input: 
            @HWUSI-EAS300R_0005_FC62TL2AAXX:8:30:9097:12305#0/1
            GAGTGATAGCGGCACCAAGATGATCCATCTCGGCAGCAAC
            +HWUSI-EAS300R_0005_FC62TL2AAXX:8:30:9097:12305#0/1
            hhhhhhhhhhhghhhhghhhhhhhhhhhhhhhhhhhhhhc
    
    Output:
            @HWUSI-EAS300R_0005_FC62TL2AAXX:8:30:18447:12115#0/1
            CGTAGCTGTGTGTACAAGGCCCGGGAACGTATTCACCGTG
            +HWUSI-EAS300R_0005_FC62TL2AAXX:8:30:18447:12115#0/1
            BDEE?BB@;?E?EEDA?@2@BBBA@EED=EGEGGGG=GGG

'''
import sys
import math
class FastQreader :
    '''
    Takes in a FastQ file and yields the four lines of a FastQ file: @'header, sequence, '+' header/sequence 
    identifier, and quality value sequence
    
    
    Author: Christie Kwok and Priyanka Khaware
    Date: June 5, 2020
    Input: a FastQ file 
    Initalized: fname
    Class Attributes: N/A
    Instance Attributes: self.fname
    Methods: __init__(), doOpen(), readFasta()
    
    FastQreader frame provided by David Bernick in the form of FastAreader.
    '''
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        '''
        Input: FastQ file
        Output: @'header, DNA sequence, '+' header/sequence identifier, and quality value sequence
        
        Some things for the logic:
            I denote the following lines in the FastQ files as such:
            Line1: is the FastA header - this line starts with '@'
            Line2: is the DNA sequence
            Line3: is the qualityHeader - this line starts with '+'
            Line4: is the quality values of the sequence
            
        Logic: This first checks for an '@' character. When it find the '@', then
                it yields the header, sequence, qualityHeader, and qualitySequence. 
                The '@' is removed and the rest of the string as the header(Line1).
                It also resets sequence, qualityHeader, and qualitySequence to 
                nothing, and resets qualityCheck to False.
                
                Next, this moves down to the next line in the file, it checks if 
                it has a '+' in the first index. If it does, than it proceeds
                and removes the '+'. This is the qualityHeader(Line3). It also sets 
                qualityCheck to True. 
                    
                For the current line, if qualityCheck is True, then it strips the 
                line of whitespace and makes it uppercase. 
                (This is the quality values of the sequence: Line4.)
                    
                For the current line, if qualityCheck if False, then it strips
                the line of whitespace and makes it uppercase 
                (This is the DNA sequence: Line2.)
            
            '''
        
        header = ''
        sequence = ''
        qualityHeader = ''
        qualitySequence = ''
        qualityCheck = False
        #used to check if the FastA header and the quality header are the same
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            qualityHeader = ''
            qualitySequence = ''
            qualityCheck = False
            #used to check if the current line is sequence or quality score
            
            # skip to first fasta header
            line = fileH.readline()
            #checks if first line has a FastQ header
            while not line.startswith('@'):
                if not line: # we are at EOF
                    return header, sequence, qualityHeader, qualitySequence
                line = fileH.readline()
            header = line.rstrip()

        
            for line in fileH:
                if line.startswith ('@'):
                    #checks if its is a FastQ header
                    yield header,sequence, qualityHeader, qualitySequence
                    header = line.rstrip()
                    #resets the following to nothing, qualityCheck to False
                    sequence = ''
                    qualityHeader = ''
                    qualitySequence = ''
                    qualityCheck = False
                    
                else :
                    if line.startswith('+'):
                        #checks for quality sequence header 
                        qualityHeader = line.rstrip()
                        #removes the '+' from the header
                        qualityCheck = True
                        #confirms if the the current line is a quality header
                    else :
                        if qualityCheck:
                            #if qualityCheck is True
                            #this is the quality sequence values
                            #replaces extraneous characters with N
                            qualitySequence += ''.join(line.rstrip().split()).replace("*","N").replace(".","N").replace("n","N")
                        else :
                            sequence += ''.join(line.rstrip().split()).upper()
             

        yield header, sequence, qualityHeader, qualitySequence
        
class ConvertPHRED():
    '''
    Takes in the quality value sequence from FastQreader, and the parameters from CommandLine. Returns 
    the converted quality value sequence in the user wanted format.
    
    
    Author: Christie Kwok and Priyanka Khaware
    Date: June 5, 2020
    Input: qualitySequence, PHRED33input, PHRED64input, PHRED64Bin, PHRED64SOLin, PHRED33output, PHRED64output
    Initalized: qualitySequence, PHRED33input, PHRED64input, PHRED64Bin, PHRED64SOLin, PHRED33output, PHRED64output
    Class Attributes: N/A
    Instance Attributes: self.qualitySequence, self.PHRED33input, self.PHRED64input, self.PHRED64Bin, 
        self.PHRED64SOLin, self.PHRED33output, self.PHRED64output
    Methods: __init__(), qScore(), convertPHRED), convert()
    
    '''
    
    def __init__(self, sequence, qualitySequence, PHRED33input, PHRED64input, PHRED64Bin, PHRED64SOLin, PHRED33output, PHRED64output):
        '''Initialize Components'''
        
        self.sequence = sequence #initalizes the DNA sequence
        self.qualitySequence = qualitySequence #initializes the quality value sequence
        #initalizes the parameter booleans
        self.PHRED33input = PHRED33input 
        self.PHRED64input = PHRED64input
        self.PHRED64Bin = PHRED64Bin
        self.PHRED64SOLin = PHRED64SOLin
        self.PHRED33output = PHRED33output
        self.PHRED64output = PHRED64output
        
        
    
        
    def qScore(self, letter):
        '''
        Input: each letter from the quality value sequence
        Output: The converted  value 
        
        Logic: I check for each input parameter to determine if its True/False. If the parameter is true,
        then it will do the calculations for that parameter. For all inputs, the ASCII value of the character 
        of the quality sequence are substracted from an offset. For P33 input, offset is 33. For P64 input, 
        offset is 64. For P64 'B' input, offset is 64. Additionally, every B character is set to a P(error) of 1,
        equivalent to a qScore of 0. For P64 Solexa, a formula must be used to determine the correct P(error) and 
        qScore. 
        '''
                

        if self.PHRED33input: #if the input is True
            qScore = ord(letter) - 33 #offset -33 to get qScore
            #ord gets the ASCII value 
            return qScore
        
        elif self.PHRED64input: #if the input is True
            qScore = ord(letter) - 64 #offset -64 to get qScore
            return qScore
        
        elif self.PHRED64Bin: #if the input is True
            qScore = ord(letter)
            if qScore == 66:
                #check for B
                #if B, then its error = 1
                #so set this to equivalent ASCII of 0
                qScore = 0
                return qScore
            else:
                #if its not a 'B', proceed as normal
                qScore = qScore - 64
                return qScore
            
        elif self.PHRED64SOLin:
            # the forumla to determine pError 
            qSolexa = ord(letter) - 64 #qScore of te Solexa value
            pError = (10**(qSolexa/-10)) / (1 + 10**(qSolexa/-10)) # the forumla to determine pError of PHRED using Solexa qScore
            qScore = int(-10 * math.log(pError, 10)) #calculates actual qScore using the Perror of PHRED
            return qScore

        
    def convertPHRED(self,qScore):
        '''
        Input: qScore, output parameter
        Output: converted ASCII with correct offset
        
        Logic: Once the correct output parameter and qscore are determined, the correct offset is added back 
        converted back to ASCII characters.
        '''
        if self.PHRED64output:
            converted = chr(qScore + 64) #converts to P64
            return converted
        else: 
            #bypasses PHRED33output True by forcing it to print tihs if user wants P33 output
            # also works by forcing it to be default.
            converted = chr(qScore + 33)
            return converted

       
        
    def convert(self):
        '''
        Input: input/ output parameter, qualitySequence
        Output: a full converted quality sequence
        
        Logic: This determines if a the input and output pararmeters are the same. If they are, then 
        it would return the quality sequence as it was. Otherwise, it appends each ASCII character back
        into a longer string.
        '''
        
        qualityResultSequence = "" #setting up string
        replacedSequence = ""
        if (self.PHRED33input and self.PHRED33output) or (self.PHRED64input and self.PHRED64output):
            #if the commandline parameters are the same, then it returns the sequence as it was
            return self.qualitySequence
    
        for letter in self.qualitySequence:
            qualityResultSequence += self.convertPHRED(self.qScore(letter)) #appends individual characters into a string
            
                
        return qualityResultSequence
                        
        
    def convertDNA(self):
        '''
        Input: input/ output parameter, qualitySequence, (DNA)sequence
        Output: a full converted DNA sequence or non-converted DNA sequence
        
        Logic: This determines if a the input pararmeters is P64Bin. If it is, then 
        it would return the DNA sequence with all characters replaced with 'N's where 
        there are 'B's in the qualitySequence. Otherwise if not P64Bin, it returns 
        the DNA sequence as it was.
        
        First, an empty string is created; this will contain the replaced DNA sequence.
        While the index is less than the length of the DNA sequence, it finds the character
        at that index(letter), for the DNA sequence. At that same index in the qualitySequence, if 
        the character at that index is a 'B' then the letter is set to 'N'. Each character is then
        appended to the (empty) string. 
        '''
        if self.PHRED64Bin:
            DNAsequence = "" #setting up an empty string
            sequenceLength = len(self.sequence) #length of the DNA sequence
            index = 0 

            while index < sequenceLength: 
                letter = self.sequence[index] #sets letter to be the character at that index
                if self.qualitySequence[index] == "B": #if that index character is a 'B'
                    letter = "N" #then that letter is changed to 'N'
                DNAsequence += letter #add each character to the string
                index += 1 #increment index
        else: #the DNA sequence should be the same for everything else
            DNAsequence = self.sequence
        return DNAsequence
                
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse. Taken from Lab5: findORFS. 
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        #self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False, help='Illumina 1.8/Sanger: PHRED+33 input format.')
        self.parser.add_argument('-P64in', '--PHRED64input', action = 'store', nargs='?', const=True, default=False, help='Illumina 1.3: PHRED+64 input format.')
        self.parser.add_argument('-P64Bin', '--PHRED64Bin', action = 'store', nargs='?', const=True, default=False, help='Illumina 1.3: PHRED+64 with B sequence input')
        self.parser.add_argument('-P64SOLin', '--PHRED64SOLin', action = 'store', nargs='?', const=True, default=False, help='Solexa: Solexa+64 input format.')
        self.parser.add_argument('-P33out', '--PHRED33output', action = 'store', nargs='?', const=True, default=False, help='PHRED33 output format.')
        self.parser.add_argument('-P64out', '--PHRED64output', action = 'store', nargs='?', const=True, default=False, help='PHRED64 output format.')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
            
def main(inCL = None):
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    
    fileReader = FastQreader(myCommandLine.args.inFile)#takes input from Command Line
    
    #setting up values
    header = ''
    sequence = ''
    qualityHeader = ''
    qualitySequence = ''
    
    #setting up some lists 
    headerList = list()
    sequenceList = list()
    qualityHeaderList = list()
    qualitySequenceList = list()
    convertedList = list()
    convertedSequenceList = list()
    
    for header, sequence, qualityHeader, qualitySequence in fileReader.readFasta():
        #parses through the headers and sequnces in the file
        headerList.append(header) #appends the '@' header
        
        sequenceList.append(sequence)#appends the DNA sequence header
        
        qualityHeaderList.append(qualityHeader) #appends the '+' header
        
        qualitySequenceList.append(qualitySequence) #appends te unconveted value sequence
        
        #calls ConvertPHRED class
        convertP = ConvertPHRED(sequence, qualitySequence, myCommandLine.args.PHRED33input, myCommandLine.args.PHRED64input,
                                             myCommandLine.args.PHRED64Bin, myCommandLine.args.PHRED64SOLin,
                                             myCommandLine.args.PHRED33output, myCommandLine.args.PHRED64output)
        
        convertedQualitySequence = (convertP.convert()) #appends each converted sequence
        convertedSequence = (convertP.convertDNA())
        #printing each thing
        print(header)
        print(convertedSequence)
        print(qualityHeader)
        print(convertedQualitySequence + "\n") #add a blank line too make it easier to see
    
main()

