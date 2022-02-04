#!/usr/bin/env python3
# Name: Christie Kwok(cdkwok)
# Group Members: Priyanka Khaware (pkhaware), 
#I got help from LSS (Tutor: Mohammad Abdulqader; worked with Natalie Siquerios, Aishwarya Basude) for ORFfinder
# check extra credit

########################################################################
# CommandLine
########################################################################
'''
Reads a nucleotide string from file input(FASTA .fa format) and outputs a file containing the header of the FASTA file,
frame, start codon position, stop codon position, and length of the (longest) gene. 

This program can take command line arguments to adjust for the output. Can change the start codons, stop codons, 
minimun gene length, and/or all genes.
the protein.

Example: 
 Input: tass2.fa

 Output: (this is shortened from the total output)
    tass2 NODE_159_length_75728_cov_97.549133
    +1 57166..61908  4743
    -1  8192..11422  3231
    +2 65963..69004  3042
    -3 14589..16862  2274
    -2  2968.. 4872  1905
    +1 64093..65952  1860
    -3    30.. 1694  1665
    +1 69475..71052  1578
    +1 48805..50223  1419
'''
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
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
        


########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   

def main(inCL=None):
    '''
    Find some genes.  
    '''
   
    
    
    if inCL is None:
        myCommandLine = CommandLine()
    else :
        myCommandLine = CommandLine(inCL)
    
###### replace the code between comments.
    import sequenceAnalysis
    
    ORFreader = sequenceAnalysis.FastAreader(myCommandLine.args.inFile)
        #FastAreader class to read input
    
    outFile = open(myCommandLine.args.outFile, 'w')
        #setting up write to outFile
        #cite: https://docs.python.org/2/tutorial/inputoutput.html
        
    header = ' '
    sequence = ' '
    
    
    for header, sequence in ORFreader.readFasta():
        findORF = sequenceAnalysis.ORFfinder(sequence, myCommandLine.args.start, myCommandLine.args.stop,
                                             myCommandLine.args.minGene, myCommandLine.args.longestGene)
            #findORF class with parameters
        orfData = findORF.orf()
            #calling .orf()
        orfData.sort(key = lambda orf: (orf[3]), reverse = True)
            #sorting .orf() data by length in descending order while sorting start position
            # and stop position in ascending order
        
        outFile.write(header + "\n")
         #write header to outFile
            
        for data in orfData:
            if data[0] > 0:
                #if the frame is positive
                outFile.write('{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(data[0], data[1], data[2], data[3]))
                    #write results to outFile
            elif data[0] < 0:
                #if the frame is negative
                 outFile.write('{:d} {:>5d}..{:>5d} {:>5d}\n'.format(data[0], data[1], data[2], data[3]))
    print(myCommandLine.args.outFile)
    print
    
    #myCommandLine.args.inFile #has the input file name
    #myCommandLine.args.outFile # has the output file name
    #myCommandLine.args.longestGene # is True if only the longest Gene is desired
    #myCommandLine.args.start # is a list of start codons
    #myCommandLine.args.minGene # is the minimum Gene length to include
                            
#######
    
if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN
