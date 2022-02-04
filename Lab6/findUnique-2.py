#!/usr/bin/env python3
# Name: Christie Kwok
# Group Members: Priyanka Khaware (pkhaware)
# got help from: LSS(Tutor: Mohammad Abdulqader; worked with Natalie Siquerios, Aishwarya Basude) for findUnique class
#                and some of main()
#                MSI(Tutor: Stephen Hwang) helped with main()

'''
Receicves a FastA (extension .fa) file of tRNA headers and sequences from stdin, and returns 
a text file (or not, user's choice) from stdout containing the tRNA headers and their essential substrings. 

Example: Note: this is only a partial output because this is gonna be really long even for a single tRNA.
    Input: > tRNA | Glu | 7UC | Bos taurus | mitochondrial
    
    Output:
    tRNA|Lys|∃UU|Bostaurus|mitochondrial
    CACUAAGA"LCUAUAUAGCACPAACCU∃UU6AGUUAGAGAUUGAGAGCCAU"UACUCUCCUUGGUGACCA
    CACU
    .ACUA
    ...UAA
    ....AAG
    .......A"
    .........L
    ..........CUAU
    
Note: this program utilizes the FastAreader class to parse through the input. 

Another Note: this replaces the E-hat characters with S~ and A-hat characters for some reason. (E -> SA)
    So the Leu and Trp tRNAs, and some sequences have replacements. (when I tested this). 
    It appended the headers with the E-hat replaced. So maybe a Windows issue?
    My output is not 100% as a result. used bos-tRNA-7.fa file
'''

import sys
class FastAreader :
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
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

class findUnique:
    '''
    Determines the powerset, unique, powerset, and the essential powerset of a tRNA sequence. Returns the essential 
    substring of a tRNA sequence.
    
    Author: Christie Kwok
    Date: May 25, 2020
    Input: a FastA file that contains tRNA header and sequences.
    Initalized: a list of sequences from the FastA file, lists for the powersets, unique powerset, and essential set.
    Class Attributes: N/A
    Instance Attributes: self.sequence, self.powerSetList, self.uniqueSetList, self.essentials
    Methods: __init__(), cleanSequence(), powerSets(), uniqueSets(), essentialSets().
    '''
    def __init__(self, sequence):
        ''' Initializes some lists and the input.'''
        self.sequence = sequence 
            #initialized input of tRNA sequences
        self.powerSetList = list()
            #initalzed a list for powersets
        self.uniqueSetList =list()
            #initialized a list for the unique powersets
        self.essentials = list()
            #initialized a list for the essential substrings
    

        
    def powerSets(self):
        '''
        Input: the cleaned sequence list
        Output: a list containing the powersets of each tRNA
        
        Logic:The first for loop removes all of the extraneous characters, ".", "_", "-", and replaces them 
            with nothing. It then appeneds the cleaned sequence to a list.
            
            I then used a three for loops. the first for loop goes through the cleaned sequences list for each tRNA,
            and generates a new list(eachPowerSet) for each tRNA. The second for loop iterates through the 
            length of each tRNA; this is the start of each substrings of the tRNA(startBase). The third for loop 
            dictates where the end of each substring is(endBase). Then, substrings are generated for each start 
            and end and appended to the eachPowerSet list for the specific tRNA. The lists are converted 
            to sets, and appended to the powerSetList.
        '''
        cleanSequenceList = list()
        for sequence in self.sequence:
            cleanSequence = sequence.replace(".", "")
            cleanSequence = cleanSequence.replace("_","")
            cleanSequence = cleanSequence.replace("-", "")
            cleanSequence = cleanSequence.upper()
            cleanSequenceList.append(cleanSequence)
        
        
        for cleanSequence in cleanSequenceList:
            #iterates through self.sequence for each tRNA
            eachPowerSet = list()
                # makes a new list for each tRNA
            for startBase in range(0, len(cleanSequence)):
                # the start of eah substring
                for endBase in range(1, len(cleanSequence) - startBase + 1):
                    #the end of each substring (is the smallest size a substring can be)
                    substring = cleanSequence[startBase : startBase+endBase] 
                        #how far to go from base (base + end)
                    eachPowerSet.append(substring)
                        #add each substring to the specific tRNA list
            eachPowerSet = set(eachPowerSet)
                #converts each tRNA powerset to a set to be used in set operations, later on
            self.powerSetList.append(eachPowerSet)
                #adds each tRNA set to a list

        

    def uniqueSets(self):
        '''
        Input: tRNA powerset list of sets
        Output: a list of unique powersets for each tRNA
        
        Logic: I used a single for loop. It iterates through each of the sets in the powerset list. It then slices
            setOfInterest for a single powerset. All other powersets are unioned together. The difference of the
            single set and the unioned sets are taken, to get the unique powersets for that tRNA. It is then appened
            to a list (uniqueSetList).
        '''
        for powerSet in range(len(self.powerSetList)):
            #iterates for the items in self.powerSetList
            originalSet = self.powerSetList[powerSet] 
                #get the powerset for a single tRNA by slicing 
            others = self.powerSetList[:powerSet] + self.powerSetList[powerSet + 1:] 
                #all of the other sets, by slicing for all things before and after the original set
            unionSet = set().union(*others)
                #unions the other sets
            uniqueSet = originalSet.difference(unionSet)
                #subtract the original set from unioned sets 
            self.uniqueSetList.append(uniqueSet) 
                #append the uniqueSet to a list
            
            
    def essentialSets(self):
        '''
        Input: the list of unique sets for each tRNA
        Output: all of the essential sets for each tRNA
        
        Logic: I used two for loops and an if statement. The first for loop iterates through the each 
        of the uniqueSets in  uniqueSetList, to get a set. I then sort the set by length; longest 
        string to shortest(sortedSet) of the current set. The current set is sorted by length because 
        the larger sequences are more likely to have equivalent shorter substrings. The second for 
        loop iterates through sortedSet for each unique string. In the if statement, the first character 
        is removed and checked against the rest of the substrings. If an equivalent substring is found, 
        it is added to nonEssentials. The same is done for removal of the last character.These unique 
        substrings are considered non-essential, as a truncated version shows up within it. After all
        non-essential substrings are found, the difference of the current/original set and the 
        non-essential set is done.This removes all non-essential sets, and are left with essential 
        substrings. the essential substrings are appendedto a list (for each tRNA).
        
        '''
        essentialList = list()
            #a list that will contain essential substrings
        nonEssentials = set()
            #a list that will contain non-essential substrings
            
        for sets in self.uniqueSetList:
            #the current set
            sortedSet = sorted(sets, key = lambda entry: len(entry), reverse = True) 
            #sorts the current set from big to small, by length of the substring
        
            for item in sortedSet:
                #the current substring
                if item[1:] in sortedSet or item[:-1] in sortedSet:
                    #truncate the current substring
                    #checks if the truncated substring is essential or not
                    nonEssentials.add(item)
                    #not essential substrings added to this set
                
            differences = sets.difference(nonEssentials)
                # essentials = unique set - nonEssential set
            essentialList.append(differences)
                #add essentials to a list (for each tRNA)
        
          
        return essentialList
    
def main(inCL=None):
    '''Calls the findUnique class and its methods. This also sorts the output of the findUnique class and outputs it 
    in the format:
    tRNA header
    tRNA sequence
    .substrings
    
    Got help in LSS and MSI for this part as well.
    '''
    fileReader = FastAreader()
    
    header = ''
    sequence = ''
    
    #lists is used because it is ordered, so headers + sequences will match
    headerList = list()
        # a list that will contain tRNA headers
    sequenceList = list()
        # alist that will contain tRNA sequences
    
    
    for header, sequence in fileReader.readFasta():
        #parses through the headers and sequnces in the file
        headerList.append(header)
        sequenceList.append(sequence)
    
    
    #calling class and methods
    findPS = findUnique(sequenceList)
    findPS.powerSets()
    findPS.uniqueSets()
    essentialSet = findPS.essentialSets()
    
    
    # got help from MSI(tutor: Stephen Hwang)
    tupleList = list()
        # a list to hold tuples
    
    for item in range(0,len(essentialSet)):
        #for items in the range of
        header = headerList[item].replace(' ', '')
            #removes whitespace from header of each tRNA
        sequence = sequenceList[item].replace('_', '').replace('.', '').replace('-','')
            #the sequnces of each tRNA
            #remove extraneous characters
        essential = essentialSet[item]
            #essential substring set of each tRNA
        essential = list(essential)
            #turns the set into a list
        
        triplet = (header, sequence, essential)
            #tuple containing tRNA header, sequence, and its essential substrings
        tupleList.append(triplet)
            #add the tuple to a list
        
    
    for triplet in sorted(tupleList):
        #for each tuple, header sequence and essential substrings in respective tuple indexes
        header = triplet[0]
        sequence = triplet[1]
        essential = triplet[2]
    
        #essential = sorted(essential, key = lambda entry:len(entry))
        print(header) 
            #prints header for each tRNA
        print(sequence) 
            #prints sequences for each tRNA
        
        align = "."
            #for use in printing "."
        
        outputList = list()
            # a list that contains the substrings  with their "."
        
        for substring in essential:
            #find the amount of indexes to count for and replace with "."
            #for n in range(0, len(substring)):
            position = sequence.find(substring)
                #finds the index the substring first occurs at
            essentialOutput = str(align * position + substring)
                #multiplies the "." by how the index, and adds that many "."
                #concatonates the substring to the amount of "."
            outputList.append(essentialOutput)
                #adds it to the output list
        
        
        for item in sorted(outputList, key = lambda entry: len(entry)):
            #sorts the substrings for by length of each tRNA (shortest to longest)
            print(item)

        
        
            
    
    
    
if __name__ == "__main__":
    main() 