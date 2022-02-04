#!/usr/bin/env python3
# Name: Christie Kwok (cdkwok)
# Group Members: Priyanka Khawere

from sequenceAnalysis import NucParams, FastAreader 

def main ():
    '''Input: a file (.fa)
        Output: Sequence Length, GC content, Codon Usage
        
        Logic: I called the dictionaries from NucParams(), and set them to variable names; these are used
        throughout for calculations. 
        
            SequenceLength: The output is in Mb, and Mb is 10**6. The total nucleotide
                count was divided by Mb to get sequence length.
                
            GC Content: got the repspective counts from the nucComp dictionary, summed them, and divded it by 
                total nuc count. (multiplied by 100 for percentage)
                
            Codon Usage: I sorted rnaCodonTable by values, for correct alphabetical order of amino acids, while keeping
            the key:value pair. I then iterated through it. I used the key to get the codon count from codonComp, 
            and the value to get the animo acid count from aaComp, to calculate the frequency.
    '''
    
    
    myReader = FastAreader() # make sure to change this to use stdin
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    
        
    nucComp = myNuc.nucComposition()
        
    nucCount = myNuc.nucCount()
        
    codonComp = myNuc.codonComposition()
        
    aaComp = myNuc.aaComposition()
    

    # Sequence Length
    Mb = 10**6 #size of a megabase
    seqLength = (nucCount / Mb)
    
    print("sequence length = {:.2f}Mb\n".format(seqLength))
    
    
    # GC Content
    gCount = nucComp['G'] #G count
    cCount = nucComp['C'] #C count
    gcContent = (gCount + cCount) / nucCount * 100
        #total gc content in percentage
        
    print("GC content = {:.1f}%\n".format(gcContent))
     
        
    # sort codons in alpha order, by Amino Acid
    
    aaSorted = {k: v for k, v in sorted(myNuc.rnaCodonTable.items(), key=lambda item: item[1])} 
        #this sorts the rnaCodonTable by values, while also mainintaing the correct key:value pair
        #cite: https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
    for codon, aa in aaSorted.items():
        #iterating through the sorted dictionary

        codonCount = codonComp[codon]
        aaCount = aaComp[aa]
        codonFrequency = (codonCount / aaCount) * 100
            #calulating the codon frequency by taking the codon count and
            #amino acid count from their respective dictionaries
        
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, codonFrequency, codonComp[codon]))

if __name__ == "__main__":
    main()
    