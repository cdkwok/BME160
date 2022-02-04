#!/usr/bin/env python3 
# Name: Christie Kwok (cdkwok)
# Group Members: none

'''
Input: A RNA or DNA codon/ 3-letter amino acid code/ single-letter amino acid code in the form of a string.
Ouput: The complement of the codon or amino acid code.

    Example: 
        Input: CUU
        Output: Leu
        
I used an if-elif statement to check for the key:value pair in the dictionaries. I checked to see if the input was in both the 
key and values, and printed the corresponding value . I wasn't sure if the tables were supposed to be in my class, so I left it
outside of it. 
'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

RNA_codon_table = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C 
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in RNA_codon_table.items()}

class convert(str):
    '''Takes the input(code) from user and returns the complement of the input.'''
    
    def __init__(self, str):
        self.code = str
        
    def converter(self):
        '''Returns the single letter or three-letter amino acid code or the amino acid based on what the input is'''
        
        want = self.code.upper() #converts the input to uppercase
        default = short_AA.setdefault("unknown","unknown") #creates a default value in short_AA

        if want in short_AA.keys(): #checks if the wanted value is in the dictionary short_AA keys
            return print (want + " = " + short_AA.get(want)) #gets the value and prints it

        elif want in short_AA.values(): #checks if the value is in short_AA values
            return print (want + " = " + long_AA.get(want)) #gets the value from the converted long_AA dictionary

        elif want in RNA_codon_table.keys(): #checks if the input corresponds to an RNA_codon_table key
            wanted = RNA_codon_table.get(want) #sets the value to a wanted(variable)
            return print(want + " = " + wanted.upper()) #capitalizes and prints the value

        elif want in dnaCodonTable.keys(): #checks if the input corresponds to a dnaCodonTable key
            wanted2 = dnaCodonTable.get(want) #sets the value to wanted2
            return print(want + " = " + wanted2.upper()) #capitalizes and prints the value

        else:
            return print(short_AA.get("unknown")) #if input is not in the tables, then it prints 

            
def main():
    ''' Takes in user input and passes it to the class. calls the method from the class. '''
   
    code = input("What is the code? ") #user input
    
    converted = convert(code)
    convertedCode = converted.converter()
        
    

main()