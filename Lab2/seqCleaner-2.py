#!/usr/bin/env python3 
# Name: Christie Kwok (cdkwok)
# Group Members: none
'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example: 
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        '''returns length of the string'''
        return (length(self))
    
    
    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        
        upperInput = self.upper() #makes everything in self uppercase
        
        for letter in upperInput: #iterates through each index in upperInput to do the following
            countN = upperInput.count("N") #counts "N" in self
            replacementN = upperInput.replace("N", str({countN}), 1) #replaces the first found "N" with countN
            removeN = replacementN.replace("N", "") #replaces the other instances of "N" with nothing
        return removeN
    
def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)
    
main()