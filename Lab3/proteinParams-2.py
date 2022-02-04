#!/usr/bin/env python3
# Name: Christie Kwok (cdkwok)
# Group Members: None, got help from LSS in a couple parts(Mohammad Abdulqader)

'''
Reads a protein string from user input and returns the total amino acid count, molecular weight,
Molar extinction coefficient, Mass Extinction Coefficient, pI, and amino acid composition of 
the protein.

Example: 
 Input: VLSPADKTNVKAAW

 Output:
    Number of Amino Acids: 14
    Molecular Weight: 1499.7
    molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67
    Theoretical pI: 9.88
    Amino acid composition:
    A = 21.43%
    C = 0.00%
    D = 7.14%
    E = 0.00%
    F = 0.00%
    G = 0.00%
    H = 0.00%
    I = 0.00%
    K = 14.29%
    L = 7.14%
    M = 0.00%
    N = 7.14%
    P = 7.14%
    Q = 0.00%
    R = 0.00%
    S = 7.14%
    T = 7.14%
    V = 14.29%
    W = 7.14%
    Y = 0.00%
'''
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
    

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
            
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()