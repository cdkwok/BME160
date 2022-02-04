#!/usr/bin/env python3 
# Name: Christie Kwok (cdkwok) 
# Group Members: Priyanka Khawere

'''
Recieves a FASTQ seqname from user input and returns a the fields of the input.

Example:
    Input: @EAS139:136:FC706VJ:2:2104:15343:197393 
    Output: Instrument = EAS139
            Run ID = 136
            Flow Cell ID = FC706VJ
            Flow Cell Lane = 2
            Tile Number = 2104
            X-coord = 15343
            Y-coord = 197393

I chose to check if the string was a FASTQ seqname, so I impletement an if statement to check if the first character is '@'. If it is not, 
it would return a statement saying that it was not a FASTQ seqname. If there was confirmation of a FASTQ seqname, then it would proceed with 
sepearating each field.
'''

class FastqString (str):
    ''' Reads a FASTQ seqname from input and returns each field from the input..'''
    def parse(self):
        ''' Returns a printed statement of the fields and their corresonding places ochecks if the input is a FASTQ seqname 
        by checking if the first letter is '@''. It replaces the @ with  nothing. It then splits the string at each ":". 
        Then it prints the various fields based on what index each field is in.'''
        
        if self[0] == "@": #checks if the string is FASTQ by checking if the first index is '@'
            withoutAt = self.replace("@", "")  #replaces '@' with nothing
            fields = withoutAt.split(':',6) #withoutAt split at each ':'

            print("Instrument = " + fields[0] + "\nRun ID = " + fields[1] + "\nFlow Cell ID = " + fields[2] +
                  "\nFlow Cell Lane " + fields[3] + "\nTile Number = " + fields[4] + "\nX-coord = " + fields[5] 
                  +"\nY-coord = " + fields[6])  
                #this prints each field with the required field retrieved from its concurrent index
        else:
            print("This is not a FASTQ file.")
def main():
    ''' Function docstring goes here.'''
    FASTQ = input("What is the seqname of the FASTQ file? ")
    seqname = FastqString(FASTQ)
    parsedSeqname = seqname.parse()
    

main()