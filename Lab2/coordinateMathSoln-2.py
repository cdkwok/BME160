#!/usr/bin/env python3 
# Name: Christie Kwok (cdkwok)
# Group Members: Priyanka Khawere

'''
This program takes a single line of input from the user and strips all it of certain characters. Then it splits the string based
on a charcter and puts the wanted indexes into  tuple(s). The tuple(s) are passed into the Triad class, and its methods are called
to calculate the wanted values.
'''

import math
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /   math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r)))
    
        
        
def getTuple(coords):
    '''Takes in input, converts it into a tuple, calls Triad class methods, and return calculated numbers.'''
    
    #got help in LSS for the replacements below
    coords = coords.replace("C", "")
    coords = coords.replace("a", "")
    coords = coords.replace("N", "")
    coords = coords.replace(" = ", "")
    coords = coords. replace("(", "")
    coords2 = coords.replace(")", ",") 
    coords3 = coords2.split(",") #splits the string into a list
        
    p = float(coords3[0]), float(coords3[1]), float(coords3[2]) #makes tuple p
    q = float(coords3[3]), float(coords3[4]), float(coords3[5]) #makes tuple q
    r = float(coords3[6]), float(coords3[7]), float(coords3[8]) #makes tuple r
    tupleList = [p,q,r] #putting the tuple into a list
    triadCoords = Triad(p = tupleList[0], q = tupleList[1], r = tupleList[2]) #setting each parameter to the tuple
    ncBondLength = triadCoords.dPQ() 
    nCaBondLength = triadCoords.dQR()
    cNCaBondAngle = math.degrees(triadCoords.angleQ())
    
    
    
    return [ncBondLength, nCaBondLength, cNCaBondAngle]
    
def main():
    ''' Takes input and returns a printe statement'''
    
    coords = input("What are the coordinates? ")
    
    process = getTuple(coords)
    
    print("N-C bond length = {0:0.2f}".format(process[0]) 
         + "\nN-Ca bond length = {0:0.2f}".format(process[1]) 
         + "\nC-N-Ca bond angle = {0:0.1f}".format(process[2]))
     
    

main()