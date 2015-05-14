'''
Created on Oct 30, 2014

@author: Patrick
'''
from numpy import *

me = 'Patrick'
print 'Jason Sucks, A Poem by '+me
lines = [line.strip() for line in open('json.txt')]
print lines

slines = []

for i in lines:
    slines.append(i.split());
    
#print slines

for i in slines:
    if i[len(i)-1] == 'COMMENT':
        print 'IGNORED'
    elif i[len(i)-1] == 'TITLE':
        title = i[0]
    elif i[len(i)-1] == 'FIVE':
        print i[len(i)-1]
    elif i[len(i)-1] == 'SEVEN':
        print i[len(i)-1]
    elif i[len(i)-1] == 'AUTHOR':
        print i[0]
    elif i[len(i)-1] == 'HAIKU':
        if i[len(i)-2] == 'END':
            print"THE END"
            
b = array( [ (1.5,2,3), (4,5,6) ] )
print b


    