'''
Created on Oct 30, 2014
Reads in Observation and Navigation Receiver Independent Exchange Format 
(RINEX) files into dictionaries to be used in other functions.

@author: Patrick And Sonal
'''
import math
from pprint import pprint

'''

'''   
def RINEXREADER(file):
    head = True
    header = {}
    obs = {}
    obs['LIST']=[]
    nl = 0
    nlr = 0
    obshead = True
    sw = True
    epoch = 's'
    it2 = 0
    
    
    with open(file) as f:
        
        for line in f:
            if not ((line[0:36].strip()== 'other post-header comments skipped') or (line[28:34].strip()== '4  1')):
                if head:
                    lines=(line[0:60],line[60:])
    
                    #print lines
                    HT= HTYPER((lines[1]))
                    header = ASSIGNDIC(header,lines[0],HT)
                    if HT == 20:
                        head = False
                else:
                    if  not ((line[60:].strip() == 'COMMENT') or (line[28:34].strip()== '4 18')):
                        if obshead: #Observation Header
                            if sw:        
                                epoch= line[0:3].strip()+':'+line[3:6].strip()+':'+line[6:9].strip()+':'+line[9:12].strip()+':'+line[12:15].strip()+':'+line[15:26].strip()
                                # print epoch
                                obs['LIST'].append(epoch)
                                sw = False
                                obs[epoch]={}
                                obs[epoch]['EFLAG']=line[26:29].strip()
                                obs[epoch]['NUMSAT']=int(line[29:32])
                                nl = obs[epoch]['NUMSAT']
                                #nl = obs[epoch]['NUMSAT']*math.ceil(header['OBSTYP']['NUM']/5.0)+math.ceil(obs[epoch]['NUMSAT']/12.0)
                                S=line[68:80].strip() 
                                if  bool(S):
                                    obs[epoch]['OFFSET']=float(S)
                                else:
                                    obs[epoch]['OFFSET']=float('nan')
                                    
                                obs[epoch]['LIST']=[]
                                
                                if (obs[epoch]['NUMSAT']>=12):
                                    nlr = math.ceil(obs[epoch]['NUMSAT']/12.0)-1
                                    
                                    for x in xrange(0, 12):
                                        S=line[(32+ x*3):(35+x*3)]
                                        S=S.replace(' ','0')
                                        obs[epoch]['LIST'].append(S)
                                        obs[epoch][S] = {}
                                else:
                                    for x in xrange(0, obs[epoch]['NUMSAT']):
                                        S=line[(32+ x*3):(35+x*3)]
                                        S=S.replace(' ','0')
                                        obs[epoch]['LIST'].append(S)
                                        obs[epoch][S] = {} 
                                    obshead = False
                                    sw = True                                                                                                                                                                                                                                                                                                                                                                                                                               
                            else:
                                if (nlr == 1):
                                    for x in xrange(0, obs[epoch]['NUMSAT']%12):
                                        S=line[(32+ x*3):(35+x*3)]
                                        S=S.replace(' ','0')
                                        obs[epoch]['LIST'].append(S)
                                        obs[epoch][S] = {} 
                                    obshead = False
                                    sw = True 
                                else:
                                    nlr = nlr -1
                                    for x in xrange(0, 12):
                                        S=line[(32+ x*3):(35+x*3)]
                                        S=S.replace(' ','0')
                                        obs[epoch]['LIST'].append(S)
                                        obs[epoch][S] = {}                
                        else: 
                            #Observations
                            
                            if sw:
                                it2 = 0
                                sw = False
                                sat = obs[epoch]['LIST'][obs[epoch]['NUMSAT']-nl]
                                if (header['OBSTYP']['NUM']>=5):
                                    nlr = math.ceil(header['OBSTYP']['NUM']/5.0) - 1
                                    it = 5
                                else:
                                    it = header['OBSTYP']['NUM']
                                    sw = True
                                    nl = nl -1
                                    
                                obs[epoch][sat]={}
                                
                                for x in xrange(0, it):
                                    S=line[(0+x*16):(14+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]]=float(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]]= float('nan')
                                    S=line[(14+x*16):(15+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]+'LL']= int(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]+'LL']= float('nan')
                                    S=line[(15+x*16):(16+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]+'STR']= int(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x]+'STR']= float('nan')
                            else:
                                sat = obs[epoch]['LIST'][obs[epoch]['NUMSAT']-nl]
                                it2 = it2 +1
                                if (nlr == 1):
                                    it = header['OBSTYP']['NUM']%5
                                    sw = True
                                    nl = nl -1
                                else:    
                                    nlr = nlr -1
                                    it = 5
                                 
                                for x in xrange(0, it):
                                    S=line[(0+x*16):(14+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]]=float(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]]= float('nan')
                                    S=line[(14+x*16):(15+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]+'LL']= int(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]+'LL']= float('nan')
                                    S=line[(15+x*16):(16+x*16)].strip()
                                    if S:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]+'STR']= int(S)
                                    else:
                                        obs[epoch][sat][header['OBSTYP']['OBS'][x+5*it2]+'STR']= float('nan')   
                            
                            if nl == 0:
                                obshead = True
    f.close()
    '''fielddict_file = open("obs_data.txt","w")
    pprint(obs, fielddict_file)
    fielddict_file.close() 
    fielddict_file = open("obsHead_data.txt","w")
    pprint(header, fielddict_file)
    fielddict_file.close() '''
    #print obs['LIST']
    r = {'HEAD':header, 'OBS': obs}
    return r
'''
    fielddict_file = open("obs.txt","w")
    pprint(obs, fielddict_file)
    fielddict_file.close()
'''
def HTYPER(ss): 
    ss = ss.strip().split()   
    if ss[len(ss)-1] == 'TYPE':
        if ss[len(ss)-3] == 'VERSION':
            return 1
        else:
            return 8
    elif ss[len(ss)-1] == 'COMMENT':
        return 2
    elif ss[len(ss)-1] == 'DATE':
        return 3
    elif ss[len(ss)-1] == 'NAME':
        return 4
    elif ss[len(ss)-1] == 'NUMBER':
        return 5
    elif ss[len(ss)-1] == 'AGENCY':
        return 6
    elif ss[len(ss)-1] == 'VERS':
        return 7
    elif ss[len(ss)-1] == 'XYZ':
        return 9
    elif ss[len(ss)-1] == 'H/E/N':
        return 10
    elif ss[len(ss)-1] == 'L1/2':
        return 11
    elif ss[len(ss)-1] == 'OBSERV':
        return 12
    elif ss[len(ss)-1] == 'INTERVAL':
        return 13
    elif ss[len(ss)-1] == 'OBS':
        if ss[len(ss)-2] == 'FIRST':
            return 14
        elif ss[len(ss)-2] == 'LAST':
            return 15
        else:
            return 19
    elif ss[len(ss)-1] == 'APPL':
        return 16
    elif ss[len(ss)-1] == 'SECONDS':
        return 17
    elif ss[len(ss)-1] == 'SATELLITES':
        return 18
    elif ss[len(ss)-1] == 'HEADER':
        return 20
    else:
        print "ERROR, DOES NOT EXIST"

"""
This function ASSIGNDIC() takes in 3 values the main dictionary HE, the current string line being,
and the Header type HT being currently read.
"""    
def ASSIGNDIC(HE, S, HT):
    if HT == 1:
        RVDT = {}
        RVDT['VER'] = S[0:20].strip()
        RVDT['OBSTYP'] = S[20:40].strip()
        RVDT['SATSYS'] = S[40:60].strip()
        #print RVDT
        HE['RVDT']=RVDT
        return HE
    elif HT == 2:
        if not 'COMMENT' in HE.keys():
            HE['COMMENT']= S.strip()
        else: 
            t = HE['COMMENT']
            
            HE['COMMENT']=t+'\n'+S
        return HE
    elif HT == 3:
        PRBD = {}
        PRBD['PGEN']=S[0:20].strip()
        PRBD['RUNBY']=S[20:40].strip()
        PRBD['DATE']=S[40:60].strip()
        HE['PRBD'] = PRBD
        #print PRBD #get rid later
        return HE
    elif HT == 4:
        HE['MRKR'] = S.strip()
        return HE
    elif HT == 5:
        HE['MKNUM'] = S.strip()
        return HE
    elif HT == 6:
        DONEBY = {}
        DONEBY['OBSV'] = S[0:20].strip()
        DONEBY['AGEN'] = S[20:60].strip()
        HE['DONEBY'] = DONEBY
        #print DONEBY
        return HE    
    elif HT == 7:
        RECV = {}
        RECV['NUM'] = S[0:20].strip()
        RECV['TYP'] = S[20:40].strip()
        RECV['VERS'] = S[40:60].strip()
        #print RECV
        HE['RECV'] = RECV
        return HE
    elif HT == 8:
        ANT = {} 
        ANT['NUM'] = S[0:20].strip()
        ANT['TYP'] = S[20:40].strip()
        #print ANT
        HE['ANT'] = ANT   
        return HE
    elif HT == 9:
        POS = {}
        POS['X'] = float(S[0:15].strip())
        POS['Y'] = float(S[15:30].strip())
        POS['Z'] = float(S[30:45].strip())
        #print POS
        HE['POS'] = POS
        return HE
    elif HT == 10:
        ANTDEL = {}
        ANTDEL['HT'] = S[0:15].strip()
        ANTDEL['EAEC'] = S[15:30].strip()
        ANTDEL['NOEC'] = S[30:45].strip()
        #print ANTDEL
        HE['ANTDEL'] = ANTDEL
        return HE
    elif HT == 12:
        if not 'OBSTYP' in HE.keys():
            OBSTYP = {}
            OBSTYP['NUM']=int(S[0:6])            
            lt=[]
            
            if(OBSTYP['NUM']>=8):
                for x in xrange(0, 8):
                    lt.append(S[6+x*6:12+x*6].strip())
            else:
                for x in xrange(0, OBSTYP['NUM']):
                    lt.append(S[6+x*6:12+x*6].strip())
            OBSTYP['OBS']=lt
            HE['OBSTYP']=OBSTYP
        else: 
            t = HE['OBSTYP']
            for x in xrange(0, t['NUM']-8):
                    t['OBS'].append(S[6+x*6:12+x*6].strip())        
            HE['OBSTYP']=t
        return HE
    elif HT == 14:
        TFIRST = {}
        TFIRST['YEAR'] = S[0:6].strip()
        TFIRST['MON'] = S[6:12].strip()
        TFIRST['DAY'] = S[12:18].strip()
        TFIRST['HR'] = S[18:24].strip()
        TFIRST['MIN'] = S[24:30].strip()
        TFIRST['SEC'] = S[30:43].strip()
        TFIRST['TS'] = S[43:51].strip()
        # print TFIRST
        HE['TFIRST'] = TFIRST
        return HE   
    else:
        return HE

'''
NavReader, Reads in and stores the Ephemeris of GPS constellation from the RINEX 
navigation file.
NavReader() is the main function where the whole navigation file is being read. 
The sub functions HTYPE() reads the current stream of string line form the RINEX 
file and determines what kind of header line it is. 
Based on the Headertype From HTYPE NAV_ASSIGN() stores the header information based on
the type of header information it is. 
'''
def NAVREADER(nav_file):
    head = True
    nav_header = {}
    nav_data = {}
    count = 1
    with open(nav_file) as f:
        
        for line in f:
            
            if head:
                line = line
                lines=(line[0:60],line[60:])
                HT= HTYPEN((lines[1]))
                nav_header = NAV_ASSIGN(nav_header,lines[0],HT)
                if HT == 8:
                    head = False
            else:

                emp = ""
                pn = line[0:2].strip()
                if pn != emp and count == 1:
                    if line[0] == " ":
                        prn = 'G0'+line[1]
                    else:   
                        prn = 'G'+line[0:2]
                    if not prn in nav_data.keys():
                        nav_data[prn] = {}
                        
                    epoch = line[2:5].strip()+':'+line[5:8].strip()+':'+line[8:11].strip()+':'+line[11:14].strip()+':'+line[14:17].strip()+':'+line[17:22].strip()
                    EPOCH = {}        
                    EPOCH['SV_CLB'] = float(line[22:41].strip().replace('D','E'))
                    EPOCH['SV_CLD'] = float(line[41:60].strip().replace('D','E'))
                    EPOCH['SV_CLR'] = float(line[60:80].strip().replace('D','E')) 
                    nav_data[prn][epoch] = EPOCH
                    count = 2
                elif pn == emp and count == 2:
                    LINE2 = {}
                    LINE2['IODE'] = float(line[0:22].strip().replace('D','E'))
                    LINE2['Crs'] = float(line[22:41].strip().replace('D','E'))
                    LINE2['Deln'] = float(line[41:60].strip().replace('D','E'))
                    LINE2['Mo'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE2)
                    count = 3
                elif pn == emp and count == 3:
                    LINE3 = {}
                    LINE3['Cuc'] = float(line[0:22].strip().replace('D','E'))
                    LINE3['Ecc'] = float(line[22:41].strip().replace('D','E'))
                    LINE3['Cus'] = float(line[41:60].strip().replace('D','E'))
                    LINE3['SqrtA'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE3)
                    count = 4
                elif pn == emp and count == 4:
                    LINE4 = {}
                    LINE4['TOE'] = float(line[0:22].strip().replace('D','E'))
                    LINE4['Cic'] = float(line[22:41].strip().replace('D','E'))
                    LINE4['OMEGA'] = float(line[41:60].strip().replace('D','E'))
                    LINE4['Cis'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE4)                    
                    count = 5
                elif pn == emp and count == 5:
                    LINE5 = {}
                    LINE5['Io'] = float(line[0:22].strip().replace('D','E'))
                    LINE5['Crc'] = float(line[22:41].strip().replace('D','E'))
                    LINE5['Omega'] = float(line[41:60].strip().replace('D','E'))
                    LINE5['OMEGA_DOT'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE5)                    
                    count = 6
                elif pn == emp and count == 6:
                    LINE6 = {}
                    LINE6['IDOT'] = float(line[0:22].strip().replace('D','E'))
                    LINE6['L2_CC'] = float(line[22:41].strip().replace('D','E'))
                    LINE6['GPS_W'] = float(line[41:60].strip().replace('D','E'))
                    LINE6['L2_P'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE6)                    
                    count = 7
                elif pn == emp and count == 7:
                    LINE7 = {}
                    LINE7['SV_Acc'] = float(line[0:22].strip().replace('D','E'))
                    LINE7['SV_Health'] = float(line[22:41].strip().replace('D','E'))
                    LINE7['TGD'] = float(line[41:60].strip().replace('D','E'))
                    LINE7['IODC'] = float(line[60:80].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE7)                    
                    count = 8
                elif pn == emp and count == 8:
                    LINE8 = {}
                    LINE8['Trans_time'] = float(line[0:22].strip().replace('D','E'))
                    LINE8['Fit_int'] = float(line[22:41].strip().replace('D','E'))
                    nav_data[prn][epoch].update(LINE8)                    
                    count = 1
                else:
                    break
                #print count
    ''''f.close()
    fielddict_file = open("nav_data.txt","w")
    pprint(nav_data, fielddict_file)
    fielddict_file.close()        '''
    return nav_data

def HTYPEN(ss): 

    ss = ss.strip().split()   
    if ss[len(ss)-1] == 'TYPE':
            return 1
    elif ss[len(ss)-1] == 'DATE':
        return 2
    elif ss[len(ss)-1] == 'COMMENT':
        return 3
    elif ss[len(ss)-1] == 'ALPHA':
        return 4
    elif ss[len(ss)-1] == 'BETA':
        return 5
    elif ss[len(ss)-1] == 'W':
        return 6
    elif ss[len(ss)-1] == 'SECONDS':
        return 7
    elif ss[len(ss)-1] == 'HEADER':
        return 8
    else:
        print "ERROR: DOES NOT EXIST"    
def NAV_ASSIGN(HE, S, HT):
    if HT == 1:
        RVDT = {}
        RVDT['VERS'] = S[0:20].strip()
        RVDT['OBSTYP'] = S[20:40].strip()
        RVDT['SATSYS'] = S[40:60].strip()
        HE['RVDT'] = RVDT
        return HE
        #print RVDT
    elif HT == 2:
        PRBD = {}
        PRBD['PGEN'] = S[0:20].strip()
        PRBD['RUNBY'] = S[20:40].strip()
        PRBD['DATE'] = S[40:60].strip()
        #print PRBD
        HE['PRBD'] = PRBD
        return HE
    elif HT == 3:    
        if not 'COMMENT' in HE.keys():
            HE['COMMENT']= S.strip()
        else: 
            t = HE['COMMENT']
            HE['COMMENT']=t+'\n'+S
        return HE
    elif HT == 4:
        ION_ALPHA = {}
        ION_ALPHA['A0'] = S[0:15].strip()
        ION_ALPHA['A1'] = S[15:30].strip()
        ION_ALPHA['A2'] = S[30:45].strip()
        ION_ALPHA['A3'] = S[45:60].strip()
        HE['ION_ALPHA'] = ION_ALPHA
        return HE
    elif HT == 5:
        ION_BETA = {}
        ION_BETA['A0'] = S[0:15].strip()
        ION_BETA['A1'] = S[15:30].strip()
        ION_BETA['A2'] = S[30:45].strip()
        ION_BETA['A3'] = S[45:60].strip()
        HE['ION_BETA'] = ION_BETA
        return HE        
    elif HT == 6:
        UTC = {}
        UTC['A0'] = S[0:21].strip()
        UTC['A1'] = S[21:42].strip()
        UTC['T'] = S[42:51].strip()
        UTC['W'] = S[52:60].strip()
        HE['UTC'] = UTC
        return HE
    elif HT == 7:
        HE['LEAP'] = S[0:6].strip()
        return HE     
    else:
        return HE    

    
    
