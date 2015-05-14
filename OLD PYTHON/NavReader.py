'''
NavReader, Reads in and stores the Ephemeris of GPS constellation from the RINEX 
navigation file.
NavReader() is the main function where the whole navigation file is being read. 
The sub functions HTYPE() reads the current stream of string line form the RINEX 
file and determines what kind of header line it is. 
Based on the Headertype From HTYPE NAV_ASSIGN() stores the header information based on
the type of header information it is. 
'''
from pprint import pprint
def NavReader(nav_file):
    head = True
    nav_header = {}
    nav_data = {}
    count = 1
    with open(nav_file) as f:
        
        for line in f:
            
            if head:
                line = line
                lines=(line[0:60],line[60:])
                HT= HTYPE((lines[1]))
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
                    EPOCH['SV_CLB'] = line[22:41].strip()
                    EPOCH['SV_CLD'] = line[41:60].strip()
                    EPOCH['SV_CLR'] = line[60:80].strip() 
                    nav_data[prn][epoch] = EPOCH
                    count = 2
                    #nav_data,count,prn,epoch = NAV_DATA_head(nav_data,line)
                elif pn == emp and count == 2:
                    LINE2 = {}
                    LINE2['IODE'] = line[0:22].strip()
                    LINE2['Crs'] = line[22:41].strip()
                    LINE2['Deln'] = line[41:60].strip()
                    LINE2['Mo'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE2)
                    count = 3
                elif pn == emp and count == 3:
                    LINE3 = {}
                    LINE3['Cuc'] = line[0:22].strip()
                    LINE3['Ecc'] = line[22:41].strip()
                    LINE3['Cus'] = line[41:60].strip()
                    LINE3['SqrtA'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE3)
                    count = 4
                elif pn == emp and count == 4:
                    LINE4 = {}
                    LINE4['TOE'] = line[0:22].strip()
                    LINE4['Cic'] = line[22:41].strip()
                    LINE4['OMEGA'] = line[41:60].strip()
                    LINE4['CIS'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE4)                    
                    count = 5
                elif pn == emp and count == 5:
                    LINE5 = {}
                    LINE5['Io'] = line[0:22].strip()
                    LINE5['Crc'] = line[22:41].strip()
                    LINE5['Omega'] = line[41:60].strip()
                    LINE5['OMEGA_DOT'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE5)                    
                    count = 6
                elif pn == emp and count == 6:
                    LINE6 = {}
                    LINE6['IDOT'] = line[0:22].strip()
                    LINE6['L2_CC'] = line[22:41].strip()
                    LINE6['GPS_W'] = line[41:60].strip()
                    LINE6['L2_P'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE6)                    
                    count = 7
                elif pn == emp and count == 7:
                    LINE7 = {}
                    LINE7['SV_Acc'] = line[0:22].strip()
                    LINE7['SV_Health'] = line[22:41].strip()
                    LINE7['TGD'] = line[41:60].strip()
                    LINE7['IODC'] = line[60:80].strip()
                    nav_data[prn][epoch].update(LINE7)                    
                    count = 8
                elif pn == emp and count == 8:
                    LINE8 = {}
                    LINE8['Trans_time'] = line[0:22].strip()
                    LINE8['Fit_int'] = line[22:41].strip()
                    nav_data[prn][epoch].update(LINE8)                    
                    count = 1
                else:
                    break
                #print count
    f.close()
    fielddict_file = open("nav_data.txt","w")
    pprint(nav_data, fielddict_file)
    fielddict_file.close()        
    return nav_data

def HTYPE(ss): 

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
