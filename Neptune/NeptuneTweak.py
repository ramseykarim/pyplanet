def modify(gas,cloud,C,Cl):

    comment = "'Standard' Neptune atmospheric tweaking from Imke's code..."
    
    sol = 0.0
    zSol = 0.0
    nAtm = len(gas[C['P']])
    for i in range(nAtm):
        Plyr = gas[C['P']][i]
        Tlyr = gas[C['T']][i]
        
        ### Process H2S
        depleteh2s =  True
        option='A'  #the NH4SH "normal" option
        constAmth2s = 4e-9
        if depleteh2s:
            if option=='A':
                frach2s = [1.0,0.05]
                Pfrh2s = [43.0, 25.0]
                if Plyr > Pfrh2s[0]:
                    fr = frach2s[0]
                elif Plyr < Pfrh2s[1]:
                    fr = frach2s[1]
                    comment+='\nReduce H2S by %.2f%%' % (100.0*(1.0 - fr))
                else:
                    m = (frach2s[1]-frach2s[0]) / (Pfrh2s[1]-Pfrh2s[0])
                    fr = frach2s[0] + m*(Plyr - Pfrh2s[0])
                    comment+='\nReduce H2S by %.2f%%' % (100.0*(1.0 - fr))
                #print 'depleting H2S', Plyr, fr
                gas[C['H2S']][i] = fr*gas[C['H2S']][i]
                if gas[C['H2S']][i] < constAmth2s:
                    gas[C['H2S']][i] = constAmth2s
            if option =='B':
                if Plyr < 200.0:
                    comment+='\nReduce H2S to .001 at %f'%(Plyr)
                    gas[C['H2S']][i]*=0.001
                    comment+='===>>>and NH3!!!'
                    gas[C['NH3']][i]*=0.001

        ### Process NH3 
        enhancenh3 = True
        constAmtnh3 = 12e-9
        #if Tlyr > 400.0:
        #    gas[C['NH3']][i] = 1.94E-4
        if depleteh2s and option=='B' and Plyr < 200.0:
            comment+='\nNH3 was modified above under H2S'
        if enhancenh3:
            fracnh3 = [1.0,1.3]
            Pfrnh3 = [43.0, 25.0]
            if Plyr > Pfrnh3[0]:
                fr = fracnh3[0]
            elif Plyr < Pfrnh3[1]:
                fr = fracnh3[1]
                comment+='\nEnhance NH3 by %.2f%%' % (100.0*(fr - 1.0))
            else:
                m = (fracnh3[1]-fracnh3[0]) / (Pfrnh3[1]-Pfrnh3[0])
                fr = fracnh3[0] + m*(Plyr - Pfrnh3[0])
                comment+='\nReduce NH3 by %.2f%%' % (100.0*(fr - 1.0))
            #print 'enhancing NH3', Plyr, fr
            gas[C['NH3']][i] = fr*gas[C['NH3']][i]
            if gas[C['NH3']][i] < constAmtnh3:
                comment+='\nEnhance NH3 to '+str(constAmtnh3)
                #print 'Enhance NH3 to ',constAmtnh3
                gas[C['NH3']][i] = constAmtnh3

        ### Process CO
        if gas[C['P']][i] > 0.1585:
            gas[C['CO']][i] = 0.0
        else:
            gas[C['CO']][i] = 1.0E-6
        ### Process CO13
        gas[C['CO13']][i] = 1.0E-2*gas[C['CO']][i]
        ### Process HCN, SOLN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['SOLN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

        ### compute dz and integrated solution cloud and solution cloud height
        if i==0:
            gas[C['DZ']][i] = 0.0
        else:
            gas[C['DZ']][i] = abs(gas[C['Z']][i]-gas[C['Z']][i-1])*1.0E5
            sol+=gas[C['SOLN']][i]*gas[C['DZ']][i]
            if gas[C['SOLN']][i] > 0.0:
                zSol+=gas[C['DZ']][i]

    return comment, gas, cloud
