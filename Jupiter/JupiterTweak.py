def modify(gas,cloud,C,Cl):

    comment = "Jupiter tweaking"
    
    sol = 0.0
    zSol = 0.0
    nAtm = len(gas[C['P']])
    for i in range(nAtm):
        Plyr = gas[C['P']][i]
        Tlyr = gas[C['T']][i]
        
        ### Process H2S


        ### Process NH3 


        ### Process CO
        gas[C['CO']][i] = 0.0
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
