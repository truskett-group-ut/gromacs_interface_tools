import re

def GromacsTime():
    f = open('./grompp.mdp', 'r')
    gromacs_data = f.read()
    
    #initial check on zero time
    search = re.search(r'tinit\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)', gromacs_data)
    if search.group(1):
        tinit = float(search.group(1))
        if tinit > 0.00000001:
            raise Exception('"tinit" in grompp file is non-zero.')
        
    #read in the time data
    try:
        dt = float(re.search(r'dt\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)', gromacs_data).group(1))
        nsteps = float(re.search(r'nsteps\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)', gromacs_data).group(1))
    except:
        raise Exception('Could not extract "dt" and/or "nsteps" from the grompp file.')
        
    f.close()
        
    return (dt*nsteps)