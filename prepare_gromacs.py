import json
import sys
import re

with open('../settings.xml', 'r') as outfile:  
    path = re.search(r'<code_directory>(.+)</code_directory>', outfile.read()).group(1)

sys.path.insert(0, path)

from xml_extractor import xml_extractor 
from potentials import potential_combiner
from gromacs_interface_tools import gromacs_potential_maker

def MakeBool(text):
    text = text.lower()
    if re.match(r'\s*true\s*', text):
        return True
    elif re.match(r'\s*false\s*', text):
        return False

###########################################################################################

#read in the initialization data
with open('../potential_specs__params_state.json') as data_file:    
    potential_specs__params_state = json.load(data_file)
    potential_specs = potential_specs__params_state['specs']
    params_state = potential_specs__params_state['state']
    
#create the potential
potential = potential_combiner.Potential(potential_specs)
potential.SetParamsState(params_state)

#read in potential parameter values
with open('./params_val.json') as data_file:    
    params_val = json.load(data_file)

#load in xml settings
xml_extractor = xml_extractor.XMLExtractor()
xml_extractor.Parse('../settings.xml')
xmlParser = lambda x1, x2: xml_extractor.GetText('simulation', 'gromacs', x1, x2)

#generate gromacs data
gromacs_potential_maker = gromacs_potential_maker.SimulationPotentialConverter()
gromacs_potential_maker.LoadPotential(potential, params_val)
gromacs_potential_maker.TabulatePotential(r_max=float(xmlParser('table', 'r_max')), 
                                          dr=float(xmlParser('table', 'dr'))
                                         )
gromacs_potential_maker.CutShiftTabulated(e_max=float(xmlParser('potential', 'e_max')), 
                                          f_max=float(xmlParser('potential', 'f_max')), 
                                          r_max=float(xmlParser('potential', 'r_max')), 
                                          shift=MakeBool(xmlParser('potential', 'shift'))
                                         )
gromacs_potential_maker.MakeTable(filename='./table.xvg')
gromacs_potential_maker.InsertGromppCutoff(r_buffer=float(xmlParser('grompp', 'r_buffer')), 
                                           filename='./grompp.mdp')