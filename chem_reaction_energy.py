#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 00:26:05 2021

@author: lyshen
"""

from pymatgen.ext.matproj import MPRester
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.units import FloatWithUnit
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.io.vasp import Vasprun
from pymatgen.core.periodic_table import Element

#This initializes the REST adaptor. Put your own API key in if necessary.
a = MPRester('tLf90ZQdzmWJelG4')
compatibility = MaterialsProjectCompatibility()
#%%
#This gets all entries belonging to the system.
vasprun_Li2OHCl = Vasprun("vasprun_pd.xml") 
vasp_entry_Li2OHCl = vasprun_Li2OHCl.get_computed_entry(inc_structure=True) 
Li2OHCl_entry = compatibility.process_entry(vasp_entry_Li2OHCl)
#%%
mp_entries = a.get_entries_in_chemsys(['Li', 'Cl', 'O', 'H'])
#mp_entries = a.get_entries_in_chemsys(['Li', 'Cl', 'O'])

#%%
all_entries = compatibility.process_entries([Li2OHCl_entry] + mp_entries)
#%%
#This method simply gets the lowest energy entry for all entry with the same composition.
def get_most_stable_entry(formula):
    relevant_entries = [entry for entry in all_entries if entry.composition.reduced_formula == Composition(formula).reduced_formula]
    relevant_entries = sorted(relevant_entries, key=lambda e: e.energy_per_atom)
    return relevant_entries[0]
#def get_decom_energy(reactants,products)
#%% Precursors | Reactants
LiOH = get_most_stable_entry("LiOH")
LiCl = get_most_stable_entry("LiCl")
Li2O = get_most_stable_entry("Li2O")
H2O = get_most_stable_entry("H2O")
#%% (all possible) products

#Li2OHCl = get_most_stable_entry("Li2OHCl")
Li3OCl = get_most_stable_entry("Li3OCl")
#%% Reaction calculation
#reaction = ComputedReaction([LiCl,LiOH], [Li2OHCl])  #if I just considered Li2OHCl, less stable and stoichometry is not balanced.
reaction = ComputedReaction([LiCl, LiOH], [Li3OCl, H2O])  #if I just considered Li2OHCl, less stable and stoichometry is not balanced.
                   
energy = FloatWithUnit(reaction.calculated_reaction_energy, "eV atom^-1")

print("Caculated")
print(reaction)
print("Reaction energy = {}".format(energy.to("kJ mol^-1")))
print(energy)