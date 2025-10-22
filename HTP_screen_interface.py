#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 16:43:47 2023

@author: lyshen
"""

'''from pymatgen.ext.matproj import MPRester

with MPRester("WBDGk31i0QdCZVmz") as m:

    # Structure for material id
    structure = m.get_structure_by_material_id("mp-1")
    mp_entries = m.get_entries_in_chemsys(["Li"])
    docs = m.summary.search(elements=["Si", "O"],  band_gap=(0.5, 1.0), fields=["material_id",  "band_gap",  "volume", "formula_pretty"])
'''    
#OLD API KEY: WBDGk31i0QdCZVmz
#NEW API KEY: 5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj
import numpy as np
from collections import OrderedDict
from pymatgen.io.vasp.sets import Vasprun
from pymatgen.core import Structure
#from m3gnet.models import Relaxer
from mp_api.client import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.core import Element, Composition
from pymatgen.core.units import FloatWithUnit
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.analysis.interface_reactions import GrandPotentialInterfacialReactivity
## PS: 与SSE的反应能计算必须避开Li单质，GrandPotentialInterfacialReactivity类中不能用开放元素作为反应物之一
dirname = 'D:\OneDrive - HKUST Connect\M3GNET\Li_stable_interface\Li6PS5Cl\\'
with open(dirname + 'data.csv', 'a') as f1:
    f1.write('mp_id' + '|' + 'formula' + '|' + 'E_hull' + '|' + 'E_form' + '|' + 'reaction_Li' + '|' + 'E_react_Li'  + '|' + 'reaction_LiSiPSCl' + '|' + 'E_react_LiSiPSCl' + '|' + '\n')#)

vasprun_Li3InCl6 = Vasprun("vasprun-LPSCl.xml") 
vasp_entry_Li3InCl6 = vasprun_Li3InCl6.get_computed_entry(inc_structure=True)
compatibility = MaterialsProjectCompatibility() 
Li3InCl6_entry = compatibility.process_entry(vasp_entry_Li3InCl6)
with MPRester("5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj") as rester:
    LIC_entries = rester.get_entries_in_chemsys(['Li', 'P', 'S', 'Cl'])

LIC_struc = Structure.from_file("LPSCl.cif")
LIC_comp = LIC_struc.composition.reduced_formula
#%% Download the structures that contain Li element in Materials Project and relax the structures
with MPRester("5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj") as mpr: 
    docs = mpr.summary.search(elements=["Li"], energy_above_hull=(0, 0.035), formation_energy=(-1000, 0), fields=["material_id", "volume", "formula_pretty", "energy_above_hull", "formation_energy_per_atom", "chemsys", "composition"])
    #docs = mpr.summary.search(elements=["Li"], energy_above_hull=(0, 0.035), fields=["material_id", "volume", "formula_pretty", "energy_above_hull", "formation_energy_per_atom", "chemsys", "composition"])
num = len(docs)

for n in range(num):
    example_doc = docs[n]
    mpid = example_doc.material_id # a Materials Project ID
    formula = example_doc.formula_pretty # a formula
    volume = example_doc.volume # a volume
    ehull = example_doc.energy_above_hull
    eformation = example_doc.formation_energy_per_atom
    chemical_system = example_doc.chemsys
    comp = example_doc.composition
    
    if formula == str('Li'):
        continue
    else:    
        with MPRester("5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj") as mpr: 
            structure = mpr.get_structure_by_material_id(str(example_doc.material_id))
        
        #example_doc.fields_not_requested # list of unrequested fields
        structure.to(str(example_doc.material_id) + '-' + formula + '.cif')
        '''
        print(mpid)
        print(formula)
        print('volume:', volume, "A^3")
        print('energy above the hull:', ehull, "eV/atom")
        print('formation energy:', eformation, "eV/atom")
        print('chemical system:', chemical_system)'''
        #print(structure)
        
        # prepare parameters for reaction energy calculation
        chem_sys = []
        for i in range(len(chemical_system.split('-'))):
            element = chemical_system.split('-')[i]
            chem_sys.append(element)
        
        with MPRester("5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj") as rester:
            mp_entries = rester.get_entries_in_chemsys(chem_sys)
        
        compatibility = MaterialsProjectCompatibility()
        entries = compatibility.process_entries(mp_entries)
        pd = PhaseDiagram(entries)
        
        entry = rester.get_entry_by_material_id(str(example_doc.material_id))
        #%% Electrochemical Stability with LITHIUM metal
        #identify a reference for lithium chemical potential using the bulk Li energy uli0
        li_entries = [e for e in entries if e.composition.reduced_formula == "Li"]
        uli0 = min(li_entries, key=lambda e: e.energy_per_atom).energy_per_atom
        
        #el_profile = pd.get_element_profile(Element("Li"), entry.composition)
        el_profile = pd.get_element_profile(Element("Li"), comp)
        for i, d in enumerate(el_profile):
            voltage = -(d["chempot"] - uli0)
            '''print("Voltage: %s V" % voltage)
            print(d["reaction"])
            print("")'''
            break
        
        def get_most_stable_entry(formula):
            relevant_entries = [entry for entry in entries if entry.composition.reduced_formula == Composition(formula).reduced_formula]
            relevant_entries = sorted(relevant_entries, key=lambda e: e.energy_per_atom)
            return relevant_entries[0]
        
        
        Li = get_most_stable_entry("Li")
        reactant = get_most_stable_entry(formula)
        
        product = []
        prods = []
        #print('products:')
        for j in range(len(d["reaction"].products)):
            product.append(d["reaction"].products[j].reduced_formula)
            locals()['product'+str(j)] = get_most_stable_entry(product[j])
            prods.append(locals()['product'+str(j)])
        
        
        # Reaction calculation
        reaction = ComputedReaction([Li, reactant], prods)#[locals()['product'+str(k)] for k in range(len(d["reaction"].products))])
        energy = FloatWithUnit(reaction.calculated_reaction_energy, "eV atom^-1")
        if energy >= -0.2:#若不与Li反应
            with open(dirname + 'data.csv', 'a') as f1:
                f1.write(str(mpid) + '|' + str(formula) + '|' + str(ehull) + '|' + str(eformation) + '|' + str(reaction) + '|' + str(energy) + '|')# + 'reaction_LIC' + '|' + 'E_react_LIC')
            print("Caculated")
            print(reaction)
            print("Reaction energy = {}".format(energy.to("kJ mol^-1")))
            print(energy)
            
            #%%计算与LIC反应能
            mpr = MPRester('5we9cnJa3XVOHuXzwqgXqVbZAdxCXlLj')

            reactant1 = 'Li6PS5Cl'
            reactant2 = formula
            grand = True
            if grand:
               open_el = 'P' 
               relative_mu = -1 #不要轻易改这个数值，反应会出偏差

            comp1 = Composition(reactant1)
            comp2 = Composition(reactant2)
            elements = [e.symbol for e in comp1.elements + comp2.elements]
            if grand:
                elements.append(open_el)
            elements = list(set(elements)) 

            entries = mpr.get_entries_in_chemsys(elements)
            pd = PhaseDiagram(entries)


            if grand:
                mu = pd.get_transition_chempots(Element(open_el))[0]
                #mu = pd.el_refs[Element("Li")]
                chempots = {open_el: relative_mu + mu}
                gpd = GrandPotentialPhaseDiagram(entries, chempots)
                interface = GrandPotentialInterfacialReactivity(
                    comp1, comp2, gpd, norm=True, include_no_mixing_energy=True, pd_non_grand=pd, use_hull_energy=False)
            else:
                interface = GrandPotentialInterfacialReactivity(
                    comp1, comp2, pd, norm=True, include_no_mixing_energy=False, pd_non_grand=None, use_hull_energy=False)

            plt = interface.plot()

            reaction_energy = np.zeros(100)
            for _, ratio, reactivity, rxn, rxn_energy in interface.get_kinks():
                critical_rxns = [
                    OrderedDict([
                                ("Atomic fraction", round(ratio, 3)),
                                ("Reaction equation", rxn),
                                ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 1)),
                                ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 3)),
                                ])
                                ]
                print(_, rxn, reactivity)
                reaction_energy[_] = str(reactivity)
                with open(dirname + 'data.csv', 'a') as f1:
                    f1.write(str(_) +'  ' + str(rxn) +'   ' + 'E=' + str(reactivity) + 'eV/atom;')
            
            with open(dirname + 'data.csv', 'a') as f1:
                f1.write('|' + str(reaction_energy[2]) + '|' + '\n')
            
            #calculate the Wulff shape and obtain the most exposed surface for the interfacial formation energy calculation
            
            


#%%relax Li3InCl6
'''
LIC = Structure.from_file('Li3InCl6_mp-676109_computed.cif')
atom_num = len(LIC.atomic_numbers)

relaxer = Relaxer()  # This loads the default pre-trained model
relax_results = relaxer.relax(LIC)
LIC_structure = relax_results['final_structure']
LIC_energy = relax_results['trajectory'].energies[-1]# / atom_num

print(f"Relaxed lattice parameter is {LIC_structure.lattice.abc[0]: .3f},{LIC_structure.lattice.abc[1]: .3f},{LIC_structure.lattice.abc[2]: .3f} Å")
print(f"Final energy is {LIC_energy.item(): .3f} eV")
print('--------------------------------------------------------')
LIC_structure.to('LIC.cif', 'cif')



#%% combine LIC with the downloaded structure, need to consider different exposure surface 计算界面形成能





'''

