#!/usr/bin/env python
# coding: utf-8

# In[25]:


import math as m


# In[26]:


def graphene(x_rings, y_rings):
    
    bond_length = 1.4 # angstrom
    theta = (30*m.pi)/180.0
    
    x_length = bond_length*m.cos(theta)
    y_length = bond_length*m.sin(theta)
    
    total_atoms = []
    x_displace = []
    for i in range(y_rings):
        if i%2 == 0:
            atoms = x_rings + 1
            total_atoms.append(atoms)
            total_atoms.append(atoms)
            x_displace.append(0.0000)
            x_displace.append(0.0000)
        elif i%2 == 1:
            atoms = x_rings
            total_atoms.append(atoms)
            total_atoms.append(atoms)
            x_displace.append(x_length)
            x_displace.append(x_length)
    
    if y_rings%2 == 0:
        total_atoms.append(x_rings - 1)
        x_displace.append(2*x_length)
    elif y_rings%2 == 1:
        total_atoms.append(x_rings)
        x_displace.append(x_length)
    
    total_atoms.insert(0, x_rings)
    x_displace.insert(0, x_length)
    
    x_move = 2*x_length
    
    x_coords = []
    for i in range(x_rings + 1):
        x_coords.append(x_move*i)
        
    y_layers = y_rings*2 + 2
    
    y_coords = [0.0]
    next_coord = 0.0
    for i in range(y_layers - 1):
        if i%2 == 0:
            next_coord += 0.7
            y_coords.append(next_coord)
        elif i%2 == 1:
            next_coord += 1.4
            y_coords.append(next_coord)
    
    dump_coords = []
    for i in y_coords:
        partial_dump = []
        for j in x_coords:
            dump_coord = []
            dump_coord.append(j)
            dump_coord.append(i)
            partial_dump.append(dump_coord)
        dump_coords.append(partial_dump)
    
    for i in range(len(dump_coords)):
        for j in range(len(dump_coords[i])):
            dump_coords[i][j][0] += x_displace[i]
    
    finalized_coords = []  
    for i in range(len(dump_coords)):
        for j in range(total_atoms[i]):
            finalized_coords.append(dump_coords[i][j])
    
    cup_base = []
    ids = 1
    if y_rings%2 == 1:
        for i in range(y_rings//2 + 1):
            cup_base.append(ids)
            ids += 2*x_rings + 2*(x_rings + 1)
    elif y_rings%2 == 0:
        for i in range(y_rings//2):
            cup_base.append(ids)
            ids += 2*x_rings + 2*(x_rings + 1)
    
    cup_bases = []
    for i in cup_base:
            for j in range(x_rings):
                cup_bases.append(i + j)
    
    cap_peak = []
    ids = 3*(x_rings + 1)
    for i in range(y_rings//2):
        cap_peak.append(ids)
        ids += 2*(x_rings + 1) + 2*x_rings
    
    cap_peaks = []
    for i in cap_peak:
        for j in range(x_rings):
            cap_peaks.append(i + j)
    
    bonds = []
    for i in cup_bases:
        bond = []
        bond.append(i)
        bond.append(i + x_rings)
        bonds.append(bond)
        bond = []
        bond.append(i)
        bond.append(i + x_rings + 1)
        bonds.append(bond)
    
    for i in cap_peaks:
        bond = []
        bond.append(i)
        bond.append(i - x_rings - 1)
        bonds.append(bond)
        bond = []
        bond.append(i)
        bond.append(i - x_rings)
        bonds.append(bond)
        
    modified_cup_bases = cup_bases[x_rings:]
    
    for i in range(len(modified_cup_bases)):
        bond = []
        bond.append(modified_cup_bases[i])
        bond.append(cap_peaks[i])
        bonds.append(bond)
    
    jump_start = 1
    add_start = 1
    cup_cap_base = []
    for i in range(y_rings//2):
        cup_cap_base.append(x_rings*jump_start + add_start)
        jump_start += 4
        add_start += 2
    
    cup_cap_bases = []
    for i in cup_cap_base:
        for j in range(x_rings + 1):
            cup_cap_bases.append(i + j)
    
    for i in cup_cap_bases:
        bond = []
        bond.append(i)
        bond.append(i + x_rings + 1)
        bonds.append(bond)
    
    if y_rings%2 == 0:
        ceiling_bases = []
        for i in range(sum(total_atoms) - x_rings + 2, sum(total_atoms) + 1):
            ceiling_bases.append(i)
    
        for i in ceiling_bases:
            bond = []
            bond.append(i)
            bond.append(i - x_rings)
            if i == ceiling_bases[0]:
                hammer_start = i - x_rings
            bonds.append(bond)
            bond = []
            bond.append(i)
            bond.append(i - x_rings + 1)
            if i == ceiling_bases[-1]:
                hammer_stop = i - x_rings + 1
            bonds.append(bond)
    
        for i in range(hammer_start, hammer_stop + 1):
            bond = []
            bond.append(i)
            bond.append(i - x_rings)
            bonds.append(bond)
    
    elif y_rings%2 == 1:
        ceiling_bases = []
        for i in range(sum(total_atoms) - x_rings + 1, sum(total_atoms) + 1):
            ceiling_bases.append(i)
        print(ceiling_bases)    
        for i in ceiling_bases:
            bond = []
            bond.append(i)
            bond.append(i - x_rings - 1)
            if i == ceiling_bases[0]:
                hammer_start = i - x_rings - 1
            bonds.append(bond)
            bond = []
            bond.append(i)
            bond.append(i - x_rings)
            if i == ceiling_bases[-1]:
                hammer_stop = i - x_rings
            bonds.append(bond)
        
        for i in range(hammer_start, hammer_stop + 1):
            bond = []
            bond.append(i)
            bond.append(i - x_rings - 1)
            bonds.append(bond)
    
    f = open("graphene.mol2", "w")
    f.writelines("# Name: graphene.mol2")
    f.writelines("\n")
    f.writelines("# Written by: graphene.ipynb developed by Aamir Alaud Din")
    f.writelines("\n")
    f.writelines("# Date: 2023.04.23")
    f.writelines("\n\n")
    f.writelines("@<TRIPOS>MOLECULE")
    f.writelines("\n")
    f.writelines("Graphene")
    f.writelines("\n")
    f.writelines(f'  {sum(total_atoms)}  {len(bonds)}  0  0  0')
    f.writelines("\n")
    f.writelines("SMALL")
    f.writelines("\n")
    f.writelines("NO CHARGES")
    f.writelines("\n")
    f.writelines("*****")
    f.writelines("\n")
    f.writelines("@<TRIPOS>ATOM")
    for i in range(len(finalized_coords)):
        f.writelines("\n")
        f.writelines(f'    {i + 1:<6d}    C  {finalized_coords[i][0]:<9.4f}  {finalized_coords[i][1]:<9.4f}  0.0000  C.ar  1  GRAPHENE')
    f.writelines("\n")
    f.writelines("@<TRIPOS>BOND")
    for i in range(len(bonds)):
        f.writelines("\n")
        f.writelines(f'  {i + 1:<10d}    {bonds[i][0]:<8d}    {bonds[i][1]:<8d}    1')
    f.close()
    
    f = open("graphene.lt", "w")
    f.writelines("# Name: graphene.lt")
    f.writelines("\n")
    f.writelines("# Written by: graphene.ipynb developed by Aamir Alaud Din")
    f.writelines("\n")
    f.writelines("# Date: 2023.04.26")
    f.writelines("\n\n")
    f.writelines("import \"gaff.lt\"")
    f.writelines("\n\n")
    f.writelines("Graphene inherits GAFF {")
    f.writelines("\n")
    f.writelines("  write('Data Atoms') {")
    f.writelines("\n")
    for i in range(len(finalized_coords)):
        f.writelines("\n")
        f.writelines(f'    $atom:c{i + 1:<6d}    $mol:.    @atom:ca    0.0000    {finalized_coords[i][0]:<9.4f}    {finalized_coords[i][1]:<9.4f}    0.0000')
    f.writelines("\n  }\n")
    f.writelines("  write('Data Bond List')  {")
    for i in range(len(bonds)):
        f.writelines("\n")
        f.writelines(f'    $bond:b{i + 1:<6d}    $atom:c{bonds[i][0]:<8d}    $atom:c{bonds[i][1]:<8d}')
    f.writelines("\n  }")
    f.writelines("\n} # End of Graphene")
    f.close()
    f = open("system.lt", "w")
    f.writelines("# Name: system.lt")
    f.writelines("\n")
    f.writelines("# Written by graphene.ipynb developed by Aamir Alaud Din")
    f.writelines("\n")
    f.writelines("# Date: 2023.04.26")
    f.writelines("\n")
    f.writelines("import \"graphene.lt\"")
    f.writelines("\n\n")
    f.writelines("write_once('Data Boundary')  {")
    f.writelines("\n")
    f.writelines(f'  0.0  {finalized_coords[2*x_rings][0]:<9.4f}  xlo  xhi')
    f.writelines("\n")
    f.writelines(f'  0.0  {finalized_coords[-1][1]:<9.4f}  ylo  yhi')
    f.writelines("\n")
    f.writelines(f'  0.0  0.0  zlo  zhi')
    f.writelines("\n}")
    f.writelines("\n\n")
    f.writelines("GRAPHENE = new Graphene")
    f.close()


# In[27]:


graphene(15, 22)

