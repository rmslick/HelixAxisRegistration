from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
import numpy as np
# Open 3-d points of the file that was
# transformed
with open("ATrans_1.txt") as f:
    coordinates = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
coordinates = [x.strip() for x in coordinates] 
coordArray = []
for coord in coordinates:
    # break up here
    coords = coord.split(",")
    coordArray.append(coords)

#for c in coordArray:
#    print(c)
# read
def CreatePDB(coordArray):
    sloppyparser = PDBParser(PERMISSIVE = True, QUIET = True) 
    structure = sloppyparser.get_structure("MD_system", "helix1-axisCopy.pdb")

    sb = StructureBuilder()
    sb.set_header(structure.header)
    # Iterate through models
    for i in range(len(list(structure.get_models()))):
        # Iterate through chains
        models = list(structure.get_models())
        counter = 0
        for j in range(len(list(models[i].get_chains()))):
            chains = list(models[i].get_chains())
            #Iterate thgouth residues
            for k in range(len(list(chains[j].get_residues()))):
                #Iterate through 
                residues = list(chains[j].get_residues())
                for l in range(len(list(residues[k].get_atoms()))):
                    #Set coord for each
                    for atom in structure[i][chains[j].id][residues[k].id].get_atoms():
                        structure[i][chains[j].id][residues[k].id][atom.id].set_coord(np.array((float(coordArray[counter][0]),float(coordArray[counter][1]),float(coordArray[counter][2]))))
                        print(structure[i][chains[j].id][residues[k].id][atom.id].get_vector())
                    counter += 1


    io=PDBIO()
    io.set_structure(structure)
    io.save("bio-pdb-pdbio-out.pdb")

'''
io=PDBIO()
model = structure.get_models()
models = list(model)
for model in models:
    chains = list(model.get_chains()) 
    counter = 0
    for chain in chains:
        residues = list(chain.get_residues())
        for residue in residues:
            atoms = list(residue.get_atoms())
            for atom in atoms:
                atomVex = atom.get_vector()
                for i in range(3):
                    atomVex.__setitem__(i,float(coordArray[counter][i]))
                print(atomVex)
            counter += 1
'''


