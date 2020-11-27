import chimera
from numpy import asarray
from numpy import save
from Bio import PDB

def PDBToNPY(fpathin, fpathout):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    struct = parser.get_structure('1ABZ',fpathin)

    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x,y,z = atom.get_coord()
                    print(x,y,z)
# define data
data = asarray([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]])
# save to npy file
#save('data.npy', data)
# Open each protein database file
helixAxis1=chimera.openModels.open('helix1-axis.pdb',type="PDB")
helixAxis2=chimera.openModels.open('helix2-axis.pdb',type="PDB")
hAx1 = helixAxis1[0]
hAx2 = helixAxis2[0]
count1 = 0
residues1 = hAx1.residues
residues2 = hAx2.residues
atoms = []
# Has all the coordinates for corresponding atom
coord = []
allcoords1 = []
for i in range(len(residues1)):
    resAtoms = residues1[i].atoms
    for i in range(len(resAtoms)):
        atoms.append(resAtoms[i])
        coord.append(resAtoms[i].coord())
        someCord = resAtoms[i].coord()
        cSet = []
        cSet.append(someCord[0])
        cSet.append(someCord[1])
        cSet.append(someCord[2])
        allcoords1.append(cSet)       
print("Done")
atoms2 = []
# Has all the coordinates for corresponding atom
coord2 = []
allcoords2 = []
for i in range(len(residues2)):
    resAtoms = residues2[i].atoms
    for i in range(len(resAtoms)):
        atoms2.append(resAtoms[i])
        coord2.append(resAtoms[i].coord())
        coord = resAtoms[i].coord()
        cSet = []
        cSet.append(coord[0])
        cSet.append(coord[1])
        cSet.append(coord[2])
        allcoords2.append(cSet)
print("Done2")
#the x-coordinate use: my_atomCoords[0]
#the y-coordinate use: my_atomCoords[1]
#the z-coordinate use: my_atomCoords[2]
print(len(coord))
print(len(coord2))
#print(allcoords2)
#print(allcoords1)
coordArray1 = asarray(allcoords1)
coordArray2 = asarray(allcoords2)
save('coord1_t.npy',coordArray1)
save('coord2_t.npy',coordArray2)
