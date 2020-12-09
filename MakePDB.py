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
# Create a pdb
with open("helix1-axisCopy.pdb") as f:
    pdbLines = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
pdbs = [x.strip() for x in pdbLines] 
for pdb in pdbs:
    print(pdb)
transformedRows = []
for i in range(len(coordArray)):
    rowSplit = pdbs[i].split(' ')
    newCoords = coordArray[i]
    k = 0
    for count in range(len(rowSplit)) :
        if rowSplit[count].find('.') != -1:
            rowSplit[count] = str('{:.3f}'.format(round(float(newCoords[0]),3)))
            break
        k = count+1
    
    z = k+1
    for j in range(z, len(rowSplit)):
        if rowSplit[z].find('.') != -1:
            rowSplit[z] = str('{:.3f}'.format(round(float(newCoords[1]),3)))
            break
        z += 1
    z +=1
    
    for j in range(z, len(rowSplit)):
        if rowSplit[z].find('.') != -1:
            rowSplit[z] = str('{:.3f}'.format(round(float(newCoords[2]),3)))
            break
    newRow = ""
    counter = 0
    for ele in rowSplit:
        if counter == len(rowSplit) -1:
            newRow += ele
        elif len(ele) > 0:
            newRow += ele
            newRow += " "
        else:
            newRow += " "
        counter += 1 
    #print(newRow)
    transformedRows.append(newRow)
pdbFileNew = []
for i in range(len(transformedRows)):
    pdbs[i] = transformedRows[i]

f = open("helix1-axisCopy2.pdb","w")
for pdb in pdbs:
    print(pdb)
    f.write(pdb+"\n")

#Write out transformed file
