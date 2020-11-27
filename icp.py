import numpy as np
from sklearn.neighbors import NearestNeighbors
from numpy import load
from random import randint
from Bio import PDB
from numpy import asarray

deletedIndex = []
# Accepts pdb, returns numpy array
def PDBToNPY(fpathin):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    struct = parser.get_structure('1ABZ',fpathin)
    allcoords1 = []
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x,y,z = atom.get_coord()
                    cSet = []
                    cSet.append(x)
                    cSet.append(y)
                    cSet.append(z)
                    allcoords1.append(cSet)
    #coordArray1 = asarray(allcoords1)
    return allcoords1

def best_fit_transform(A, B):
    '''
    Calculates the least-squares best-fit transform that maps corresponding points A to B in m spatial dimensions
    Input:
      A: Nxm numpy array of corresponding points
      B: Nxm numpy array of corresponding points
    Returns:
      T: (m+1)x(m+1) homogeneous transformation matrix that maps A on to B
      R: mxm rotation matrix
      t: mx1 translation vector
    '''

    assert A.shape == B.shape

    # get number of dimensions
    m = A.shape[1]

    # translate points to their centroids
    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    AA = A - centroid_A
    BB = B - centroid_B

    # rotation matrix
    H = np.dot(AA.T, BB)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)

    # special reflection case
    if np.linalg.det(R) < 0:
       Vt[m-1,:] *= -1
       R = np.dot(Vt.T, U.T)

    # translation
    t = centroid_B.T - np.dot(R,centroid_A.T)

    # homogeneous transformation
    T = np.identity(m+1)
    T[:m, :m] = R
    T[:m, m] = t

    return T, R, t


def nearest_neighbor(src, dst):
    '''
    Find the nearest (Euclidean) neighbor in dst for each point in src
    Input:
        src: Nxm array of points
        dst: Nxm array of points
    Output:
        distances: Euclidean distances of the nearest neighbor
        indices: dst indices of the nearest neighbor
    '''

    assert src.shape == dst.shape

    neigh = NearestNeighbors(n_neighbors=1)
    neigh.fit(dst)
    distances, indices = neigh.kneighbors(src, return_distance=True)
    return distances.ravel(), indices.ravel()


def icp(A, B, init_pose=None, max_iterations=20, tolerance=0.001):
    '''
    The Iterative Closest Point method: finds best-fit transform that maps points A on to points B
    Input:
        A: Nxm numpy array of source mD points
        B: Nxm numpy array of destination mD point
        init_pose: (m+1)x(m+1) homogeneous transformation
        max_iterations: exit algorithm after max_iterations
        tolerance: convergence criteria
    Output:
        T: final homogeneous transformation that maps A on to B
        distances: Euclidean distances (errors) of the nearest neighbor
        i: number of iterations to converge
    '''

    assert A.shape == B.shape

    # get number of dimensions
    m = A.shape[1]

    # make points homogeneous, copy them to maintain the originals
    src = np.ones((m+1,A.shape[0]))
    dst = np.ones((m+1,B.shape[0]))
    src[:m,:] = np.copy(A.T)
    dst[:m,:] = np.copy(B.T)

    # apply the initial pose estimation
    if init_pose is not None:
        src = np.dot(init_pose, src)

    prev_error = 0

    for i in range(max_iterations):
        # find the nearest neighbors between the current source and destination points
        distances, indices = nearest_neighbor(src[:m,:].T, dst[:m,:].T)

        # compute the transformation between the current source and nearest destination points
        T,_,_ = best_fit_transform(src[:m,:].T, dst[:m,indices].T)

        # update the current source
        src = np.dot(T, src)

        # check error
        mean_error = np.mean(distances)
        if np.abs(prev_error - mean_error) < tolerance:
            break
        prev_error = mean_error

    # calculate final transformation
    T,_,_ = best_fit_transform(A, src[:m,:].T)

    #return T, distances, i
    return T
# Accepts two numpy arrrays
def TransformAllPoints(coord1,coord2):
    global deletedIndex
    # Ensure they are the same size
    #coord1 = list(load(fname1))
    #coord2 = list(load(fname2))
    coord1Copy =[] #coord1
    coord2Copy =[]# coord2
    for c in coord1:
        coord1Copy.append(c)
    for c in coord2:
        coord2Copy.append(c)        
    #print("Point total of pdb1: " +str(len(coord1)))
    #print("Point total of pdb2: " +str(len(coord2)))
    if len(coord1Copy) != len(coord2):
        print("Resolving PDB size difference...")
        if len(coord1) > len(coord2):
            diff = len(coord1) - len(coord2)
            # assume shape is same, density is diff
            for i in range(diff):
                randInd = randint(0,len(coord1))
                del coord1[randInd]
                deletedIndex.append(randInd)
        else:
            diff = len(coord2) - len(coord1)
            # assume shape is same, density is diff
            for i in range(diff):
                randInd = randint(0,len(coord2))
                del coord2[randInd] 
                deletedIndex.append(randInd)
      
        #print("Point total of pdb1: " +str(len(coord1Copy)))
        #print("Point total of pdb2: " +str(len(coord2Copy)))
    
    # (A,B)
    T = icp(np.array(coord1),np.array(coord2))

    coord1Adjust = []
    for coord in coord1Copy:
        tempList = coord
        tempList.append(1)
        coord1Adjust.append(np.array(tempList))
    coord2Adjust = []
    for coord in coord2Copy:
        tempList = coord
        tempList.append(1)
        coord2Adjust.append(np.array(tempList))

    # Print side by side

    #print("Point total of pdb1: " +str(len(coord1Adjust)))
    #print("Point total of pdb2: " +str(len(coord2Adjust)))
    #https://github.com/ClayFlannigan/icp/blob/master/test.py
    #np.dot(T,coord1[0])

    print("Using transformation matrix:\n"+str(T))
    #print("Length of coord1 "+str(len(coord1Adjust)))
    #print("Length of coord1 "+str(len(coord2)))
    print("Length of coord 1 adjust: "+str(len(coord1Adjust)))
    print("Length of coord 2 adjust: "+str(len(coord2Adjust)))
    coord1 = np.array(coord1Adjust)
    coord1T = []
    for coord in coord1:
        coord1T.append(np.dot(T,coord))

    coord1T.reverse()
    with open('ATrans_1.txt', 'w') as f:
        for item in coord1T:
            f.write("%s\n" % str(item.tolist()))
    B = np.array(coord2).tolist()
    with open('B_1.txt', 'w') as f:
        for item in B:
            f.write("%s\n" % str(item))
def TransformPDBFile(fpath):
    with open(fpath) as f:
        lines = f.readlines()
    overHundred = False
    makeOverHundred = False
    foundIt = True
    for line in lines:
        pdbLine = line.split(' ')
        #print(pdbLine)
        xCoord = 0
        yCoord = 0
        zCoord = 0
        totalNums = 0
        
        for i in range(len(pdbLine)):
            #print(totalNums)
            try:
                coordVal = float(pdbLine[i])
                totalNums += 1
                #print(coordVal)
                if not overHundred:
                    if totalNums == 3:
                        #print(coordVal)
                        xCoord = coordVal
                    elif totalNums == 4:
                        #print(coordVal)
                        yCoord = coordVal
                    elif totalNums == 5:
                        zCoord = coordVal
                        #print(coordVal)
                        #print(pdbLine)
                        print(xCoord,yCoord,zCoord)
                else:
                    if totalNums == 4:
                        #print(coordVal)
                        xCoord = coordVal
                    elif totalNums == 5:
                        #print(coordVal)
                        yCoord = coordVal
                    elif totalNums == 6:
                        zCoord = coordVal
                        #print(coordVal)
                        #print(pdbLine)
                        print(xCoord,yCoord,zCoord) 
                if totalNums == 2 and int(coordVal) == 99 and foundIt:
                    makeOverHundred = True
                    print("Switching")
            except:
                pass
        if makeOverHundred:
            overHundred = True
            foundIt = False
def MenuDriver():
    while True:
        userChoice = input("Enter '1' for simple transformation and '2' for full-point registration:")
        if userChoice == '2':
            #fpath1 = input("Enter filepath to pdb1: ")
            #fpath2 = input("Enter filepath to pdb2: ")
            #fpath1 = 'coord1.npy'
            #fpath2 = 'coord2.npy'
            coord1 = PDBToNPY("helix1-axis.pdb")
            coord2 = PDBToNPY("helix2-axis.pdb")
            TransformAllPoints(coord1,coord2)
        elif userChoice == "1":
            print("Performing simple, single-point alignment...")
        else:
            print('Whoops, bad input.  Enter 1 or 2 for single point or multiple point alignment, respectively.')
#TransformPDBFile("helix1-axis.pdb")
MenuDriver()