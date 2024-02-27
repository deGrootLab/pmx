import sys, os
import copy as xcopy
#import networkx as nx
import matplotlib
#matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pmx import *
import pmx.parser as pmxparser
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

###################################################################
################## graph plotting options #########################
###################################################################
class graphOptions:
    def __init__(self, nodeDeg, bMST=True, filename=None):
        self.opt = {}

        ### by default, MST graph options ###
        self.opt['dpi'] = 300
        self.opt['layout'] = "spring" # spring, random, shell, spectral, circular, fruchterman_reingold
        # node options
        self.opt['nodeSizeScale'] = 100
        self.opt['nodeSizeSlope'] = 0.2 # scaling for the node size
        self.opt['nodeShape'] = "o" # s o ^ > v < d p h 8
        self.opt['nodeColor'] = nodeDeg # nodeDeg, b: blue, g: green, r: red, c: cyan, m: magenta, y: yellow, k: black, w: white
        self.opt['nodeCMAP'] = "hot" # http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
        self.opt['nodeVmin'] = 0.1 #
        self.opt['nodeVmax'] = 1.5 #
        self.opt['nodeTransparency'] = 1.0
        self.opt['nodeLineWidth'] = 0.0 # scalar, array
        # node label options
        self.opt['nodeLabels'] = True
        self.opt['nodeFontColor'] = 'w'
        self.opt['nodeFontSize'] = 3
        self.opt['nodeFontFamily'] = 'sans-serif'
        self.opt['nodeFontWeight'] = 'bold' # http://matplotlib.org/api/font_manager_api.html
        self.opt['nodeLabelOffsetX'] = 0.0
        self.opt['nodeLabelOffsetY'] = 0.0
        # edge options
        self.opt['edgeWidth'] = 0.1
        self.opt['edgeColor'] = 'k'
        self.opt['edgeColorExtra'] = None
        self.opt['edgeStyle'] = 'solid' # solid|dashed|dotted,dashdot
        self.opt['edgeTransparency'] = 0.5
        self.opt['edgeCMAP'] = None
        self.opt['edgeVmin'] = None
        self.opt['edgeVmax'] = None
        # edge label options
        self.opt['edgeLabels'] = None
        self.opt['edgeFontColor'] = 'r'
        self.opt['edgeFontSize'] = 8
        self.opt['edgeFontFamily'] = 'sans-serif'
        self.opt['edgeFontWeight'] =  'normal' #'bold'
        self.opt['edgeLabelPos'] = 0.5 # 0=head, 0.5=center, 1=tail)
        self.opt['edgeLabelTransp']= 1.0

        ### default for the RefLig
        if bMST==False:
            self.opt['layout'] = "circular" # spring, random, shell, spectral, circular, fruchterman_reingold

        ### read from file
        if filename is not None:
            self.read_from_file(filename,nodeDeg)

        self.opt['nodeVmin'] = min(nodeDeg)*self.opt['nodeVmin'] #
        self.opt['nodeVmax'] = max(nodeDeg)*self.opt['nodeVmax'] #

    def read_from_file(self,filename,nodeDeg):
        fp = open(filename,'r')
        lines = fp.readlines()
        fp.close()
        lines = pmxparser.kickOutComments( lines, comment = '#')
        for l in lines:
            l = l.rstrip()
            l = l.lstrip()
            foo = l.split()
            foo[-1] = foo[-1].replace("'","")
            foo[-1] = foo[-1].replace("\"","")
            if 'None' == foo[-1]:
                self.opt[foo[0]] = None
                continue
            if 'False' == foo[-1]:
                self.opt[foo[0]] = None
                continue
            if 'True' == foo[-1]:
                self.opt[foo[0]] = True
                continue
            self.opt[foo[0]] = foo[-1]
            if ('Size' in foo[0]) or ('min' in foo[0]) or ('max' in foo[0]) or ('Transp' in foo[0]) or ('Width' in foo[0]) or ('Pos' in foo[0]) or ('Offset' in foo[0]):
                self.opt[foo[0]] = float(foo[-1])
        # some options need adjustment
        if self.opt['nodeColor']=='var':
            self.opt['nodeColor']=nodeDeg
###################################################################
###################################################################
###################################################################

def plotImages(cmdl,layout,G,nsize,nodeNames,ligImg):
    ax = plt.gca()
    fig = plt.gcf()
    trans=ax.transData.transform
    trans2=fig.transFigure.inverted().transform
    const_imsize=0.00015 # this is the image size
    impos=0.0 # image position

    for n in G.nodes():
        namePNG = namePathFromPNG(ligImg)
        nodeName = nodeNames[n]
        imgFilename = namePNG+"/"+nodeName+".png"
        if os.path.isfile(imgFilename)==False:
            print("Image not found: ",imgFilename)
            continue
#        img=mpimg.imread('/home/vgapsys/project/make_hybrid/thrombin_test_set/o3a_mapping/foo/lig_2ZC9.png')
        img=mpimg.imread(imgFilename)
        G.node[n]['image'] = img
        imsize = const_imsize*nsize[n]
        imsize = 0.1
        xx,yy=trans(layout[n]) # figure (graph) coordinates
        xa,ya=trans2((xx,yy)) # axes coordinates
        a = plt.axes([xa-imsize/2.0,ya-imsize/2.0,imsize,imsize])
        a.imshow(G.node[n]['image'],zorder=1)
        a.set_aspect('equal')
        a.axis('off')

    ax.axis('off')

def edgesRefLig(refLig,nodeNames,mat,cycle):
    edges = []
    # first edges are special
    if cycle==0:
        for i in range(0,shape(mat)[0]):
            if nodeNames[i]==refLig:
                for j in range(0,shape(mat)[1]):
                    if i==j:
                        continue
                    newEdge = [i,j]
                    edges.append(newEdge)
                    mat[i,j] = np.inf
                    mat[j,i] = np.inf
                break
    # other edges are different: constructing cycles
    else:
        for i in range(0,shape(mat)[0]):
            if nodeNames[i]==refLig:
                continue
            else:
                j = np.argsort(mat[i])[0,1]
                newEdge = [i,j]
                edges.append(newEdge)
                mat[i,j] = np.inf
                mat[j,i] = np.inf
    return(edges)

def identifyRefLig(mat,cmdlopt,nodeNames):
    refName = ''
    if cmdlopt.is_set:
        refName = cmdlopt.value
        if refName not in nodeNames:
            print("The name provided for the reference ligand is not found among the molecules: ",refName)
            sys.exit(1)
    else:
        minRowSum = 999.99
        for i in range(0,shape(mat)[0]):
            rowSum = np.sum(mat[i])
            if rowSum < minRowSum:
                minRowSum = rowSum
                refName = nodeNames[i]
    return(refName)

def plotMol(pdb,filename):
    pdbName = reformatPDB(pdb,1)
    mol = Chem.MolFromPDBFile(pdbName,removeHs=False,sanitize=False)
    os.remove(pdbName)
#    AllChem.GenerateDepictionMatching2DStructure(mol,refMol)
    opt = Draw.DrawingOptions()
    opt.bgColor=None
    Draw.MolToFile(mol,filename,kekulize=False,options=opt)

def plotMolTransparent(pdb,filename):
    pdbName = reformatPDB(pdb,1)
    mol = Chem.MolFromPDBFile(pdbName,removeHs=False,sanitize=False)
    os.remove(pdbName)
    newmol = Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol,forceCoords=True,kekulize=True)
    mc = Chem.Mol(newmol.ToBinary())
    Chem.Kekulize(mc)
    img = Draw.MolToImage(mc,size=(450, 450),kekulize=True)

    # transparent background
    img = img.convert("RGBA")
    datas = img.getdata()
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    img.putdata(newData)

    img.save(filename, "PNG",dpi=(300.0, 300.0))

def writeFormatPDB(fname,m,title="",nr=1):
    fp = open(fname,'w')
    for atom in m.atoms:
        foo = cp.deepcopy(atom)
        # chlorine
        if( 'CL' in atom.name or 'Cl' in atom.name or 'cl' in atom.name ):
            foo.name = "CL"+"  "
            print(foo, file=fp)
        # bromine
        elif( 'BR' in atom.name or 'Br' in atom.name or 'br' in atom.name ):
            foo.name = "BR"+"  "
            print(foo, file=fp)
        elif( len(atom.name) > 4): # too long atom name
            foo = cp.deepcopy(atom)
            foo.name = foo.name[:4]
            print(foo, file=fp)
        else:
            print(atom, file=fp)
    print('ENDMDL', file=fp)
    fp.close()

def reformatPDB(filename,num):
#    newname = filename.split(".")[0]+"_tempFormat.pdb"
    newname = "tempFormat"+str(num)+".pdb"
    m = Model().read(filename)
    rem = []
    for a in m.atoms:
        if 'EP' in a.name:
            rem.append(a)
    # remove
    for a in rem:
        m.remove_atom(a)
    writeFormatPDB(newname,m)
    return(newname)

def findInPythonpath(filename):
    path = ""
    pythonpath = os.environ.get('PYTHONPATH')
    foo = pythonpath.split(":")
    for folder in foo:
        script = folder+"/"+filename
        if(os.path.isfile(script)==True):
            path = script
            break
    return(path)

def getNodeNames(mat,fileName=None):
    nodeNames = []
    nodeNum = len(mat)
    # node names are simply numbers
    if(fileName==None):
        count = 0
        for n in range(0,nodeNum):
            nodeNames.append(count)
            count = count + 1
    # node names are read from a file
    else:
        with open(fileName) as fp:
            for line in fp:
                line = line.rstrip()
                nodeNames.append(line.split()[0])
        fp.close()
    return(nodeNames)

def getNodeDegrees(nodeNames,edgeList):
    nodeConn = [0 for n in range(0,len(nodeNames))]
    for edge in edgeList:
        nodeConn[edge[0]] = nodeConn[edge[0]] + 1
        nodeConn[edge[1]] = nodeConn[edge[1]] + 1
    return(nodeConn)

def writeNodes(fileName,nodeNames,nodeDeg):
    fp = open(fileName,"w")
    for n,conn in zip(nodeNames,nodeDeg):
        fp.write("%s    %s\n" % (n,conn))
    fp.close()

def writeEdges(filename,nodeNames,edgeList,mat):
    fp = open(filename,"w")
    for edge in edgeList:
        n1 = nodeNames[edge[0]]
        n2 = nodeNames[edge[1]]
        fp.write("%s    cr      %s      %s\n" % (n1,n2,mat[edge[0],edge[1]]))
    fp.close()

def modifyMat(edges,mat):
    for edge in edges:
        e1 = edge[0]
        e2 = edge[1]
        mat[e1,e2] = np.inf
        mat[e2,e1] = np.inf
    return mat

def genGraph(mat,edgeList,nx):
    G=nx.Graph()
    nodes = list(range(0,mat.shape[0]))
    G.add_nodes_from(nodes)
    edgesG = []
    for edges in edgeList:
        edge = (edges[0],edges[1],mat[edges[0],edges[1]])
        edgesG.append(edge)
    G.add_weighted_edges_from(edgesG)
    return(G)

def genNodeLabels(mat,nodeNames,layout,opt):
    labels = {}
    nodes = list(range(0,mat.shape[0]))
    i = 0
    for node,name in zip(nodes,nodeNames):
        labels[i] = name
        i = i + 1

    # adjust label positions
    outlayout = cp.deepcopy(layout)
    for key in list(layout.keys()):
        outlayout[key][0] += opt['nodeLabelOffsetX']
        outlayout[key][1] += opt['nodeLabelOffsetY']

    return(labels,outlayout)

def genEdgeLabels(G):
    labels = {}
    i = 0
    for edge in G.edges():
        ind = (edge[0],edge[1])
        labels[ind] = round(G[edge[0]][edge[1]]['weight'],2)
        i = i + 1
    return(labels)

# courtesy Andreas Mueller
def mstPrim(matRMSD):
    mat = matRMSD.copy()
    vertexNum = mat.shape[0]
    spanningEdges = []
    vertexVisited = [0]
    visitedNum = 1
    diagInd = np.arange(vertexNum)
    mat[diagInd, diagInd] = np.inf
    while visitedNum != vertexNum:
        newEdge = np.argmin(mat[vertexVisited], axis=None)
        newEdge = divmod(newEdge, vertexNum)
        newEdge = [vertexVisited[newEdge[0]], newEdge[1]]
        spanningEdges.append(newEdge)
        vertexVisited.append(newEdge[1])
        mat[vertexVisited, newEdge[1]] = np.inf
        mat[newEdge[1], vertexVisited] = np.inf
        visitedNum += 1
    return(np.vstack(spanningEdges))

def readMat(fileName):
    mat = []
    rowNum = 0
    with open(fileName) as fp:
        for line in fp:
            line = line.rstrip()
            line = line.lstrip()
            foo = line.split(" ")
            dim = len(foo)
            if(rowNum == 0):
                mat = initMat(dim)
            for i in range(0,dim):
                mat[rowNum][i] = float(foo[i])
            rowNum = rowNum + 1
            #print(line)
    fp.close()
    return(mat)

def initMat(dim):
    mat = [[0 for x in range(dim)] for x in range(dim)]
    return(mat)

def symMat(mat,row):
    for i in range(0,row):
        mat[row][i] = mat[i][row]
    mat[row][row] = 0
    return(mat)

def nameFromPDB(pdb):
    name = xcopy.deepcopy(pdb)
    foo = name.split("/")[-1]
    bar = foo.split(".")[0]
    return(bar)

def namePathFromPNG(png):
    name = xcopy.deepcopy(png)
    bar = name.split(".")[0]
    return(bar)

def doMorphes(pdbs,callString,outMatFile,ligImg=None):
    visited = []
    nodeNames = []
    dim = len(pdbs)
    rmsdMat = initMat(dim)
    rowCount = 0

    for pdb1 in pdbs:
        if(ligImg!=None):
            namePDB = nameFromPDB(pdb1)
            namePNG = namePathFromPNG(ligImg)
            imgFilename = namePNG+"/"+namePDB+".png"
            # check if png exists
            if os.path.isfile(imgFilename)==False:
                plotMolTransparent(pdb1,imgFilename)
            #sys.exit(0)
        rmsdMat = symMat(rmsdMat,rowCount)
        visited.append(pdb1)
        colCount = rowCount + 1
        nodeNames.append( nameFromPDB(pdb1) )
        for pdb2 in pdbs:
            if pdb2 in visited:
                continue
            callS = callString+" -i1 "+pdb1+" -i2 "+pdb2+" -score "+outMatFile
            print(callS)
            os.system(callS)
            fp = open(outMatFile,"r")
            foo = fp.read().rstrip()
            rmsdMat[rowCount][colCount] = float(foo.split(" ")[-1])
            fp.close()
            colCount = colCount + 1
        rowCount = rowCount + 1

    fp = open(outMatFile,"w")
    for row in rmsdMat:
        l = ''
        for col in row:
            l = l+" "+str(col)
        fp.write("%s\n" % l)
    fp.close()

    return(rmsdMat,nodeNames)

def main(argv):

    desc=('Build a Minimum Spanning Tree or a Tree based on a reference ligand suggesting the alchemical pairs of ligands to morphe.',
        '',
        '-ligImg: provide a path to the ligand images. The script will search for /path/NodeName.png files.',
        'If -ligImg path has no images (as .png files), the images will be created from the .pdb files (if -pdb is provided).',
        )

# define input/output files

    files= [
        FileOption("-rmsd", "r/o",["dat"],"rmsd_mat.dat", "an rmsd matrix"),
        FileOption("-pdb", "r/o/m",["pdb"],"path/*.pdb", "a set of pdbs"),
        FileOption("-ormsd", "w/o",["dat"],"outRMSD.dat", "output dissimilarity matrix"),
        FileOption("-oNodes", "w/o",["dat"],"nodeNames.dat", "output names of the graph nodes"),
        FileOption("-oEdges", "w",["dat"],"edges.dat", "output names of the graph edges"),
        FileOption("-iNodes", "r/o",["dat"],"nodeNames.dat", "read names of the graph nodes"),
        FileOption("-graph", "w/o",["png","svg","pdf"],"graph.png", "plot a graph"),
        FileOption("-graphOpt", "r/o",["dat"],"graphOpt.dat", "graph options file (see example)"),
        FileOption("-ligImg", "r/w/o",["png"],"path_to_images.png", "path to images of ligands"),
        ]

# define options

    options=[
        Option( "-mst", "bool", "True", "build a minimum spanning tree, otherwise use one reference ligand"),
        Option( "-refName", "string", "none", "optionally provide a name of the reference ligand"),
        Option( "-alignment", "bool", "False", "if -pdb==True method 1: 3D alignment"),
        Option( "-mcs", "bool", "False", "if -pdb==True method 2: maximum common substructure"),
        Option( "-H2H", "bool", "False", "should hydrogen be morphed into hydrogen?"),
        Option( "-H2heavy", "bool", "False", "should hydrogen be morphed into a heavy atom (also sets -H2H to true)"),
        Option( "-dMCS", "bool", "False", "find MCS, superimpose the structures based on the MCS, apply distance in\
                Cartesian coordinate space to define the morphes"),
        Option( "-d", "float", "0.05", "distance (nm) between atoms to consider them morphable"),
        Option( "-cycleNum", "int", "0", "number of cycles to build"),
        Option( "-timeout", "float", "60.0", "maximum time (s) for an MCS search"),
        ]

    help_text = ("Build a Minimum Spanning Tree suggesting the alchemical pairs of ligands to morphe.",
                 "Provide an RMSD matrix of distances between ligands, a path to the ligand pdbs or a path to the distance scores.",
                 "If the pdbs are provided, the atoms_to_morph.py will be called to generate an RMSD matrix."
                )

    cmdl = Commandline( argv, options = options, fileoptions = files, program_desc = help_text, check_for_existing_files = False, version = "0.0" )

    # deal with flags
    outMatFile = "outRMSD.dat"
    cycleNum = 0
    ligImg = None
    if cmdl.opt['-ormsd'].is_set:
        outMatFile = cmdl['-ormsd']
    if cmdl.opt['-cycleNum'].is_set:
        cycleNum = cmdl['-cycleNum']
    if cmdl.opt['-ligImg'].is_set:
        ligImg = cmdl['-ligImg']
    bMST = True
    if cmdl.opt['-mst'].is_set:
        bMST = cmdl.opt['-mst'].value

    # pdb_list or rmsd_mat
    if( cmdl.opt['-pdb'].is_set==True ):
        print("Calculate dissimilarities between the ligands.")
#       if( (cmdl.opt['-alignment'].is_set==False) and (cmdl.opt['-mcs'].is_set==False) ):
#           print "For ligand mapping need to select -alignment or -mcs"
#           sys.exit(1)
        if( len(findInPythonpath("atoms_to_morph.py"))==0 ):
            print("Place atoms_to_morph.py in the PYTHONPATH")
            sys.exit(1)
        callString = "python "+findInPythonpath("atoms_to_morph.py")
        if cmdl.opt['-alignment'].is_set:
            callString = callString+" -alignment"
        if cmdl.opt['-mcs'].is_set:
            callString = callString+" -mcs"
        if cmdl.opt['-H2H'].is_set:
            callString = callString+" -H2H"
        if cmdl.opt['-H2heavy'].is_set:
            callString = callString+" -H2heavy"
        if cmdl.opt['-dMCS'].is_set:
            callString = callString+" -dMCS"
        if cmdl.opt['-timeout'].is_set:
            callString = callString+" -timeout "+str(cmdl['-timeout'])
        if cmdl.opt['-d'].is_set:
            callString = callString+" -d "+str(cmdl['-d'])
        rmsdMat,nodeNames = doMorphes(cmdl['-pdb'],callString,outMatFile,ligImg)
    elif( cmdl.opt['-rmsd'].is_set==True ):
        rmsdMat = readMat(cmdl['-rmsd'])
        if( cmdl.opt['-iNodes'].is_set==True ):
            nodeNames = getNodeNames(rmsdMat,cmdl['-iNodes'])
        else:
            nodeNames = getNodeNames(rmsdMat)
    else:
        print("Need to provide either an RMSD matrix or a set of ligand pdbs.")
        sys.exit(0)

    # mst or one reference ligand
    rmsdMat = np.matrix(rmsdMat)
    nodeNum = shape(rmsdMat)[0]
    mat = xcopy.deepcopy(rmsdMat)
    edgeList = []
    edgeListMain = []
    edgeListExtra = []
    if bMST: # mst
        for mst in range(0,cycleNum+1):
            edges = mstPrim(mat)
            mat = modifyMat(edges,mat)
            edgeList.extend(edges)
            if mst==0:
                edgeListMain.extend(edges)
            else:
                edgeListExtra.extend(edges)
    else: # reference
        refLig = identifyRefLig(mat,cmdl.opt['-refName'],nodeNames)
        for cycle in range(0,cycleNum+1):
            edges = edgesRefLig(refLig,nodeNames,mat,cycle)
            edgeList.extend(edges)
            if cycle==0:
                edgeListMain.extend(edges)
            else:
                edgeListExtra.extend(edges)
        print("Using reference ligand: ",refLig)

    # output nodes
    nodeDeg = getNodeDegrees(nodeNames,edgeList)
    if( cmdl.opt['-oNodes'].is_set==True ):
        writeNodes(cmdl['-oNodes'],nodeNames,nodeDeg)
    # output edges
    writeEdges(cmdl['-oEdges'],nodeNames,edgeList,rmsdMat)

    # plot graph
    if( cmdl.opt['-graph'].is_set ):
        import networkx as nx

        fnameGraphOpt = None
        if cmdl.opt['-graphOpt'].is_set:
            fnameGraphOpt = cmdl['-graphOpt']

        G = genGraph(rmsdMat,edgeList,nx)
        GO = graphOptions(nodeDeg,bMST=bMST,filename=fnameGraphOpt)
        # graph layout
        if( GO.opt['layout'] == "spring" ):
            layout=nx.spring_layout(G,weight='weight')
        elif( GO.opt['layout'] == "random" ):
            layout=nx.random_layout(G)
        elif( GO.opt['layout'] == "shell" ):
            layout=nx.shell_layout(G)
        elif( GO.opt['layout'] == "spectral" ):
            layout=nx.spectral_layout(G,weight='weight')
        elif( GO.opt['layout'] == "circular" ):
            layout=nx.circular_layout(G)
        elif( GO.opt['layout'] == "fruchterman_reingold" ):
            layout=nx.fruchterman_reingold_layout(G)

        # node size
        nsize = [ ( ( (i-min(nodeDeg))*GO.opt['nodeSizeSlope']+min(nodeDeg) ) ) * GO.opt['nodeSizeScale'] for i in nodeDeg]

        # draw nodes
        if(cmdl.opt['-ligImg'].is_set==False):
            nx.draw_networkx_nodes(G,pos=layout,node_size=nsize,node_color=GO.opt['nodeColor'],alpha=GO.opt['nodeTransparency'], vmin=GO.opt['nodeVmin'], vmax=GO.opt['nodeVmax'], cmap=GO.opt['nodeCMAP'], linewidths=GO.opt['nodeLineWidth'],node_shape=GO.opt['nodeShape'])

        # node labels
        if( (GO.opt['nodeLabels'] != None) and (GO.opt['nodeLabels'] != False) ):
            nodeLabels,nodeLabelLayout=genNodeLabels(rmsdMat,nodeNames,layout,GO.opt)
            nx.draw_networkx_labels(G,pos=nodeLabelLayout,labels=nodeLabels,font_color=GO.opt['nodeFontColor'],font_weight=GO.opt['nodeFontWeight'],font_size=GO.opt['nodeFontSize'],font_family=GO.opt['nodeFontFamily'])

        # draw edges
        if (GO.opt['edgeColorExtra']!=None) and (GO.opt['edgeColorExtra']!=False) and (GO.opt['edgeColorExtra']!=GO.opt['edgeColor']):
            G1 = genGraph(rmsdMat,edgeListMain,nx)
            nx.draw_networkx_edges(G1,pos=layout,width=GO.opt['edgeWidth'],edge_color=GO.opt['edgeColor'],style=GO.opt['edgeStyle'],alpha=GO.opt['edgeTransparency'],edge_cmap=GO.opt['edgeCMAP'],edge_vmin=GO.opt['edgeVmin'],edge_vmax=GO.opt['edgeVmax'])
            # extra edges
            if cycleNum>0:
                G2 = genGraph(rmsdMat,edgeListExtra,nx)
                nx.draw_networkx_edges(G2,pos=layout,width=GO.opt['edgeWidth'],edge_color=GO.opt['edgeColorExtra'],style=GO.opt['edgeStyle'],alpha=GO.opt['edgeTransparency'],edge_cmap=GO.opt['edgeCMAP'],edge_vmin=GO.opt['edgeVmin'],edge_vmax=GO.opt['edgeVmax'])
        else:
            # draw edges
            nx.draw_networkx_edges(G,pos=layout,width=GO.opt['edgeWidth'],edge_color=GO.opt['edgeColor'],style=GO.opt['edgeStyle'],alpha=GO.opt['edgeTransparency'],edge_cmap=GO.opt['edgeCMAP'],edge_vmin=GO.opt['edgeVmin'],edge_vmax=GO.opt['edgeVmax'])

        # edge lables
        if( (GO.opt['edgeLabels'] != None) and (GO.opt['edgeLabels'] != False) ):
            edgeLabels=genEdgeLabels(G)
            nx.draw_networkx_edge_labels(G,pos=layout,edge_labels=edgeLabels,font_color=GO.opt['edgeFontColor'],font_size=GO.opt['edgeFontSize'],font_family=GO.opt['edgeFontFamily'],font_weight=GO.opt['edgeFontWeight'],label_pos=GO.opt['edgeLabelPos'],alpha=GO.opt['edgeLabelTransp'])

        ### if needed, plot images ###
        if(cmdl.opt['-ligImg'].is_set==True):
#           nodeLabels,nodeLabelLayout=genNodeLabels(rmsdMat,nodeNames,layout,GO.opt)
            plotImages(cmdl,layout,G,nsize,nodeNames,ligImg)
#            nx.draw_networkx_labels(G,pos=nodeLabelLayout,labels=nodeLabels,font_color=GO.opt['nodeFontColor'],font_weight=GO.opt['nodeFontWeight'],font_size=GO.opt['nodeFontSize'],font_family=GO.opt['nodeFontFamily'])

        plt.axis('off')
        plt.savefig(cmdl['-graph'],dpi=GO.opt['dpi'],transparent=True)

main( sys.argv )
