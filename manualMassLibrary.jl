push!(LOAD_PATH, pwd())
using MasslistFunctions

CO2 = createCompound(C=1, O=2)
NH3 = createCompound(N=1, H=3)
PYRIDINE = createCompound(N=1,H=5,C=5)
ACETONE = createCompound(C=3,H=6,O=1)
APINENE = createCompound(C=10, H=16)
BCARY = createCompound(C=15, H=24)
ISOPRENE = createCompound(C=5, H=8)
HEXANONE =  createCompound(C=6, H=12, O=1)
PINONALDEHYDE = createCompound(C=10, H=16, O=2)
PINONALDEHYDEPAN = createCompound(C=10, H=15, O=6, N=1)
PINONICACID = createCompound(C=10, H=16, O=3) #C10H16O3
PINICACID = createCompound(C=9, H=14, O=4) #C9H14O4
ACETIC = createCompound(C=2, H=4, O=2)
ACETICFRAG = createCompound(C=2, H=2, O=1)
NH3 = createCompound(N=1, H=3)
DMA = createCompound(C=2, H=7, N=1)
TMA = createCompound(C=3, H=9, N=1)
TMB = createCompound(C=9, H=12)

ACETONITRILE = createCompound(C=2, H=3, N=1)

OrgNitratNO = createCompound(C=10, H=15, O=5, N=1)
NORPINONALDEHYDE = createCompound(C=9, H=14, O=2)
NORPINONALDEHYDEPAN = createCompound(C=9, H=13, O=6, N=1)
H3O = createCompound(H=2, O=1)
H3OH2O = createCompound(H=4, O=2)
H3OH2OH2O = createCompound(H=6, O=3)

########## Sulphur ###############
DMS = createCompound(C=2, H=6, S=1)
DMSO = createCompound(C=2, H=6, S=1, O=1)
DMSO2 = createCompound(C=2, H=6, S=1, O=2)
MSIA = createCompound(C=1, H=4, S=1, O=2)
MSA = createCompound(C=1, H=4, S=1, O=3)

########## BCARY Specific ########
#Fast
C15H22O2 = createCompound(C=15, H=22, O=2)
C15H24O2 = createCompound(C=15, H=24, O=2)
C15H24O3 = createCompound(C=15, H=24, O=3)
C15H26O4 = createCompound(C=15, H=26, O=4)
#Slow / Sticky
C14H20O3 = createCompound(C=14, H=20, O=3)
C13H20O4 = createCompound(C=13, H=20, O=4)
C14H22O4 = createCompound(C=14, H=22, O=4)
C15H22O4 = createCompound(C=15, H=22, O=4)
C15H24O4 = createCompound(C=15, H=24, O=4)
C14H25NO4 = createCompound(C=14, H=25, O=4, N=1)
C15H29NO6 = createCompound(C=15, H=29, O=6, N=1)
#Nitrates
C15H23NO4 = createCompound(C=15, H=23, O=4, N=1)
C13H19NO6 = createCompound(C=13, H=19, O=6, N=1)
C15H25NO4 = createCompound(C=15, H=25, O=4, N=1)
C15H21NO5 = createCompound(C=15, H=21, O=5, N=1)
C15H23NO5 = createCompound(C=15, H=23, O=5, N=1)
C15H23NO6 = createCompound(C=15, H=23, O=6, N=1)
C15H17NO7 = createCompound(C=15, H=17, O=7, N=1)

###### END BCARY ###############


###### NAPHTHA #################

NAPHTHA = createCompound(C=10,H=8)


###### END NAPHTHA #############


compound2 = createCompound(C=11, H=12, O=1)
compound3 = createCompound(C=4, H=6, O=1)
compound4 = createCompound(C=5, H=8, O=1)
compound5 = createCompound(C=4, H=6, O=2)
compound6 = createCompound(C=5, H=10, O=1)
compound7 = createCompound(C=6, H=10, O=1)
compound8 = createCompound(C=7, H=12, O=1)
compound9 = createCompound(C=6, H=10, O=2)
compound10 = createCompound(C=5, H=8, O=3)

# Dimers

dimer1 = createCompound(C=19, H=28, O=9)
dimer2 = createCompound(C=18, H=28, O=10)
dimer3 = createCompound(C=20, H=32, O=9)
dimer4 = createCompound(C=19, H=28, O=11)
dimer5 = createCompound(C=19, H=26, O=19)

dimer6 = createCompound(C=20, H=31, N=1, O=7)
dimer7 = createCompound(C=19, H=29, N=1, O=8)
dimer8 = createCompound(C=20, H=31, N=1, O=9)
dimer9 = createCompound(C=20, H=31, N=1, O=13)
