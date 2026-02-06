# -*- coding: utf-8 -*-
"""
Created on March 4th 2025

This code takes the input SMILES, kOH, and saturation vapor concentration and speciates based on RACM2B-VCP3
"""

import warnings
import rdkit
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem

def get_racm_roc(smiles_input,koh,log10cstar,cook,phase=None):
    '''
    Function maps input VOCs to RACM2B-VCP3.
    Uses functional group and molecule info from RDKit http://www.rdkit.org/
    Function inputs, for ONE compound (make loop outside this function):
        smiles_input: smiles string (should be canonical for explicit species, some alt added) 
        koh: kOH (in cm3/molec-s)
        log10cstar: (Cstar in micrograms/m3)
        cook: 0 for non-cooking, 1 for cooking
        phase (optional; default value is None; options are 'gas', 'particle', or None)
              This is only used to label species that can exist in both gas and particle phases.
              It does not do any calculations on what phase the species should be in. Semivolatile
              partitioning should be calculated external to this function.
    kOH and C* may not be used if the compound can be mapped without it. 
    RACM2 reference: https://ars.els-cdn.com/content/image/1-s2.0-S1352231012011065-mmc1.pdf
    '''
    
    # R2SMH
    
    # Prep inputs
    if smiles_input == '-':
        unksmiles_msg = (
            f'SMILES {smiles_input} not recognized. Mapping to UNKSMILES.'
        )
        warnings.warn(unksmiles_msg)
        return 'UCM'
    smiles       = Chem.CanonSmiles(smiles_input) # standardize input SMILES string
    m            = Chem.MolFromSmiles(smiles)
    smiles_upper = smiles_input.upper()

    # Count C=C and atoms
    nCdblC  = smiles_upper.count('=C')-smiles_upper.count('O=C')
    if nCdblC < 0:
        nCdblC = 0
    nC      = smiles_upper.count('C')-smiles_upper.count('CL')
    nO      = smiles_upper.count('O')
    nN      = smiles_upper.count('N')
    nSi     = smiles_upper.count('SI')
    nH      = 0
    OVOC    = 0
    for atom in m.GetAtoms():
        nH += atom.GetTotalNumHs()
    # O:C ratio
    if nC > 0:
        OtoC = nO/nC
    
    # Count functional groups (http://rdkit.org/docs/source/rdkit.Chem.Fragments.html)
    nacid     = rdkit.Chem.Fragments.fr_COO(m,countUnique=True)     # carboxylic acid
    nketone   = rdkit.Chem.Fragments.fr_ketone(m,countUnique=True)
    naldehyde = rdkit.Chem.Fragments.fr_aldehyde(m,countUnique=True)
    ncarbonyl = nketone + naldehyde
    nbenzene  = rdkit.Chem.Fragments.fr_benzene(m,countUnique=True)
    nalcohol  = rdkit.Chem.Fragments.fr_Al_OH(m,countUnique=True) + \
                  rdkit.Chem.Fragments.fr_Ar_OH(m,countUnique=True)      # aliphatic and aromatic
    nfuran    = rdkit.Chem.Fragments.fr_furan(m,countUnique=True) # number of furan rings
    nester =rdkit.Chem.Fragments.fr_ester(m,countUnique=True) # number of esters
    nether =rdkit.Chem.Fragments.fr_ether(m,countUnique=True) # number of ethers
    nether =nether-nester
    nepox=rdkit.Chem.Fragments.fr_epoxide(m,countUnique=True) # number of epoxide
    ch3 = Chem.MolFromSmarts('[CH3]')
    ntertH=len(m.GetSubstructMatches(ch3))
    ch2oh = Chem.MolFromSmarts('[CH2O]')
    npriOH=len(m.GetSubstructMatches(ch2oh))
    choh = Chem.MolFromSmarts('[CHO]')
    nsecOH=len(m.GetSubstructMatches(choh))

    # gotring variable is never used so this could be removed
    #     it may be useful to keep it here commented out in case it is needed in the future
    #for atom in range(len(m.GetAtoms())):
        #if m.GetAtomWithIdx(atom).IsInRing():
        #    gotring = 1 # 0 = no ring, 1 = ring
        #    break
        #else: gotring = 0
    nnitrate = smiles_upper.count('O[N+](=O)[O-]') + smiles_upper.count('O[N+]([O-])=O')
    nsilox = smiles_upper.count('O[SI]')
    nperoxide = smiles_upper.count('COO') +  smiles_upper.count('OOC') - smiles_upper.count('COOC')
    #namine = smiles_upper.count('CN') + smiles_upper.count('NC') + smiles_upper.count('C(N') + smiles_upper.count('N(C')
    tfmonoterpene = (nC == 10 and nH == 18 and nO == 1 and smiles != 'CCCCCCCC=CC=O') or (nC == 10 and nH == 16 and smiles != 'CCCCCC=CC=CC=O')

    # Mapper is for ROC only and not elemental carbon
    if   ( nC <= 0 ):                 mechspecies = ['UCM','UCM']
    elif ( smiles == '[C]' ):         mechspecies = ['UCM','UCM']

     # Map CO to UNKCRACMM; CO will be mapped to SLOWROC if not handled explicitly
    elif ( smiles == '[C-]#[O+]' ):   mechspecies = ['UCM','UCM']
    # The same applies to CO2
    elif ( smiles == 'O=C=O' ):       mechspecies = ['UCM','UCM']

    # Explicit species
    elif ( smiles == 'CC=O' ):        mechspecies=  ['HC15', 'CCHO']   # acetaldehyde
    elif ( smiles == 'C#C' ):         mechspecies = ['HC37', 'ACETYLENE']   # acetylene
    elif ( smiles == 'CC(C)=O' ):     mechspecies = ['HC18','ACET']   # acetone
    elif ( nC==6 and nH==6 and nO==0 and nbenzene==1 ):
                                      mechspecies = ['HC38','BENZENE']   # benzene
    elif ( smiles == 'C'  ):          mechspecies = ['HC01','ECH4']  # methane
    elif ( smiles == 'CCO'):          mechspecies = ['HC48','ETHANOL']   # ethanol
    elif ( smiles == 'C=C'):          mechspecies = ['HC07','ETHENE']   # ethene aka ethylene
    elif ( smiles == 'OCCO'):         mechspecies = ['HC49','ETEG']  # ethylene glycol
    elif ( smiles == 'CC' ):          mechspecies = ['HC02','ALK1']   # ethane
    elif ( smiles == 'C=O'):          mechspecies = ['HC14','HCHO']  # formaldehyde
    elif ( smiles == 'C=CC(=C)C' ):   mechspecies = ['HC10','ISOPRENE']   # isoprene (output of Chem.CanonSmiles)
    elif ( smiles == 'CO'):           mechspecies = ['HC21','MEOH'] # methanol
    elif ( smiles == 'O=CO'):         mechspecies = ['HC30','HCOOH']  # formic acid
    elif ( smiles == 'CC(=O)O'):      mechspecies = ['HC31','ORA2']  # acetic acid / HC31
    elif ( smiles == 'C=Cc1ccccc1'):  mechspecies = ['HC47','STYRENES']   # styrene (added in CRACMM2)
    #elif ( smiles == 'CCC(C)=O' ):   mechspecies = ['HC19','MEK']   # methyl ethyl ketone
    #elif ( smiles == 'C=CC(C)=O' ):  mechspecies = ['HC28','MVK']   # methly vinyl ketone
    elif ( smiles == 'Cc1ccccc1' ):   mechspecies = ['HC41','TOLUENE']   # toluene
    elif ( smiles == 'C=CC=C' ):      mechspecies = ['HC46','DIENES'] # 1,3 butadiene 
    elif ( smiles == 'Cc1cccc(C)c1' ):         mechspecies = ['HC42','M-XYLENE'] # m-xylene
    elif ( smiles == 'Cc1ccccc1C' ):           mechspecies = ['HC43','O-XYLENE'] # o-xylene
    elif ( smiles == 'Cc1ccc(C)cc1' ):         mechspecies = ['HC44','P-XYLENE'] # p-xylene
    elif ( smiles == 'C(C(CO)O)O'):            mechspecies = ['HC53','GLYC'] #glycerol
    elif ( smiles =='C=C(C)C=O'):              mechspecies = ['HC27','MACR'] #Methacrolein
    elif (cook==1 and smiles == 'CCC=O'  ):        mechspecies = ['HC60','CALD']      # Propanal
    elif (cook==1 and smiles == 'CCCC=O'  ):       mechspecies = ['HC61','CALD']       # Butanal
    elif (cook==1 and smiles == 'CCCCC=O'  ):      mechspecies = ['HC62','CALD']       # Pentanal
    elif (cook==1 and smiles == 'CCCCCC=O'  ):     mechspecies = ['HC63','CALD']       # Hexanal
    elif (cook==1 and smiles == 'CCCCCCC=O'  ):    mechspecies = ['HC64','CALD']       # Heptanal
    elif (cook==1 and smiles == 'CCCCCCCC=O'  ):   mechspecies = ['HC65','CALD']       # Octanal
    elif (cook==1 and smiles == 'CCCCCCCCC=O'  ):  mechspecies = ['HC66','CALD']       # Nonanal
    elif ( smiles == 'CC(C)O'       ):         mechspecies = ['HC51','IPOH']      # Isopropyl alcohol
    elif ( smiles == 'OCc1ccccc1'):            mechspecies = ['HC83','BENZOH']      # Benzyl Alcohol
    elif ( smiles == 'FC(F)(F)c1ccc(Cl)cc1'):  mechspecies = ['HC58','PCBTF'] # parachlorobenzotrifluoride PCBTF
    elif ( smiles == 'Clc1ccc(Cl)cc1'):        mechspecies = ['HC59','PDCBZ'] # p-dichlorobenzene PDCBZ
    elif ( smiles == 'CC=C'):                  mechspecies= ['HC36','PROPYLENE'] # propylene 
    elif ( smiles == 'C(=O)C=O'):              mechspecies= ['HC22','GLY'] # glyoxal
    elif ( smiles == 'CC(=O)C=O'):             mechspecies= ['HC23','MGLY'] # methyl glyoxal 
    elif ( nC == 3 and nH == 8 and nO == 0 and nSi == 0 and nN == 0):  mechspecies= ['HC45','PROPANE'] # propane
    elif ( nC == 4 and nH == 10 and nO == 0 and nSi == 0 and nN == 0): mechspecies= ['HC39','BUTANES'] # butanes
    elif ( nC == 5 and nH == 12 and nO == 0 and nSi == 0 and nN == 0): mechspecies= ['HC40','PENTANES'] # pentanes
    # Glyoxal and glycoaldehyde (here due to solubility alt mappings: ACD, ETEG)
    elif ( nC==2 and nO==2 and naldehyde>=1 ): mechspecies = ['HC22','GLY']   # glyoxal
    # Propylene glycol
    #elif ( nC==3 and nO==2 and nalcohol==2):  mechspecies = ['HC52','PROG'] # propylene glycol and other 3 carbon dialcohols                                
    # Low reactivity species (new to v0.1)
    elif ( koh < 3.5e-13 ):                    mechspecies = ['HC57','NROG'] # low reactivity gas

    # SVOC species binned
    elif ( log10cstar < -0.5): # C* bin centered on 0.1 ug/m3
            mechspecies = ['HC88','SVOC1']
    elif ( log10cstar < 0.5 ): # C* bin centered on 1
            mechspecies = ['HC87','SVOC2']
    elif ( log10cstar < 1.5): # C* bin centered on 10^1
            mechspecies = ['HC86','SVOC3']
    elif ( log10cstar < 2.5): # C* bin centered on 10^2
            mechspecies = ['HC85','SVOC4']

    # Lumped terpene species 
##########################################################
    elif ( nC == 15 and nH == 24 ) and ( nCdblC>=1 ): mechspecies = ['HC77','SESQ'] # sesquiterpenes
    elif ( tfmonoterpene and nCdblC==1 ): mechspecies = ['HC78','API'] # a-pinene monoterpenes
    elif ( tfmonoterpene and nCdblC>=2 ): mechspecies = ['HC79','LIM'] # limonene monoterpenes
    elif ( tfmonoterpene and nCdblC==0):  mechspecies = ['HC80','EUC']      # Eucalyptol

    # Multi-ring aromatics (PAH and NAPH can be collapsed together if necessary)  
    # elif ( nbenzene >= 1 and log10cstar < 3.5 and (nO/nC) == 0 ): mechspecies = 'PAH'  # PAH and other lower-volatility aromatics (v0.1)
    #elif ( nbenzene > 1 and nO/nC == 0 ): mechspecies = 'HC12' # Naphthalene-like, PAH with 2 rings
    # Single-ring aromatics (excluding explicit species)
    elif ( nbenzene > 0 ): # Single-ring aromatics
        if ( log10cstar < 5.5):             mechspecies = ['HC76','IVOCP5-ARO'] # C* bin centered on 10^5 (v0.1)
        elif ( log10cstar < 6.5):           mechspecies = ['HC75','IVOCP6-ARO'] # C* bin centered on 10^6 (v0.1)
        elif ( naldehyde > 0 ):              mechspecies = ['HC17','BALD']     # Benzaldehyde and arom. aldehydes
        elif ( nC>=6 and nalcohol>=2 ):      mechspecies = ['HC34','CATECHOLS']      # catechols
        elif ( nC>7 and nalcohol==1 ):       mechspecies = ['HC33','XYLENEOLS']      # xylenols
        elif ( nC==7 and nalcohol==1 ):      mechspecies = ['HC26','CSL' ]      # cresol
        elif ( nC==6 and nalcohol==1 ):      mechspecies = ['HC25','PHEN']     # phenol   
        # any single-ring aromatics that have not been mapped by rules above
        elif ( koh < 1.3E-11 ):              mechspecies = ['HC12','ARO1']      # other aromatics (CRACMM2)
        else:                                mechspecies = ['HC13','ARO2']                   
  
    elif ( log10cstar < 3.5 ): # C* bin centered on 1000 ug/m3
        if OtoC > 0.1: OVOC=1
        else: mechspecies = ['HC74','IVOCP3-ALK']
    elif ( log10cstar < 4.5): # C* bin centered on 10^4
        if OtoC > 0.1: OVOC=1
        else: mechspecies = ['HC73','IVOCP4-ALK']
    elif ( log10cstar < 5.5 ): # C* bin centered on 10^5
        if OtoC > 0.05: OVOC=1
        else: mechspecies = ['HC72','IVOCP5-ALK']
    elif ( log10cstar < 6.5 ): # C* bin centered on 10^6
        if OtoC > 0.05: OVOC=1
        else: mechspecies = ['HC71','IVOCP6-ALK']

    # Species with double bonds, not aromatic
    elif ( nCdblC>=2 ):                                       mechspecies = ['HC46','DIENES']    #  dienes (nC>=4 gauranteed)    
    elif ( nCdblC>=1 and nC>=3 and naldehyde>=1 and cook==1 ):mechspecies = ['HC67','CUALD']  # unsaturated cooking aldehydes 
    elif ( nCdblC>=1 and nC>=3 and naldehyde==1 ):            mechspecies = ['HC29','UALD']     # unsaturated aldehydes 
    elif ( nCdblC==1 and ncarbonyl>=2 ):                      mechspecies = ['HC70','DCB1']     # unsaturated dicarbonyls         
    elif ( nCdblC>=1 and nketone>=1 ):                        mechspecies = ['HC28','MVK']     # MVK and unsaturated ketones

    # map OLI and OLT
    elif ( nCdblC==1 ):
        atoms = [a for a in m.GetAtoms()]
        bonds = [b for b in m.GetBonds()]
        bondtype = [b.GetBondType() for b in bonds]
        #        double bond at end                         bond isn't part of a ring    the double bond is to a carbon
        if bondtype[0]==rdkit.Chem.rdchem.BondType.DOUBLE and not bonds[0].IsInRing() and atoms[0].GetSymbol()=='C':
                                               mechspecies = ['HC08','OLE1']
        elif bondtype[-1]==rdkit.Chem.rdchem.BondType.DOUBLE and not bonds[-1].IsInRing() and atoms[-1].GetSymbol()=='C':
                                               mechspecies = ['HC08','OLE1']
        else:                                  mechspecies = ['HC09','OLE2']
    
    elif (nsilox > 0):  
        if (smiles == 'C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1' ): mechspecies = ['HC54','D4SILOXANE']  #d4 siloxane
        elif (smiles == 'C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1' ): mechspecies = ['HC55','D5SILOXANE']  #d5 siloxane
        else : mechspecies = ['HC56','OTHSILOXANE']  #other siloxanes

    
    # Oxygenated species without double bonds (mapped in order of decreasing koh)
    elif ( naldehyde>=1 and nketone>=1 ):               mechspecies = ['HC24','MGLY'] # methylglyoxal and similar like C4H6O2 (1.5e-11)
    elif (nC>9 and naldehyde>=1 and cook==1):           mechspecies = ['HC68','CALD'] 
    elif ( naldehyde>=1 ):                              mechspecies = ['HC16','ALD']  # higher aldehydes (C>3) (1.98e-11)   
    elif ( nketone>=1 and koh >= 5E-13 and koh <5E-12): mechspecies = ['HC19','MEK']  # ketones 
    elif ( nketone>=1 and koh >= 5E-12):                mechspecies = ['HC20','KET']  # ketones 
    elif ( npriOH==2):                                  mechspecies = ['HC49','ETEG']
    elif ( nalcohol==2):                                mechspecies = ['HC52','PROPG']
    elif ( nalcohol>=3):                                mechspecies = ['HC53','GLYC']
    elif ( npriOH>=1):                                  mechspecies = ['HC69','ROH']  # C3 and higher primary alcohols 
    elif ( nalcohol==1):                                mechspecies = ['HC51','IPOH']  # C3 and higher secondary alcohols 
    elif ( nacid>=1 ):                                  mechspecies = ['HC32','ORA2'] # acetic acid and higher acids (C>=2) (6.5e-13) 

  
    # HC Series, koh in cm3/s, 298 K, 1 atm
    elif ( koh <= 3.4E-13 ):                    mechspecies = ['HC02','ALK1']
    elif ( koh >= 3.4E-13 and koh < 1.7E-12 ):  mechspecies = ['HC03','ALK2']  
    elif ( koh >= 1.7E-12 and koh < 3.4E-12 ):  mechspecies = ['HC04','ALK3'] 
    elif ( koh >= 3.4E-12 and koh < 6.8E-12 ):  mechspecies = ['HC05','ALK4'] 
    elif ( koh >= 6.8E-12 ):                    mechspecies = ['HC06','ALK5']

    
    if (OVOC == 1):
      if (nC < 5 and npriOH==2): mechspecies = ['HC49','ETEG']
      elif (nC < 5 and nalcohol==2): mechspecies = ['HC52','PROPG']
      elif (nC < 5 and nalcohol>=3): mechspecies = ['HC53','GLYC']
      elif (ntertH < 3 and nalcohol>0 and nC > 4 and nC < 8 and nether+nester>0): mechspecies = ['HC81','CARB'] 
      elif (ntertH < 3 and nalcohol>0 and nC > 7 and nether+nester>0): mechspecies = ['HC82','BCARB']
      elif (ntertH < 3 and nalcohol<1 and nC > 4 and nether+nester+nepox>1): mechspecies = ['HC84','DPGMEA']     
      elif ( log10cstar < 3.5 ): mechspecies = ['HC74','IVOCP3-ALK'] # C* bin centered on 1000 ug/m3       
      elif ( log10cstar < 4.5 ): mechspecies = ['HC73','IVOCP4-ALK'] # C* bin centered on 10^4        
      elif ( log10cstar < 5.5 ): mechspecies = ['HC72','IVOCP5-ALK'] # C* bin centered on 10^5       
      elif ( log10cstar < 6.5 ): mechspecies = ['HC71','IVOCP6-ALK']# C* bin centered on 10^6
         

    bin=mechspecies[0]
    sprac=mechspecies[1]
    return bin,sprac
    # end of function
