#all the modules used will have to pip install rdkit, csv, pysmiles
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import AllChem
import random
import csv
from pysmiles import read_smiles
from pysmiles import write_smiles
from rdkit.Chem import BondType

sums=0
#list of balanced complexes calculated

list_balanced_cat=[]
#ligands from csv file stored in this list
ligands2=[]
charge_list=[]

#for name of unsorted catalyst list

stringnum = str(random.randint(1,200000))
passed=False

#defining the class that has all the code and functions could not be fucked to rename it all
class catalyst_create:


 def init__(self,number_cycles):
  self.number_cycles=number_cycles
  #number of times/number of unsorted complexes generated. The number of sorted complexes will be way less than this

 def Catalyst_generation(number_cycles):
  catlist=[]
  #list for unsorted complex
  for x in range(number_cycles):
   #getting metals from text file
   met_file= open('metals.txt','r')
   met=met_file.read()
   met=met.split()
   #where ligands are stored in unsorted complex phase
   ligands=[]
   #gets ligands from csv file
   with open('ligands.csv', 'r', newline='') as csvfile:
           reader = csv.reader(csvfile)
           data = list(reader)
           for row in data:
            charge_got = row[0].split(' ')
            ligands.append(charge_got[0])

   string_lig=''
   lig=[]
   lig2=[]
   list_lines=[]
   #the maximum number of ligands sorounding the metal center
   max_coordination=7
   #random number of ligands added
   for l in range(random.randint(1,max_coordination)):
    lig.append(ligands[random.randint(1,len(ligands))-1])
    string_lig=string_lig+str(lig[l-1])
    lig2.append(string_lig)



   for x in range(len(lig2)):
    #random complex is created with a random metal and random ligand enviroment again could not be fucked to rename anything

    cat = str(met[random.randint(1,len(met)-1)]+lig2[random.randint(1,len(lig2))-1])


   catlist.append(cat)
#writes unsorted complexes to a text file not the most efficient method but it was useful for what i was originally doing
  with open('catcreatefile'+stringnum+'.txt','w+') as cat_create_file:
   for catn in range(len(catlist)):
    cat_create_file.write("{}\n".format(catlist[catn]))

#makes sure the complexes are charge balanced and have a sensible coordination environment aka most annoying part
 def Catalyst_sort():
  sum_l=[]
  sum_ml=[]
  with open('catcreatefile'+stringnum+'.txt','r') as cat_read_file:
   cat_read= (cat_read_file.read()).split()
#again reads ligands from csv file but also gets the charges as well
  with open('ligands.csv', 'r', newline='') as csvfile:
   reader = csv.reader(csvfile)
   data = list(reader)
   for row in data:
    charge_got = row[0].split(' ')
    print(charge_got)
    print(charge_got[0])
    ligands2.append(charge_got[0])
    charge_list.append(int(charge_got[1]))
    print(charge_list)
  with open('metals.txt','r') as metals_file:
      metals = metals_file.read()
      metals = metals.split()
  for c in range(len(cat_read)):
#had a bug rdkit didnt like the way i wrote smiles so had to give it to another module then use the output they gave for rdkit, stupidest solution of my life
# this part translates the smiles into a graph yady yadda makes chemistry happen
   mol = read_smiles(cat_read[c])
   cat_smiles = write_smiles(mol)

   mol = Chem.MolFromSmiles(cat_smiles)
   total_charge_list=[]
   negative_sum=[]
   found_ligands=[]
#checks if a ligand from the csv file is in the randomly generated complex
   for item in ligands2:
           count = cat_read[c].count(item)
           found_ligands.extend([item] * count)

#calculates the overall charge of the ligands from charge data provided by the csv file
   for item in found_ligands:

         negative_sum.append(charge_list[ligands2.index(item)]*found_ligands.count(item))
#only allows for coordination environments of 2,4,6 you can change it to suit your needs
   acceptable_coordination=[2,4,6]

   if len(found_ligands) in acceptable_coordination:
    passed=True
#this bit uses smiles data of the metal to get the metal charge
   for atom in mol.GetAtoms():
             total_charge_list.append(atom.GetFormalCharge())

   total_charge=sum(total_charge_list)+sum(negative_sum)
#if the coordination environment is good and the complex is balanced then will write complex to a txt file
   if total_charge==0 and passed==True:
    # Convert covalent bonds to dative bonds for non-halogen ligands
    mol_modified = convert_to_dative_bonds(mol, found_ligands)
    if mol_modified:
        # Convert back to SMILES with dative bonds
        cat_smiles_dative = Chem.MolToSmiles(mol_modified)
        list_balanced_cat.append(cat_smiles_dative)

    with open('balancedcatalyst.txt','w') as balanced_file:
     for bal in range(len(list_balanced_cat)):
      balanced = balanced_file.write("{}\n".format(list_balanced_cat[bal]))

 #displays complex as an image dosnt work very well the connectivities are all fucked as well as the number of hydrogens wasnt really part of my original project but thought you might wanna know how to do it

 def display_complex():
     with open('balancedcatalyst.txt','r') as balanced_file:
         read_bal=balanced_file.read()
         read_bal=read_bal.split()
     for c in range(len(read_bal)):
      mol_bal = Chem.MolFromSmiles(read_bal[c])
      Chem.AddHs(mol_bal)
      img = Draw.MolsToGridImage([mol_bal])
      img.save('complex_image'+str(c)+'.png')

# New function to convert covalent bonds to dative bonds for non-halogen ligands
def convert_to_dative_bonds(mol, found_ligands):
    """
    Convert covalent bonds between metal and non-halogen ligands to dative bonds
    """
    # Halogens that should keep covalent bonds
    halogens = ['Cl', 'F', 'Br', 'I']

    # Create editable molecule
    rw_mol = Chem.RWMol(mol)

    # Find metal atoms (typically transition metals with positive charge or low electronegativity)
    metal_atoms = []
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() in [3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,39,40,41,42,43,44,45,46,47,48,49,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103]:
            metal_atoms.append(atom.GetIdx())

    if not metal_atoms:
        return None

    # Convert bonds to dative for non-halogen ligands
    bonds_to_modify = []
    for bond in rw_mol.GetBonds():
        begin_idx = bond.GetBeginAtom().GetIdx()
        end_idx = bond.GetEndAtom().GetIdx()

        # Check if this is a metal-ligand bond
        if (begin_idx in metal_atoms) or (end_idx in metal_atoms):
            # Determine which atom is the metal and which is the ligand
            if begin_idx in metal_atoms:
                metal_idx = begin_idx
                ligand_idx = end_idx
            else:
                metal_idx = end_idx
                ligand_idx = begin_idx

            # Get the ligand symbol
            ligand_atom = rw_mol.GetAtomWithIdx(ligand_idx)
            ligand_symbol = ligand_atom.GetSymbol()

            # Check if ligand is a non-halogen
            is_halogen = False
            for halogen in halogens:
                if ligand_symbol == halogen or any(halogen in lig for lig in found_ligands if ligand_symbol in lig):
                    is_halogen = True
                    break

            # Convert to dative bond if not halogen
            if not is_halogen:
                bonds_to_modify.append((metal_idx, ligand_idx, bond.GetIdx()))

    # Remove old bonds and add dative bonds
    for metal_idx, ligand_idx, bond_idx in bonds_to_modify:
        # Remove the existing bond
        rw_mol.RemoveBond(metal_idx, ligand_idx)
        # Add dative bond (pointing from ligand to metal)
        rw_mol.AddBond(ligand_idx, metal_idx, BondType.DATIVE)

    return rw_mol.GetMol()


catalyst_create.Catalyst_generation(4000)
catalyst_create.Catalyst_sort()

check_display = input('do you want to display png images of all balanced complexes? y/n')
if check_display=='y':
 catalyst_create.display_complex()