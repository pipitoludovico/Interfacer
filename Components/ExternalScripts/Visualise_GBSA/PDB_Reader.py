import os


class Parser:
    ligand_residues_list = {}
    receptor_residues_list = {}

    def __init__(self):
        self.changeLigdHis()
        self.changeRecHis()
        self.readLigand()
        self.readReceptor()

    @staticmethod
    def changeLigdHis():
        os.system(f"sed -i 's/HSP/HIS/g' ligand.pdb")
        os.system(f"sed -i 's/HSE/HIS/g' ligand.pdb")
        os.system(f"sed -i 's/HSD/HIS/g' ligand.pdb")

    @staticmethod
    def changeRecHis():
        os.system(f"sed -i 's/HSP/HIS/g' receptor.pdb")
        os.system(f"sed -i 's/HSE/HIS/g' receptor.pdb")
        os.system(f"sed -i 's/HSD/HIS/g' receptor.pdb")

    def readLigand(self):
        with open('ligand.pdb', 'r') as ligand:
            for line in ligand.readlines():
                if line.startswith('ATOM'):
                    segid = str(line[72:75].strip())
                    value = line[17:20].strip() + " " + line[22:26].strip()
                    if segid not in self.ligand_residues_list:
                        self.ligand_residues_list[segid] = []
                    if value not in self.ligand_residues_list[segid]:
                        self.ligand_residues_list[segid].append(value)

    def readReceptor(self):
        with open('receptor.pdb', 'r') as ligand:
            for line in ligand.readlines():
                if line.startswith('ATOM'):
                    segid = str(line[72:75].strip())
                    value = line[17:20].strip() + " " + line[22:26].strip()
                    if segid not in self.receptor_residues_list:
                        self.receptor_residues_list[segid] = []
                    if value not in self.receptor_residues_list[segid]:
                        self.receptor_residues_list[segid].append(value)

    def getLists(self):
        return self.ligand_residues_list, self.receptor_residues_list
