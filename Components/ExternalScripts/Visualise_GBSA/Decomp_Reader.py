import os


class DecompReader:

    def __init__(self):
        self.rec_decomp_results = []
        self.lig_decomp_results = []
        root = os.getcwd()
        self.datFile = f'{root}/FINAL_DECOMP_MMPBSA.dat'
        self.CleanDecompHis()
        self.WriteTotalFromDecomp()
        self.TakeLigandResults()
        self.TakeReceptorResults()

    def CleanDecompHis(self):
        os.system(f"sed -i 's/HSP/HIS/g' {self.datFile}")
        os.system(f"sed -i 's/HSE/HIS/g' {self.datFile}")
        os.system(f"sed -i 's/HSD/HIS/g' {self.datFile}")

    def WriteTotalFromDecomp(self):
        decomp_purged = open("total_purged.csv", 'w')
        with open(f'{self.datFile}', 'r') as f:
            for purgedLine in f.readlines():
                decomp_purged.writelines(purgedLine)
                if purgedLine == "\n" or purgedLine == " ":
                    break

    def TakeReceptorResults(self):
        with open('total_purged.csv', 'r') as f:
            after_header = f.readlines()[7:]
            for lines in after_header:
                if lines == "" or lines == " " or lines == "\n":
                    break
                if lines.split()[1].split(',')[1] == "R":
                    total_energy = float(lines.split(",")[-3])
                    self.rec_decomp_results.append(total_energy)

    def TakeLigandResults(self):
        with open('total_purged.csv', 'r') as f:
            after_header = f.readlines()[7:]
            for lines in after_header:
                if lines == "" or lines == " " or lines == "\n":
                    break
                if lines.split()[1].split(',')[1] == "L":
                    total_energy = float(lines.split(",")[-3])
                    self.lig_decomp_results.append(total_energy)

    def GetDecompResults(self):
        return self.lig_decomp_results, self.rec_decomp_results
