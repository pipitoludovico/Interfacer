import os
from pathlib import Path
import concurrent.futures
from subprocess import run


class Refiner:
    def __init__(self, IDs: dict, ROOT_: str):
        self.IDS: dict = IDs
        self.ROOT_: Path = Path(ROOT_)
        self.CleanPDB()

    def Grepping(self, ID_):
        INNERPATH = f'{self.ROOT_ / "structures" / ID_}'
        os.system(f'grep -w "ATOM" {INNERPATH}/{ID_}.pdb > {INNERPATH}/{ID_}_clean.pdb')
        run(f"sed -i 's/HSP/HIS/g' {INNERPATH}/{ID_}_clean.pdb", shell=True)
        run(f"sed -i 's/HSD/HIS/g' {INNERPATH}/{ID_}_clean.pdb", shell=True)
        run(f"sed -i 's/HSE/HIS/g' {INNERPATH}/{ID_}_clean.pdb", shell=True)

    def CleanPDB(self):
        try:
            processes = []
            with concurrent.futures.ProcessPoolExecutor() as executor:
                for ID in self.IDS.keys():
                    future = executor.submit(self.Grepping, ID)
                    processes.append(future)
        except Exception as e:
            print("Grepping failed: ", e)
            exit()

    def GetChains(self):
        for ID in self.IDS.keys():
            with open(str(self.ROOT_ / "structures" / ID / ID) + "_clean.pdb", 'r') as PDBfile:
                for line in PDBfile.readlines():
                    if line.split()[0] == 'ATOM':
                        chainID = line[21]
                        if chainID not in self.IDS[ID]:
                            self.IDS[ID].append(chainID)
        return self.IDS
