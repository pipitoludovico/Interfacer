import os.path
from os import cpu_count, chdir
from multiprocessing import Pool, Manager
import pandas as pd
from ..Interfacer.VMD import GetVDWcontacts

features = ['DES', 'SUM', "AVG_S", "STD_S", "AVG_GBSA", "SD_GBSA", "AVG_RMSD", "SD_RMSD",
            'GBSA', 'NumContacts']


class FeaturizerClass:
    def __init__(self, dbDict: dict, root: str):
        self.dbDict = dbDict
        self.root = root
        os.chdir(self.root)
        self.datFile = "FINAL_DECOMP_MMPBSA.dat"
        self.dataframe = pd.DataFrame(columns=[*features], dtype=float)

    def Featurize(self, pdbID: str, chainsID: list, shared_list: list) -> None:
        agChain = ''
        abChain = ''
        if len(chainsID) == 2:
            abChain = chainsID[0]
            agChain = chainsID[1]
        if len(chainsID) > 2:
            abChain = chainsID[0] + " " + chainsID[1]
            agChain = chainsID[2]
        chdir("structures/" + pdbID)
        if not os.path.exists('contacts.int'):
            GetVDWcontacts(abChain, agChain)
        rec_results = self.GetDecomp()
        df: pd.DataFrame = self.dataframe.copy()
        self.update_dataframe(df, pdbID, rec_results)
        shared_list.append(df)  # Append modified dataframe to the shared list
        chdir(self.root)

    def update_dataframe(self, df_, pdbID: str, rec_results: dict) -> None:
        """Updated the dataframe which will be then transferred to the shared memory
        and concatenated to make column-wise observation for conserved destibilizing interactions"""
        interfaceRes = []

        with open('minimal_interface/selection_i0.pdb', 'r') as interfaceFile:
            for line in interfaceFile.readlines():
                if len(line.split()) > 2:
                    if float(line.split()[10]) > 0.15:
                        residueID = line.split()[3] + "-" + line[23:27].strip()
                        if residueID not in interfaceRes:
                            interfaceRes.append(residueID)
        for resname, energy in rec_results.items():
            if float(energy) > 0.005:
                if resname in interfaceRes:
                    if resname not in df_.columns:
                        df_[resname] = None
                        df_.at[pdbID, resname] = energy

        with open('gbsa/results_mmgbsa.dat', 'r') as GBSA:
            for line in GBSA.readlines():
                if "DELTA TOTAL" in line:
                    GBSA = float(line.split()[2])
                    df_.at[pdbID, 'GBSA'] = GBSA
                    break

        with open('contacts.int', 'r') as NumContacts:
            for line in NumContacts.readlines():
                numContact = int(line.split(",")[0])
                df_.at[pdbID, 'NumContacts'] = numContact
        with open('../../DynamicScores_amber/DynamicScores.csv', 'r') as DESscores:
            for line in DESscores.readlines():
                if pdbID in line:
                    df_.at[pdbID, 'DES'] = line.split(",")[1]
                    df_.at[pdbID, 'SUM'] = line.split(",")[2]
                    df_.at[pdbID, 'AVG_S'] = line.split(",")[3]
                    df_.at[pdbID, 'STD_S'] = line.split(",")[4]
                    df_.at[pdbID, 'AVG_GBSA'] = line.split(",")[5]
                    df_.at[pdbID, 'SD_GBSA'] = line.split(",")[6]
                    df_.at[pdbID, 'AVG_RMSD'] = line.split(",")[7]
                    df_.at[pdbID, 'SD_RMSD'] = line.split(",")[8]
        print(df_.loc[pdbID])

    @staticmethod
    def GetDecomp() -> dict:
        """We are interested in finding what are the destabilising contributions
        in our Ab-Ag binding patch. The data will be stored in a dictionary to be comfortably
        transferred to a dataframe laster on for ML"""
        rec_decomp_results = {}

        def GetResults():
            with open('gbsa/total_purged.csv', 'r') as f:
                after_header = f.readlines()[7:]
                for lines in after_header:
                    if lines == "" or lines == " " or lines == "\n":
                        break
                    if lines.split()[1].split(',')[1] == "R":
                        resname = lines.split()[2]
                        resnum = lines[14:17].strip()
                        total_energy = float(lines.split(",")[-3])
                        if lines.split()[3] not in rec_decomp_results:
                            rec_decomp_results[resname + "-" + resnum] = total_energy

        GetResults()
        return rec_decomp_results

    def ParallelFeaturize(self) -> None:
        """Everything is conceptualized to be scalable. The stronger the machine, the more jobs
        it'll be able to run on parallel"""
        cpuUnits = int(cpu_count() // 2)
        manager = Manager()
        shared_list = manager.list()  # Create a shared list
        with Pool(processes=cpuUnits) as p:
            for pdbID, chains in self.dbDict.items():
                p.apply_async(self.Featurize, args=(pdbID, chains, shared_list,))
            p.close()
            p.join()
        for df in shared_list:
            self.dataframe = pd.concat([self.dataframe, df])

    def GetDataAsDataframe(self) -> pd.DataFrame:
        return self.dataframe
