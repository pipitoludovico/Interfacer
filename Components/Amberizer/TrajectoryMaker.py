import os.path

from multiprocessing import Pool, cpu_count
from subprocess import DEVNULL
from .ComplexMaker import *
from ..Interfacer.VMD import GetCloseSelection, RestoreChains
from ..ReceptorTools.PrepareReceptor_AMBER import PrepareProtein
from ..LigandTools.Tleap import RunTleap, MakePrmtop

external_scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'ExternalScripts')
OPENMM_PIPELINE = os.path.join(external_scripts_dir, 'openmm_pipeline.py')
DYNAMIC_SCORER = os.path.join(external_scripts_dir, 'DynamicScorer.py')
VISUALIZER = os.path.join(external_scripts_dir, 'Visualise_GBSA/Visualise_GBSA.py')


def RunMMPBSA():
    run("$AMBERHOME/bin/MMPBSA.py -i mmgbsa.in -o results_mmgbsa.dat -cp gbsa/complex.prmtop -rp gbsa/receptor.prmtop -lp gbsa/ligand.prmtop -y gbsa/complex.dcd -eo gbsa.csv > logs/GBSA.log 2>&1",
        check=True, shell=True, stdout=DEVNULL)


def GetCrystalCoords(filePath):
    with open(f'{filePath}', 'r') as f:
        lines = f.readlines()
    x_coords, y_coords, z_coords = [], [], []
    for line in lines:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            x_coords.append(x)
            y_coords.append(y)
            z_coords.append(z)

    x_range = np.max(x_coords) - np.min(x_coords)
    y_range = np.max(y_coords) - np.min(y_coords)
    z_range = np.max(z_coords) - np.min(z_coords)
    return x_range, y_range, z_range


class TrajectoryMaker:
    def __init__(self, dbDict):
        self.dbDict = dbDict
        self.ROOT = os.getcwd()

    @staticmethod
    def TrajectorizePDB(chains) -> None:
        agChain = ''
        abChain = ''
        if len(chains) == 2:
            abChain = chains[0]
            agChain = chains[1]
        if len(chains) > 2:
            abChain = chains[0] + " " + chains[1]
            agChain = chains[2]

        if not os.path.exists('./system') or not os.path.exists('solvated.pdb'):
            FILES = ["complex_final.pdb", "receptor_final.pdb", "ligand_final.pdb"]
            for splitFile in FILES:
                MakePrmtop(splitFile)
            RunTleap(ionize=False, conc=None)
            numberWaters = run('grep "WAT" solvated.pdb | wc -l', shell=True, capture_output=True, text=True)
            output_string = int(numberWaters.stdout.strip())
            # formula taken from https://computecanada.github.io/molmodsim-amber-md-lesson/12-Adding_Ions/index.html
            ionConc = (0.0028798 * output_string) // 2
            # calculating the waters in the solvated complex
            RunTleap(ionize=True, conc=ionConc)
            # running with OpenMM
        if not os.path.exists('./system/run_0/Traj_0.dcd'):
            os.chdir("./system")
            run(f'python -u {OPENMM_PIPELINE} -r 3 -e 2 3 -pt 2 -et 2 -prt 10 -in 4 -eq E -run R > openmm.log 2>&1',
                shell=True)
            os.chdir("../")
        if not os.path.exists("DynamicScores_amber"):
            # Compute a weighted ensemble dynamic scores
            RestoreChains()
            run(f"python {DYNAMIC_SCORER} True {agChain}", shell=True)
        if not os.path.exists("gbsa/GBSA_on_PDB.pdb"):
            os.chdir("./gbsa")
            run("cp ../final/complex_final.pdb complex.pdb", shell=True)
            run("cp ../final/receptor_final.pdb receptor.pdb", shell=True)
            run("cp ../final/ligand_final.pdb ligand.pdb", shell=True)
            run(f'python {VISUALIZER}', shell=True)
            os.chdir("../")
        if not os.path.exists("minimal_interface"):
            GetCloseSelection(abChain, agChain)

    def MakeTrajectoryFromPDB(self, pdb, chains) -> None:
        chdir('./structures')
        chdir(pdb)
        os.makedirs('./logs', exist_ok=True)

        if not os.path.exists('gbsa/') and not os.path.exists('initial/'):
            PrepareProtein(pdb=pdb)
            x_range, y_range, z_range = GetCrystalCoords(f"./final/complex_final.pdb")
            CopyAndSplitSystem(self.ROOT, pdb, chains, False, x_range, y_range, z_range)

        if not os.path.exists('gbsa/complex.dcd') or not os.path.exists("DynamicScores_amber"):
            self.TrajectorizePDB(chains)
        chdir(self.ROOT)

    def ParallelPipeline(self) -> None:
        cpuUnits = int(cpu_count() // 2)
        processes = []
        with Pool(processes=cpuUnits) as p:
            for pdb, chains in self.dbDict.items():
                processes.append(p.apply_async(self.MakeTrajectoryFromPDB, args=(pdb, chains,)))
            for proc in processes:
                proc.get()
            p.close()
            p.join()
