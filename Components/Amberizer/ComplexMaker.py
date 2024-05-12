from os import chdir
from subprocess import run
import MDAnalysis as Mda
import numpy as np

import os

external_scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'ExternalScripts')
SCRIPTPATH = os.path.join(external_scripts_dir, 'prepareProtein.py')


def CopyAndSplitSystem(ROOT: str, pdb: str, chains: list, htmd: bool, x_range, y_range, z_range) -> None:
    if not htmd:
        u = Mda.Universe(f"./final/complex_final.pdb")
        u.dimensions = np.array([[x_range, y_range, z_range, 90, 90, 90]])
        SeparateComplex(u, chains)
    else:
        try:
            run(f"{SCRIPTPATH} ./final/complex_final.pdb cytosol", shell=True)
            u = Mda.Universe(f'protein_H.pdb')
            SeparateComplex(u, chains)
        except Exception as e:
            print(repr(e))
            chdir(ROOT)
            with open('failed.txt', 'a') as failFile:
                failFile.write(pdb + " failed.\n")
            return


def SeparateComplex(universe, chains: list):
    abChain = ''
    agChain = ''
    if len(chains) == 2:
        abChain = chains[0]
        agChain = chains[1]
    if len(chains) > 2:
        abChain = chains[0] + " " + chains[1]
        agChain = chains[2]

    sel = universe.select_atoms(f'protein or nucleic and chainID {abChain} and not name H*')
    sel.write('final/receptor_final.pdb')
    # run('pdb4amber -i final/receptor_pre.pdb -o ./final/receptor_final.pdb -y -a', shell=True)

    sel = universe.select_atoms(f'protein or nucleic and chainID {agChain} and not name H*')
    sel.write('final/ligand_final.pdb')
    # run('pdb4amber -i final/ligand_pre.pdb -o ./final/ligand_final.pdb -y -a', shell=True)
