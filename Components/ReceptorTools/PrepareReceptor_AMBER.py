import os.path

from subprocess import run


def PrepareProtein(pdb: str) -> None:
    """Prepares the system for AMBER using htmd."""
    os.makedirs('final', exist_ok=True)
    os.makedirs('temp', exist_ok=True)

    if not os.path.exists('./.check'):
        import htmd.ui as htmdmodule
        tick = None
        with open(f"./{pdb}_clean.pdb", 'r') as originalREC:
            for line in originalREC:
                if 'MEMB' in line:
                    tick = 32.0
        try:
            protein = htmdmodule.Molecule(f'./{pdb}_clean.pdb')
            receptor = htmdmodule.systemPrepare(protein, hydrophobic_thickness=tick, ignore_ns_errors=True,
                                                hold_nonpeptidic_bonds=True, titration=True)
            receptor.write('./temp/system_H.pdb')
            # we need to remove hydrogens as systemPrepare adds CHARMM-like Hs atomtypes to the system...
            run('pdb4amber -i ./temp/system_H.pdb -o ./final/complex_final.pdb -y -a --noter; rm -r ./temp;touch ./.check',
                shell=True)
            del protein, receptor  # clearing memory as we don't need those anymore
        except Exception as e:
            print("Building the system raised this exception:\n\n", e)
            print("\nTrying pdb4amber first.\n")
            run(f'pdb4amber -i {pdb}_clean.pdb -o ./temp/complex_pdb4amber.pdb -y -a', shell=True)
            protein = htmdmodule.Molecule('./temp/complex_pdb4amber.pdb')
            receptor = htmdmodule.systemPrepare(protein, hydrophobic_thickness=tick, ignore_ns_errors=True,
                                                hold_nonpeptidic_bonds=True, titration=True)
            receptor.write('./temp/receptor_H.pdb')
            run('pdb4amber -i ./temp/receptor_H.pdb -o ./final/complex_final.pdb -y -a --noter; rm -r ./temp;touch ./.check',
                shell=True)
    else:
        print(
            "Found hidden .check file from a previous pdb4amber run. Remove this file if you want to run pdb4amber again")
