from os import makedirs
from subprocess import run, CalledProcessError, DEVNULL


def RunTleap(ionize=None, conc=1) -> None:
    makedirs('system', exist_ok=True)
    makedirs('gbsa', exist_ok=True)

    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",

         "complex = loadpdb ./final/complex_final.pdb" if not ionize else '',
         'setBox complex "vdw"' if not ionize else '',
         "savepdb complex ./gbsa/complex.pdb" if not ionize else '',
         "solvateBox complex TIP3PBOX 15 0.75" if not ionize else '',
         "savepdb complex solvated.pdb" if not ionize else '',
         "complex = loadpdb solvated.pdb" if ionize else "",
         "addIons complex Na+ 0" if ionize else "",
         "addIons complex Cl- 0" if ionize else "",
         f"addIons complex Na+ {conc}" if ionize else "",
         f"addIons complex Cl- {conc}" if ionize else "",
         'setBox complex "vdw"' if ionize else "",
         "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd" if ionize else "",

         "quit"]

    with open('inleap', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run("tleap -f inleap", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()


def MakePrmtop(structureFile) -> None:
    makedirs('gbsa', exist_ok=True)
    output = structureFile.split("_")[0]
    _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
         "source leaprc.lipid21",
         "set default PBRadii mbondi3",

         f"complex = loadpdb ./final/{structureFile}",
         'setBox complex "vdw"',
         f"saveamberparm complex ./gbsa/{output}.prmtop ./gbsa/{output}.inpcrd",
         f"savepdb complex ./gbsa/{output}_H.pdb",
         "quit"]

    with open('inleap_prmtop', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run("tleap -f inleap_prmtop", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()
