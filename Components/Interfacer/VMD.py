import os
from subprocess import run, DEVNULL


def RestoreChains():
    chains_: list = []
    with open("./final/complex_final.pdb", 'r') as pdb4amberFile:
        for pdb4amberLine in pdb4amberFile.readlines():
            if len(pdb4amberLine.split()) > 2:
                if "ATOM" in pdb4amberLine:
                    chainID = pdb4amberLine[21]
                    resNum = pdb4amberLine[23:26].strip()
                    if (chainID, resNum) not in chains_:
                        chains_.append((chainID, resNum))
    with open("./gbsa/complex_H.pdb", 'r') as minimizedPDB:
        minimizedLines = minimizedPDB.readlines()
        with open('./gbsa/complex_chains.pdb', 'w') as test:
            for line in minimizedLines:
                if len(line.split()) > 2:
                    resNum_new = line[23:26].strip()
                    for pair in chains_:
                        if resNum_new == pair[1]:
                            newline = line[:21] + str(pair[0]) + line[22:]
                            test.write(newline)
                else:
                    test.write(line)


def GetVDWcontacts(abChains, agChains):
    tcl_script_lines = [
        'mol load pdb final/complex_final.pdb',
        '',
        'proc contactFreq {sel1 sel2 outFile mol} {',
        '  set allCounts {}',
        '  set numberOfFrames [molinfo $mol get numframes]',
        '',
        '  if { $outFile != "stdout" } {',
        '     set outFile [open $outFile w]',
        '  }',
        '',
        '  for {set i 0} {$i < $numberOfFrames} {incr i} {',
        '    molinfo $mol set frame $i',
        '',
        '    set frameCount1 [atomselect $mol "$sel1 and noh and within 3.5 of ($sel2 and noh)"]',
        '    set frameCount2 [atomselect $mol "$sel2 and noh and within 3.5 of ($sel1 and noh)"]',
        '',
        '    set uniqueContacts [list]',
        '',
        '    foreach a [$frameCount1 get {chain resname resid}] {',
        '      foreach b [$frameCount2 get {chain resname resid}] {',
        '        set contact [concat $a "-" $b]',
        '        if {[lsearch -exact $uniqueContacts $contact] == -1} {',
        '          lappend uniqueContacts $contact',
        '        }',
        '      }',
        '    }',
        '',
        '    set numContacts [llength $uniqueContacts]',
        '',
        '    lappend allCounts $numContacts',
        '',
        '    if { $outFile == "stdout" } {',
        '      puts "Frame $i $numContacts contacts"',
        '    } else {',
        '      set contactInfo {}',
        '      foreach contact $uniqueContacts {',
        '        lappend contactInfo [join $contact ""]',
        '      }',
        '      puts $outFile "$numContacts,[join $contactInfo ","]"',
        '    }',
        '',
        '    $frameCount1 delete',
        '    $frameCount2 delete',
        '  }',
        '',
        '  if { $outFile != "stdout" } {',
        '    close $outFile',
        '  }',
        '}',
        '',
        f'set sel1 "chain {" ".join(abChains)}"',
        f'set sel2 "chain {agChains}"',
        'set outName "contacts.int"',
        '',
        '',
        'puts "GETTING CONTACTS"',
        'contactFreq $sel1 $sel2 $outName top',
        '',
        'exit'
    ]

    with open('getContacts.tcl', 'w') as vmdContactsFile:
        for line in tcl_script_lines:
            vmdContactsFile.write(line + "\n")
    run(f'vmd -dispdev text -e getContacts.tcl > logs/getContacts.log', shell=True, check=True)


def GetCloseSelection(abChain, agChain):
    os.makedirs('minimal_interface', exist_ok=True)
    selectVMD = ("package require psfgen",
                 "resetpsf",
                 f"mol load pdb gbsa/complex.pdb",
                 f'set sel [atomselect top "(not chain {" ".join(abChain)} and same residue as protein within 10 of chain {" ".join(abChain)}) or'
                 f' (not chain {agChain} and same residue as protein within 10 of chain {agChain})"]',
                 "$sel writepdb minimal_interface/selection.pdb",
                 "quit")
    with open('writeSel.tcl', 'w') as writer:
        for line in selectVMD:
            writer.write(line + "\n")
    run('vmd -dispdev text -e writeSel.tcl > logs/writeSel.log 2>&1; cp gbsa/complex.pdb complex_minimized.pdb',
        shell=True,
        stdout=DEVNULL)
