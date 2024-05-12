import warnings
import argparse
from argparse import Namespace
from os import getpid, system, path
from typing import Any

warnings.filterwarnings(action='ignore')


class ArgParser:
    def __init__(self):
        pass

    @staticmethod
    def argumentParser() -> tuple[Namespace, dict[str, list[Any]]]:
        """
        Parse command-line arguments for the program.
        Returns:
            argparse.Namespace: Parsed command-line arguments.
        """
        ap = argparse.ArgumentParser()
        ap.add_argument('infile', type=str, help="Path to your input file with a column of PDB ids.")

        ap.add_argument("-k", '--kill', required=False, action='store_true', default=False,
                        help="Kills the current pid process.")

        args = ap.parse_args()

        # adopting a dictionary to avoid duplicates
        STRUCTURES = {}
        if not path.exists(args.infile):
            print("The specified input file does not exist.")
            exit()
        else:
            with open(args.infile, 'r') as structureIDsFile:
                for line in structureIDsFile:
                    STRUCTURES[line.replace("\n", "")] = []

        if args.kill:
            try:
                with open(".mypid", 'r') as previousPID:
                    oldPID = previousPID.read()
                system(f'kill -9 -{oldPID}')
                exit()
            except Exception as e:
                print("No existing process found to kill.")
                print(e)
                exit()

        newPID = getpid()
        with open(".mypid", "w") as PIDFILE:
            PIDFILE.write(str(newPID))

        return args, STRUCTURES
