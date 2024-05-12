import os
from os import getcwd
import pandas as pd

from Components.Parser.Parser import ArgParser
from Components.Fetcher.Fetcher import Fetcher
from Components.StructureRefiner.StructureRefiner import Refiner
from Components.Amberizer.TrajectoryMaker import TrajectoryMaker
from Components.Interfacer.PeStO import ParallelPesto
from Components.Featurizer.FeatureMaker import FeaturizerClass

ROOT = getcwd()
args, structures = ArgParser.argumentParser()


def main():
    Fetcher().FetchPDB(structures)
    refiner = Refiner(structures, ROOT)
    dbDict = refiner.GetChains()

    TrajMaker = TrajectoryMaker(dbDict)
    TrajMaker.ParallelPipeline()
    os.chdir(ROOT)
    ParallelPesto(dbDict, ROOT)
    featurizer = FeaturizerClass(dbDict, ROOT)
    featurizer.ParallelFeaturize()
    df = featurizer.GetDataAsDataframe()
    pd.DataFrame.to_csv(df, 'output.csv', index=True)


if __name__ == '__main__':
    main()
