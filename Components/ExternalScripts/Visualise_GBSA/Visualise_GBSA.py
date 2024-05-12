import Decomp_Reader
import PDB_Reader
import PDB_Writer
import Plotter


def main():
    fileParser = PDB_Reader.Parser()
    ligandList, receptorList = fileParser.getLists()
    decomp = Decomp_Reader.DecompReader()
    ligandDecompResults, recDecompResults = decomp.GetDecompResults()
    PDB_Writer.DataWriter(ligandList, receptorList, ligandDecompResults, recDecompResults)
    Plotter.Plotter()


if __name__ == '__main__':
    main()
