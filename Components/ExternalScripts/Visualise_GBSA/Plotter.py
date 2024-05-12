import matplotlib.pyplot as plt
import pandas as pd


class Plotter:

    def __init__(self):
        self.Plot()

    @staticmethod
    def Plot():
        files = ['ligand.dat', 'receptor.dat']
        for file in files:
            df = pd.read_csv(file, sep='\t', usecols=['Residue', 'GBSA'])
            df_filtered = df[(df['GBSA'] < -0.5) | (df['GBSA'] > 0.5)]
            plt.style.use('ggplot')  # uses gnuplot style
            df_filtered.plot(x='Residue', y='GBSA', kind='bar', legend=False)
            plt.xlabel('Residues', labelpad=30)
            plt.ylabel('Kcal/mol', labelpad=15)
            plt.tight_layout()
            plt.savefig(f"{file.split('.')[0]}.jpg", dpi=300)
            plt.show()
            plt.close()


