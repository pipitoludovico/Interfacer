import os

import requests
from urllib import request


class Fetcher:
    def __init__(self, url="https://files.rcsb.org/download/"):
        self.url = url
        self.response = None
        self.entryID = None

    def FetchPDB(self, entryID):
        self.entryID = [str(ids).upper() for ids in entryID]
        if not os.path.exists('./structures'):
            os.makedirs("./structures", exist_ok=True)
            for id_to_fetch in self.entryID:
                try:
                    self.response = requests.get(self.url + id_to_fetch + ".pdb")
                    os.makedirs(f'structures/{id_to_fetch}', exist_ok=True)
                    request.urlretrieve(f'http://files.rcsb.org/download/{id_to_fetch}.pdb',
                                        f'./structures/{id_to_fetch}/{id_to_fetch}.pdb')
                except ConnectionRefusedError:
                    with open('FailedDownloads.txt', 'a') as Failed:
                        Failed.write(str(id_to_fetch) + "\n")
