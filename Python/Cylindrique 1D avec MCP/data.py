import numpy as np


# On crée ici une classe DataDict qui hérite du comportement du dictionaire classique (dict) mais
# pour lequel on souhaite étendre les fonctionnalitées.
class DataDict(dict):

    # Ces deux fonctions permettent de d'utiliser la notation pointée.
    # Par exemple,  pour une instance "param" de DataDict, on pourra faire un data.M = 2 au lieu de data['M']
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as k:
            raise AttributeError(k)

    def __setattr__(self, key, value):
        self[key] = value

    # cette méthode est appelée lorqu'on a besoin de récupérer une représentation String de l'objet,
    # par exemple lorsqu'on fait un print()
    def __repr__(self):
        str = ''
        for k in self.keys():
            str += f"{k} = {self[k]}\n"
        return str


    # On en profite pour ajouter la méthode pour charger les données depuis un fichier
    def load(self,filename):
        try:
            file = open(filename, "r")
        except OSError:
            print(f"Oops: problème lors de l'ouverture du fichier {filename}.")
            exit()

        # Si j'arrive ici c'est que le fichier existe (et est lisible)

        lines = file.readlines() # Je récupère un tableau de String contenant chaque ligne du fichier
        for str in lines: # Pour chaque ligne
            str = str.strip()        # je commence par nettoyer un peu (ie j'enlève les espaces)
            if str.startswith("#") or len(str)==0:
                # C'est une ligne de commentaire (ou ligne vide), je ne traite pas
                continue
            # Je découpe autour du signe '='
            tmp = str.split("=")
            if len(tmp)<2:
                # on vérifie qu'on a ce qu'on attend
                print(f"Oops: problème lors de l'analyse de la ligne {str}.")
                exit()
            # Ici la "clé" est dans tmp[0] et la valeur dans tmp[1] --> à mettre dans le dictionnaire "GLOBALS"
            key = tmp[0].strip()

            # pour la valeur, on se donne la possibilité de considérer tout ce qui se trouve après une espace comme un commentaire

            try:
                val = float(tmp[1].strip().split(" ")[0]) # Je vous laisse décomposer cette ligne pour comprendre...
            except ValueError:
                print(f"Oops: problème lors lors du traitement du paramètres {key}: '{tmp[1]}' n'est pas un nombre valide.")
                exit()
            # ici on a un nombre valide dans la variable val
            if val == int(val):
                # petit "trick" pour convertir les nombres entier en variables de type INT
                # -> Important pour éviter les erreurs ensuite lors de la création des tableau par exemple car 10.0 (float)
                # n'est pas considéré par Python comme 10 (entier)
                val=int(val)
            # Finalement, ajout de la variable dans le dictionnaire !!
            self[key]=val


