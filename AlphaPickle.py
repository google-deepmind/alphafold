#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# ##AlphaPickle: making AlphaFold2 outputs interpretable
# 
# AlphaPickle is multipurpose Python script for producing plots and user-legible files from the output of [AlphaFold2](https://www.nature.com/articles/s41586-021-03819-2) \([notebook](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)\) and [Colabfold](https://www.biorxiv.org/content/10.1101/2021.08.15.456425v2) \([notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=UGUBLzB3C6WN)\).
# 
# The functions provided here aim to tackle the problem of output metadata from these programs being difficult to process for users who don't have Python training. Currently, PAE outputs from AlphaFold 2.1.0 are in the form of .pkl files, which can only be read using a Python script. For data from AlphaFold DB and ColabFold, these data are in .json format (not typically easy to process for non-code writing users). 
# 
# For all of the above software, pLDDT values are outputted in the B-factor field of the PDB file for each prediction. This is very useful for visualisation (e.g. using the ChimeraX command `color bfactor palette alphafold`), but may be difficult in terms of customisable visualisation for non-code writing users.
# 
# This collection of code will take any of the above output files and provide a .csv file (which can be opened and used for plotting in Excel, Numbers, Google Sheets) as well as a downloadable plot.
# 
# To use, please run every cell in order (or select Runtime > Run all). When prompted, select as many files of the given type (Please upload PAE files and pLDDT files separately) as you wish in the file browser.


# @title Import modules and define classes and functions
# @markdown Please run this cell first to prepare the python environment.
import argparse
import glob
import pathlib
import pickle as pkl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import json
from sys import exit
import os
from Bio import PDB as pdb
# @title Import modules and define classes and functions
# @markdown Please run this cell first to prepare the python environment.
import argparse
import glob
import json
import os
import pathlib
import pickle as pkl
from sys import exit

import numpy as np
import pandas as pd
from Bio import PDB as pdb
from matplotlib import pyplot as plt


# Define class for AlphaFold metadata file and class methods


class AlphaFoldMetaData(object):
  def __init__(self, PathToFile, SavePath=None, FastaSequence=None, ranking=None):
    # Define attributes
    self.PathToFile = PathToFile
    self.FastaSequence = FastaSequence
    self.saving_filename = self.PathToFile.stem
    self.saving_pathname = self.PathToFile.parent if SavePath == None else SavePath
    if ranking:
      self.saving_filename = "ranked_{}".format(ranking)

  # Generate a plot of pLDDT value
  def plot_pLDDT(self, size_in_inches=12, axis_label_increment=100):
    x = list(range(0, len(self.pLDDT), 1))
    y = list(self.pLDDT)

    plt.figure(figsize=(size_in_inches, (size_in_inches / 2)))
    ticks = np.arange(0, len(self.pLDDT), axis_label_increment)
    plt.xticks(ticks, fontname="Helvetica")
    plt.yticks(fontname="Helvetica")
    plt.title(
      f"pLDDT for {self.saving_filename}, ranked {ranking.index(str(self.saving_filename).replace('result_', ''))}")
    plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
    plt.ylabel("Predicted LDDT", size=14, fontweight="bold", fontname="Helvetica")
    plt.plot(x, y)

    plt.savefig(
      f"{self.saving_pathname}/plddt_ranked_{ranking.index(str(self.saving_filename).replace('result_', ''))}_{self.saving_filename}_pLDDT.png",
      dpi=300)
    # plt.show()

    # Use standard AlphaFold colors
    """
    cmap = cols.LinearSegmentedColormap.from_list("", ["red","orange","yellow","cornflowerblue","blue"])

    plt.figure(figsize=(size_in_inches,(size_in_inches/2)))
    ticks = np.arange(0, len(self.pLDDT), axis_label_increment)
    plt.xticks(ticks, fontname="Helvetica")
    plt.yticks(fontname="Helvetica")
    plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
    plt.ylabel("Predicted LDDT",size=14, fontweight="bold", fontname="Helvetica")
    plt.scatter(x,y, c=y, cmap=cmap, s=5) 
    plt.clim(0,100)
    scale = plt.colorbar(shrink=0.5)
    scale.set_label(label="Predicted LDDT",size=12, fontweight="bold", fontname="Helvetica")
    # Save to directory with pickle file in
    plt.savefig('{}/{}_pLDDT.png'.format(self.saving_pathname, self.saving_filename), dpi=300)
    """

    # Generate a plot from PAE measurements

  def plot_PAE(self, size_in_inches=12, axis_label_increment=100):
    ticks = np.arange(0, self.PAE[1].size, axis_label_increment)
    plt.figure(figsize=(size_in_inches, size_in_inches))
    PAE = plt.imshow(self.PAE, cmap="bwr")
    plt.title(
      f"PAE for {self.saving_filename}, ranked {ranking.index(str(self.saving_filename).replace('result_', ''))}")
    plt.xticks(ticks, fontname="Helvetica")
    plt.yticks(ticks, fontname="Helvetica")
    plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
    plt.ylabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
    scale = plt.colorbar(PAE, shrink=0.5)
    scale.set_label(label="Predicted error (Å)", size=12, fontweight="bold", fontname="Helvetica")

    # Save plot
    plt.savefig(
      f"{self.saving_pathname}/pae_ranked_{ranking.index(str(self.saving_filename).replace('result_', ''))}_{self.saving_filename}_PAE.png",
      dpi=300)

    # Generate dataframe from PAE data and save to csv
    pd_PAE = pd.DataFrame(self.PAE)
    pd_PAE.to_csv(
      f"{self.saving_pathname}/pae_ranked_{ranking.index(str(self.saving_filename).replace('result_', ''))}_{self.saving_filename}_PAE.csv")


class AlphaFoldPickle(AlphaFoldMetaData):

  def __init__(self, PathToFile, FastaSequence=None, ranking=None):
    super().__init__(PathToFile, FastaSequence, ranking)  # Define attributes
    if ranking:
      self.saving_filename = "ranked_{}".format(ranking)
    self.data = []
    self.PAE = None

    # Extract pickled data
    with (open("{}".format(self.PathToFile), "rb")) as openfile:
      while True:
        try:
          self.data.append(pkl.load(openfile))
        except EOFError:
          break

    # Try statement accounts for data run using non-pTM models, with no PAE output
    try:
      self.PAE = self.data[0]['predicted_aligned_error']
    except:
      print("PAE model data not present. To access this performance metric, run AlphaFold"
            "using pTM-enabled models.")

    # Define pLDDT
    self.pLDDT = self.data[0]['plddt']

  # Generate a ChimeraX attribute file from pLDDT measurements
  def write_pLDDT_file(self):
    seqMismatch = False
    pd_lDDT = pd.DataFrame(self.pLDDT)
    # Name dataframe column
    pd_lDDT.columns = ["pLDDT"]

    # If the fasta file was provided:
    if self.FastaSequence != None:

      # Open the fasta file in read mode
      with (open("{}".format(self.FastaSequence), "r")) as openfile:
        fasta = openfile.read()

      # Delete header line and remove line-breaks
      sequence = fasta.split("\n", 1)[1].replace("\n", "")

      # Check that the lengths of the two sequences match
      if len(sequence) != len(pd_lDDT):

        # If not, ignore the fasta file
        print(
          "Length of sequence in fasta file provided ({}) does not match length of sequence used in AlphaFold prediction ({}). Ignoring fasta file.".format(
            len(sequence), len(pd_lDDT)))
        seqMismatch = True
      # If they do,
      else:
        # Convert the fasta sequence into a residue list
        list_sequence = []
        for item in sequence:
          list_sequence.append(item)

        # Convert the list into a pandas series
        pd_sequence = pd.Series(list_sequence)

        # Insert the series into the dataframe at column 1 to act as labels for the data
        pd_lDDT.insert(0, "Residue", pd_sequence)

    # Otherwise, remind user to check that they have used corret input files
    else:
      print("Number of residues for which pLDDT is provided: ", len(pd_lDDT),
            "If this does not match the length of your sequence, please double check the input file.")

    # Tell python not to elide middle rows of dataframe when printing to std.out
    pd.set_option("display.max_rows", None, "display.max_columns", None)

    # Save dataframe to ./outputfiles with same name as original pickle and .csv extension
    pd_lDDT.to_csv('{}/{}_pLDDT.csv'.format(self.saving_pathname, self.saving_filename))
    # Delete residue ID
    if self.FastaSequence != None and seqMismatch == False:
      lDDT_table = pd_lDDT.drop('Residue', axis=1)
    else:
      lDDT_table = pd_lDDT

    # Initialise list to store Chimera-style residue identifiers (":x", where x = residue number)
    residue_list = []

    # Populate this list
    for residue in range(0, len(lDDT_table)):
      residue_list.append(":{}".format(residue + 1))

    # Save to pandas format
    chimerax_numbering = pd.Series(residue_list)

    # Insert in the first column of the dataframe, to satisfy ChimeraX formatting
    lDDT_table.insert(0, 'Numbering', chimerax_numbering)

    # Tidy indices so the first label is 1 not 0
    pd_lDDT.index += 1

    # Create a file to save the Chimera attribute output into

    with (open('{}/{}_lDDT.txt'.format(self.saving_pathname, self.saving_filename), 'w+')) as openfile:

      # Write file header in correct format
      openfile.write('attribute: pLDDTvalue\nmatch mode: 1-to-1\nrecipient: residues\n')

      # Iterate over rows of dataframe, writing residue ID and lDDT value to file with correct formatting
      for i, row in lDDT_table.iterrows():
        openfile.write("\t{}\t{}\n".format(row['Numbering'], row['pLDDT']))

    return pd_lDDT


class AlphaFoldJson:
  def __init__(self, PathToDirectory):
    self.PathToDirectory = PathToDirectory
    self.RankingDebug = []
    try:
      with open("{}/ranking_debug.json".format(self.PathToDirectory)) as jsonfile:
        self.RankingDebugRaw = json.load(jsonfile)
      for index in enumerate(self.RankingDebugRaw['order']):
        self.RankingDebug.append(index)
    except:
      exit(
        "To use batch processing, please ensure that the ranking_debug.json file and the result_model_n.pkl files are present in the directory issued in the command. Exiting AlphaPickle now...")


class AlphaFoldPDB(AlphaFoldMetaData):
  def loadCleanStructure(self, id, filePath):
    standardResidues = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    parser = pdb.PDBParser()
    parsedStructure = parser.get_structure(id, filePath)
    for chain in parsedStructure.get_chains():
      removeResidues = list()
      for i, residue in enumerate(chain.get_residues()):
        if residue.resname not in standardResidues:
          removeResidues.append(residue.id)
          print(residue.id)
      [chain.detach_child(id) for id in removeResidues]

    return parsedStructure

  def extractPLDDT(self, PDBobject):
    pLDDT = []
    for residue in PDBobject.get_residues():
      i = 0
      for atom in residue.get_atoms():
        while i < 1:
          pLDDT.append(atom.bfactor)
          i += 1
    pLDDT_series = pd.Series(pLDDT)
    return pLDDT_series

  def __init__(self, PathToFile, FastaSequence=None, ranking=None):
    super().__init__(PathToFile, FastaSequence, ranking)
    # Define attributes
    if ranking:
      self.saving_filename = "ranked_{}".format(ranking)
    self.structure = self.loadCleanStructure("test", PathToFile)
    self.pLDDT = self.extractPLDDT(self.structure)
    self.data = []
    self.PAE = None

  def PDB_write_pLDDT(self):
    residueNumbers = pd.Series(range(1, len(self.pLDDT) + 1))
    if len(residueNumbers) != len(self.pLDDT):
      print("Oops")
    else:
      pd_lDDT = pd.DataFrame(self.pLDDT)
      pd_lDDT.columns = ["pLDDT"]
      pd_lDDT.insert(0, "Residue", residueNumbers)
      pd_lDDT.to_csv('{}/{}_pLDDT.csv'.format(self.saving_pathname, self.saving_filename))


class AlphaFoldPAEJson(AlphaFoldMetaData):
  def extractPAEfromJson(self, PathToFile):

    with open(PathToFile, 'r') as file:
      jsonstring = json.load(file)

      residue1 = jsonstring[0]['residue1']
      residue2 = jsonstring[0]['residue2']
      pae = jsonstring[0]['distance']

    paeArray = np.ones((max(residue1), (max(residue2))))

    for i, j, n in zip(residue1, residue2, pae):
      paeArray[int(i - 1), int(j - 1)] = n

    return paeArray

  def __init__(self, PathToFile, FastaSequence=None, ranking=None):
    super().__init__(PathToFile, FastaSequence, ranking)
    if ranking:
      self.saving_filename = "ranked_{}".format(ranking)

    self.PAE = self.extractPAEfromJson(PathToFile)
    self.pLDDT = None


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Run alphapickle to Alphafold results")
  parser.add_argument('--res_dir', dest='res_dir', help='PDB file for scan', required=True, type=str)
  parser.add_argument('--save_dir', dest='save_dir', action="store", type=str, default="",
                      help="Path to save plots and csvs")
  args_dict = vars(parser.parse_args())

  # read args
  res_dir = pathlib.Path(args_dict['res_dir']).resolve()
  save_dir = pathlib.Path(args_dict['save_dir'] if args_dict['save_dir'] != "" else res_dir.joinpath('plot')).resolve()

  # check args
  assert os.path.exists(res_dir)

  cwd = str(res_dir)

  os.makedirs(save_dir, exist_ok=True)

  ranking_record = pathlib.Path(f"{cwd}/ranking_debug.json").resolve()
  ranking = json.load(open(ranking_record, 'r'))['order']

  # @title Upload files
  # @markdown Upload as many metadata files as you need to process here. The program will prompt for two separate uploads:
  # @markdown - PAE files: select multiple files containing PAE data of either (or a mixture of) .pkl or .json format. If this function is not desired, press cancel upload.
  # @markdown - pLDDT files: select multiple files containing pLDDT data of either (or a mixture of) .pkl or .pdb format. If this function is not desired, press cancel upload.

  # print("Select PAE files for upload")
  PAEfiles = [pathlib.Path(i).resolve() for i in list(glob.glob(f'{cwd}/*.pkl')) if "features.pkl" not in i]
  # print(PAEfiles)
  # print("Select pLDDT files for upload")
  pLDDTfiles = [pathlib.Path(i).resolve() for i in list(glob.glob(f'{cwd}/*.pkl')) if "features.pkl" not in i]
  # print(pLDDTfiles)

  # @title Process files
  # @markdown Input parameters for plotting here. The program will plot all of the previously uploaded files and save the metadata values to a .csv file, for easy legibility and downstream processing.

  # @markdown Don't worry about copying and pasting plots; the outputs can all be downloaded below.

  # @markdown Input desired plot size, in inches.
  plot_size = 12  # @param{type: "integer"}

  # @markdown Input value to increment plot axes by (this may need finetuning based on output)
  plot_increment = "50"  # @param[10,25,50,100,250,500]
  plot_increment = int(plot_increment)

  filesForDownload = []

  for PAEfile in PAEfiles:
    print(PAEfile)
    fileType = PAEfile.suffix
    fileName = PAEfile.stem
    if fileType == ".pkl":
      results = AlphaFoldPickle(PAEfile, None)
      results.saving_pathname = save_dir
      results.saving_filename = fileName
      if type(results.PAE) == np.ndarray:
        print("Plotting PAE for {} and saving to csv".format(PAEfile))
        results.plot_PAE(size_in_inches=plot_size, axis_label_increment=plot_increment)
      filesForDownload.extend(
        ["{}_PAE_.png".format(results.saving_filename), "{}_PAE.csv".format(results.saving_filename)])
    else:
      raise TypeError(
        "Expected file of type .pkl or .json. Check the extensions of uploaded PAE files match one of these and rerun the upload step.")

  for pLDDTfile in pLDDTfiles:

    fileType = pLDDTfile.suffix
    fileName = pLDDTfile.stem
    if fileType == ".pkl":
      results = AlphaFoldPickle(pLDDTfile, None)
      results.saving_filename = fileName
      results.saving_pathname = save_dir
      # results.PDB_write_pLDDT()
      print("Plotting pLDDT for {} and saving to csv".format(pLDDTfile))
      results.plot_pLDDT(size_in_inches=plot_size, axis_label_increment=plot_increment)
      filesForDownload.extend(
        ["{}_pLDDT.png".format(results.saving_filename), "{}_pLDDT.csv".format(results.saving_filename)])
    elif fileType == ".pdb":
      results = AlphaFoldPDB(pLDDTfile)
      results.saving_filename = fileName
      results.saving_pathname = save_dir
      print("Plotting pLDDT for {} and saving to csv".format(pLDDTfile))
      # results.PDB_write_pLDDT()
      results.plot_pLDDT(size_in_inches=plot_size, axis_label_increment=plot_increment)
      filesForDownload.extend(
        ["{}_pLDDT.png".format(results.saving_filename), "{}_pLDDT.csv".format(results.saving_filename)])
    else:
      raise TypeError(
        "Expected file of type .pkl or .pdb. Check the extensions of uploaded pLDDT files match one of these and rerun the upload step.")
