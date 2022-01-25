import pandas as pd
import numpy as np
import yaml
from zipfile import ZipFile
from biom.table import Table
from biom.util import biom_open
import copy

import plotly.graph_objects as go

"""
import argparse

parser = argparse.ArgumentParser(description='Preparing a barplot from phyloseq barplot data.')

parser.add_arglument('--taxonomy','-t', dest='TAX', type=str, help='Taxonomy file.', required=True)
parser.add_argument('--feature','-f', dest='FEAT', type=str, help='Feature (OTU) file.', required=True)
parser.add_argument('--mapping','-m', dest='MAP', type=str, help='Mapping file.', required=True)
parser.add_argument('--taxonomySep','-ts', dest='TAX_SEP', type=str, help='Taxonomy file separator [class default is TAB].')

args = parser.parse_args()

########################################################################################################################
"""

class seqObject():
    """
    This class represents a sequencing file object managing a mapping file, feature/otu file and taxonomy with some associated functions
    """
    def __init__(self,
                 mappingFile=None,
                 taxonomyFile=None,
                 featureFile=None,
                 mappingSep = '\t',
                 taxonomySep = '\t',
                 featureSep = '\t',
                 sampleNamesColumn="Sample",
                 taxonomyDatabase=None,
                 representativeSeqFile=None,
                 featureFormat="automatic",
                 featureColumnName="None",
                 featureRowName="",
                 taxonomyFormat=None,
                 featureFilePickle=False,
                 matchFeaturesToTaxonomy=True,
                 threading=10):
        """
        Creates and organizes class variables from the input files.

        :param mappingFile: A file containing mapping file
            :type mappingFile: str
        :param taxonomyFile: A file containing the taxonomy with separate columns for taxonomic levels
            :type taxonomyFile: str
        :param featureFile: A file containing the feature/otu table with samples in columns and features/otus in rows
            :type featureFile: str
        :param mappingSep: Separator in the mapping file, default is '\t'
            :type mappingSep: char
        :param taxonomySep: Separator in the taxonomy file, default is '\t'
            :type taxonomySep: char
        :param featureSep: Separator in the feature/otu file, default is '\t'
            :type featureSep: char
        :param sampleNamesColumn: Name of the column containing sample names.
            :type sampleNamesColumn: str
        :param taxonomyDatabase: Selection of a database for a Silva taxonomy output (SILVA, RDP, GTDB, LTP, EMBL)
            :type taxonomyDatabase: str
        :param featureFilePickle: A flag determining whether the loaded feature file is in a pickle format for faster loading (default is False)
            :type featureFilePickle: bool

        """

        print("Initializing data loading.")
        self.mapSampleNames = sampleNamesColumn

        self.load_mapping(file=mappingFile, sep=mappingSep)
        self.load_taxonomy(file=taxonomyFile, sep=taxonomySep, database=taxonomyDatabase, taxonomyFormat=taxonomyFormat)
        self.load_feature(file=featureFile,
                          sep=featureSep,
                          format=featureFormat,
                          featureColumnName=featureColumnName,
                          representativeSeqFile=representativeSeqFile,
                          featureFilePickle=featureFilePickle,
                          threading=threading,
                          matchFeaturesToTaxonomy=matchFeaturesToTaxonomy)

    def load_mapping(self, file, sep):
        """
        Loading the mapping file and formatting it.

        :param file: input mapping file
            :type file: str
        :param sep: separator character used in the mapping file
            :type sep: char
        :return: mapping file based dictionary
        """
        print("Loading mapping file")
        # Loading mapping file
        #self.mapDict = {}
        #with open(file, 'r') as f_map:
            # Loading a header of the file
        #    self.mapHeader = f_map.readline().split(sep)
        #    self.mapHeader[-1] = self.mapHeader[-1].rstrip()

            # Looking up a sample name column
            #if self.sampleName in self.mapHeader:
            #    name_pos = self.mapHeader.index(self.sampleName)
            #else:
            #    raise ValueError("Mapping file : Column '{0}' missing in mapping file. {0} is a mandatory field to identify sample names.".format(self.sampleName))

            # Parsing the file into a dictionary
            #for line in f_map:
            #    line = line.rstrip()
            #    line_list = line.split(sep)
            #    key = line_list[name_pos]
            #    line_list.pop(name_pos)

            #    self.mapDict[key] = line_list

            # Removing the sample name from the mapping
            #self.mapHeader.pop(name_pos)

        self.df_map = pd.read_csv(file, sep=sep)

    def load_taxonomy(self, file, sep, database, taxonomyFormat):
        """
        Loading the taxonomy file and formatting it
        :param file: Input taxonomy file.
            :type file: str
        :param sep: Separator character used in the taxonomy file.
            :type sep: char
        :return: taxonomy based dictionary
        """
        print("Loading taxonomy")

        # If file ending is .qza, extract taxonomic table from the file and proceed to treat it regularly
        if file[-4:] == ".qza":
            file = self._load_qiime_zip(file, "taxonomy.tsv")

        self.taxDict = {}

        # Switch cases functions for different inputs
        def name_fun(header):
            print("Recognizing a 'name' column.")
            name_pos = header.index("Name")
            format = "custom"
            return name_pos, format

        def qiime2_fun(header):
            print("Recognizing a Qiime2 format.")
            name_pos = header.index("Feature ID")
            format = "Qiime2"
            return name_pos, format

        def mothur_fun(header):
            print("Recognizing a mothur format.")
            name_pos = 0
            format = "mothur"
            return name_pos, format

        def silva_fun(header):
            print("Recognizing a silva format.")
            name_pos = 0
            format = "silva"
            return name_pos, format

        def dada2_fun(header):
            print("Recognizing a DADA2 format.")
            name_pos = 0
            format = "dada2"
            return name_pos, format

        # Guessing the format of the first cell that should contain the id od the taxonomy columns
        def format_cases(argument):
            switcher = {
                "Name": name_fun,
                "Feature ID": qiime2_fun,
                "Otu1": mothur_fun,
                '"job_id"': silva_fun,
                '"name"' : silva_fun,
                "name" : silva_fun,
                '""' : dada2_fun
            }
            func = switcher.get(argument, lambda argument: (0, 0))
            output = func(argument)
            return output

        with open(file, 'r') as f_tax:
            # Loading a header of the file
            self.taxHeader = f_tax.readline().split(sep)
            self.taxHeader[-1] = self.taxHeader[-1].rstrip()

            # Retrieving format of the data file and the position of the name column (just in case this is random)
            name_pos, format = format_cases(self.taxHeader[0])
            if taxonomyFormat:
                format = taxonomyFormat

            if format == 0:
                raise ValueError("Taxonomy file : Unknown format of a taxonomy file. Please specify the taxonomy file format in 'taxonomyFormat'."
                                 "Supported formats are: Qiime2, mothur, silva, dada2\n "
                                 "Automatically recognized formats are: \n"
                                 "Qiime2 format with a 'Feature ID' column and a 'Taxon' column. \n"
                                 "Silva format with a 'name' column and taxonomy split in separate columns. \n"
                                 "DADA2 format with an unnamed sample column and separate columns for each taxonomic category.")

            ## This should be in the new custom function
            if format == "custom":
                # Parsing the file into a dictionary
                for line in f_tax:
                    line = line.rstrip()
                    line_list = line.split(sep)
                    key = line_list[name_pos]
                    line_list.pop(name_pos)
                    # Adding the OTU name at the end instead of the beginning of the list for the multiple index later
                    line_list.append(key)

                    self.taxDict[key] = line_list

                # Removing the sample name from the mapping
                self.taxHeader.pop(name_pos)
            elif format == "seqDataClass":
                self.df_tax = pd.read_csv(file, sep=sep)

            def line_formatting(line):#, name_pos, taxonomy_pos):

                line[1] = line[1].replace('"', '')

                line_list = []
                try:
                    for taxon in line[1].split(";"):
                            # Removing crazy irregular taxon levels that RDP got creative with. Thanks Obama.
                            if taxon[-4:] == "idae" or taxon[-5:] == "ineae":
                                #print(f"skipping {entry}")
                                continue
                            taxon = taxon.split('(')[0]
                            if len(taxon) != 0:
                                line_list.append(taxon)
                except:
                    raise ValueError("It was not possible to process taxonomy line : {}. Exiting program. "
                                     "Check the taxonomy file.".format(line))

                if len(line_list) < 6:
                    filler = ["Unclassified"]*(6 - len(line_list))
                    line_list = line_list + filler

                if len(line_list) > 7:
                    print(line_list)

                self.taxDict[line[0]] = line_list
                return 1

            if format == "mothur":

                line = self.taxHeader
                line = sep.join(line)
                line_formatting(line)

                self.taxHeader = ["feature-name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

                #for line in f_tax:
                #    line_formatting(line, name_pos=0, taxonomy_pos=1)

            # Silva csv format
            if format == "silva":
                # Choosing the right database
                if database in ["EMBL", "GTDB", "LTP", "RDP", "SILVA"]:
                    if database == "EMBL":
                        tax_column = "lca_tax_embl_ebi_ena"
                        print("Using the EMBL database.")
                    if database == "GTDB":
                        tax_column = "lca_tax_gtdb"
                        print("Using the GTDB database.")
                    if database == "LTP":
                        tax_column = "lca_tax_ltp"
                        print("Using the LTP database.")
                    if database == "RDP":
                        tax_column = "lca_tax_rdp"
                        print("Using the RDP database.")
                    if database == "SILVA":
                        tax_column = "lca_tax_slv"
                        print("Using the SILVA database.")
                else:
                    raise Exception("Wrong database for Silva output was defined. Please select one from : GTDB, RDP, SILVA, "
                                    "LTP, EMBL")

                # Correcting redundant set of quotation marks
                newHeader = []
                for entry in self.taxHeader:
                    newHeader.append(entry.replace('"', ''))
                self.taxHeader = newHeader

                # Retrieving sequence identifier name
                if 'sequence_identifier' in self.taxHeader:
                    name_pos = self.taxHeader.index("sequence_identifier")
                elif "name" in self.taxHeader:
                    name_pos = self.taxHeader.index("name")
                else:
                    raise Exception("Could not find the Feature identifier column. Should be either 'sequence_identifier"
                                    " or 'name'. Correct that, please. Thank you.")
                database_pos = self.taxHeader.index(tax_column)

                #for line in f_tax:
                #    line_formatting(line, name_pos=name_pos, taxonomy_pos=database_pos)

                silva_tax = pd.read_csv(file, sep=sep)
                tax_cols = list(silva_tax.columns)
                tax_cols.pop(name_pos)
                tax_cols.pop(tax_cols.index(tax_column))

                silva_tax.drop(tax_cols, axis=1, inplace=True)

                import math

                for i,row in silva_tax.iterrows():
                    line_formatting(row)

                self.taxHeader = ["feature-id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]#, "Species"]

            if format == "dada2":
                self.df_tax = pd.read_csv(file)
                #name_list = []
                #for i,entry in enumerate(df_tax.iloc[:,0]):
                #    name_list.append("ASV_" + str(i))

                columns = list(self.df_tax.columns)
                columns[0] = "feature-id"
                self.df_tax.columns = columns
                # Renaming ASVs for convenience
                #df_tax.iloc[:,0] = name_list
                self.df_tax.fillna("Unclassified", inplace=True)

            # We will re-open the file later in the appropriate function
        if format == "Qiime2":
            self.qiime2_taxonomy(file, name_pos, sep)

    def qiime2_taxonomy(self, file, name_pos, sep):
        """
        Loading the taxonomy file in the Qiime2 format and re-formatting it
        :param file: Input taxonomy file.
            :type file: str
        :param name_pos: Position of the name column (column with the taxon codes - hash codes in Qiime2)
            :type file: int
        :param sep: Separator character used in the taxonomy file.
            :type sep: char
        :return: taxonomy based dictionary
        """
        # Removing confidence column from the tax header.
        #self.taxHeader.pop(self.taxHeader.index("Confidence"))
        #self.taxHeader.pop(self.taxHeader.index("Taxon"))
        self.taxHeader = ["feature-id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

        # Loading tables directly
        with open(file, 'r') as f_in:
            header = f_in.readline().split(sep)
            if "Taxon" in header:
                tax_val_index = header.index("Taxon")
            else:
                raise ValueError("Column 'Taxon' not found in the taxonomy file. This may be due to different version "
                                 "of Qiime2 used.")

            # This is to get an idea of distribution of taxonomic lengths (depth) within the tax file. For debugging.
            #tax_len_dict = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0}

            for line in f_in:
                line = line.split('\t')
                # Taxon identifier code (name)
                key = line[name_pos]
                tax_value_raw = line[tax_val_index].split(';')
                # Cleaning the names
                tax_value = []
                for value in tax_value_raw:
                    tax_value.append(value.split('_')[-1])

                # Filling in the empty unidentified levels
                # For now the taxonomic depth is hard-coded at 6
                tax_len_level = 6
                # This is again distribution of taxonomic lengths (depths) within the file. For debugging.
                #tax_len_dict[len(tax_value)] = tax_len_dict[len(tax_value)] + 1

                # First trimming the lengths
                if len(tax_value) > tax_len_level:
                    tax_value = tax_value[0:tax_len_level]
                # Now filling in the empty spaces
                if len(tax_value) < tax_len_level:
                    how_much_shorter = tax_len_level - len(tax_value)
                    tax_value = tax_value + ["Unidentified"]*how_much_shorter

                self.taxDict[key] = tax_value

            #import sys
            #sys.exit()

    def load_feature(self, file, sep, format, featureColumnName, representativeSeqFile, featureFilePickle,  threading, matchFeaturesToTaxonomy):
        """
        Loading the feature/otu file and creating the final multi indexed table

        :param format: Specification of an input format (Qiime2, dada2, custom)
            :type format: str
        :param file: Input mapping file.
            :type file: str
        :param sep: Separator character used in the feature/otu file.
            :type sep: char
        :param representativeSeqFile: Representative sequence file used to generate taxonomy. If taxonomy is generated
         externally using DADA2, ASVs should be renamed  to match that file exactly.
            :type file: str
        :param featureFilePickle: A flag that determines whether the loaded feature file is a pickle format.
            :type featureFilePickle: bool
        :param threading: Number of threads to be spawned by the process in moments where it is possible.
            :type threading: int
        :param matchFeaturesToTaxonomy: Boolean that turns off matching feature index  to the taxonomy index. This can
        be done to increase speed (default = True)
            :type matchFeaturesToTaxonomy: bool
        :return: Finalized multi indexed table containing both taxonomy and mapping file information.
        """

        # If file ending is .qza, extract table from a biom hdf5 file
        if file[-4:] == ".qza":
            table_hdf5  = self._load_qiime_zip(file, "feature-table.biom")
            with biom_open(table_hdf5) as f:
                table_biom = Table.from_hdf5(f)

                self.data = table_biom.to_dataframe(dense=True)
        elif format == "automatic":
            if featureFilePickle == True:
                raise ValueError("Automatic feature table parsing unavailable for pickle file. Please specify the format (Qiime2, dada2, custom)")
            with open(file, 'r') as f_in:
                line = f_in.readline()
                if line == "# Constructed from biom file\n":
                    format = "Qiime2"
                if len(line.split(',')[0].replace('"', '')) == 0: # or len(line.split(",")[1]) > 100:
                    format = "dada2"
                else:
                    format = "custom"

        ######################################
        # Custom format
        ######################################
        if format == "custom":
            self.data = pd.read_csv(file, sep=sep)
            if featureColumnName:
                try:
                    self.data.index = self.data[featureColumnName]
                    self.data = self.data.drop([featureColumnName], axis=1)
                except:
                    "Error: featureColumName does not match a colum in an input file"
            else:
                try:
                    self.data.index = self.data["#OTU ID"]
                    self.data = self.data.drop(["#OTU ID"], axis=1)
                    self.data.index = self.data.index.rename("ID")
                except:
                    try:
                        self.data.index = self.data["Unnamed: 0"]
                        self.data = self.data.drop(["Unnamed: 0"], axis=1)
                        self.data.index = self.data.index.rename("ID")
                    except:
                        raise ValueError(
                            "ID/OTU column not found. Check the OTU name format.")

            self.data.columns = self.data.columns.rename(self.mapSampleNames)

            '''
            # Renaming otus from OTU1 to OTU0001 etc.
            new_rownames = []
            # Retrieving number of characters in the longest feature name
            feature_name_max_len = self.data.index.str.len().max()
            # extracting a feature/otu name
            name_len = 0
            for letter in self.data.index[0]:
                if letter.isdigit():
                    break
                name_len = name_len + 1
            feature_name = self.data.index[0][0:name_len]
            for row in self.data.index:
                if len(row) < feature_name_max_len:
                    # Stands for "Otu..."
                    row_number = row[name_len:]
                    zeros = (feature_name_max_len - name_len - len(row_number)) * '0'
                    new_row = feature_name + zeros + row_number
                new_rownames.append(new_row)

            self.data.index = new_rownames
            
            for entry in self.taxDict:
                if len(entry) < feature_name_max_len:
                    row = entry[3:]
                    zeros = (feature_name_max_len - 3 - len(row)) * '0'
                    row = "Otu" + zeros + row
                    self.taxDict[row] = self.taxDict.pop(entry)
            '''

            # Sorting the otu table
            #self.data = self.data.sort_index()



        ######################################
        # Qiime2 format
        ######################################
        if format == "Qiime2":
            # Skipping the first line of the table and using the second as a header
            if featureFilePickle == True:
                # Loading a pickle format file
                self.data = pd.read_pickle(file)
            else:
                self.data = pd.read_csv(file, sep=sep, header=1)
            self.data.index = self.data["#OTU ID"]
            self.data = self.data.drop(["#OTU ID"], axis=1)

            self.data.columns = self.data.columns.rename(self.mapSampleNames)


        ######################################
        # DADA2 format
        ######################################
        if format == "dada2":
            if featureFilePickle == True:
                self.data = pd.read_pickle(file)
            else:
                self.data = pd.read_csv(file, sep=sep, index_col=0)
            self.data = self.data.transpose()
            # If we have some representative sequences to check
            if representativeSeqFile:
                rep_seq_dict = {}
                with open(representativeSeqFile, 'r') as f_rep:
                    for i, line in enumerate(f_rep):
                        if line[0] == '>':
                            feature_name = line[1:].rstrip()
                            sequence = f_rep.__next__()
                            sequence = sequence.rstrip()
                            if not sequence.isalpha():
                                raise ValueError(f"Representative sequence file seems "
                                                 f"corrupted on the line: {i}, sequence: {feature_name}.")
                        else:
                            # If fasta is spread over several lines, they are merged together.
                            # Performing a reverse lookup of a key (sequence) belonging to a specific value (seq name)
                            key = next(key for key, value in rep_seq_dict.items() if value == feature_name)
                            # Adding the new sequence at the end of the old one, saving the new key and deleting the old one
                            new_key = key + line.rstrip()
                            rep_seq_dict[new_key] = rep_seq_dict[key]
                            del rep_seq_dict[key]
                        # Saving sequences into the representative sequence dictionary
                        rep_seq_dict[sequence] = feature_name
                new_index = []
                for sequence in self.data.index:
                    try:
                        new_index.append(rep_seq_dict[sequence])
                    except:
                        raise ValueError(f"Sequence '{sequence}' not found in the representative set.")
                self.data.index = new_index
                print("Representative sequences checked successfuly and features renamed using  thir values.")
            # Sorting the rows to match the feature table

        ######################################
        # seqDataClass format
        ######################################
        if format == "seqDataClass":
            self.data = pd.read_csv(file, sep=sep, index_col=0)

        if self.taxDict:
            for entry in self.taxDict:
                if len(self.taxDict[entry]) != 6:
                    print(entry)
                    print("Tax Dictionary entry wrong legnth! It should be 6 and is : {}".format(len(self.taxDict[entry])))
                break
            # Adding a multi index
            self.df_tax = pd.DataFrame.from_dict(self.taxDict, orient='index')
            # Inserting column in the first place of the table
            self.df_tax.insert(0, "feature-id", self.df_tax.index)
            # This is because some taxnomic classifiers add even species level and we don't want that kind of stuff here
            if 6 in self.df_tax.columns:
                self.df_tax.drop(labels=6, axis=1, inplace=True)
            self.df_tax.columns = ["feature-id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

        sorter = list(self.data.index)
        sorterIndex = dict(zip(sorter, range(len(sorter))))
        self.df_tax["tax_rank"] = self.df_tax["feature-id"].map(sorterIndex)
        self.df_tax.sort_values(["tax_rank"], ascending=True, inplace=True)
        #del self.df_tax["Rank"]
        #df1 = df1.transpose()

        if len(sorter) > len(self.df_tax.index):
            #TODO: This bit is probably broken at the moment
            for entry in self.df_tax.index:
                sorter.pop(sorter.index(entry))
            raise ValueError("The following features were not present in the taxonomy : {}. Please add them.".format(sorter))
        elif len(sorter) < len(self.df_tax.index):
            print("Taxonomy has some redundant features. Attempting to prune these.")
            for entry in self.df_tax["feature-id"]:
                if entry not in sorter:
                    self.df_tax = self.df_tax[self.df_tax['feature-id'] != entry]


        # Sanity check that all features in sorter (feature file) are also in self.df_tax_index (taxonomy file)
        #for entry in self.df_tax.index:
        if matchFeaturesToTaxonomy:
            print("Matching features between feature file and taxonomy file.")
            tot_len = len(self.df_tax)
            import concurrent.futures

            def feature_match(entry_list):
                for entry in entry_list:
                    # if entry is just a number
                    try:
                        if int(entry) not in sorter:
                            raise ValueError("Feature file and taxonomy file features mismatch at {}".format(entry))
                    except:
                        if entry not in sorter:
                            raise ValueError("Feature file and taxonomy file features mismatch at {}".format(entry))
                    # Printing progress percentage, but only every 10th entry to not overwhelm the jupyter lab
                #print("finished entries 50")
                return 2500

            def chunks(lst, n):
                """Yield successive n-sized chunks from lst."""
                for i in range(0, len(lst), n):
                    yield lst[i:i + n]

            i = 0
            with concurrent.futures.ThreadPoolExecutor(max_workers=threading) as executor:
                futures = []
                for chunk in chunks(self.df_tax["feature-id"], n = 2500):
                    futures.append(executor.submit(feature_match, entry_list=chunk))
                    #print(chunk)
                for future in concurrent.futures.as_completed(futures):
                    i = i + future.result()
                    print(f"\r Progress {100 * (i / tot_len):.2f}%", end="\r", flush=True)
            print("finished matching")


            for i, entry in enumerate(self.df_tax["feature-id"]):
                # if entry is just a number
                try:
                    if int(entry) not in sorter:
                        raise ValueError("Feature file and taxonomy file features mismatch at {}".format(entry))
                except:
                    if entry not in sorter:
                        raise ValueError("Feature file and taxonomy file features mismatch at {}".format(entry))
                # Printing progress percentage, but only every 10th entry to not overwhelm the jupyter lab
                if i%10 == 0:
                    print(f"\r Progress {100*(i/tot_len):.2f}%", end="\r", flush=True)
            print("Features match.")

        # Renaming the data frame columns names (but not for DADA2, which doesn't produce taxDict
        #if self.taxDict:
        #    new_columns = self.taxHeader[0:len(self.df_tax.columns)]
        #    #new_columns.append("tax_rank")
        #    self.df_tax.columns = new_columns


        tax_index = pd.MultiIndex.from_frame(self.df_tax)

        # Setting the taxonomy indices
        self.data.index = tax_index

        # Arranging the values of the Multiindex according to the data frame

        #self.df_map = pd.DataFrame.from_dict(self.mapDict)
        #self.df_map = self.df_map.transpose()
        #self.df_map["sample-id"] = list(self.df_map.index)
        # Sorting based on the columns present in the feature table
        sorter = list(self.data.columns)
        sorterIndex = dict(zip(sorter, range(len(sorter))))
        for name in self.df_map[self.mapSampleNames]:
            if name not in sorter:
                raise ValueError(f" Sample names in mapping file are not matching sample names in feature table! {name}.")
        self.df_map["map_rank"] = self.df_map[self.mapSampleNames].map(sorterIndex)
        self.df_map.sort_values(["map_rank"], ascending=True, inplace=True)

        #del self.df_map["Rank"]
        # Creating a data frame from the map ditionary
        #df2 = df2.transpose()
        # Setting the names of columns
        #self.mapHeader.append("sample-id")
        #self.mapHeader.append("map_rank")
        #self.df_map.columns = self.mapHeader

        del self.df_map["map_rank"]

        map_index = pd.MultiIndex.from_frame(self.df_map)

        self.data.columns = map_index

    # Loading a file from qiime zip archive and returning its position
    def _load_qiime_zip(self, archive, file):
        with ZipFile(archive, 'r') as zip:
            for entry in zip.filelist:
                #print(entry.filename.split('/')[-1])
                if entry.filename.split('/')[-1] == file:
                    #print(entry.filename.split('/')[1])
                    table = zip.extract(entry.filename, path="/tmp/")
        return table

    def rarefy_to_even_depth(self, seqDepth, seed):
        """
        Rarefaction to even depth is a normalization method in which an equal number of sequences is randomlz retrieved
         from each sample to form a new feature/OTU table with equal sequencing depths. If some samples fall bellow a defined
         depth of sequencing treshold, they are removed from the dataset.

        :param seqDepth: A number of sequences that are randomly sampled from the dataset.
            :type seqDepth: int
        :param seed: A random seed value to maintain repeatability.
            :type seed: int
        :return: A rarefied feature/OTU table is saved in place of the original table.
        """
        rndSeed = np.random.RandomState(seed)

        sampleSums = np.sum(self.data, axis=0)
        nSpecies = self.data.shape[0]

        for i, column in self.data.iteritems():

            # If total number of sequences is lower than sequencing depth, remove the sample
            if column.sum() < seqDepth:
                sampleName = column.name[self.data.columns.names.index(self.mapSampleNames)]
                print("Sample {} was removed from the dataset because it contains insufficient amount of sequences ({}).".format(sampleName, column.sum()))

                self.remove_sample(category=self.mapSampleNames, sample=sampleName)
            else:
                relativeFreq = column.values/sampleSums[i]

                randomChoice = rndSeed.choice(a=nSpecies, size=seqDepth, p=relativeFreq)
                self.data[i] = np.bincount(randomChoice, minlength=nSpecies)

        # Removing rows that no longer contain any features/OTUs
        print("Removing {} features/OTUs that no longer appear in any sample.".format(sum(self.data.sum(axis=1) == 0)))
        self.data = self.data[self.data.sum(axis=1) != 0]

    def remove_sample(self, sample, category):
        catPosition = self.data.columns.names.index(category)

        if sample in self.data.columns.get_level_values(level=catPosition):

            self.data = self.data.drop(level=catPosition, columns=sample)
            self.data.columns = self.data.columns.remove_unused_levels()
        else:
            print("Error: Trying to remove {}, which is not available among samples.".format(sample))
            quit()

    def add_multiindex_index(self, dict, label):
        # Constructing a data frame from the supplied dictionary
        df2 = pd.DataFrame.from_dict(dict, orient="index", columns=[label])
        # Transforming multiindex into a data frame
        df1 = self.data.index.to_frame()
        df1.index = df1.ID
        # Joining the two data frames (adding the new column)
        df1 = df1.join(df2)
        # Creating a new MultiIndex from the new data frame
        newIndex = pd.MultiIndex.from_frame(df1)

        self.data.index = newIndex

    def add_multiindex_column(self, dict, label):
        df1 = self.data.columns.to_frame()

        nameColNumber = self.data.columns.names.index(self.mapSampleNames)
        df1.index = df1.iloc[:,nameColNumber]

        names = list(df1.Name)
        newDict = {}
        i = 0
        for name in names:
            newDict[name] = "Something" + str(i)
            i = i + 1

        df2 = pd.DataFrame.from_dict(dict, orient="index", columns=[label])

        df1 = df1.join(df2)

        newIndex = pd.MultiIndex.from_frame(df1)

        self.data.columns = newIndex

    def _add_multiindex_column(self, input_df, dict):
        # More general version of add_multiindex_column
        df1 = input_df.columns.to_frame()
        df1.reset_index(drop=True, inplace=True)

        for key in dict:
            df2 = pd.DataFrame.from_dict(dict)
            df1 = df1.join(df2)

        newIndex = pd.MultiIndex.from_frame(df1)

        input_df.columns = newIndex

        return input_df

    def add_feature_parameter(self, yamlParamDictFile, paramName):
        with open(yamlParamDictFile, 'r') as stream:
            try:
                paramDict = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        self.add_multiindex_index(dict=paramDict, label=paramName)

    def add_column_from_dict(self, yamlDictFile, label):
        with open(yamlDictFile, 'r') as stream:
            try:
                valueDict = yaml.load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        otuValues = self.data.index.get_level_values("OTU")

        output = []
        for value in otuValues:
            output.append(valueDict[value])

        idx = pd.IndexSlice

        #self.data.loc[:,idx[label, label]] = output
        self.data[label, label] = output

    def extract_features_csv(self, fileName, sep = ',', tax_level="feature-id"):
        #TODO: Remove this function, since it is replace by extracct table
        # Flattening the data frame's multiindex for easier downstream analysis
        output = self.data.copy()
        output.columns = output.columns.get_level_values(self.mapSampleNames)
        output.index = output.index.get_level_values(tax_level)


        output.to_csv(path_or_buf=fileName, sep=sep)

    def extract_rep_seq(self, filename, renameFeatures=True):
        # This function extract representative sequences and renames the
        with open(filename, 'w') as f_out:
            name = "feature"
            name_list = []
            for i,entry in enumerate(self.data.index.get_level_values("feature-id")):
                f_out.write('>' + name + str(i) + '\n')
                f_out.write(entry + '\n')
                name_list.append(name + str(i))

        if renameFeatures:
            # Renaming the oridinal dataset as well
            tax_index = self.data.index.to_frame(index = False)
            tax_index["feature-id"] = name_list
            tax_idx = pd.MultiIndex.from_frame(tax_index)
            self.data.index = tax_idx

    def extract_taxonomy(self, filename, sep=','):
        output = self.data.copy()
        output = output.index.to_frame(index=False)
        output.to_csv(filename, sep=sep, index=False)

    def extract_features(self, tax_level, map_level=None, file_name=False, sep='\t', show_output=True):
        """
        Simplifying a table by summarizing it to a defined taxonomic level and one or more levels of metadata.

        :param tax_level: Desired taxonomic level of the outcome. All the lower levels will be summarized in each
        category.
            :type tax_level: str
        :param map_level: Desired level or a list of levels to which the data should be subset. All finer categories
        will be added together
            :type map_level: char
        :return: simplified table
        """
        # If map level isn't defined, coerce automatically to sample name
        if not map_level:
            map_level = self.mapSampleNames

        if isinstance(map_level, str):
            df1 = self.data.copy()
            df1 = df1.sum(level=tax_level, axis=0)
            # Summarising over taxonomic groups and mapping categories provided
            df1 = df1.sum(level=map_level, axis=1)
        elif isinstance(map_level, list):
            '''
            Table simplification
            Group and sum each of the desired columns/categories effectively removing all other levels of metadata.
            Resulting data frame is left with only the columns levels in the map_level list.
            '''
            df1 = self.data.copy()
            df1 = df1.groupby(level=map_level, axis=1).sum()
            df1 = df1.sum(level=tax_level, axis=0)
            # Groupping all rows from the desired taxonomic level and summing them up.
            #df1 = df1.groupby(level=tax_level, axis=0).sum()

        if file_name:
            df1.to_csv(path_or_buf=file_name, sep=sep)
        if show_output:
            return df1
        else:
            return f"Features extracted in the file {file_name}"

    def colour_library(self, code, length=None):

        bright_palette = ["#ffd37f","#ffa700","#7f5300", # Yellow
                          "#7fc3a1","#008744","#004322", # Green
                          "#7fabf3","#0057e7","#002b73", # Blue
                          "#ea968f","#d62d20","#6b1610", # Red
                          "#c3a17f","#874400","#432200", # Brown
                          "#c37fa1","#870043","#430021", # Purple
                          "#b18db3","#641B67","#320d33", # Violet
			 ]

        def retrieve_longer_palettes(palette, length):
            output = []
            n = len(palette)
            for i in range(length):
                j = i % n
                output.append(palette[j])

            return output

        if length is not None:
            bright_palette = retrieve_longer_palettes(palette=bright_palette, length=length)

        output = {"bright_palette": bright_palette}

        if isinstance(code, list):
            return code
        else:
            if code not in output:
                raise ValueError("Colour code not found. Currently available options: {}".format(output))
            return output[code]

    def copy(self):
        return copy.deepcopy(self)

    def stacked_barplot(self,
                        tax_level,
                        map_level,
                        treshold_mean=0.001,
                        plotter="plotly",
                        x_axis_order=None,
                        show_x_axis_ordering = False,
                        subplot_titles=None,
                        colors="bright_palette",
                        plot_theme: object = None,
                        width=1200,
                        height=600):

        """
        Plotting a stacked barplot using plotly (for now) plotting engine and streamlined subsettig options.

        :param tax_level: Desired taxonomic level of the outcome. All the lower levels will be summarized in each
        category.
            :type tax_level: str
        :param map_level: Desired level or a list of levels to which the data should be subset. All finer categories
        will be added together
            :type map_level: char
        :param treshold_mean: Abundance treshold for considering a group 'Low abundance'. This helps to declutter
        the resulting plots. Should be a number lower then 1, default value is 0.001.
            :type treshold_mean: num
        :param x_axis_order: Forces order of the categorical elements of the X axis.
            :type x_axis_order: none or list of str
        :param subplot_titles: Names of the subplots.
            :type subplot_titles: none or list of str
        :param colors: Name of a color palette from the colour_library.
            :type colors: str
        :param plot_theme: Plotly theme for the plot. One of the following :'ggplot2', 'seaborn', 'simple_white',
        'plotly', 'plotly_white', 'plotly_dark', 'presentation', 'xgridoff', 'ygridoff', 'gridon', 'none']
            :type plot_theme: str
        :param width: Width of the final plot in pixels.
            :type width: int
        :param height: Height of the final plot in pixels.
            :type height: int
        :return: simplified table
        """

        df1 = self.extract_features(tax_level=tax_level, map_level=map_level)

        colSum = df1.sum(axis=0)

        df1 = df1.div(colSum)

        # If threshold mean is defined, summarise low abundance taxa in a 'low abundance' category
        if treshold_mean:
            lowAbundance = df1.loc[df1.mean(axis=1) < treshold_mean, :].sum()
            la = lowAbundance.to_frame(name="Low abundance")
            la = la.transpose()

            df1 = df1.loc[df1.mean(axis=1) > treshold_mean, :]
            df1 = df1.append(la)

            colSum = df1.sum(axis=0)
            df1 = df1.div(colSum)
        # If map level is just a string, use basic plotting with no subplots
        if isinstance(map_level, str):
            '''
            df1 = self.data.copy()
            # Summarising over taxonomic groups and mapping categories provided
            df1 = df1.sum(level=tax_level, axis=0)
            df1 = df1.sum(level=df1 = self.data.copy()map_level, axis=1)
            '''

            df1 = df1.sort_index(axis=0, ascending=False)

            if plotter == "plotly":

                #df1 = df1.reset_index()
                #df1 = df1.melt(id_vars="index")
                if x_axis_order:
                    df1 = df1[x_axis_order]
                    #cat_dtype = pd.api.types.CategoricalDtype(categories=x_axis_order, ordered=True)
                    #df1.columns = df1.columns.astype(cat_dtype)
                    #df1 = df1.reindex(sorted(df1.columns), axis=1)
                    #df1 = df1.sort_values(df1.columns)

                #import plotly.express as px
                #fig = px.bar(df1, x=map_level, y="value", color="index")
                import plotly.graph_objs as go

                color_list = self.colour_library(colors, length=df1.shape[0])
                fig = go.Figure()

                for row, colour in zip(df1.index, color_list):

                    try:
                        x = df1.columns.categories
                    except:
                        x = df1.columns.to_list()

                    fig.add_trace(go.Bar(x=x,
                                         y=df1.loc[row, :],
                                         name=row,
                                         marker_color=colour
                                         ),
                                  )

                fig.update_layout(title=dict(text=f'Relative abundance per {map_level}'),
                                  barmode='stack',
                                  width=width,
                                  height=height,
                                  template=plot_theme,
                                  xaxis_type="category")
            return fig
        # If map level is a list, proceed to a plotting with subplots
        if isinstance(map_level, list):

            if not ((x_axis_order is None) or isinstance(x_axis_order, dict)):
                raise Exception("Error: Type of x_axis_order must be either None (not included) or dictionary. Provided value is {}".format(type(x_axis_order)))

            if x_axis_order:
                # In order to make the column names ordered, we need to create a new column index to replace the old one
                # This is done by copying the old column
                df_cols = df1.columns.to_frame()
                df_cols.reset_index(drop=True, inplace=True)
                if show_x_axis_ordering:
                    return df_cols
                for key in x_axis_order:
                    cat_dtype = pd.api.types.CategoricalDtype(categories=x_axis_order[key], ordered=True)
                    df_cols[key] = df_cols[key].astype(cat_dtype)

                newIndex = pd.MultiIndex.from_frame(df_cols)

                df1.columns = newIndex

            # Loading a necessary plotly library for construction of subplots.
            from plotly.subplots import make_subplots
            import plotly.graph_objs as go
            #print(df1)

            map_dict = {}
            for level in map_level:
                map_dict[level] = {}
                map_level_vals = df1.columns.get_level_values(level).unique()
                for value in map_level_vals:
                    map_dict[level][value] = df1.xs(value, level=level, axis=1)

            ncols = len(df1.columns.get_level_values(map_level[0]).unique())

            if subplot_titles:
                if len(subplot_titles) != ncols:
                    raise ValueError("Wrong amount of provided supblot titles (subplot_titles)."
                                     " Provided {} titles, but attempting to plot {} plots.".format(len(subplot_titles), ncols))

                fig = make_subplots(rows=1, cols=ncols, subplot_titles=subplot_titles)
            else:
                block_names = df1.columns.get_level_values(map_level[0]).unique()
                block_names = block_names.to_list()
                fig = make_subplots(rows=1, cols=ncols, subplot_titles=block_names)

            color_list = self.colour_library(colors, length=df1.shape[0])

            # Looping over the first user supplied level
            for level in map_level:
                # Looping through the values of the first level supplied
                # If a provided set is categorical, use categories
                try:
                    categories = df1.columns.get_level_values(level).categories
                # Otherwise use just unique values
                except:
                    categories = df1.columns.get_level_values(level).unique()
                for i, category in enumerate(categories):
                    #print("Category: {}".format(category))
                    # Looping through taxonomic categories (map_dict...index is taxonomy informations)
                    # Adding colour information from the provided color_list

                    #df_temp = map_dict[level][category]
                    df_temp = df1.xs(key=category, axis=1, level=level)
                    df_temp.sort_index(axis=0, inplace=True, ascending=False)
                    for row, color in zip(df_temp.index, color_list):
                        try:
                            x = df_temp.columns.categories
                        except:
                            x = df_temp.columns.to_list()
                        #print("i: {}, item: {}, color: {}".format(i, item, color))
                        fig.add_trace(go.Bar(x=x,
                                             y=df_temp.loc[row,:],
                                             name=row,
                                             marker_color=color
                                             ),
                                      row=1,
                                      col=i+1
                                      )
                    fig.update_layout(title=dict(text='Relative abundance {} per {}'.format(map_level[0], map_level[1])),
                                      barmode='stack',
                                      width=width,
                                      height=height,
                                      template=plot_theme,
                                      xaxis_type="category",
                                      legend_traceorder="reversed")
                # Legend is being added for each subplot separately so we need to manually turn it off.
                # This is done by scanning for duplicate entries and turning their legends off.
                unique_legend = []
                for entry in fig.data:
                    if entry.name not in unique_legend:
                        unique_legend.append(entry.name)
                    else:
                        entry["showlegend"] = False
                print("end")

                return fig
            #print(map_dict)


    #def NMDS(self, plotter="plotly"):
    def lefse_data_export(self, filename, taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'], suppress_output = True):

        """
        Plotting a stacked barplot using plotly (for now) plotting engine and streamlined subsettig options.

        :param filename: Name of the output tab separated file.
            :type filename: str
        :param taxonomy_levels: List of taxonomy levels in order from general to specific.
            :type taxonomy_levels: list
        :return: lefse formatted table
        """

        print(f"Exporting data in lefse format into a file : {filename}")

        from skbio.stats.composition import multiplicative_replacement

        df1 = self.data.copy()
        df2 = pd.DataFrame(multiplicative_replacement(df1))
        df2.columns = df1.columns
        df2.index = df1.index

        cols = df2.columns
        col_sums = df2.sum(axis=0)
        output_df = pd.DataFrame(col_sums)
        output_df = output_df.transpose()
        output_df.index = ["column sums"]

        previous_levels = []
        output_index = []
        i = 0
        for level in taxonomy_levels:

            index_frame = df2.groupby(by=level, axis=0).sum()
            if previous_levels:
                output_index = []
                for entry in index_frame.index:
                    higher_levels = []
                    for previous_level in previous_levels:
                        higher_levels.append(df2.xs(key=entry, level=level, axis=0).index.to_frame(index=0).iloc[0][previous_level])

                    higher_levels.append(entry)
                    output_index.append("|".join(higher_levels))

            previous_levels.append(level)

            temp_df = df2.groupby(by=level, axis=0).sum() / col_sums
            if output_index:
                temp_df.index = output_index
            output_df = output_df.append(temp_df, ignore_index=False)




        print("done")

        output_df.to_csv(filename, sep = "\t")

        if not suppress_output:
            return(output_df)

if __name__ == "__main__":

    #plotly_template = "plotly_white"
    #tax_database = "RDP"

    filepath = "/home/user/qiime/thalium/meta_analysis/meyerhofer/"

    seq1 = seqObject(mappingFile=filepath + "mayerhofer_metadata.csv",
                     taxonomyFile=filepath + "meyerhofer_DNA_dada2_taxonomy.csv",
                     featureFile=filepath + "mayerhofer_feature_table_pickle.bz2",
                     mappingSep=',',
                     taxonomySep=',',
                     featureSep=',',
                     sampleNamesColumn="ENA-RUN",
                     featureFilePickle=True,
                     featureFormat="dada2",
                     threading=30,
                     matchFeaturesToTaxonomy=False
                     # taxonomyDatabase="RDP",
                     # representativeSeqFile="/home/user/qiime/Zygophilum/Taxonomy/rep_seq_cDNA.fasta"
                     )

    # fig.write_image("composition_DNA.pdf")

    #data.rarefy_to_even_depth(seqDepth=5241, seed=123)
