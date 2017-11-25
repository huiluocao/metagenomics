
import itertools as it
import os
import re
import sys
from hashlib import md5
from os import path
import pandas as pd
import csv

hit_row_regex = re.compile("^\s*\d\s*((\?)|(\!))\s*")
# ========
# Classes:
# ========

class HMMHit(object):
    def __init__(self, target_protein, hmm_name, score, e_value, hmm_from, hmm_to, ali_from, ali_to, hmm_length):
            self.target_protein = str(target_protein)
            self.hmm_name = str(hmm_name)
            self.score = float(score)
            self.e_value = float(e_value)
            self.hmm_from = int(hmm_from)
            self.hmm_to = int(hmm_to)
            self.hmm_overlap = self.hmm_to - self.hmm_from
            self.ali_from = int(ali_from)
            self.ali_to = int(ali_to)
            self.ali_length = self.ali_to - self.ali_from
            self.hmm_coverage = float(self.hmm_overlap) / float(hmm_length)

    def __repr__(self):
            out_string = """
===========================
Target Protein HMM Name Score e_value hmm_from hmm_to hmm_overlap ali_from ali_to ali_length hmm_coverage
%s %s %f %f %d %d %d %d %d %d %f
""" % (self.target_protein, self.hmm_name, self.score, self.e_value, self.hmm_from, self.hmm_to,
                    self.hmm_overlap, self.ali_from, self.ali_to, self.ali_length, self.hmm_coverage)

            return out_string

    def __str__(self):
            return self.__repr__()

    def get_md5(self):
            """
            Concatenates all object properties as strings and hashes them using an MD5 checksum.

            :return: MD5 checksum of all object properties.
            """
            hash_string = "".join([str(x) for x in self.__dict__.values()])  # Join all attributes into a single string.
            hash_md5 = md5(hash_string.encode('utf-8')).hexdigest()  # Create md5 hash.
            return hash_md5

# ------------------------------------------------------------------------------------------------------------
def get_hmm_length(file_name):
    with open("C:/Users/Huiluo/Desktop/HKU-microbiology/microbiome/tools/Resfam/Resfams-full.hmm/length.txt", "rU") as inFile:
        for line in inFile.readlines():
            if file_name in line:
                hmm_length = int(line.split("\t")[1])
    """
    Gets HMM Length.

    :param hmm_path: path to the HMM file.
    :return: The length of the hmm.
    """
    return hmm_length

# ------------------------------------------------------------------------------------------------------------
def filter_hmm_hit_list(hmm_hit_list, e_value_cutoff="1e-10", hmm_coverage=0.7, max_align_overlap=0.5):
	"""
	Filters HMM gits by E-Value, Coverage and Overlap between hits.

	:param hmm_hit_list: List of HMM hit objects.
	:param e_value_cutoff: The E-Value cutoff for hits.
	:param hmm_coverage: The HMM coverage cutoff for hits.
	:param max_align_overlap: The maximum overlap percentage between overlapping HMM hits.
	:return: List of filtered HMM hit objects.
	"""
	hmm_hit_list = [hit for hit in hmm_hit_list if hit.e_value < float(e_value_cutoff)]  # Filters hits by E-value.

	i = 0
	while i < (len(hmm_hit_list) - 1):
		hit_one = hmm_hit_list[i]  # Current Row in hit table.
		hit_two = hmm_hit_list[i + 1]  # Row below.
		if hit_one.target_protein == hit_two.target_protein:
			overlap_between_hits = hit_one.ali_to - hit_two.ali_from
			if overlap_between_hits > 0:
				# If the overlap is greater than 50% of either alignment.
				if ((float(overlap_between_hits) / float(hit_one.ali_length)) > max_align_overlap) or (
							(float(overlap_between_hits) / float(hit_two.ali_length)) > max_align_overlap):

					if hit_one.e_value < hit_two.e_value:
						hmm_hit_list.remove(hit_two)
					else:
						hmm_hit_list.remove(hit_one)
					i -= 1  # Resets list index.
		i += 1

	hmm_hit_list = [hit for hit in hmm_hit_list if hit.hmm_coverage > hmm_coverage]  # Filters by Query Coverage.
	return hmm_hit_list

    
# ------------------------------------------------------------------------------------------------------------
def parse_hmmsearch_results(hmm_results_string, hmm_name, hmm_length):
    """
    Parses HMM searches text output and generates a two dimensional array of the domain alignments results.

    :param hmm_results_string: hmmsearch text output as a string.
    :param hmm_name: The name of the HMM file.
    :param hmm_length: The length of the HMM file.
    :return: List of HMM hit objects.
    """

    hmm_hit_list = hmm_results_string.split(">>")  # Splits output at domain alignments.
    del hmm_hit_list[0]  # Removed string content above first >> which does not need to be iterated through.

    # Removes alignment visualization by spiting on the visualization header
    # and keeping all that is above it (first element) which is the HMM domain table.
    hmm_hit_list_cleaned = [x.split("Alignments")[0] for x in hmm_hit_list]

    hmm_object_list = []
    for protein in hmm_hit_list_cleaned:
        domain_table = protein.splitlines()
        target_protein_name = domain_table[0].split()[0]  # Gets target protein accession from domain table header row.

        for row in domain_table:
            if hit_row_regex.match(row):  # If row is a domain table line.
                column = row.split()

                hmm_hit_object = HMMHit(
                                    target_protein=target_protein_name,
                                    hmm_name=hmm_name,
                                    score=column[2],
                                    e_value=column[5],
                                    hmm_from=column[6],
                                    hmm_to=column[7],
                                    ali_from=column[9],
                                    ali_to=column[10],
                                    hmm_length=hmm_length)

                hmm_object_list.append(hmm_hit_object)
                
    return hmm_object_list


def get_hmm_name(input_item):
    global hmm_name
    for line in input_item.split("\n"):
        if line.startswith("Accession"):
            hmm_name = line.split(":")[1].strip()
    return hmm_name

def get_filtered_results(entry):
    hmm_name = get_hmm_name(entry)
    hmm_length = get_hmm_length(hmm_name)
    hmm_hit_list = parse_hmmsearch_results(entry, hmm_name, hmm_length)
    return hmm_hit_list

def de_duplicate(all_hits):
    df = pd.DataFrame(all_hits)
    print(df)
    ndf= pd.DataFrame(df[0].str.split("\s").tolist())
    print(ndf)
    ndf1 = ndf.loc[ndf.groupby([0], as_index=False)[2].idxmax()]
    return ndf1            
# ------------------------------------------------------------------------------------------------------------
def main():
    with open ("Andes-1_arg_hmmr1.txt") as input_handel:
        output_file= open("Andes-1_hmmr_filtered.txt","a")
        all_hits = []
        contents = input_handel.read()
        for entry in contents.split("//\n"):  
            hmm_hit_list = get_filtered_results(entry)                      
            if len(hmm_hit_list) != 0:
                filtered_hmm_hit_list = filter_hmm_hit_list(hmm_hit_list)
                if len(filtered_hmm_hit_list) !=0:
                    print(filtered_hmm_hit_list )
                    for i in filtered_hmm_hit_list:
                        all_hits.append(str(i).split("\n")[3])
                        output_file.write(str(i).split("\n")[3]+"\n")
        output_file.close()
        de_duplicate_hits = de_duplicate(all_hits)
        print(de_duplicate_hits)
        with open ("Andes-1_arg_hmmr1_cleaned.txt","w") as output_handel:
            de_duplicate_hits.to_csv(output_handel, index =False, sep ="\t")
            

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

        
	
            
