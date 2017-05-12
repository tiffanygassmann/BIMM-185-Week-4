#Tiffany Gassmann
#BIMM 185
#week 4 - exercise 4 - Parsing SwissProt files


from Bio import SwissProt
import gzip,wget, itertools

url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.dat.gz"
file_name = wget.download(url)

#formats file for tab seperation
def format (list1):
    return str(list1).replace('[','').replace(']','\t').replace("'",'').replace(',','').replace('(','').replace(')','\t')

def file_parse():
    file = gzip.open("uniprot.gz")

    #Declaration of arrays which check for repitions
    non_rep_id = []
    non_rep_org = []
    non_rep_tax = []

    swiss_records = SwissProt.parse(file)

    for swiss_record in swiss_records:

        #NCBI ID
        id = swiss_record.taxonomy_id
        if id not in non_rep_id:
            non_rep_id.append(id)

        #ORGANISM NAME
        organism = (swiss_record.organism.strip('.'))
        if organism not in non_rep_org:
            non_rep_org.append(organism)

        #TAXONOMY
        taxonomy= (swiss_record.organism_classification)
        if taxonomy not in non_rep_tax:
            non_rep_tax.append(taxonomy)

    #ZIP arrays to column/tab seperated output
    for i in zip(non_rep_id, non_rep_org, non_rep_tax):
        print ("".join(map((str), list(format(i)))))

print file_parse()

