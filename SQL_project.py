#Assignment 4 BIMM 185
#Tiffany Gassmann
#A10684434

from Bio import SeqIO
import gzip, re

#---------------------------------E-COLI----------------------------------------------------------------------

#Creates File to test each sectiona and to view
def file_handle_ecoli():

    file = gzip.open("ecoli_genomic.gbff.gz")

    gb_record = SeqIO.read(file,"genbank")


    return gb_record

#Creates file which we use BioPython Methods on
def file_parse_ecoli():
    file = gzip.open("ecoli_genomic.gbff.gz")

    record = SeqIO.parse(file, "genbank")

    return record

#---------------------------------AGRO-BACTERIA----------------------------------------------------------------------
#Creates file which we use biopy6thon methods on
def file_handle_agro():

    file = gzip.open("agrobac_genomic.gbff.gz")
    records = list(SeqIO.parse(file, "genbank"))
    gb_record = records[0]

    return gb_record

#Creates file which we use BioPython Methods on
def file_parse_agro():
    file = gzip.open("agrobac_genomic.gbff.gz")
    records = list(SeqIO.parse(file, "genbank"))

    return records

#----------------------------END FILE OPEN/PARSE-----------------------------------------------------------------------

#Method to remove ectra characters after converted to string

def format (list1):
    return str(list1).replace('[','').replace(']','').replace("'",'').replace(',','').replace('(','').replace(')','')

#------------------------------------------------------------------------------------------------------------------------
#GENOMIC TABLE SQL FORMAT:
# GENOME_ID  TAX_ID  SHORT_NAME  LONG_NAME   SIZE_BP DOMAIN ACCESSION   RELEASE_DATE

#Extract source info -  TAX_ID, NAME
def source_info(record):

    rec = next(record)
    for f in rec.features:
        if f.type == 'source':
            tax_ID = (f.qualifiers['db_xref'])
            long_name = f.qualifiers['organism']

    num_only= "[0-9]+"
    for ID in tax_ID:
        tax = re.findall(num_only, ID)

    return format(tax),format(long_name)

#function to format and create tab seperated genome table for my sql
def create_genomic_table(read,parse):

    genome_id = str(read.id)
    tax_id,long_name = source_info(parse)
    short_name = ""
    seq_length = str(len(read))
    domain = "Bacteria"
    accession =  read.annotations['accessions'][0]
    release_date = read.annotations['date']

    print "\t".join([genome_id, tax_id,short_name, long_name, seq_length,domain, accession,str(release_date)])

#variables for diferent file reads
read_record_ecoli = file_handle_ecoli()
parse_record_ecoli = file_parse_ecoli()

read_record_agro = file_handle_agro()
parse_record_agro = file_parse_agro()


#Ecoli_genome = create_genomic_table(read_record_ecoli,parse_record_ecoli)
#Agro_genome = create_genomic_table(read_record_agro,parse_record_agro)

#-------------------------------------END GENOME-TABLE----------------------------------------------------------------

#REPLICON TABLE SQL FORMAT:
# REPLICON_ID   GENOME_ID NAME  GENE_NUMBER TYPE    STRUCTURE


def create_replicon_table(records):
    #number_replicons = len(records)
    names = []
    types = []
    structures = []
    counts = []

    #iterate over all replicons - agrobactera
    for record in records:

        #name AGRO
        name = record.description

        #name e. Coli
        #name = record.features
        names.append(name)

        #reg ex to find structure/type
        structure_type = record.description

        #type = chromosome/plasmid
        type_ex = "(?:plasmid|chromosome)"
        #structure = linear/circular
        structure_ex = "(?:linear|circular)"

        types.append(re.findall(type_ex,structure_type))
        structures.append(re.findall(structure_ex,structure_type))

        #Check for "Complete Genome" -- means chromosome for e. Coli

        if record.description.contains("complete genome"):
            type_ex = "Chromosome"

        #count number of genes per replicon
        CDS_count = 0
        for f in record.features:

            if f.type == 'CDS':
                CDS_count += 1
        counts.append(CDS_count)

        #print "\t".join([format([names, gene_numbers, types, structures])])

    #variable tests
    print "Names:", names
    #print "Gene_Numbers", gene_numbers
    print "Types:" ,types
    print "structures",structures
    print "CDS_count", counts


#Create Replicon Table command
file = gzip.open("agrobac_genomic.gbff.gz")
#file = gzip.open("ecoli_genomic.gbff.gz")

records = list(SeqIO.parse(file, "genbank"))
#record = file_parse_ecoli()

#create_replicon_table(records)
#create_replicon_table(records)

#----------------------------END-REPLICONS----------------------------------------------------------------------------

#REAL_GENES TABLE SQL FORMAT:
#GENE_ID GENOMIC_ID REPLICON_ID LOCUS_TAG NAME  STRAND  EXON_NUMBER LENGTH_BP   PRODUCT
def cds_info_genetab(records):

    #ecoli
    genome_id = 1
    replicon_id = 1
    gene_id = 0

    #agrobac
    #genome_id = 2
    #gene_id = 4319

    for record in records:
        #increment replicon ID for each record in Agrobacteria
        replicon_id += 1
        #go through each entry
        #rec = next(record)
        for f in record.features:
            #ALL other info
            if f.type == 'CDS':

                #increment Gene_ID
                gene_id += 1

                #Location and Strand Direction taken from Location NOT qualifiers
                locations_strands = (str(f.location))
                #regular expression to seperate the start, stop, direction
                reg = "\[([0-9]+):([0-9]+)\]\((.)\)"
                loc_str = re.findall(reg, locations_strands)

                for start, stop, strand in loc_str:
                    lengths =''.join(str(int(stop) - int(start)))
                    strands = ''.join(str(strand))

                #Genes, Locus Tags, Synonyms taken from Qualifiers
                #Note: Gene Name for AgroBacteria = Locus Tag
                gene_name = "\t".join(f.qualifiers['gene'])
                #gene_name = "\t".join(f.qualifiers['locus_tag'])
                locus_tags = '\t'.join(f.qualifiers['locus_tag'])


                #check for key error if no product exists: return "pseduo"

                if 'product' in f.qualifiers:
                    #prod_ex = "^.*(?=(\:))"
                    product = (f.qualifiers['product'])

                else: product = "Pseudo"

                products = ''.join(product)

               #exon number

                exon = f.location
                exon_count = str(exon).count("[")


                print '\t'.join([str(gene_id),str(genome_id), str(replicon_id), locus_tags,gene_name,strand, str(exon_count), lengths, products])

#ECOLI PARSE
#records = list(file_parse_ecoli())
#cds_info_genetab(records)

#AGRO PARSE
#file = gzip.open("agrobac_genomic.gbff.gz")
#records = list(SeqIO.parse(file, "genbank"))
#cds_info_genetab(records)


#------------------------------END-GENE-TABLE---------------------------------------------------------------------------

#FUNCTIONS SQL FORMAT:
#GENE ID    FUNCTION

def cds_info_fuctiontab(records):
    #ECOLI START
    #gene_id = 0

    #AGRO Start
    gene_id = 4319

    for record in records:
        for f in record.features:
            # ALL other info
            if f.type == 'CDS':
                gene_id+=1
                # Function taken from Qualifiers
                if 'function' in f.qualifiers:
                    function = (format(f.qualifiers['function']))
                else:
                    function = "Pseudo"

                functions = ''.join(function)
                print '\t'.join([str(gene_id),functions])

#ECOLI PARSE
#records = list(file_parse_ecoli())

#AGRO PARSE
#file = gzip.open("agrobac_genomic.gbff.gz")
#records = list(SeqIO.parse(file, "genbank"))


#cds_info_fuctiontab(records)
#-----------------------------END-FUNCTION-TABLE------------------------------------------------------------------------

#EXTERNAL REFERENCES: SQL FORMAT:
#GENE_ID        EXTERNAL_DATABASE   EXTERNAL_ID
def cds_info_exreftab(records):

    # ECOLI START
    #gene_id = 0

    # AGRO Start
    gene_id = 4319
    for record in records:
        for f in record.features:

            if f.type == 'CDS':
                gene_id +=1
                if 'db_xref' in f.qualifiers:
                    # GENE ID taken from Qualifiers

                    ex_ref = f.qualifiers['db_xref']

                    data_bases = ex_ref[0:3]
                    # data_bases[0] = GI - gene identifier number
                    # data_bases[1] = ASAP -A Systematic annotation Package for Community Analysis of Genomes
                    # data_bases[2] = Uniport/Swiss-Prot -	A Systematic Annotation Package for Community Analysis of Genomes

                    for db in data_bases:

                        ex_db_rex = "^[^\:]*"
                        ex_db = format(re.findall(ex_db_rex,db))
                        ex_ID_rex = '\:(.*)'
                        ex_ID = format(re.findall(ex_ID_rex,db))


                        print '\t'.join([str(gene_id),ex_db,ex_ID])


#ECOLI PARSE
#records = list(file_parse_ecoli())

#AGRO PARSE
#file = gzip.open("agrobac_genomic.gbff.gz")
#records = list(SeqIO.parse(file, "genbank"))

#cds_info_exreftab(records)

#-----------------------------END-EXTERNAL-REF-TABLE--------------------------------------------------------------------

#SYNOYMS: SQL FORMAT:
#GENE_ID        SYNONYM

def cds_info_syntab(records):
    # ECOLI START
    #gene_id = 0

    # AGRO Start
    gene_id = 4319

    for record in records:
        for f in record.features:

            if f.type == 'CDS':
                gene_id +=1


                #ECOLI
                if 'gene_synonym' in f.qualifiers:
                    gene_synonyms = str(f.qualifiers['gene_synonym']).split(';')

                #AGRO
                elif 'old_locus_tag' in f.qualifiers:
                    gene_synonyms = str(f.qualifiers['old_locus_tag']).split(';')

                for syn in gene_synonyms:
                    format(syn)
                    print '\t'.join([str(gene_id), format(syn)])


#ECOLI PARSE
#records = list(file_parse_ecoli())

#AGRO PARSE
#file = gzip.open("agrobac_genomic.gbff.gz")
#records = list(SeqIO.parse(file, "genbank"))

#cds_info_syntab(records)

#-----------------------------END-SYNONYM-TABLE--------------------------------------------------------------------

#EXONS: SQL FORMAT:
#GENE_ID       EXON     LEFT_POSITION   RIGHT_POSITION  LENGTH_BP


def cds_info_exontab(records):

    # ECOLI START
    #gene_id = 0
    #exon = 0

    # AGRO Start
    gene_id = 4319
    exon = 4319

    exons = []

    for record in records:
        for f in record.features:

                if f.type == 'CDS':
                    exon +=1
                    gene_id +=1

                    # exon number


                    left_right = (str(f.location))
                    # regular expression to seperate the start, stop, direction
                    reg = "\[([0-9]+):([0-9]+)\]\((.)\)"
                    loc_str = re.findall(reg, left_right)

                    for start, stop, strand in loc_str:

                        lengths = ''.join(str(int(stop) - int(start)))
                        lefts = ''.join(str(start))
                        rights = ''.join(str(stop))


                    print '\t'.join([str(gene_id),str(exon), lefts,rights,lengths])


#ECOLI PARSE
#records = list(file_parse_ecoli())

#AGRO PARSE
file = gzip.open("agrobac_genomic.gbff.gz")
records = list(SeqIO.parse(file, "genbank"))


cds_info_exontab(records)





#-----------------------------END-EXON-TABLE--------------------------------------------------------------------



#Tests to view each section in the larger file and export each to new text file


def extract_source(gb_record):
    my_source = gb_record.features[0]
    return my_source

def extract_gene(gb_record):
    my_gene = gb_record.features[1]
    return my_gene

def extract_cds(gb_record):
    my_CDS = gb_record.features[2]
    return my_CDS


#Gene Section
#gene = extract_gene(file_handle_ecoli())
#with open('gene.txt', 'w') as file:
#    print >> file, gene
#CDS Section
cds = extract_cds(file_handle_agro())
#with open('cds_agro.txt', 'w') as file:
#    print >> file, cds
#Source Section
#source = extract_source(file_handle_ecoli())
#with open('source.txt', 'w') as file:
 #   print >> file, source