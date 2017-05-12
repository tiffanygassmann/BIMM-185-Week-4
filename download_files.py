import os

#Download Genomes from addresses
def findGenAdd(genome_ID,type):
    file = open("assembly_summary_refseq.txt",'r')
    address = ''

    #search in file line
    for line in file:
        if genome_ID in line:
            if type == "reference":
                if "reference" in line:
                    address = line.strip().split('\t')[-1]
                else:
                    address = line.strip().split('\t')[-1]
    return address

#use rsync to download addresses
def downloadReference(address, genome_ID):
    address = address.replace('ftp','rysnc')
    os.system("rsync -azvL {} ./{}".format(address,genome_ID))

#Call
if __name__ == '__main__':

    #Accession numbers ecoli, agro bactera
    ref1 = str(511145)
    ref2 = str(1435057)

    address1 = findGenAdd(ref1,'')
    downloadReference(address1,ref1)

