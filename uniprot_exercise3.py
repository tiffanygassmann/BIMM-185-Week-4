import wget, re, urllib3.request, urllib3


#url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README_.txt"

#file_name = wget.download(url)


readme_open = open("README")
readme_read = readme_open.read()


def find_ID(README):

    genome_name = '(UP[0-9]+)'
    IDS =  re.findall(genome_name, README)

    return IDS[0:3]

print find_ID(readme_read)


#HEYY


IDS = find_ID(readme_read)

urls = []
for ID in IDS:
    urls.append("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/{}_*".format(ID))



print urls
