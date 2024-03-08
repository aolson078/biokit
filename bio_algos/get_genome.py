from Bio import Entrez, SeqIO

def fetch_genome(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        return record
    except Exception as e:
        print(f"Error fetching genome: {e}")
        return None

def get_assembly_summary(id):
    try:
        esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
        esummary_record = Entrez.read(esummary_handle, validate=False)
        return esummary_record
    except Exception as e:
        print(f"Error getting assembly summary: {e}")
        return None



def get_assemblies(term, download=True, path='./summaries'):
    """Download genbank assemblies for a given search term."""
    Entrez.email = "aolson078@gmail.com"
    Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
    handle = Entrez.esearch(db="assembly", term=term, retmax='10')
    record = Entrez.read(handle)
    ids = record['IdList']

    links = []
    print(ids)
    for id in ids:
        summary = get_assembly_summary(id)
        summary = summary['DocumentSummarySet']['DocumentSummary'][0]
        if 'AssemblyAccession' in summary:
            accession = summary['AssemblyAccession']
            genome_record = fetch_genome(id)
            link = f"{path}/{accession}.gbk"
            try:
                SeqIO.write(genome_record, link, "genbank")
                links.append(link)
                if download:
                    print(f"Downloaded: {link}")
            except Exception as e:
                print(f"Invalid accession ID")
        else:
            print(f"Assembly accession not found for ID: {id}")

    return links


# Call the function with your search term
#links = get_assemblies("mycobacterium tuberculosis", download=True)
Entrez.email = "aolson078@gmail.com"
Entrez.apikey = "4a1d5a80996f3691a67335ded0b79299c708"
#get_assemblies("Human")


get_assemblies("Human")

n = ['47010088', '21654361', '21654351', '21654321', '21653681', '21653661', '21145841', '21145831', '21145821', '21145811']
#handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")

# genome_record = fetch_genome('21145811')
# print(genome_record)