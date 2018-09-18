from Bio import Entrez, SeqIO
from beautifultable import BeautifulTable as BT
import re, sys, json

# function to index details in genbank features. note this function return a list
def __index_genbank_features__(gb_record, feature_type, qualifier) :
    answer = []
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                for value in feature.qualifiers[qualifier] :
                    answer.append(value)
    return answer

def ezsearch(enzyme, organism, email, rettype='json', pathprefix='./ezsearch_ret'):

    # provide email for entrez
    Entrez.email = email

    # look up for the gene name and annotation from NCBI
    term_str = "{}[orgn] AND {}".format(organism, enzyme)
    search_ret = Entrez.esearch(db='gene', term=term_str, sort='relevance')
    try:
        id_first = Entrez.read(search_ret)['IdList'][0]
    except IndexError:
        print("No match for that!")
        return None
    gene_info = Entrez.efetch(db='gene', id=id_first, retmode='text').read()
    # following is a piece of example of return for a gene query
    """

    1. lacZ
    beta-D-galactosidase [Escherichia coli str. K-12 substr. MG1655]
    Other Aliases: b0344, ECK0341, JW0335
    Annotation:  NC_000913.3 (363231..366305, complement)
    ID: 945006
    """
    gene_info_lines = gene_info.split('\n')
    gene_name = gene_info_lines[1][3:]
    try:
        gene_enzyme_name = re.findall(r"[ a-zA-Z0-9-]+(?= \[)", gene_info_lines[2])[0]
        gene_species = re.findall(r"(?<=\[)[ a-zA-Z0-9-\.\(\):]+(?=\])", gene_info_lines[2])[0]
        gene_annotation = re.findall(r"(?<=Annotation: )[ A-Za-z0-9\._]+", gene_info)[0]
        gene_start = re.findall(r"(?<=\()[\d]+(?=\.\.)", gene_info)[0]
        gene_end = re.findall(r"(?<=\.\.)[\d]+", gene_info)[0]
    except IndexError:
        print("Error in parsing the gene informatin. Here is the info from NCBI:{}".format(gene_info))
        return None

    # look up for whole sequence and other stuff according to the chromosome coordinate
    search_ret = Entrez.esearch(db='nucleotide', term=gene_annotation, sort='relevance')
    try:
        id_first = Entrez.read(search_ret)['IdList'][0]
    except IndexError:
        print("Matches for gene, but no match for the genome position!")
        return None
    gb_ret = Entrez.efetch(db='nuccore', id=id_first,\
                           seq_start=gene_start, seq_stop=gene_end, rettype='gb', retmode='text').read()

    # store genbank file
    gb_file = open(pathprefix+'.gb', 'w')
    gb_file.writelines(gb_ret)
    gb_file.close()

    # get other detailed info from genbank file
    gb_parse = SeqIO.parse(pathprefix+'.gb', 'gb')
    gb_record = list(gb_parse)[0]
    gene_chromosome = __index_genbank_features__(gb_record, 'source', 'chromosome')
    if len(gene_chromosome) == 0:
        gene_chromosome = None
    else:
        gene_chromosome = gene_chromosome[0]    #only one chromosome possible
    gene_note = __index_genbank_features__(gb_record, 'gene', 'note')
    if len(gene_note) == 0:
        gene_note = None
    else:
        gene_note = gene_note[0]
    gene_description = __index_genbank_features__(gb_record, 'CDS', 'function')
    if len(gene_description) == 0:
        gene_description = gene_note
    else:
        gene_description = gene_description[0]
    gene_EC = __index_genbank_features__(gb_record, 'CDS', 'EC_number')
    if len(gene_EC) == 0:
        gene_EC = None
    else:
        gene_EC = gene_EC[0]    #only one chromosome possible
    gene_nt_seq = gb_record.seq
    gene_tran_seq = __index_genbank_features__(gb_record, 'CDS', 'translation')
    #due to alternative splicing, there could be multiple translations

    # store all query data to dict
    gene = {}
    gene['name'] = gene_name
    gene['enzyme_name'] = gene_enzyme_name
    gene['species'] = gene_species
    gene['start'] = gene_start
    gene['end'] = gene_end
    gene['description'] = gene_description
    gene['chromosome'] = gene_chromosome
    gene['EC'] = gene_EC
    gene['nt_seq'] = str(gene_nt_seq)
    gene['tran_seq'] = gene_tran_seq

    # return required format
    if rettype == 'json':
        return json.dumps(gene)
    elif rettype == 'dict':
        return gene
    else:
        print('Wrong returning format! Returning json instead!')
        return json.dumps(gene)

def eztable(enzyme_info, width=100) :
    enzyme_table = BT(max_width=width)
    enzyme_table.append_row(['name', 'enzyme_name', 'species', 'start', 'end'])
    enzyme_table.append_row([enzyme_info['name'],\
                            enzyme_info['enzyme_name'],\
                            enzyme_info['species'],\
                            enzyme_info['start'],\
                            enzyme_info['end']])
    enzyme_table.append_row(['description', 'chromosome', 'EC', 'nt_seq', 'tran_seq'])
    try:
        description = enzyme_info['description']
        enzyme_table.append_row([description if description is None or len(description)<30 else description[:30]+'...',\
                               enzyme_info['chromosome'],\
                               enzyme_info['EC'],\
                               enzyme_info['nt_seq'][:30]+'...',\
                               str(len(enzyme_info['tran_seq']))+' types; '+enzyme_info['tran_seq'][0][:25]+'...'])
    except:
        print('Unexpected data structure! May lack in import field!')
        return None
    return enzyme_table
