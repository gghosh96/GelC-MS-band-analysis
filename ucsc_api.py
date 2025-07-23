import requests
import mysql.connector

"""Connect to the UCSC public MySQL genome database."""
def connect_ucsc():
    return mysql.connector.connect(user='genome', host='genome-mysql.soe.ucsc.edu', database='hg38')

"""
Fetch transcript info for a gene from UCSC refGene table.

Args:
    gene_name (str): The gene symbol

Returns:
    dict: One transcript record with highest exonCount, or None if not found.

"""
def fetch_transcript_info(gene_name):
    conn = connect_ucsc()
    cur = conn.cursor(dictionary=True)
    query = """
        SELECT * FROM refGene
        WHERE name2 = %s AND cdsStart < cdsEnd
        ORDER BY exonCount DESC LIMIT 1
    """
    cur.execute(query, (gene_name,))
    row = cur.fetchone()
    conn.close()
    return row
    
"""
Fetch genomic DNA sequence from UCSC API.

Args:
    chrom (str): Chromosome name, e.g. 'chr1'
    start (int): Start coordinate (0-based)
    end (int): End coordinate (exclusive)

Returns:
    str: DNA sequence in uppercase letters (A,T,C,G)

"""
def fetch_genomic_sequence(chrom, start, end):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    res = requests.get(url)
    res.raise_for_status()
    return res.json()["dna"].upper()