import requests
import mysql.connector

def connect_ucsc():
    return mysql.connector.connect(
        user='genome',
        host='genome-mysql.soe.ucsc.edu',
        database='hg38'
    )

def fetch_all_transcripts(gene_name):
    conn = connect_ucsc()
    cur = conn.cursor(dictionary=True)
    query = """
        SELECT * FROM refGene
        WHERE name2 = %s AND cdsStart < cdsEnd
    """
    cur.execute(query, (gene_name,))
    transcripts = cur.fetchall()
    conn.close()
    return transcripts

def select_longest_cds_transcript(transcripts):
    if not transcripts:
        return None
    return max(transcripts, key=lambda tx: tx['cdsEnd'] - tx['cdsStart'])

def fetch_genomic_sequence(chrom, start, end):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    res = requests.get(url)
    res.raise_for_status()
    return res.json()["dna"].upper()
