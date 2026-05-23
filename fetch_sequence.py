from Bio import Entrez, SeqIO
import json
import os

Entrez.email = "maylayaffa@gmail.com"

# daftar subspesies yang dicari
# format: (nama_display, query_ncbi)
TARGETS = [
  ("P_t_sumatrae",   "Panthera tigris sumatrae cytochrome b"),
  ("P_t_sondaica",   "Panthera tigris sondaica cytochrome b"),
  ("P_t_balica",     "Panthera tigris balica cytochrome b"), # nggak ada data, tapi tetap dimasukin biar lengkap
  ("P_t_tigris",     "Panthera tigris tigris cytochrome b"),
  ("P_t_altaica",    "Panthera tigris altaica cytochrome b"),
  ("P_t_corbetti",   "Panthera tigris corbetti cytochrome b"),
  ("P_t_jacksoni",   "Panthera tigris jacksoni cytochrome b"),
]

os.makedirs("data", exist_ok=True)


def fetch_sequence(label, query, max_results=3):
  """
  Cari sekuens di NCBI berdasarkan query, ambil accession pertama yang ketemu, return SeqRecord atau None kalau gagal.
  """
  print(f"Mencari: {label} ...")
  
  # Step 2a: Cari accession ID yang cocok
  handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
  record = Entrez.read(handle)
  handle.close()
  
  ids = record["IdList"]
  if not ids:
    print(f"Tidak ditemukan data untuk: {label}")
    return None, None
    
  accession_id = ids[0]  
  
  handle = Entrez.efetch(
    db="nucleotide",
    id=accession_id,
    rettype="fasta",       
    retmode="text"
  )
  seq_record = SeqIO.read(handle, "fasta")
  handle.close()
  
  if not seq_record.seq or len(seq_record.seq) == 0:
    print(f"Sekuens kosong untuk: {label}, skip")
    return None, None
  
  seq_record.id = label
  seq_record.description = f"{label} | cytochrome b | {accession_id}"
  
  print(f"[SUCCESS] {label}: {accession_id} ({len(seq_record.seq)} bp)")
  return seq_record, accession_id
  
def main():
  all_records = []
  accession_map = {} 
  
  for label, query in TARGETS:
    record, acc_id = fetch_sequence(label, query)
    if record:
      all_records.append(record)
      accession_map[label] = acc_id
  
  fasta_path = "data/raw_sequences.fasta"
  SeqIO.write(all_records, fasta_path, "fasta")
  print(f"\n{len(all_records)} sekuens disimpan ke {fasta_path}")
  
  json_path = "data/accession_ids.json"
  with open(json_path, "w") as f:
    json.dump(accession_map, f, indent=2)
  print(f"Accession IDs disimpan ke {json_path}")

if __name__ == "__main__":
  main()
