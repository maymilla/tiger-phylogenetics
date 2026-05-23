from Bio import Entrez, SeqIO
import json
import os

Entrez.email = "maylayaffa@gmail.com"

TARGETS = [
    ("P_t_sumatrae", "AF053054.1"),   # Harimau Sumatera
    ("P_t_sondaica", "OQ601561.1"),   # Harimau Jawa (punah)
    ("P_t_tigris",   "AF053053.1"),   # Harimau Bengal
    ("P_t_altaica",  "AF053051.1"),   # Harimau Siberian/Amur
    ("P_t_corbetti", "AF053050.1"),   # Harimau Indochinese
    ("P_t_jacksoni", "EU184702.1"),   # Harimau Malayan
]

os.makedirs("data", exist_ok=True)


def fetch_sequence(label, accession_id):
  """
  Cari sekuens di NCBI berdasarkan accession id, ambil accession pertama yang ketemu, return SeqRecord atau None kalau gagal.
  """
  print(f"Mencari: {label} ...")
  
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
  
  for label, accession_id in TARGETS:
    record, acc_id = fetch_sequence(label, accession_id)
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
