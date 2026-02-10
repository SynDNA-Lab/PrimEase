import logging
from dataclasses import dataclass, field

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



@dataclass
class Target:
    target_path: str
    record: SeqRecord = field(init=False)
    seq: str = field(init=False)
    

    def load_record(self) -> SeqRecord:
        ''' Returns the 0-th SeqRecord'''
        file_extension = self.target_path.split(".")[-1]
        print(file_extension)
        if file_extension in ["gb", "gbk"]:
            rec_type = "genbank"
        elif file_extension in ["dna"]:
            rec_type = "snapgene"
        elif file_extension in ["fa", "mpfa", "fna", "fsa", "fasta"]:
            rec_type = "fasta"
        else:
            raise ValueError(f"Unknown file extension {file_extension}")
            
        records = []
        for seq_record in SeqIO.parse(self.target_path, rec_type):
            records.append(seq_record)
        
        logging.info(f"Successfully loaded the target sequence record")
        return records[0]

    def write_fasta(self) -> None:
        SeqIO.write(self.record, "target.fasta", "fasta")

    def __post_init__(self) -> None:
        self.record = self.load_record()
        self.seq = str(self.record.seq).upper()
        self.write_fasta()
