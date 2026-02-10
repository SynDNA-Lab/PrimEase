#python3 genome_overlap_annotator.py genome.fasta ACGTACGTAC GTACTGATCG

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

class OverlappingFragmentAnnotator:
    def __init__(self,expected_structure_path,seq_frag1,seq_frag2):
        self.expected_structure_path = expected_structure_path
        self.seq_frag1 = seq_frag1
        self.seq_frag2 = seq_frag2
        self.output_path = "annotated_sequence_expected_structure.gb"
        self.expected_structure_record = None
        self.sequence_of_interest_seq = None

    def __post_init__(self): 
        self.run()

    def load_expected_structure(self): 
        file_name = self.expected_structure_path.lower()
    
        if any(ext in file_name for ext in ['fa', 'fasta']):
            with open(self.expected_structure_path, "r") as handle: 
                self.expected_structure_record = SeqIO.read(handle, "fasta") 
        elif any(ext in file_name for ext in ['gb', 'genbank']):
            with open(self.expected_structure_path, "r") as handle: 
                self.expected_structure = SeqIO.read(handle, "genbank") 

    def find_overlap(self, seq1, seq2):
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        max_overlap = ""
        
        #seq1 to seq2
        for i in range(1, min(len(seq1), len(seq2)) + 1):
            if seq1[-i:] == seq2[:i]:
                max_overlap = seq1[-i:]
        
        #seq2 to seq1
        for i in range(1, min(len(seq1), len(seq2)) + 1):
            if seq2[-i:] == seq1[:i]:
                if len(seq2[-i:]) > len(max_overlap):
                    max_overlap = seq2[-i:]
        
        return max_overlap if max_overlap else None

    def annotate_sequence_of_interest(self):
        overlap = self.find_overlap(self.seq_frag1, self.seq_frag2)

        if not overlap:
            print("No overlap found between the two fragments.")
            raise ValueError(
            "No overlap found between the two fragments. "
            "Please check that the sequences are correct and overlapping."
        )
        if len(overlap) < 10:
            raise ValueError(
                f"Overlap found is too short ({len(overlap)} bp). "
                "Overlap must be at least 10 base pairs."
            )
    
        expected_structure_seq = str(self.expected_structure_record.seq).upper()
        overlap_start = expected_structure_seq.find(overlap)

        if overlap_start == -1:
            print(f"Overlap '{overlap}' not found in genome.")
            raise ValueError ("The overlap between the two fragments is not in the DNA file given above. Check that both fragment are correct ")


        overlap_end = overlap_start + len(overlap)
        print(f"Overlap found in genome: {overlap} at positions {overlap_start} - {overlap_end}")

        self.sequence_of_interest_seq = overlap

        # Annotate the region in the genome
        sequence_of_interest_feature = SeqFeature(
            FeatureLocation(overlap_start, overlap_end),
            type="fragment",
            qualifiers={
                "note": "Overlapping region between input sequences",
                "label": "SequenceOfInterestOverlap"
            }
        )

        self.expected_structure_record.features.append(sequence_of_interest_feature)


    def write_genbank(self):
        gb_record = SeqRecord(
            Seq(str(self.expected_structure_record.seq)),
            id=self.expected_structure_record.id,
            name=self.expected_structure_record.name,
            description=self.expected_structure_record.description,
            annotations={"molecule_type": "DNA"},
            features=self.expected_structure_record.features
        )
        with open(self.output_path, "w") as handle:
            SeqIO.write(gb_record, handle, "genbank")
        print(f"Annotated GenBank written to {self.output_path}")

    def run(self):
        self.load_expected_structure()
        self.annotate_sequence_of_interest()
        self.write_genbank()
        
class SequenceAnnotator: 

    def __init__(self, expected_structure_path, search_seq, output_path="annotated_sequence_expected_structure.gb"): 
        self.expected_structure_path = expected_structure_path 
        self.search_seq = search_seq.upper() 
        self.output_path = output_path 
        self.expected_structure_record = None 
    
    def __post_init__(self): 
        self.run()

    def load_expected_structure(self): 
        file_name = self.expected_structure_path.lower()
    
        if any(ext in file_name for ext in ['fa', 'fasta']):
            with open(self.expected_structure_path, "r") as handle: 
                self.expected_structure_record = SeqIO.read(handle, "fasta") 
        elif any(ext in file_name for ext in ['gb', 'genbank']):
            with open(self.expected_structure_path, "r") as handle: 
                self.expected_structure_record = SeqIO.read(handle, "genbank") 

    def find_sequence(self): 
        expected_structure_seq = str(self.expected_structure_record.seq).upper() 
        index = expected_structure_seq.find(self.search_seq) 
        return index if index != -1 else None 

    def annotate_sequence(self, start_index): 
        end_index = start_index + len(self.search_seq) 
        feature = SeqFeature(FeatureLocation(start_index, end_index),type="fragment",qualifiers={"note": "Matched search sequence", "label": "FoundSequence"} ) 
        self.expected_structure_record.features.append(feature) 

    def write_genbank(self): 
        gb_record = SeqRecord(
            Seq(str(self.expected_structure_record.seq)),id=self.expected_structure_record.id,name=self.expected_structure_record.name, 
            description=self.expected_structure_record.description, 
            annotations={"molecule_type": "DNA"}, 
            features=self.expected_structure_record.features) 

        with open(self.output_path, "w") as output_handle: 
            SeqIO.write(gb_record, output_handle, "genbank") 

    def run(self): 
        self.load_expected_structure() 
        pos = self.find_sequence() 
        if pos is not None: 
            self.annotate_sequence(pos) 
            self.write_genbank() 
            return True, pos 
        else: 
            raise ValueError ("Pasted sequence of interest is not in the file given above")

 