import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from target import Target
from roi import ROI
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


@dataclass
class Selection(ABC):
    target: Target
    regions: list[ROI] = field(init=False)

    @abstractmethod
    def selection(self) -> list[ROI]:
        pass

@dataclass
class SequenceOfInterest(Selection):
    target: Target
    regions: list[ROI] = field(init=False)
    def __post_init__(self) -> None:
        self.selection()

    def selection(self) -> list[ROI]:
        targets = []
        for feature in self.target.record.features : 
            if (feature.type == "misc_feature" and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == "fragment") or feature.type == 'fragment' :
                targets.append([feature.location.start,feature.location.end])
        rois = []
        for idx, trg in enumerate(targets):
            #breakpoint()
            rois.append(ROI(
                name=f"target_{trg[0]}",
                target_start = trg[0],
                target_end = trg[1],
                expected_structure_sequence= self.target.seq
            ))
        print(rois[0])
        self.regions = rois

@dataclass
class GenbankfileHandling(Selection) : 
    target : Target
    genbank_path : str
    expected_structure_record : SeqRecord = field(init=False)
    regions: list[ROI] = field(init=False)
    list_sequence_of_interest : list[str] = field(init=False)
    def __post_init__(self) : 
        self.list_sequence_of_interest = []
        self.load_expected_structure()
        self.selection()
    
    def load_expected_structure(self):
        with open(self.genbank_path, "r") as handle:
            self.expected_structure_record = SeqIO.read(handle, "genbank")
        print(f"Expected structure loaded: {self.expected_structure_record.id} ({len(self.expected_structure_record.seq)} bp)")

    def find_overlap(self, seq1, seq2,min_overlap=10):
        seq1, seq2 = seq1.upper(), seq2.upper()
        max_len = min(len(seq1), len(seq2))
        
        for i in range(max_len, 0, -1):  # Start from the longest possible overlap
            if seq1[-i:] == seq2[:i] and i > min_overlap:
                return seq1[-i:]
        
        return None


    def selection(self):

        sequences_of_interest = []
        overlapping_features = set()

        features = [feature for feature in self.expected_structure_record.features if ((feature.type == "misc_feature" and 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] == "fragment") or feature.type == 'fragment')]

        if not features:
            raise ValueError("GenBank file isn't annotated with 'misc_feature' elements.")

        seq = self.expected_structure_record.seq

        for i in range(len(features)):
            feature_i_seq = features[i].extract(seq)

            overlap_found = False
            for j in range(i + 1, len(features)):
                feature_j_seq = features[j].extract(seq)
                
                max_overlap = self.find_overlap(feature_i_seq, feature_j_seq)

                if max_overlap:
                    overlap_found = True
                    start = seq.find(max_overlap)
                    end = start + len(max_overlap)

                    if start == -1:
                        print(f"Overlap '{max_overlap}' not found in genome.")
                        continue

                    print(f"Overlap found in genome: {max_overlap} at positions {start} - {end}")
                    sequences_of_interest.append([start, end])
                    overlapping_features.add(i)
                    overlapping_features.add(j)
            
        for idx, feature in enumerate(features):
            if idx not in overlapping_features:
                start = int(feature.location.start)
                end = int(feature.location.end)
                sequences_of_interest.append([start, end])
                print(f"No overlap for feature {idx}, adding its range: {start}-{end}")

        rois = []

        if sequences_of_interest :
            for bound in sequences_of_interest:
                rois.append(ROI(
                    name=f"target_{bound[0]}",
                    target_start = bound[0],
                    target_end = bound[1],
                    expected_structure_sequence = str(self.expected_structure_record.seq)
                ))
                self.list_sequence_of_interest.append(str(self.expected_structure_record.seq)[bound[0]:bound[1]])
        self.regions = rois
        print(rois)