import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
import matplotlib.pyplot as plt


class VisualizeGenbank(BiopythonTranslator):

    def __init__(self,path): 
        super().__init__()
        self.path = path 
        self.color_map = {
            "CDS": "#66c2a5",
            "gene": "#fc8d62",
            "misc_feature": "#8da0cb",
            "rRNA": "#e78ac3",
            "tRNA": "#a6d854",
            "source": "#ffd92f",
        }
        self.output_path = "Genbank_vizualisation.png"

    def draw(self) :
        record = SeqIO.read(self.path, "genbank")
        graphic_record = self.translate_record(record)
        ax, _ = graphic_record.plot(figure_width=12)
        ax.set_title("Position of the Sequence Of Interest (SOI)", fontsize=14)


        handles = [
            patches.Patch(color=color, label=ftype)
            for ftype, color in self.color_map.items()
        ]
        ax.legend(handles=handles, title="Feature Type", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()
        plt.savefig(self.output_path, dpi=300)
        plt.show()


    
    def compute_feature_color(self, feature):
        # Pick color based on type, default gray if unknown
        return self.color_map.get(feature.type, "#cccccc")
    
    def compute_feature_label(self, feature):
        # Label only those with locus_tag == "fragment"
        if ("locus_tag" in feature.qualifiers and feature.qualifiers["locus_tag"][0] == "fragment") or (feature.type == "fragment"):
            # Use gene name, label, or locus_tag as text
            return ("SOI"
                # feature.qualifiers.get("gene", [""])[0]
                # or feature.qualifiers.get("label", [""])[0]
                # or feature.qualifiers.get("locus_tag", ["fragment"])[0]
            )
        else:
            return None  # no label for other features



class Visualization:
    def __init__(self, sequence, primers_df, sequence_of_interest_seq):
        self.sequence = str(sequence.upper())
        self.primers = primers_df
        self.primers = self.primers.reset_index(drop=True)
        self.sequence_of_interest = sequence_of_interest_seq.upper()
        self.processed_features = []

        self.primer_colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())

        self.draw_diagram()
        

    def _find_primer_positions(self):
        """Process primers and store start/end/strand info."""
        self.processed_features = []
        for i, row in self.primers.iterrows():
            name = row['name']
            fw_seq = row['primer_left_sequence'].upper()
            rv_seq = row['primer_right_sequence'].upper()

            # Forward primer
            fw_start = self.sequence.find(fw_seq)
            if fw_start == -1:
                print(f"Forward primer '{name}' not found.")
                continue
            fw_end = fw_start + len(fw_seq)

            # Reverse primer (search using reverse complement)
            rv_seq_rc = str(Seq(rv_seq).reverse_complement())
            rv_end = self.sequence.find(rv_seq_rc)
            if rv_end == -1:
                print(f" Reverse primer '{name}' not found.")
                continue
            rv_start = rv_end + len(rv_seq_rc) 

            color = self.primer_colors[i % len(self.primer_colors)]

            self.processed_features.append({
                'name': name + '.fw',
                'start': fw_start,
                'end': fw_end,
                'strand': '+',
                'color': color,
                'index': i
            })
            self.processed_features.append({
                'name': name + '.rv',
                'start': rv_start,
                'end': rv_end,
                'strand': '-',
                'color': color,
                'index': i
            })

    def _find_sequence_of_interest_position(self):
        """Find and store sequence of interest position"""
        start = self.sequence.find(self.sequence_of_interest)
        if start == -1:
            print("sequence_of_interest not found.")
            return None
        end = start + len(self.sequence_of_interest)
        return (start, end)

    def draw_diagram(self, output_path="primers_plot.png"):
        self._find_primer_positions()
        sequence_of_interest_pos = self._find_sequence_of_interest_position()
        sequence_of_interest_center = (sequence_of_interest_pos[0] + sequence_of_interest_pos[1]) // 2

        # Set up plot
        fig, ax = plt.subplots(figsize=(12, 4))
        seq_length = len(self.sequence)

        # Find sequence of interest center for coordinate shift
        if sequence_of_interest_pos:
            sequence_of_interest_center = (sequence_of_interest_pos[0] + sequence_of_interest_pos[1]) // 2
        else:
            sequence_of_interest_center = seq_length // 2  # fallback if sequence of interest not found

        # Genome baseline (centered around 0)
        ax.hlines(y=0, xmin=-sequence_of_interest_center, xmax=seq_length , color='black', linewidth=2)


        # Draw primers
        primer_positions = {}  # NEW: store arrow midpoints for each primer

        for feature in reversed(self.processed_features):
            y = 1 + (len(self.primers) - 1 - feature['index']) * 0.6
            dx = feature['end'] - feature['start']
            arrow_props = dict(
                width=0.2,
                head_width=0.3,
                head_length=3,
                length_includes_head=True,
                color=feature['color']
            )

            # Draw arrow
            ax.add_patch(FancyArrow(
                feature['start'] , y,
                dx if feature['strand'] == '+' else -abs(dx),
                0,
                **arrow_props
            ))

            # Compute midpoint (for connection line)
            midpoint_x = feature['start'] + dx / 2
            primer_name = feature['name'].replace('.fw', '').replace('.rv', '')  # base name

            # Save midpoint position
            if primer_name not in primer_positions:
                primer_positions[primer_name] = []
            primer_positions[primer_name].append((midpoint_x, y))

            # Label
            ax.text(midpoint_x, y + 0.25, feature['name'], ha='center', fontsize=8)

        # Draw connecting lines between forward and reverse arrows
        for name, coords in primer_positions.items():  # NEW
            if len(coords) == 2:  # we have both fw and rv
                (x1, y1), (x2, y2) = coords
                ax.plot([x1, x2], [y1, y2], color='black', alpha = 0.4, linewidth=1.2, zorder=0)


        # Draw
        if sequence_of_interest_pos:
            max_height = 1 + len(self.primers) * 0.6
            label_offset = 0.4
            sequence_of_interest_width = sequence_of_interest_pos[1] - sequence_of_interest_pos[0]
            sequence_of_interest_x = sequence_of_interest_pos[0]
            rect = patches.Rectangle(
                (sequence_of_interest_x, 0),        # (x, y)
                sequence_of_interest_width,         # width
                max_height,            # height
                linewidth=0,
                edgecolor=None,
                facecolor='red',
                alpha=0.5,
                label="Sequence of Interest"
            )
            ax.add_patch(rect)
            # Optional: add label centered inside the box
            ax.text(
                sequence_of_interest_x + sequence_of_interest_width / 2 ,
                max_height +label_offset,
                "Sequence of interest",
                ha='center',
                va='top',
                fontsize=9,
                color='red'
            )
        # Format
        if sequence_of_interest_pos:
            sequence_of_interest_center = (sequence_of_interest_pos[0] + sequence_of_interest_pos[1]) // 2

            # Fixed ±500 bp window
            window_size = 500
            x_start = max(0, sequence_of_interest_center - window_size)
            x_end = min(seq_length, sequence_of_interest_center + window_size)

            # Set x-limits centered around the sequence of interest
            ax.set_xlim(x_start, x_end)
            ax.set_xlabel("Position (bp) in the DNA sequence")
        else:
            # Fallback: still use ±500 if the sequence of interest was not found
            window_size = 500
            ax.set_xlim(-window_size, window_size)
            ax.set_xlabel("Position (bp) in the DNA sequence")

        ax.set_ylim(-1, 1 + len(self.primers) * 0.8)
        ax.set_yticks([])
        ax.set_title("Primer Binding Sites and Sequence of interest Region")
        plt.tight_layout()

        # Save and/or show
        plt.savefig(output_path, dpi=300)
        plt.show()

        return output_path