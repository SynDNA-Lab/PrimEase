import os
import sys
import shutil
import argparse
import logging
from datetime import datetime
import glob
import streamlit as st 
import tempfile 
import pandas as pd
import uuid 

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from config import Config
from target import Target
from selection import SequenceOfInterest,GenbankfileHandling
from primer3 import Primer3Interface
from bowtie import BowtieInterface
from offtarget import OfftargetChecker
from sequence_annotator import OverlappingFragmentAnnotator,SequenceAnnotator
from Visualization import Visualization, VisualizeGenbank
#---------------------------------------------------------------Side functions-----------------------------------------
def handle_file_upload(file, filename: str) -> str:
    """Save uploaded file to disk and return the saved path."""
    path = os.path.join(os.getcwd(), filename)
    with open(path, "wb") as f:
        f.write(file.read())
    return path

def check_file_structure(path) : 
    if "fasta" in path or "fa" in path : 
        with open(path, "r") as handle:
            if any(SeqIO.parse(handle, "fasta")) == False : 
                raise TypeError ("Something is wrong in the structure of the fasta file. Check it and restart the code")
    elif "genbank" in path or "gb" in path:
        with open(path, "r") as handle:
                if any(SeqIO.parse(handle, "genbank")) == False : 
                    raise TypeError ("Something is wrong in the structure of the genbank file. Check it and restart the tool")
    return
def get_sequence_id(region_str: str) -> str:
    for line in region_str.splitlines():
        if line.startswith("SEQUENCE_ID="):
            return line.split("=")[1].strip()
    return "Unknown"

def initial_cleanup () -> None : 
    extensions = ("*.png", "*.fasta", "*.fa", "*.gb", "*.genbank")

    # Remove files in the current directory
    for ext in extensions:
        for file in glob.glob(ext):
            try:
                os.remove(file)
                print(f"Deleted: {file}")
            except Exception as e:
                print(f"Error deleting {file}: {e}")

    print("Cleanup completed!")

def cleanup(config:Config, task_dir) -> None:
    now = datetime.now()
    folder_name = f"{config.jobname}_{now.strftime('%Y%m%d_%H%M')}"
    os.makedirs(folder_name, exist_ok=True)

    files_to_move = [
        "bt_target.csv",
        "bt_host.csv",
        "potential_primers.fasta",
        "primer3_result",
        "settings",
        "qTagGer_Output.csv",
        "primers_plot.png",
        "Genbank_vizualisation.png"
    ]

    for file in files_to_move:
        if os.path.exists(file):
            try:
                shutil.move(file, os.path.join(folder_name, file))
            except Exception as e:
                logging.warning(f"Could not move {file}: {e}")
        else:
            logging.info(f"File {file} not found, skipping cleanup.")

    # Remove temporary files if they exist
    for temp_file in ["annotated_sequence_expected_structure.gb", "target.fasta"]:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except Exception as e:
                logging.warning(f"Could not remove {temp_file}: {e}")
    
    if os.path.exists(task_dir):
        try:
            print("host path removed")
            shutil.rmtree(task_dir)
        except Exception as e:
            logging.warning(f"Could not remove {task_dir}: {e}")



def copy_reference_genome(id):
    # Create a unique host directory
    task_id = str(uuid.uuid4())
    task_dir = os.path.join("bowtie_index", task_id)
    os.makedirs(task_dir, exist_ok=True)
    host_dir = os.path.join(task_dir, "host")
    os.makedirs(host_dir, exist_ok=True)
    reference_dir = "bowtie_index/reference_genomes"

    if id == "Ecoli":
        keywords = ["ecoli"]

    elif id == "Scerevisiae":
        keywords = ["scerevisiae"]

    elif id == "E.coli + S.cerevisiae":
        keywords = ["ecoli", "scerevisiae"]

    else:
        keywords = []

    for filename in os.listdir(reference_dir):
        src = os.path.join(reference_dir, filename)

        if os.path.isfile(src):
            name_lower = filename.lower()
            if any(keyword in name_lower for keyword in keywords):
                shutil.copy(src, os.path.join(host_dir, filename))

    return task_dir

#-------------------------------------------------------------------Main functions---------------------------------------------------
def main_sequence_of_interest(uploaded_expected_structure) -> None:
    st.title("Primer design from :red[_a single sequence of interest_]")
    config = Config()

    expected_structure_path = None

    if type(uploaded_expected_structure) == str :
        expected_structure_path = uploaded_expected_structure
    else : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)
    
    check_file_structure(expected_structure_path)

    st.subheader("Modification of parameters (optional) : ")
    show_uploader = st.checkbox("Set a specific primer3 settings file ")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    if show_uploader:
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")


    # --- UI ---
    if 'genome_id' not in st.session_state:
        st.session_state['genome_id'] = "None"

    options = st.selectbox(
        "Add genome to prevent off-target amplification with a host organism",
        ("None","E.coli reference genome", "S.cerevisiae reference genome", "E.coli + S.cerevisiae reference genome")
    )

    # --- predefined genomes ---
    if options == "E.coli reference genome":
        st.session_state['host_id'] = "Ecoli"

    elif options == "S.cerevisiae reference genome":
        st.session_state['host_id'] = "Scerevisiae"


    # --- custom genome ---
    elif options == "E.coli + S.cerevisiae reference genome" : 
        st.session_state['host_id'] = "E.coli + S.cerevisiae"
    
    elif options == "None" : 
        st.session_state['host_id'] = "None"


    st.divider()


    st.subheader("Paste the sequence of interest")
    seq = st.text_input("Paste Sequence of interest")
    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product",
                                              options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])

    if 'sequence_submitted' not in st.session_state:
        st.session_state['sequence_submitted'] = False

    if st.button("Submit sequence"):
        if expected_structure_path:
            if seq != "" : 
                st.session_state['sequence_submitted'] = True
                st.success(f"Sequence of interest submitted")
            else : 
                st.error("Please paste a sequence of interest before submitting.")
                st.stop()
        else:
            st.error("Please upload a valid DNA file before submitting.")
            st.stop()

    if st.session_state['sequence_submitted']:
        st.session_state['sequence_submitted'] = False
        try:
            task_dir = copy_reference_genome(id = st.session_state['host_id'])

            annotator = SequenceAnnotator(expected_structure_path=expected_structure_path, search_seq=seq)
            annotator.run()
            logging.info("Annotation completed")

            target = Target(target_path="annotated_sequence_expected_structure.gb")
            logging.info("Target loaded")

            selection = SequenceOfInterest(target=target)
            logging.info("Sequence of interest selection done")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(
                regions=selection.regions[0],
                primer3_path=config.primer3_path,
                is_specific_file=specific_file,
                size_pcr_product = min_max
            )
            p3interface.run()
            logging.info("Primer3 finished")

            logging.info("Running Bowtie...")
            bowtie = BowtieInterface(config=config,task_dir=task_dir)
            logging.info("Bowtie alignment done")

            logging.info("Checking off-targets...")
            check_offtargets = OfftargetChecker(
                primer_candidates=p3interface.primer_sites,
                bowtie_target=bowtie.result_target,
                bowtie_host=bowtie.result_host,
                config=config
            )

            st.dataframe(check_offtargets.final_primers[[
                                "name",
                                "primer_pair_product_size",
                                "primer_pair_penalty",
                                "primer_left_sequence",
                                "primer_right_sequence",
                                "primer_left_tm",
                                "primer_right_tm",
                                "primer_left_gc_percent",
                                "primer_right_gc_percent"
                                                        ]])

            visualization = Visualization(
                sequence=str(selection.target.record.seq).upper(),
                primers_df=check_offtargets.final_primers,
                sequence_of_interest_seq=annotator.search_seq
            )
            st.image("primers_plot.png")

            if task_dir is not None : 
                cleanup(config, task_dir)
                for path in [save_path_settings, expected_structure_path]:
                    if path and os.path.exists(path):
                        os.remove(path)

        except Exception as e:
            st.error(f"Failed to process DNA fragments: {e}")



def main_fragments(uploaded_expected_structure) -> None:
    st.header("Primer design from :red[_two overlapping sequences_]")
    config = Config()

    expected_structure_path = None
    if type(uploaded_expected_structure) == str :
        expected_structure_path = uploaded_expected_structure
    else : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)
    
    check_file_structure(expected_structure_path)

    st.subheader("Modification of parameters (optional) : ")
    show_uploader = st.checkbox("Set a specific primer3 settings file")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    if show_uploader:
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")

    if 'genome_id' not in st.session_state:
        st.session_state['genome_id'] = "None"

    options = st.selectbox(
        "Add genome to prevent off-target amplification with a host organism",
        ("None","E.coli reference genome", "S.cerevisiae reference genome", "E.coli + S.cerevisiae reference genome")
    )

    # --- predefined genomes ---
    if options == "E.coli reference genome":
        st.session_state['host_id'] = "Ecoli"

    elif options == "S.cerevisiae reference genome":
        st.session_state['host_id'] = "Scerevisiae"


    # --- custom genome ---
    elif options == "E.coli + S.cerevisiae reference genome" : 
        st.session_state['host_id'] = "E.coli + S.cerevisiae"
    
    elif options == "None" : 
        st.session_state['host_id'] = "None"


    st.divider()


    st.subheader("Paste the two overlapping sequences")
    seq1 = st.text_input("Enter sequence 1")
    seq2 = st.text_input("Enter sequence 2")
    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product", options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])

    if 'sequence_submitted_overlapping' not in st.session_state:
        st.session_state['sequence_submitted_overlapping'] = False

    if st.button('Submit sequence'):
        if expected_structure_path:  
                if seq1 != "" and seq2 != "" :
                    st.session_state['sequence_submitted_overlapping'] = True
                    st.success(f"Both overlapping sequences submitted")
                else : 
                    st.error("Please paste a sequence of interest before submitting.")
                    st.stop()
        else:
            st.error("Please enter a valid DNA file before submitting.")
            st.stop() 

    if st.session_state['sequence_submitted_overlapping']:
        st.session_state['sequence_submitted_overlapping'] = False
        try:
            task_dir = copy_reference_genome(id = st.session_state['host_id'])
            logging.info("Start expected DNA assembly annotation")
            annotator = OverlappingFragmentAnnotator(expected_structure_path=expected_structure_path,seq_frag1=seq1,seq_frag2=seq2)
            annotator.run()
            logging.info("Successfully performed expected DNA assembly annotation")

            target = Target(target_path="annotated_sequence_expected_structure.gb")
            logging.info("Successfully loaded target record")

            selection = SequenceOfInterest(target=target)
            logging.info("Successfully performed selection")

            logging.info("Running Primer3...")
            p3interface = Primer3Interface(
                regions=selection.regions[0],
                primer3_path=config.primer3_path,
                is_specific_file=specific_file,
                size_pcr_product = min_max
            )
            p3interface.run()
            logging.info("Primer3 finished")

            logging.info("Running Bowtie")
            bowtie = BowtieInterface(config=config,task_dir=task_dir)
            logging.info("Successfully performed bowtie alignments")

            logging.info("Offtarget checking")
            check_offtargets = OfftargetChecker(
                primer_candidates = p3interface.primer_sites,
                bowtie_target = bowtie.result_target,
                bowtie_host=bowtie.result_host,
                config=config
            )

            st.dataframe(check_offtargets.final_primers[[
                                "name",
                                "primer_pair_product_size",
                                "primer_pair_penalty",
                                "primer_left_sequence",
                                "primer_right_sequence",
                                "primer_left_tm",
                                "primer_right_tm",
                                "primer_left_gc_percent",
                                "primer_right_gc_percent"
                                                        ]])
            visualization = Visualization(
                sequence = str(selection.target.record.seq).upper(), 
                primers_df=check_offtargets.final_primers,
               sequence_of_interest_seq=annotator.sequence_of_interest_seq)
            st.image('primers_plot.png')

            try:  
                cleanup(config,task_dir)
                for path in [save_path_settings, expected_structure_path]:
                    if path and os.path.exists(path):
                        os.remove(path)
            except Exception as e : 
                ""
        except Exception as e:
            st.error(f"Failed to process DNA fragments : {e}")


def main_annotated_genbank(uploaded_expected_structure) -> None:
    st.header("Primer design from :red[_annotations in a genbank file_]")
    config = Config()

    expected_structure_path = None
    if uploaded_expected_structure : 
        expected_structure_path = handle_file_upload(uploaded_expected_structure, uploaded_expected_structure.name)
        print(expected_structure_path)
    check_file_structure(expected_structure_path)
    st.subheader("Modification of parameters (optional)?")

    show_uploader = st.checkbox("Set a specific primer3 settings file ?")

    uploaded_settings = None
    save_path_settings = None
    specific_file = False

    if show_uploader :
        uploaded_settings = st.file_uploader("Drop your file here", type=["bak"])
        st.write("Primer3 webpage for designing: https://primer3.org/manual.html")

        if uploaded_settings:
            save_path_settings = handle_file_upload(uploaded_settings, "settings_spec.bak")
            st.success("Primer3 settings file uploaded")
            specific_file = True
        else:
            st.warning("No settings file uploaded yet.")

    st.divider()
    
    if 'genome_id' not in st.session_state:
        st.session_state['genome_id'] = "None"

    options = st.selectbox(
        "Add genome to prevent off-target amplification with a host organism",
        ("None","E.coli reference genome", "S.cerevisiae reference genome", "E.coli + S.cerevisiae reference genome")
    )

    # --- predefined genomes ---
    if options == "E.coli reference genome":
        st.session_state['host_id'] = "Ecoli"

    elif options == "S.cerevisiae reference genome":
        st.session_state['host_id'] = "Scerevisiae"


    # --- custom genome ---
    elif options == "E.coli + S.cerevisiae reference genome" : 
        st.session_state['host_id'] = "E.coli + S.cerevisiae"
    
    elif options == "None" : 
        st.session_state['host_id'] = "None"


    st.divider()

    st.write("Enter minimum and maximum size for the PCR product")
    min_max = st.select_slider("Select a range of size for the PCR product", options=list(range(100, 801)), value=(300, 500),)
    st.write("You selected size between", min_max[0], "and", min_max[1])

    target = Target(target_path=expected_structure_path)
    logging.info("Successfully loaded target record")

    selection = GenbankfileHandling(genbank_path=expected_structure_path, target=target)
    logging.info("Successfully performed selection")
    sequence_labels = [region.name for region in selection.regions]

    selected_indices = st.multiselect(
        "Select one or more sequences of interest to design primers for:",
        options=list(range(len(sequence_labels))),
        format_func=lambda x: sequence_labels[x],
    )
    vizgb = VisualizeGenbank(path=expected_structure_path)
    vizgb.draw()
    st.image('Genbank_vizualisation.png', caption="Position of Sequence of Interest in the complete DNA sequence ")
    if selected_indices:
        st.markdown("### Preview of selected sequences")
        for idx in selected_indices:
            with st.expander(f"{sequence_labels[idx]}"):
                st.code(selection.list_sequence_of_interest[idx], language='text')
                
    if st.button("Run Primer Design on Selected Sequences"):
        if not selected_indices:
            st.warning("Please select at least one sequence.")
        else:
            progress_bar = st.progress(0)
            total = len(selected_indices)
            download_df = pd.DataFrame()
            task_dir = copy_reference_genome(id = st.session_state['host_id'])
            for count, i in enumerate(selected_indices, start=1):
                logging.info(f"Running Primer3")

                # Run Primer3

                p3interface = Primer3Interface(
                    regions=selection.regions[i],
                    primer3_path=config.primer3_path,
                    is_specific_file=specific_file,
                    size_pcr_product=min_max
                )
                p3interface.run()

                # Run Bowtie
                logging.info(f"Running Bowtie")

                bowtie = BowtieInterface(config=config,task_dir=task_dir)
                
                # Off-target checking
                logging.info(f"Off-target checking")
                check_offtargets = OfftargetChecker(
                    primer_candidates=p3interface.primer_sites,
                    bowtie_target=bowtie.result_target,
                    bowtie_host=bowtie.result_host,
                    config=config
                )
                check_offtargets.final_primers = (check_offtargets.final_primers.reset_index(drop=True))

                name = check_offtargets.final_primers.loc[0, 'name']
                name = "_".join(name.split("_")[:-1])

                st.markdown("---")
                st.subheader(f"Results for sequence of interest :blue[_{name}_]")

                download_df = pd.concat([download_df, check_offtargets.final_primers[[
                    "name",
                    "primer_pair_product_size",
                    "primer_pair_penalty",
                    "primer_left_sequence",
                    "primer_right_sequence",
                    "primer_left_tm",
                    "primer_right_tm",
                    "primer_left_gc_percent",
                    "primer_right_gc_percent"
                ]]])

                st.dataframe(check_offtargets.final_primers[[
                    "name",
                    "primer_pair_product_size",
                    "primer_pair_penalty",
                    "primer_left_sequence",
                    "primer_right_sequence",
                    "primer_left_tm",
                    "primer_right_tm",
                    "primer_left_gc_percent",
                    "primer_right_gc_percent"
                ]])

                # Visualization
                visualization = Visualization(
                    sequence=str(selection.target.record.seq).upper(),
                    primers_df=check_offtargets.final_primers,
                    sequence_of_interest_seq=selection.list_sequence_of_interest[i]
                )

                st.image('primers_plot.png', caption=f"Primers Visualization for {name}")

                # Update progress
                progress_bar.progress(count / total)
                st.markdown("---")

            csv_data = download_df.to_csv(index=False).encode("utf-8")
            st.download_button (label="Download all the primers in a csv file (This unloads the results)",
                        data=csv_data,
                        file_name="primease.csv",
                        icon=":material/download:")
            
            st.success("Primer design completed for all selected sequences!")
    try:  
        cleanup(config,task_dir)
        for path in [save_path_settings, expected_structure_path]:
            if path and os.path.exists(path):
                os.remove(path)
    except Exception as e:
        ""

#-------------------------------------------------------------UTILITY FUNCTION------------------------------------------------------
import re
import io

def fix_genbank_locus(uploaded_file):
    original_name = uploaded_file.name  # preserve filename

    content = uploaded_file.read()

    if isinstance(content, bytes):
        content = content.decode("utf-8")

    lines = content.splitlines()

    if lines and lines[0].startswith("LOCUS"):
        parts = lines[0].split()

        if len(parts) >= 7:
            _, name, length, bp, mol_type, topology, date = parts[:7]

            fixed_locus = (
                f"LOCUS       "
                f"{name:<16}"
                f"{int(length):>11} bp    "
                f"{mol_type:<6}"
                f"{topology:<8}"
                f"SYN "
                f"{date}"
            )

            lines[0] = fixed_locus

    fixed_content = "\n".join(lines) + "\n"

    fixed_file = io.BytesIO(fixed_content.encode("utf-8"))

    # Add name attribute so it behaves like UploadedFile
    fixed_file.name = original_name

    return fixed_file





#---------------------------------------------------------------------Script----------------------------------------------------------
logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO)

initial_cleanup()
st.title("Primease :dna: \n (Primer Design for Sequence of Interest )")

uploaded_expected_structure = st.file_uploader(
    "Select the DNA template (expected assembly, genome...)", 
    type=['fa', 'fasta', 'genbank', 'gb']
)
if uploaded_expected_structure : 
    st.success("Expected DNA assembly successfully uploaded")
file_type = None
sequence_type = None

if 'Continue' not in st.session_state:
    st.session_state['Continue'] = False
if 'Option' not in st.session_state:
    st.session_state['Option'] = None

pasted_expected_structure_construct = st.text_area("Or paste your DNA sequence here :")

if pasted_expected_structure_construct:
    sequence_type = st.selectbox(
        "Format of the sequence of interest :",
        ("Single sequence", " Two overlapping sequences")
    )
    if st.button("Continue"):
        st.session_state['Continue'] = True
        st.session_state['Option'] = sequence_type

if st.session_state.get('Continue') and st.session_state.get('Option') and pasted_expected_structure_construct:
    pasted_expected_structure_construct = pasted_expected_structure_construct.replace("\n", "").replace("\r", "").replace(" ", "")
    record = SeqRecord(Seq(pasted_expected_structure_construct), id="pasted_sequence", description="")
    
    with open("pasted_sequence.fasta", "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

    if st.session_state['Option'] == "Single sequence":
        main_sequence_of_interest("pasted_sequence.fasta")
    elif st.session_state['Option'] == "Two overlapping sequences":
        main_fragments("pasted_sequence.fasta")

elif uploaded_expected_structure:
    file_name = uploaded_expected_structure.name.lower()
    if any(ext in file_name for ext in ['fa', 'fasta']):
        file_type = "fasta"
        sequence_type = st.selectbox(
            "Format of the sequence of interest:", 
            ("Single sequence", "Two overlapping sequences")
        )

    elif any(ext in file_name for ext in ['gb', 'genbank']):
        file_type = "genbank"
        sequence_type = st.selectbox(
            "Format of the Sequence of interest:", 
            ("Single sequence", "Two overlapping sequences", "Annotation in a genbank file")
        )

    if st.button("Continue") and sequence_type:
        st.session_state['Continue'] = True
        st.session_state['Option'] = sequence_type

if st.session_state.get('Continue') and st.session_state.get('Option'):
    if file_type == "fasta":
        
        if st.session_state['Option'] == "Single sequence":
            main_sequence_of_interest(uploaded_expected_structure)
        elif st.session_state['Option'] == "Two overlapping sequences":
            main_fragments(uploaded_expected_structure)

    elif file_type == "genbank":

        uploaded_expected_structure.seek(0)
        uploaded_expected_structure = fix_genbank_locus(uploaded_expected_structure)


        if st.session_state['Option'] == "Annotation in a genbank file":
            main_annotated_genbank(uploaded_expected_structure)

        elif st.session_state['Option'] == "Single sequence":
            main_sequence_of_interest(uploaded_expected_structure)

        elif st.session_state['Option'] == "Two overlapping sequences":
            main_fragments(uploaded_expected_structure)


