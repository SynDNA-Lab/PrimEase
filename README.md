# PrimEase
This project is a variation of the qTagGer tool originally developed by RGSchindler.
It enables the automatic design of primers for regions of interest, making it easy to integrate into cloning or primer design pipelines. 
The tool simplifies the process of verifying successful recombination of different DNA fragments through PCR 
by allowing users to quickly generate the necessary primers.

## Installation 

### From Github repository
```bash
# Clone the repository
git clone https://github.com/SynDNA-Lab/PrimEase.git
cd PrimEase

#installation of the virtual environment and activation
conda env create -f environment.yml --name primease
conda activate primease

#run the main script
streamlit run tool/main_streamlit.py
```
## Use
This script opens the tool in the browser window.

It offers two main options:
  - Option 1: Drop (or paste) a FASTA sequence, then paste the sequence of interest or two overlapping sequences. The overlapping region will be treated as the sequence of interest.
  - Option 2: Drop an already annotated GenBank file with the desired sequences annotated either as type = "fragment" or as type = "misc_feature" with locus_tag = "fragment".

It is also possible to modify the settings of the Primer3 tool by dropping a modified settings file.

## In depth tutorial

After you finish the installation of the tool, open the http adress in your browser. 

<img width="750" height="403" alt="image" src="https://github.com/user-attachments/assets/af89e42f-c355-488a-9570-ae1f5dcc64cc" />

#### Uploading your DNA background (template, exepected assembly, genome..)
Here you have two options : 
1) Selecting a DNA file from your local storage and drop it in the interface
2) Pasting a DNA sequence in the tool (Only the sequence, no title). Press ctrl + enter

Depending on the type of file given during the previous step you can select different options. 

#### Single Sequence 
After pressing continue with the **single sequence** option, you will get the following interface :

<img width="626" height="553" alt="image" src="https://github.com/user-attachments/assets/a0c37e7b-6582-4de9-8d88-1aad37f675c0" />

The optional modification of the default parameters are expalined later. Let's focus on the main tool for now. 
In the box *Paste *the **sequence of interest** you can paste a DNA sequence. The tool is gonna search for this sequence in the background DNA given above and design primers on both side of the sequence. the sequence needs to be a part of the background.

Then, you can press select and you get the following table and figure : 

<img width="788" height="522" alt="image" src="https://github.com/user-attachments/assets/2090d906-da74-45cc-991d-392913005766" />

In the table you can find the most important information about the 5 best paires of primers found by the tool (penalty, PCR product size, sequence, melting temperature). This table can be downloaded as a csv file in the upper right corner.
The figure underneath shows the position of the sequence of interest and the primer in the background DNA.

#### Two overlapping sequence 
After pressing continue with the **two overlapping sequence** option, you will get the following interface :

<img width="594" height="587" alt="image" src="https://github.com/user-attachments/assets/4bc14f07-c006-41b7-a004-a096e06ed8ea" />

The optional modification of the default parameters are expalined later. Let's focus on the main tool for now. 
In the box **Paste the two overlapping sequences** you can paste two DNA sequences. The tool is gonna search for an overlap (at least 10 bp) between the two fragments and then search the overlapping sequence in the background given previously. The primers are then calculated around this overlap. Both sequences need to be part of the background.

In the table you can find the most important information about the 5 best paires of primers found by the tool (penalty, PCR product size, sequence, melting temperature). This table can be downloaded as a csv file in the upper right corner.
The figure underneath shows the position of the sequence of interest and the primer in the background DNA.

#### Annotation in genbankfile 
Only with a **genbank file**, you have access to a third option in the selection box : Annotation in a genbank file.

**Conditions !!**
For this options to work, your sequence of interest or overlapping fragment must be annotated with the following rules 
1) either "type" = "fragment"
2) or "type" = "misc_feature" AND the qualifier "locus_tag" = "fragment"

<img width="674" height="568" alt="image" src="https://github.com/user-attachments/assets/93a0b1d0-bf8a-46b3-b445-6351349da0fd" />

and the following interface 

<img width="508" height="516" alt="image" src="https://github.com/user-attachments/assets/d3a7ccc0-f43f-4d6d-919a-78c97c6efa14" />

With the **Choose options** selectionbox, you can see all the annotation that are recognized by the tool (see Conditionss above). You can also check the position in the background by looking at the figure. The other annotations are also shown with a specific color code to navigate more easily in the background.

After selecting one or several targets, you can preview their sequence to check. 

<img width="591" height="536" alt="image" src="https://github.com/user-attachments/assets/0e7142f6-c6d1-4384-af88-bf4c089eb2cb" />

After clicking **Run Primer Design on Selected Sequences**, you get a table and a figure for each target. 
In the tables you can find the most important information about the 5 best paires of primers found by the tool (penalty, PCR product size, sequence, melting temperature). This table can be downloaded as a csv file in the upper right corner.
The figures underneath shows the position of the sequence of interest and the primer in the background DNA.
You can download all the primers at the same time using the button **Download all the primers in a csv file** 
<img width="252" height="44" alt="image" src="https://github.com/user-attachments/assets/96d480cc-1eb7-4a8b-a8c7-c37663c782bc" />

#### Optional features 

1) Set specific parameters for primer3
When you check the box, you have access to a drop where you can select a settings file that should be used by primer3 instead of the default parameters. An example of the expected structure is given in the *Example *file folder of the github.

<img width="591" height="156" alt="image" src="https://github.com/user-attachments/assets/aea19a6c-3b87-4a88-9a9b-2bec57f5c2a7" />

2) Set a specific PCR product size range
You can select a range of size for the PCR produt that is gonna overwrite the value given in the parameters. It allows fast modification.

<img width="589" height="101" alt="image" src="https://github.com/user-attachments/assets/57e2f516-46fb-4a33-9b30-e412d5829d76" />

3) Test alignment of the primers on an host genome
This option enables you to quickly prevent off-target amplification in the host organism. On the selecting box you can either select
1) No genome (default)
2) E.coli genome
3) The genome of your choice (Drop the file in the box)

<img width="584" height="160" alt="image" src="https://github.com/user-attachments/assets/4f142609-188f-4da7-bed4-e52c94a4b791" />

When this option is modified, the next running might be a bit longer due to the calculation of the index.








