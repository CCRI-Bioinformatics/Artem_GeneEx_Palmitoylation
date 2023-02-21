import os
import pandas as pd
import xml.etree.ElementTree as ET

# Set the directory where the XML files are located
directory = os.path.expanduser("~/bioinf_isilon/core_bioinformatics_unit/Internal/max_vdl/coreBioinf/kalinchenko/tissueGeneEx/data")

# Get a list of the XML files in the directory
xml_files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".xml")]

# Initialize an empty dictionary to hold the RNAexpression data
rna_expression_data = {}

# Loop through each XML file
for xml_file in xml_files:

    # Parse the XML file
    tree = ET.parse(xml_file)

    # Get the root element of the XML file
    root = tree.getroot()

    # Find the tag wih the desired expression level entries
    rna_expression = root.find("./entry/rnaExpression[@assayType='consensusTissue']")
    
    # Loop through each data child tag of the rnaExpression tag
    for data in rna_expression.iter("data"):
        # Get the name of the tissue type
        tissueName = data.find('tissue').text
        # Get the nTPM value from tissue type
        expRNA_value = data.find("./level[@type='normalizedRNAExpression']").get('expRNA')
        expRNA_value = float(expRNA_value)
        # Add the tag name and value to the dictionary
        if tissueName not in rna_expression_data:
            rna_expression_data[tissueName] = expRNA_value
        else: raise Exception('a tissue type was found twice')

    # Convert the dictionary to a Pandas dataframe
    protName = root.find('entry').find('name').text
    try:
        rna_expression_df[protName] = rna_expression_data
    except NameError:
        rna_expression_df = pd.DataFrame.from_dict(rna_expression_data, orient='index', columns=[protName])

    # Empty the dict before the new file
    rna_expression_data = {}

# Sort the cols alphabetically
rna_expression_df = rna_expression_df.sort_index(axis=1)
# Save the resulting dataframe
rna_expression_df.to_csv(f'{directory}/tissueProtTable.csv')