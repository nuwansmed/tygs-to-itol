import os
import csv
from lxml import etree
import pandas as pd

# Define the XML file path and the CSV file path
xml_file = '../data/tree.xml'  # Replace with the path to your XML file
csv_file = 'extracted_data.csv'
min_max_csv_file = 'min_max_values.csv'
normalized_csv_file = 'normalized_data.csv'

# Define the namespace for phyloxml
NS = {'phyloxml': 'http://www.phyloxml.org'}
NS_TREE = {'phy': 'http://www.phyloxml.org'}

def parse_graph(xml_element):
    """Extract graph details from an XML element."""
    graph = {}
    graph['name'] = xml_element.findtext('phyloxml:name', '', namespaces=NS).strip()
    
    legend = xml_element.find('phyloxml:legend', namespaces=NS)
    if legend is not None:
        fields = legend.findall('phyloxml:field', namespaces=NS)
        graph['fields'] = [field.findtext('phyloxml:name', '', namespaces=NS).strip() for field in fields]
    else:
        graph['fields'] = []

    data_values = []
    for values in xml_element.findall('.//phyloxml:data/phyloxml:values', namespaces=NS):
        id_value = values.attrib.get('for', '')
        value_elements = values.findall('phyloxml:value', namespaces=NS)
        data_values.append((id_value, [v.text.strip() for v in value_elements]))
    
    graph['data'] = data_values
    
    return graph

def identify_field_type(fields):
    """Identify the field type based on the field names."""
    if 'Percent G+C' in fields:
        return 'Percent G+C'
    elif 'delta statistics' in fields:
        return 'delta statistics'
    elif 'Genome size (in bp)' in fields and 'Protein count' in fields:
        return 'Multibar'
    elif 'User strain?' in fields:
        return 'User strain?'
    elif 'Type species?' in fields:
        return 'Type species?'
    return None

def extract_data_from_graph(graph):
    """Extract data based on the graph type."""
    field_type = identify_field_type(graph['fields'])
    extracted_data = {}
    
    for data_id, values in graph['data']:
        if field_type == 'User strain?':
            extracted_data[data_id] = {'User strain?': values[0] if len(values) > 0 else ''}
        elif field_type == 'Type species?':
            extracted_data[data_id] = {'Type species?': values[0] if len(values) > 0 else ''}
        elif field_type == 'Percent G+C':
            extracted_data[data_id] = {'Percent G+C': values[0] if len(values) > 0 else ''}
        elif field_type == 'delta statistics':
            extracted_data[data_id] = {'delta statistics': values[0] if len(values) > 0 else ''}
        elif field_type == 'Multibar':
            if len(values) >= 2:
                extracted_data[data_id] = {
                    'Genome size (in bp)': values[0],
                    'Protein count': values[1]
                }
    
    return extracted_data

def extract_clade_data(file_path):
    """Extract clade data from a separate XML file."""
    clade_data = {}
    try:
        tree = etree.parse(file_path)
        root = tree.getroot()
        
        for clade in root.findall(".//phy:clade", namespaces=NS_TREE):
            name = clade.find('phy:name', namespaces=NS_TREE)
            clade_id = clade.find('phy:id', namespaces=NS_TREE)
            subspecies_cluster_color = clade.find('phy:subspeciesclustercolor', namespaces=NS_TREE)
            subspecies_cluster_id = clade.find('phy:subspeciesclusterid', namespaces=NS_TREE)
            species_cluster_color = clade.find('phy:speciesclustercolor', namespaces=NS_TREE)
            species_cluster_id = clade.find('phy:speciesclusterid', namespaces=NS_TREE)
            branch_length = clade.find('phy:branch_length', namespaces=NS_TREE)
            
            if (name is not None and clade_id is not None and 
                subspecies_cluster_color is not None and subspecies_cluster_id is not None and 
                species_cluster_color is not None and species_cluster_id is not None and 
                branch_length is not None):
                clade_data[clade_id.text] = {
                    'Name': name.text.strip(),
                    'Subspecies Cluster Color': subspecies_cluster_color.text,
                    'Subspecies Cluster ID': subspecies_cluster_id.text,
                    'Species Cluster Color': species_cluster_color.text,
                    'Species Cluster ID': species_cluster_id.text,
                    'Branch Length': branch_length.text
                }
    
    except etree.XMLSyntaxError as e:
        print(f"XML Syntax Error in tree file: {e}")
    except Exception as e:
        print(f"An error occurred while extracting clade data: {e}")
    
    return clade_data

def write_csv(graph_data, clade_data, file_name):
    """Write extracted data to a CSV file."""
    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file)
        
        # Write header
        writer.writerow(['ID', 'Name', 'User strain?', 'Type species?', 'Percent G+C', 'delta statistics', 'Genome size (in bp)', 'Protein count', 'Subspecies Cluster Color', 'Subspecies Cluster ID', 'Species Cluster Color', 'Species Cluster ID', 'Branch Length'])
        
        # Combine and write data rows
        for data_id in set(graph_data.keys()).union(clade_data.keys()):
            row = [data_id]
            graph_values = graph_data.get(data_id, {})
            clade_values = clade_data.get(data_id, {})
            row.extend([
                clade_values.get('Name', ''),
                graph_values.get('User strain?', ''),
                graph_values.get('Type species?', ''),
                graph_values.get('Percent G+C', ''),
                graph_values.get('delta statistics', ''),
                graph_values.get('Genome size (in bp)', ''),
                graph_values.get('Protein count', ''),
                clade_values.get('Subspecies Cluster Color', ''),
                clade_values.get('Subspecies Cluster ID', ''),
                clade_values.get('Species Cluster Color', ''),
                clade_values.get('Species Cluster ID', ''),
                clade_values.get('Branch Length', '')
            ])
            writer.writerow(row)
    
    print(f"CSV file '{file_name}' written.")

def calculate_min_max_values(csv_file, output_file):
    """Calculate and save the min and max values of selected fields."""
    df = pd.read_csv(csv_file)
    min_max_df = pd.DataFrame({
        'Field': ['Percent G+C', 'delta statistics', 'Genome size (in bp)', 'Protein count'],
        'Min': [df['Percent G+C'].min(), df['delta statistics'].min(), df['Genome size (in bp)'].min(), df['Protein count'].min()],
        'Max': [df['Percent G+C'].max(), df['delta statistics'].max(), df['Genome size (in bp)'].max(), df['Protein count'].max()]
    })
    min_max_df.to_csv(output_file, index=False)
    print(f"Min-Max values saved to '{output_file}'.")

def normalize_data(csv_file, output_file, min_max_file):
    """Normalize the selected fields using min-max normalization."""
    df = pd.read_csv(csv_file)
    min_max_df = pd.read_csv(min_max_file)

    for index, row in min_max_df.iterrows():
        field = row['Field']
        min_val = row['Min']
        max_val = row['Max']
        df[field + '_normalized'] = (df[field] - min_val) / (max_val - min_val)
    
    df.to_csv(output_file, index=False)
    print(f"Normalized data saved to '{output_file}'.")

def process_xml(file_path):
    """Process XML file and extract data for different field types."""
    try:
        tree = etree.parse(file_path)
        root = tree.getroot()

        extracted_data = {}
        
        for graph_element in root.findall('.//phyloxml:graph', namespaces=NS):
            graph = parse_graph(graph_element)
            graph_data = extract_data_from_graph(graph)
            
            # Merge the data into the main extracted_data dictionary
            for data_id, values in graph_data.items():
                if data_id not in extracted_data:
                    extracted_data[data_id] = {}
                extracted_data[data_id].update(values)
        
        clade_data = extract_clade_data(xml_file)
        
        write_csv(extracted_data, clade_data, csv_file)
        calculate_min_max_values(csv_file, min_max_csv_file)
        normalize_data(csv_file, normalized_csv_file, min_max_csv_file)

    except etree.XMLSyntaxError as e:
        print(f"XML Syntax Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

def append_to_annotation_files(csv_file, templates_dir):
    """Append data from the CSV file to each template file."""
    # Load the normalized data from the CSV file
    df = pd.read_csv(csv_file)

    # Define the mappings for template files and corresponding columns
    template_columns = {
        'annotation_delta_statistics.txt': 'delta statistics_normalized',
        'annotation_species.txt': 'Type species?',
        'annotation_genome_size.txt': 'Genome size (in bp)_normalized',
        'annotation_subspecies.txt': 'Subspecies Cluster Color',
        'annotation_percent_gc.txt': 'Percent G+C_normalized',
        'annotation_type_strain.txt': 'Type species?',
        'annotation_protein_count.txt': 'Protein count_normalized',
        'annotation_user_strain.txt': 'User strain?'
    }

    # Process each template file
    for template_file, column in template_columns.items():
        template_path = os.path.join(templates_dir, template_file)

        # Read the data for the specific column
        data = df[['Name', column]].dropna()

        # Open the template file and append data
        with open(template_path, 'a') as f:
            for _, row in data.iterrows():
                f.write(f"{row['Name']},{row[column]}\n")
        
        print(f"Data appended to '{template_file}'.")


# Example usage
process_xml(xml_file)

append_to_annotation_files(normalized_csv_file, 'templates')

