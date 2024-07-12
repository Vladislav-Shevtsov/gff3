from flask import Flask, request, send_file, jsonify
from flask_cors import CORS
import pandas as pd
import os
import pickle


app = Flask(__name__, static_folder='static')
CORS(app)

UPLOAD_FOLDER = 'uploads'
OUTPUT_FOLDER = 'outputs'

# Ensure the upload and output directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# Load the gene_product_map from the pickle file
with open('gene_product_map.pkl', 'rb') as f:
    gene_product_map = pickle.load(f)

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file:
        file_path = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(file_path)

        # Process the uploaded file
        formatted_output_path = process_gff3(file_path)

        return jsonify({'formatted_output_path': formatted_output_path})

@app.route('/download/<filename>', methods=['GET'])
def download_file(filename):
    file_path = os.path.join(OUTPUT_FOLDER, filename)
    return send_file(file_path, as_attachment=True)

def process_gff3(file_path):
    # Define the column names for the GFF3 file
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

    # Read the GFF3 file into a DataFrame
    gff3_df = pd.read_csv(file_path, sep='\t', comment='#', names=columns)

    # Drop unnecessary columns
    gff3_df = gff3_df.drop(columns=['source', 'seqid', 'phase', 'score'])

    # Filter out rows with types 'region' and 'source'
    gff3_df = gff3_df[~gff3_df['type'].isin(['region', 'source'])]

    # Skip introns
    gff3_df = gff3_df[gff3_df['type'] != 'intron']

    # Retrieve the name of the gene from the 'attributes' column
    gff3_df['name'] = gff3_df['attributes'].str.extract(r'Name=([^;]+)')

    # Drop the original 'attributes' column
    gff3_df = gff3_df.drop(columns=['attributes'])

    # Function to swap start and end positions if strand is '-'
    def adjust_positions(row):
        if row['strand'] == '-':
            row['start'], row['end'] = row['end'], row['start']
        return row

    # Apply position adjustment function
    gff3_df = gff3_df.apply(adjust_positions, axis=1)

    # Format gene and CDS entries
    formatted_list = format_gene_cds_list(gff3_df)

    # Format the output
    formatted_output = format_output(formatted_list)

    # Define output file path
    original_filename = os.path.splitext(os.path.basename(file_path))[0]
    output_filename = f"{original_filename}_formatted.txt"
    output_file_path = os.path.join(OUTPUT_FOLDER, output_filename)
    feature_line = f">feature\t{output_filename.split('_')[0]}\n"

    # Save formatted output to file
    save_to_file(formatted_output, output_file_path,feature_line )

    return output_filename

def format_gene_cds_list(df):
    formatted_list = []
    current_gene = None
    current_gene_data = {}

    for index, row in df.iterrows():
        if row['type'] == 'gene':
            # If it's a gene entry, store the current gene data if exists
            if current_gene is not None:
                formatted_list.append(current_gene_data)
            # Initialize new gene data
            current_gene = row['name']
            current_gene_data = {
                'gene': current_gene,
                'type': 'gene',
                'start': row['start'],
                'end': row['end'],
                'CDS': [],
                'tRNA': [],
                'rRNA':[],
                'repeat_region':[]
            }
        elif row['type'] == 'CDS':
            # If it's a CDS entry, check if it's contiguous with previous CDS
            if current_gene is not None:
                if (current_gene_data['CDS'] and 
                    current_gene_data['CDS'][-1]['end'] + 1 == row['start']):
                    # Merge with previous CDS if contiguous
                    current_gene_data['CDS'][-1]['end'] = row['end']
                else:
                    # Otherwise, add new CDS entry
                    current_gene_data['CDS'].append({
                        'type': 'CDS',
                        'start': row['start'],
                        'end': row['end']
                    })
        elif row['type'] == 'tRNA':
            # If it's a tRNA entry, check if it's contiguous with previous tRNA
            if current_gene is not None:
                if (current_gene_data['tRNA'] and 
                    current_gene_data['tRNA'][-1]['end'] + 1 == row['start']):
                    # Merge with previous tRNA if contiguous
                    current_gene_data['tRNA'][-1]['end'] = row['end']
                else:
                    # Otherwise, add new tRNA entry
                    current_gene_data['tRNA'].append({
                        'type': 'tRNA',
                        'start': row['start'],
                        'end': row['end'],
                        
                        
                    })
        elif row['type'] == 'rRNA':
            # If it's a tRNA entry, check if it's contiguous with previous tRNA
            if current_gene is not None:
                if (current_gene_data['rRNA'] and 
                    current_gene_data['rRNA'][-1]['end'] + 1 == row['start']):
                    # Merge with previous tRNA if contiguous
                    current_gene_data['rRNA'][-1]['end'] = row['end']
                else:
                    # Otherwise, add new tRNA entry
                    current_gene_data['rRNA'].append({
                        'type': 'rRNA',
                        'start': row['start'],
                        'end': row['end'],
                        
                        
                    })   
        elif row['type'] == 'repeat_region':     
            if current_gene is not None:
                if(current_gene_data['repeat_region'] and
                    current_gene_data['repeat_region'][-1]['end'] + 1 == row['start']):
                        # Merge with previous 
                    current_gene_data['rRNA'][-1]['end'] = row['end']
                else:
                    current_gene_data['repeat_region'].append({
                        'type':'repeat_region',
                       'start': row['start'],
                        'end': row['end']
                    })
                
        
    # Append the last gene data
    if current_gene is not None:
        formatted_list.append(current_gene_data)

    return formatted_list

def format_output(formatted_list):
    formatted_text = []
    countt = 0
    for entry in formatted_list:
        # Gene line
        gene_line = f"{entry['start']}\t{entry['end']}\tgene"
        formatted_text.append(gene_line)

        # Current gene line
        current_gene_line = f"\t\t\tgene\t{entry['gene']}"
        formatted_text.append(current_gene_line)

        # Track if 'CDS' and 'tRNA' labels have been added for the gene
        cds_label_added = False
        trna_label_added = False
        rrna_label_added = False

        # Output CDS entries
        for idx, cds_entry in enumerate(entry['CDS']):
            if not cds_label_added:
                cds_line = f"{cds_entry['start']}\t{cds_entry['end']}\tCDS"
                formatted_text.append(cds_line)
                cds_label_added = True
            else:
                cds_line = f"{cds_entry['start']}\t{cds_entry['end']}"
                formatted_text.append(cds_line)

        # Add product, codon_start, and transl_table lines only if there are CDS entries
        if entry['CDS']:
            product = f"\t\t\tproduct\t{gene_product_map[entry['gene']]}" 
            codon_start_line = "\t\t\tcodon_start\t1"
            transl_table_line = "\t\t\ttransl_table\t11"
            formatted_text.extend([product, codon_start_line, transl_table_line.rstrip()])

        # Output tRNA entries
        for idx, trna_entry in enumerate(entry['tRNA']):
            if not trna_label_added:
                trna_line = f"{trna_entry['start']}\t{trna_entry['end']}\ttRNA"
                formatted_text.append(trna_line)
                trna_label_added = True
            else:
                trna_line = f"{trna_entry['start']}\t{trna_entry['end']}"
                formatted_text.append(trna_line)
                
            # Add product line after tRNA entries
        if entry['tRNA']:
            product = f"\t\t\tproduct\t{gene_product_map[entry['gene']]}" 
            formatted_text.append(product.rstrip())

       
            
        for idx, rrna_entry in enumerate(entry['rRNA']):
            if not rrna_label_added:
                rrna_line = f"{rrna_entry['start']}\t{rrna_entry['end']}\trRNA"
                formatted_text.append(rrna_line)
                rrna_label_added = True
            else:
                rrna_line = f"{rrna_entry['start']}\t{rrna_entry['end']}"
                formatted_text.append(rrna_line)
         # Add product line after tRNA entries
        if entry['rRNA']:
            product = f"\t\t\tproduct\t{gene_product_map[entry['gene']]}" 
            formatted_text.append(product.rstrip())
            
        for idx, repeating_entry in enumerate(entry['repeat_region']):
            if entry['repeat_region']:
                countt+=1
            repeating_line = f"{repeating_entry['start']}\t{repeating_entry['end']}\trepeat_region"
            formatted_text.append(repeating_line)
            
            if countt == 1:
                product = f'\t\t\tnote\tinverted repeat A'
            elif countt ==2:
                product = f'\t\t\tnote\tinverted repeat B'
            	
            formatted_text.append(product.rstrip())
            
    return "\n".join(formatted_text)

feature_name= f'>feature\t'
def save_to_file(formatted_output, file_path,feature_line):
    #feature_line=f'>feature\t'
    with open(file_path, 'w') as file:
        
        file.write(feature_line)
        file.write(formatted_output)
    print(f"Formatted output saved to {file_path}")

if __name__ == '__main__':
    app.run(debug=True)
#some commented lines