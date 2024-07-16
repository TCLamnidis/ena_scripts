#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
from lxml import etree

# import xml.dom.minidom
# def pretty_print_xml(xml_string):
#   dom = xml.dom.minidom.parse(xml_fname) # or xml.dom.minidom.parseString(xml_string)
#   pretty_xml_as_string = dom.toprettyxml()
#   print(pretty_xml_as_string)

## Read the sample registration csv file information into a pandas dataframe
def compile_csv_info(sample_csv, runs_csv):
  with open(sample_csv, 'r') as f:
    sample_info = pd.read_csv(f, skipinitialspace=True)
  
  with open(runs_csv, 'r') as f:
    runs_info = pd.read_csv(f, skipinitialspace=True)
  
  ## Join the tw tables by the 'id' column ('sampleId' in runs_csv)
  df = pd.merge(sample_info, runs_info, left_on='id', right_on='sampleId', how='right', suffixes=['_sample', '_runs'])
  
  return df

def create_submission_xml(fn):
  submission = etree.Element('SUBMISSION')
  etree.SubElement(
    etree.SubElement(
      etree.SubElement(submission, 'ACTIONS'),
      'ACTION'
      ), 'ADD')
  xml_string=etree.tostring(submission, pretty_print=True).decode()
  
  with open(fn, 'w') as output_fn:
    output_fn.write(xml_string)

## Helper function to create multiple tags in the XML if the input is a list
def create_tag(root, tag, attribute_name, values):
  ## The name of the attribute to assign to the tag is give as a string in attribute_name
  ## If the values is a list, create multiple tags
  if isinstance(values, list):
    for value in values:
      etree.SubElement(root, tag, **{attribute_name: value})
  else:
    etree.SubElement(root, tag, **{attribute_name: values})

## Create the XML for the analysis bam
def create_xml(analysis_alias=None, analysis_title=None, analysis_description=None, study_accessions=None, sample_accessions=None, run_accessions=None, reference_name=None, file_name=None, file_md5=None):
  ## Check if all arguments are not None
  if not all([analysis_alias, analysis_title, analysis_description, study_accessions, sample_accessions, run_accessions, reference_name, file_name, file_md5]):
    raise ValueError('All arguments must be provided')
  
  ## Create the XML
  root = etree.Element('ANALYSIS_SET')
  analysis = etree.SubElement(root, 'ANALYSIS', alias=analysis_alias)
  etree.SubElement(analysis, 'TITLE').text = analysis_title
  etree.SubElement(analysis, 'DESCRIPTION').text = analysis_description
  create_tag(analysis, 'STUDY_REF', "accession", study_accessions)
  create_tag(analysis, 'SAMPLE_REF', "accession", sample_accessions)
  create_tag(analysis, 'RUN_REF', "accession", run_accessions)
  ref_alignment_analysis = etree.SubElement(etree.SubElement(analysis, 'ANALYSIS_TYPE'), 'REFERENCE_ALIGNMENT')
  create_tag(etree.SubElement(ref_alignment_analysis, 'ASSEMBLY'), 'STANDARD', "accession", reference_name)
  analysis_files = etree.SubElement(analysis, 'FILES')
  ## The fiels have multiple attributes, so they are not created with the custom function.
  etree.SubElement(analysis_files, 'FILE', checksum_method="MD5", checksum=file_md5, filename=file_name, filetype="bam")
  
  return etree.tostring(root, pretty_print=True).decode()

def main():
  parser = argparse.ArgumentParser(description='Create XML for the analysis bams for upload to the ENA.')
  parser.add_argument('sample_csv', help='Path to the ENA Sample information csv file')
  parser.add_argument('runs_csv', help='Path to the ENA Runs information csv file')
  parser.add_argument('bam_annotation', help='Path to the bam annotation csv file')
  parser.add_argument('md5sums', help='Path to the md5sum file')
  parser.add_argument('-p', '--prefix_ftp', required=False ,default='', type=str, help='The prefix of the path to the files in the ftp server. If not provided, the path will be the same as the file name.')
  parser.add_argument('-o', '--output_dir', required=False, default='.', type=str, help='The output directory for the xml files')
  args = parser.parse_args()
  
  ## Throw error if output dir is empty string or root directory
  if args.output_dir == '' or args.output_dir == '/':
    sys.exit('Error: The output directory cannot be empty or the root directory, as that would place files in the root directory. Please provide a valid directory.')
  
  ## If the desired output directory does not exist, ask if it should be created.
  if not os.path.exists(args.output_dir):
    create_dir = input('The output directory does not exist. Do you want to create it? (y/n): ')
    if create_dir.lower() == 'y':
      os.makedirs(args.output_dir)
    else:
      sys.exit('Error: The output directory does not exist.')
  
  ## Read the csv files into a pandas dataframe
  df = compile_csv_info(args.sample_csv, args.runs_csv)
  
  with open(args.bam_annotation, 'r') as f:
    bam_annotation = pd.read_csv(f, skipinitialspace=True)
  
  ## Merge bam annotation with the dataframe
  df = pd.merge(df, bam_annotation, left_on='alias_sample', right_on='sample_name', how='right', suffixes=['_sample', '_bam'])
  
  ## Drop the columns that are not needed
  df = df.drop(
    [
      'firstCreated_sample',
      'firstPublic_sample',
      'releaseStatus_sample',
      'submissionAccountId_sample',
      'secondaryId',
      'title_sample',
      'taxId',
      'scientificName',
      'commonName',
      'alias_runs',
      'instrumentModel',
      'firstCreated_runs',
      'firstPublic_runs',
      'releaseStatus_runs',
      'submissionAccountId_runs'
    ],
    axis=1
  )
  
  ## Pivot df to have the runs, experiemnts and studies in a list
  df = df.groupby(['id_sample', 'alias_sample', 'sample_name', 'bam', 'title_bam', 'description', 'reference_name']).agg({
    'id_runs': lambda x: list(set(x)),
    'experimentId' : lambda x: list(set(x)),
    'studyId' : lambda x: list(set(x))
  }).reset_index()
  
  ## Read md5sums from provided file and add them as a column to the dataframe
  with open(args.md5sums, 'r') as f:
    md5sums = pd.read_csv(f, sep=" ", names=['file_md5','file_name'], skipinitialspace=True)
  
  ## Merge md5sums with the dataframe
  df = pd.merge(df, md5sums, left_on='bam', right_on='file_name', how='left', suffixes=['_bam', '_md5'])
  
  ## Create the XML
  for line in df.itertuples():
    xml_text = create_xml(
      analysis_alias=line.alias_sample, analysis_title=line.title_bam, analysis_description=line.description,
      study_accessions=line.studyId, sample_accessions=line.id_sample, run_accessions=line.id_runs,
      reference_name=line.reference_name, file_name="{}{}".format(args.prefix_ftp, line.bam), file_md5=line.file_md5
    )
  
    with open("{}/{}.xml".format(args.output_dir, line.file_name), 'w') as output_fn:
      output_fn.write(xml_text)
  
  ## Finally create the submission action XML
  create_submission_xml("{}/submission.xml".format(args.output_dir))

if __name__ == '__main__':
  main()
