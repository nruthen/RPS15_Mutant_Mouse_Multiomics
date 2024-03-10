import numpy
import pandas
import os
import random
import shutil
import argparse

# Extracted from the following python notebook: https://app.terra.bio/#workspaces/cll-mouse/RPS15/analysis/launch/gtf_to_bed12.ipynb
# Default environemnt: (GATK 4.2.4.0, Python 3.7.12, R 4.2.1)

# List of packages in Python Notebook environment
# Package	Version
# absl-py	1.0.0
# aiohttp	3.8.1
# aiosignal	1.2.0
# ansiwrap	0.8.4
# anyio	3.5.0
# apache-beam	2.36.0
# appdirs	1.4.4
# argcomplete	2.0.0
# argon2-cffi	21.3.0
# argon2-cffi-bindings	21.2.0
# arrow	1.2.2
# arviz	0.12.1
# asn1crypto	1.4.0
# astroid	2.11.6
# astunparse	1.6.3
# async-timeout	4.0.2
# asynctest	0.13.0
# attrs	21.4.0
# Babel	2.9.1
# backcall	0.2.0
# backports.functools-lru-cache	1.6.4
# bagit	1.8.1
# beatrix-jupyterlab	3.1.7
# beautifulsoup4	4.11.1
# bgzip	0.3.5
# binaryornot	0.4.4
# biopython	1.79
# black	22.1.0
# bleach	5.0.1
# blinker	1.4
# bokeh	2.4.3
# brewer2mpl	1.4.1
# brotlipy	0.7.0
# bx-python	0.8.13
# CacheControl	0.12.11
# cached-property	1.5.2
# cachetools	5.0.0
# certifi	2022.6.15
# cffi	1.15.0
# cftime	1.6.1
# chardet	4.0.0
# charset-normalizer	2.0.12
# cli-builder	0.1.5
# click	8.0.4
# cloudpickle	2.0.0
# cmake	3.22.2
# colorama	0.4.4
# coloredlogs	15.0.1
# conda	4.11.0
# conda-package-handling	1.7.3
# cookiecutter	1.7.3
# crcmod	1.7
# cromshell	2.0.0
# cryptography	36.0.1
# cwltool	3.1.20220802125926
# cycler	0.11.0
# Cython	0.29.32
# daal	2021.6.0
# daal4py	2021.6.3
# dataclasses	0.8
# db-dtypes	1.0.3
# debugpy	1.5.1
# decorator	5.1.1
# defusedxml	0.7.1
# deprecat	2.1.1
# descartes	1.1.0
# dill	0.3.4
# dm-tree	0.1.6
# docker	5.0.3
# docker-pycreds	0.4.0
# docopt	0.6.2
# entrypoints	0.4
# explainable-ai-sdk	1.3.3
# explainers	0.1
# fastavro	1.4.9
# fasteners	0.17.3
# fastinterval	0.1.1
# fastprogress	1.0.3
# filelock	3.7.1
# firecloud	0.16.31
# flatbuffers	2.0
# flit_core	3.7.1
# fonttools	4.29.1
# frozenlist	1.3.0
# fsspec	2022.2.0
# future	0.18.2
# gast	0.4.0
# gatkpythonpackages	0.1
# gcsfs	2022.2.0
# getm	0.0.4
# ggplot	0.11.5
# gitdb	4.0.9
# GitPython	3.1.27
# google-api-core	2.8.2
# google-api-python-client	2.38.0
# google-apitools	0.5.31
# google-auth	2.6.0
# google-auth-httplib2	0.1.0
# google-auth-oauthlib	0.4.6
# google-cloud-aiplatform	1.10.0
# google-cloud-appengine-logging	1.1.0
# google-cloud-audit-log	0.2.0
# google-cloud-bigquery	2.34.4
# google-cloud-bigquery-datatransfer	3.7.0
# google-cloud-bigquery-storage	2.12.0
# google-cloud-bigtable	2.5.2
# google-cloud-core	2.2.2
# google-cloud-dataproc	3.3.0
# google-cloud-datastore	2.4.0
# google-cloud-dlp	3.6.0
# google-cloud-firestore	2.3.4
# google-cloud-kms	2.11.0
# google-cloud-language	2.3.2
# google-cloud-logging	3.0.0
# google-cloud-monitoring	2.8.0
# google-cloud-pubsub	2.9.0
# google-cloud-pubsublite	1.4.0
# google-cloud-recommendations-ai	0.2.0
# google-cloud-resource-manager	1.6.0
# google-cloud-scheduler	2.6.0
# google-cloud-spanner	3.13.0
# google-cloud-speech	2.12.0
# google-cloud-storage	1.44.0
# google-cloud-tasks	2.8.0
# google-cloud-translate	3.7.0
# google-cloud-videointelligence	2.6.0
# google-cloud-vision	2.6.3
# google-crc32c	1.1.2
# google-pasta	0.2.0
# google-resumable-media	2.3.3
# googleapis-common-protos	1.56.4
# greenlet	1.1.2
# grpc-google-iam-v1	0.12.4
# grpcio	1.44.0
# grpcio-gcp	0.2.2
# grpcio-status	1.44.0
# gs-chunked-io	0.5.2
# h5py	3.7.0
# hdfs	2.6.0
# horovod	0.23.0
# html5lib	1.1
# htmlmin	0.1.12
# httplib2	0.20.4
# humanfriendly	10.0
# idna	3.3
# ImageHash	4.2.1
# imageio	2.16.0
# importlib-metadata	4.11.1
# importlib-resources	5.4.0
# intel-openmp	2022.1.0
# ipykernel	6.9.1
# ipython	7.32.0
# ipython-genutils	0.2.0
# ipython-sql	0.3.9
# ipywidgets	7.6.5
# isodate	0.6.1
# isort	5.10.1
# jedi	0.18.1
# jeepney	0.7.1
# Jinja2	3.1.2
# jinja2-time	0.2.0
# jmespath	0.10.0
# joblib	1.1.0
# json5	0.9.5
# jsonschema	4.4.0
# jupyter	1.0.0
# jupyter-client	7.1.2
# jupyter-console	6.4.0
# jupyter-contrib-core	0.3.3
# jupyter-contrib-nbextensions	0.5.1
# jupyter-core	4.9.2
# jupyter-highlight-selected-word	0.2.0
# jupyter-http-over-ws	0.0.8
# jupyter-latex-envs	1.4.6
# jupyter-nbextensions-configurator	0.4.1
# jupyter-server	1.13.5
# jupyter-server-mathjax	0.2.5
# jupyter-server-proxy	3.2.1
# jupyterlab	3.2.9
# jupyterlab-git	0.34.2
# jupyterlab-pygments	0.1.2
# jupyterlab-server	2.10.3
# jupyterlab-widgets	1.0.2
# jupytext	1.13.7
# keras	2.7.0
# Keras-Preprocessing	1.1.2
# keras-tuner	1.1.0
# keyring	23.5.0
# keyrings.google-artifactregistry-auth	1.0.0
# kiwisolver	1.3.2
# kt-legacy	1.0.4
# kubernetes	22.6.0
# lazy-object-proxy	1.7.1
# libclang	13.0.0
# libcst	0.4.1
# llvmlite	0.38.0
# lockfile	0.12.2
# lxml	4.9.0
# Mako	1.2.1
# Markdown	3.3.6
# markdown-it-py	1.1.0
# MarkupSafe	2.0.1
# matplotlib	3.5.1
# matplotlib-inline	0.1.3
# matplotlib-venn	0.11.7
# mccabe	0.7.0
# mdit-py-plugins	0.3.0
# missingno	0.4.2
# mistune	0.8.4
# mizani	0.7.3
# mkl	2022.1.0
# msgpack	1.0.4
# multidict	6.0.2
# multimethod	1.4
# munkres	1.1.4
# mypy-extensions	0.4.3
# nb-conda	2.2.1
# nb-conda-kernels	2.3.1
# nbclassic	0.3.5
# nbclient	0.5.11
# nbconvert	6.5.0
# nbdime	3.1.1
# nbformat	5.1.3
# nest-asyncio	1.5.4
# netCDF4	1.6.0
# networkx	2.6.3
# nose	1.3.7
# notebook	6.4.8
# notebook-executor	0.2
# numba	0.55.1
# numpy	1.21.6
# oauth2client	4.1.3
# oauthlib	3.2.0
# opt-einsum	3.3.0
# orjson	3.6.7
# overrides	6.1.0
# packaging	21.3
# palettable	3.3.0
# pandas	1.3.5
# pandas-gbq	0.17.8
# pandas-profiling	3.1.0
# pandocfilters	1.5.0
# papermill	2.3.4
# parso	0.8.3
# pathspec	0.9.0
# patsy	0.5.2
# pdoc3	0.10.0
# pexpect	4.8.0
# phik	0.12.0
# pickleshare	0.7.5
# Pillow	9.0.1
# pip	22.2.2
# platformdirs	2.5.1
# plotnine	0.8.0
# pluggy	1.0.0
# poyo	0.5.0
# prettytable	3.1.1
# prometheus-client	0.13.1
# promise	2.3
# prompt-toolkit	3.0.27
# proto-plus	1.20.3
# protobuf	3.19.4
# prov	1.5.1
# psutil	5.9.0
# ptyprocess	0.7.0
# py4j	0.10.9.5
# pyarrow	7.0.0
# pyasn1	0.4.8
# pyasn1-modules	0.2.7
# pycosat	0.6.3
# pycparser	2.21
# pydantic	1.9.0
# pydata-google-auth	1.4.0
# pydot	1.4.2
# pyfasta	0.5.2
# Pygments	2.12.0
# PyJWT	2.3.0
# pylint	1.7.2
# pymc3	3.11.5
# pymongo	3.12.3
# pyOpenSSL	22.0.0
# pyparsing	2.4.7
# pyrsistent	0.18.1
# pysam	0.19.1
# PySocks	1.7.1
# python	3.7.12
# python-datauri	1.1.0
# python-dateutil	2.8.2
# python-lzo	1.14
# python-slugify	6.1.0
# pytz	2022.1
# pyu2f	0.1.5
# PyVCF3	1.0.3
# PyWavelets	1.2.0
# PyYAML	6.0
# pyzmq	22.3.0
# qtconsole	5.2.2
# QtPy	2.0.1
# rdflib	6.2.0
# requests	2.27.1
# requests-oauthlib	1.3.1
# retrying	1.3.3
# rsa	4.8
# ruamel-yaml-conda	0.15.100
# ruamel.yaml	0.17.21
# ruamel.yaml.clib	0.2.6
# schema-salad	8.3.20220801194920
# scikit-image	0.19.2
# scikit-learn	1.0.2
# scikit-learn-intelex	2021.6.3
# scipy	1.7.3
# seaborn	0.11.2
# SecretStorage	3.3.1
# semver	2.13.0
# Send2Trash	1.8.0
# setuptools	59.8.0
# shellescape	3.8.1
# simpervisor	0.4
# six	1.16.0
# smmap	3.0.5
# sniffio	1.2.0
# soupsieve	2.3.2.post1
# SQLAlchemy	1.4.31
# sqlparse	0.4.2
# statsmodels	0.13.2
# tangled-up-in-unicode	0.1.0
# tbb	2021.6.0
# tenacity	8.0.1
# tensorboard	2.8.0
# tensorboard-data-server	0.6.1
# tensorboard-plugin-wit	1.8.1
# tensorflow	2.7.0
# tensorflow-cloud	0.1.16
# tensorflow-datasets	4.4.0
# tensorflow-estimator	2.7.0
# tensorflow-hub	0.12.0
# tensorflow-io	0.23.1
# tensorflow-io-gcs-filesystem	0.24.0
# tensorflow-metadata	1.6.0
# tensorflow-probability	0.14.1
# tensorflow-serving-api	2.7.0
# tensorflow-transform	1.6.0
# termcolor	1.1.0
# terminado	0.13.1
# terra-notebook-utils	0.9.3
# testpath	0.6.0
# text-unidecode	1.3
# textwrap3	0.9.2
# tfx-bsl	1.6.0
# Theano	1.0.5
# Theano-PyMC	1.1.2
# threadpoolctl	3.1.0
# tifffile	2021.11.2
# tinycss2	1.1.1
# toml	0.10.2
# tomli	2.0.1
# tornado	6.1
# tqdm	4.64.0
# traitlets	5.1.1
# typed-ast	1.5.2
# typing-inspect	0.7.1
# typing-utils	0.1.0
# typing_extensions	4.1.1
# ujson	5.1.0
# unicodedata2	14.0.0
# Unidecode	1.3.3
# uritemplate	4.1.1
# urllib3	1.26.8
# visions	0.7.4
# wcwidth	0.2.5
# webencodings	0.5.1
# websocket-client	1.3.1
# Werkzeug	2.2.2
# wheel	0.37.1
# widgetsnbextension	3.5.2
# witwidget	1.8.0
# wrapt	1.13.3
# xai-tabular-widget	0.1.0
# xarray	0.20.2
# xarray-einstats	0.2.2
# xgboost	1.6.1
# yarl	1.7.2
# zipp	3.7.0

# Function to read in command line arguments
def getInputs():
    parser = argparse.ArgumentParser(
        description='Create stop codon readthrough reference',
        epilog='Expects a RibORF BED12 file of canonical ORFs and a subset of stop codon entries from a gtf file',
        add_help=True,
        allow_abbrev=True
    )
    parser.add_argument(
        '-bed',
        metavar='bed_file',
        help='RibORF BED12 file of canonical ORFs',
        required=True,
        type=str,
        dest='bed_file'
    )
    parser.add_argument(
        '-sub_gtf',
        metavar='stop_codon_entries',
        help='subset of stop codon entries from a gtf file',
        required=True,
        type=str,
        dest='stop_codon_entries'
    )
    parser.add_argument(
        '-out',
        metavar='Output',
        help='Path to output bed12 file.',
        default=None,
        dest='out'
    )
    return parser.parse_args()

# Function to read in BED-12 format files
def loadBed(bed):
    NAMES = ['chr','start','end','transcript_id','score','strand','tStart','tEnd','rgb','numE','lenE','startE']
    TYPES = {'chr': str,'start': int,'end': int,'transcript_id': str,'score': int,'strand': str,'tStart': int,'tEnd': int,'rgb': int,'numE': int,'lenE': str,'startE': str}
    bed = pandas.read_csv(bed,sep='\t',header=None,index_col=['transcript_id'],names=NAMES,dtype=TYPES)
    return(bed)

# Function for defining region downstream of stop codon as a bed file entry
def get_stopCodon_downstream(transcript_id):
    ### lookup the transcript in bed and stopCodon dicts
    try:
        transcript = bed_dict[transcript_id]
        stopCodon = stopCodon_dict[transcript_id]

        ### initiations 
        startE = transcript['start']+numpy.array([int(i) for i in transcript['startE'].split(",")[:-1]] )
        endE = startE+numpy.array([int(i) for i in transcript['lenE'].split(",")[:-1]] )
        span = 29
        p = stopCodon['stop_codon_end']+1
        stop_startE = []
        stop_endE = []

        
        ### point p at the correct exon to start with 
        i = None
        for j in range(len(startE)):
            if stopCodon['stop_codon_end'] == endE[j] and j+1 < len(startE):
                p = endE[j+1]
                i = j+1
                continue 
            elif startE[j] <= stopCodon['stop_codon_start']-1 <= endE[j]:
                p = stopCodon['stop_codon_start']-1
                i = j
                continue
                
        if i == None or stopCodon['stop_codon_end']+span >endE[-1]:
            print(transcript_id," : stop codon is not in transcript exon regions")
            return 
        
        ### finding the start and end of exon(s) in downstream of stop codon

        while span != 0:
                
            if startE[i] <= p + span <= endE[i]:
                error2 = ""
                stop_startE.append(p)
                stop_endE.append(p+span)
                span = 0 
            else:
                stop_startE.append(p)
                stop_endE.append(endE[i])
                span = span -(endE[i]-p) -1 
                try:
                    p = startE[i+1]
                    i += 1 
                except IndexError:
                    print(transcript_id," : stop codon downstream is not in transcript exon regions")
                    return


        new_transcript = {}
        new_transcript['chr'] = transcript['chr']
        new_transcript['score'] = transcript['score']
        new_transcript['strand'] = transcript['strand']
        new_transcript['rgb'] = transcript['rgb']
        new_transcript['start'] = stop_startE[0]
        new_transcript['end'] = stop_endE[-1]
        new_transcript['tStart'] = new_transcript['start'] 
        new_transcript['tEnd'] = new_transcript['end'] 
        new_transcript['numE'] = len(stop_startE)
        new_transcript['startE'] = numpy.array(stop_startE) - new_transcript['start']
        new_transcript['lenE'] = numpy.array(stop_endE) - numpy.array(stop_startE)
        # convert to str
        new_transcript['startE'] = ",".join(new_transcript['startE'].astype(str))
        new_transcript['lenE'] = ",".join(new_transcript['lenE'].astype(str))
    except KeyError:
        new_transcript=None
    return new_transcript

# Example of function get_stopCodon_downstream() usage 
# get_stopCodon_downstream("ENST00000382287.5_2")


# Function for defining region downstream of stop codon as a bed file entry
def get_stopCodon_upstream(transcript_id):
    ### lookup the transcript in bed and stopCodon dicts
    try:
        transcript = bed_dict[transcript_id]
        stopCodon = stopCodon_dict[transcript_id]

        ### initiations 
        startE = transcript['start']+numpy.array([int(i) for i in transcript['startE'].split(",")[:-1]] )
        endE = startE+numpy.array([int(i) for i in transcript['lenE'].split(",")[:-1]] )
        span = 29
        stop_startE = []
        stop_endE = []

        #print(startE)
        #print(endE)
        
        ### point p at the correct exon to start with 
        i = None
        for j in range(len(startE)):
            if stopCodon['stop_codon_start'] == startE[j]:
                p = endE[j-1]
                i = j-1
                continue 
            elif startE[j] <= stopCodon['stop_codon_start']-1 <= endE[j]:
                p = stopCodon['stop_codon_start']-1
                i = j
                continue
                
        if i == None or stopCodon['stop_codon_start']-span <startE[0]:
            print(transcript_id," : stop codon is not in transcript exon regions")
            return 

        ### finding the start and end of exon(s) in upstream of stop codon
        while span != 0:

            if startE[i] <= p - span <= endE[i]:
                error2 = ""
                stop_startE.append(p-span)
                stop_endE.append(p)
                span = 0 
            else:
                stop_startE.append(startE[i] )
                stop_endE.append(p)
                span = span -(p - startE[i]) -1 
                try:
                    p = endE[i-1]
                    i -= 1 
                except IndexError:
                    print(transcript_id," : stop codon upstream is not in transcript exon regions")
                    return
                
        stop_startE.reverse()
        stop_endE.reverse()

        new_transcript = {}
        new_transcript['chr'] = transcript['chr']
        new_transcript['score'] = transcript['score']
        new_transcript['strand'] = transcript['strand']
        new_transcript['rgb'] = transcript['rgb']
        new_transcript['start'] = stop_startE[0]
        new_transcript['end'] = stop_endE[-1]
        new_transcript['tStart'] = new_transcript['start'] 
        new_transcript['tEnd'] = new_transcript['end'] 
        new_transcript['numE'] = len(stop_startE)
        new_transcript['startE'] = numpy.array(stop_startE) - new_transcript['start']
        new_transcript['lenE'] = numpy.array(stop_endE) - numpy.array(stop_startE)
        # convert to str
        new_transcript['startE'] = ",".join(new_transcript['startE'].astype(str))
        new_transcript['lenE'] = ",".join(new_transcript['lenE'].astype(str))

    except KeyError:
        new_transcript=None

    return new_transcript

# Example of function get_stopCodon_upstream() usage 
# get_stopCodon_upstream("ENST00000614126.4_1")

# Determine the stop codon is closest to every ORF
def get_closest_stopCodons(bed_dict,stopCodon):
    row_list = list()
    print(len(row_list)) 
    for transcript_id in bed_dict.keys():
        sub_stopCodon = stopCodon.loc[stopCodon['transcript_id'] == transcript_id]
        transcript = bed_dict[transcript_id]
        startE = transcript['start']+numpy.array([int(i) for i in transcript['startE'].split(",")[:-1]] )
        endE = startE+numpy.array([int(i) for i in transcript['lenE'].split(",")[:-1]] )
        min_row_dist=10000000000
        min_row = sub_stopCodon.iloc[[0]]
        for index, row in sub_stopCodon.iterrows():
            distance_to_closest_stop = min(abs(endE-row['stop_codon_start']))
            if  distance_to_closest_stop < min_row_dist:
                min_row_dist = distance_to_closest_stop 
                min_row=row
        row_list.append(row)
    print(row_list[0])
    new_stopCodon = pandas.DataFrame(row_list, columns=stopCodon.columns)
    return new_stopCodon


# Main() Section

args = getInputs()
os.getcwd()
#os.chdir('stop_condon_read_through')
#os.listdir()

bed = loadBed(args.bed_file)
bed_dict = bed.to_dict(orient = 'index')

stopCodon = pandas.read_csv(args.stop_codon_entries,sep='\t')

# Drop all but first stop codon by coordinate
stopCodon = stopCodon.sort_values(by=['transcript_id','stop_codon_start'])
stopCodon = stopCodon.drop_duplicates(subset=['transcript_id'])

# Sort by transcript id and the position of the first nucleotide of each stop codon.
stopCodon = stopCodon.sort_values(by=['transcript_id','stop_codon_start'])

# Get closest stop codons to ORFs
temp_stopCodon = get_closest_stopCodons(bed_dict,stopCodon)
print(temp_stopCodon.head(n=5))
stopCodon = temp_stopCodon

# Reformating of object
stopCodon = stopCodon.rename({'chrname':'chr'}, axis='columns').set_index('transcript_id')
stopCodon_dict = stopCodon.to_dict(orient = 'index')
stopCodon.tail()

# Defines empty dictionary in BED12 format (which will be converted to "stopCodon_updownstream.bed" file)
new_bed_dict = {}

# Loops through transcript ids associated with list of stop codons
for transcript_id in stopCodon_dict:
    downstream = get_stopCodon_downstream(transcript_id)
    if downstream == None:
        pass
    else:
        # Defines dictionary key as "transcript_id + '_downstream'" and value as the modified BED12 downstream entry associated with the transcript id
        new_bed_dict[transcript_id+"_downstream"] = downstream
        ### exclude transcripts that don't have downstream region
        upstream = get_stopCodon_upstream(transcript_id)
        if upstream == None:
            pass
        else:
            # Defines dictionary key as "transcript_id + '_upstream'" and value as the modified BED12 entry upstream entry associated with the transcript id
            new_bed_dict[transcript_id+"_upstream"] = upstream



# Confirm that the majority of stop codons have defined upstream and downstream regions in the new dictionary (i.e. the following values shoulw be approximately equal, such as 37746 and 35748)
len(stopCodon_dict) *2
len(new_bed_dict)

# Convert Dictionary to Table:

new_bed = pandas.DataFrame.from_dict(new_bed_dict,orient='index')
new_bed.index.name = 'transcript'
new_bed = new_bed.reset_index()
new_bed.head()


new_bed = new_bed[['chr','start','end','transcript','score','strand','tStart','tEnd','rgb','numE','lenE','startE']]
new_bed = new_bed.sort_values(by=['chr','start'])
new_bed.head()

# output the stopCodon_updownstream.bed
new_bed.to_csv(args.out, header=None, index=False, sep='\t')


