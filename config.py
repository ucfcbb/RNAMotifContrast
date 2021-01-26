import os

##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) #####

subfamily_side_by_side_image_caption = True 	# Make it False to remove the side by side image captions for subfamilies.
family_side_by_side_image_caption = True 		# Make it False to remove the side by side image captions for families.
remove_homolog_subfamily = False 		# Make it True to remove the subfamilies, for which, the instances are homologs or nearly homologs.

allow_two_line_caption = True 			# Make it False to get single line caption in the output images.
extreme_filtering = False 				# Filter out more instaces to get more refined output.
min_align_len_threshold = 0 			# Motif instances will be filtered out if their alignment length is less than this threshold.
save_pymol_session = False 				# Make it True to save PyMOL sessions of the superimposed motif instances.
max_no_of_motifs_in_superimposition = 30	# If subfamily has more motifs than this threshold, it'll split and generate multiple superimposition

input_index_type = 'pdb'                # pdb / fasta [To define the indices type in input file. For FASTA, it is 0 indexed.]
annotation_source = 'merged'            # merged / dssr / fr3d
align_all_pair = False                  # True if pairwise alignment for all loops in all families. It will take more time but improve the zscore and the result to some extent.

whole_family_superimpose = False 		# Make it False to skip generating superimposition images for all instances of each motif family [It can save time]
draw_input_images = whole_family_superimpose	# By default, Input images are ignored when the whole family superimposed image is not generated.

collage_shape = 'rectangle'             # square / rectangle [Shape for the three parts of combined subfamily image]

# for merging components
merge_components = True                 # True: if CMGs are to be merged into subfamilies, False: to generate images for all CMGs separately.
connectivity_test_type = "percent"      # percent / count [Choose if merging or CMGs will depend on fixed connectivity or percentage]
# change following value according to connectivity_test_type
# for 'count', the value would represent at least how many acceptable edges are required to merge
connectivity_test_threshold = 50        # This value can represent percentage or fixed value based on 'connectivity_test_type'
rmsd_threshold_for_merging = 1.0        # This is one of the thresholds which is used to identify if a given alignment is good enough.
align_len_threshold_type = 'z-score'    # z-score / length [Familywise alignment length can be defined based on zscore or a fixed length]
align_len_zscore_threshold = 0          # z-score examples: -0.5, 0, 0.5, 1; length examples: 7, 8, 9, 10
connectivity_direction = "both"         # both / one-way [Connectivity need to be acceptable from both CMGs or one]
use_max_align_len_in_equation = True    # To put a upperbound to alignment length threshold (66% or maximum alignment length)

# 'root-oriented' traverses the motifs based on the similarity with the first loop of a CMG
# 'dijkstra' follows the traversal approach of Dijkstra's shortest path algorithm
traversal_algorithm = "root-oriented"   # root-oriented / dijkstra 

start_traversal_with_largest_subfamily = True 	# Force the traversal algorithm to start with the largest subfamily

# Font sizes are automatically adjusted from 'default_fontsize' to 'max_fontsize' based on text length
default_fontsize = 8                    # Default font size for the texts in combined images
max_fontsize = 25                       # Maximum font size for the texts in combined images
number_of_multiprocess = 8              # Number of multiprocessors to use for python parallel processing

# For input family names, provide a shortcode that can be used in image captions
known_motif_shortcode = {"c-loop": "CL", 
"e-loop": "EL", 
"hook-turn": "HT", 
"kink-turn": "KT",
"l1-complex": "L1C",
"reverse-kink-turn": "rKT",
"rope-sling": "RS",
"sarcin-ricin": "SR",
"tandem-shear": "TS",
"tetraloop-receptor": "TR",
"t-loop": "TL",
"GNAA": "GNAA",
"GNGA": "GNGA"}

# For input family names, provide a fullname that can be used in description of different outputs
known_motif_fullname = {"c-loop": "C-loop", 
"e-loop": "E-loop", 
"hook-turn": "Hook-turn", 
"kink-turn": "Kink-turn",
"l1-complex": "L1-complex",
"reverse-kink-turn": "reverse Kink-turn",
"rope-sling": "Rope-sling",
"sarcin-ricin": "Sarcin-ricin",
"tandem-shear": "Tandem-shear",
"tetraloop-receptor": "Tetraloop-receptor",
"t-loop": "T-loop",
"GNAA": "GNAA",
"GNGA": "GNGA"}

##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) ##### USER PARAMS (ADVANCED) #####


##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS #####

output_env = 'global'                   # global / local [global: cleaned up outputs; local: keep intermediate files]

# for analysis
is_length_adjusted_score = False
generate_bp_ann_files = True            # generate base-pair annotation file
# generate_alignment_files = True        # Making it true removes all previously generated alignments and generate new alignments

generate_similarity_graph_image = False # Make it True to generate images that represents the Similarity Graph

# for internal purpose only
generate_loop_source_info = True
show_cluster_source = True
scanx_align_to_superimposition = False
# use_pickle_file = True                  # to save time for multiple run with same input data
wait_factor = 0.1
wait_time = 10
max_wait_time = 600

is_normalized_score = False             # Need to make changes in code to provide option to make it True
global_component_organism_stat = {}
loop_cluster_source = {}
generate_multiple_orientation = False

download_attempts = 7					# if network speed is low and pdb or fasta file download becomes unsuccessfull, try increasing this value

##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS ##### DEVELOPER PARAMS #####


##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES #####
root_dir = os.getcwd()
data_dir = os.path.join(root_dir, 'data')
src_dir = os.path.join(root_dir, 'src')

scripts_dir = os.path.join(src_dir, 'scripts')
lib_dir = os.path.join(src_dir, 'my_lib')

pdbx_dir = os.path.join(data_dir, 'pdbx')
fasta_dir = os.path.join(data_dir, 'fasta')
pdb_fasta_mapping_dir = os.path.join(data_dir, 'pdb_fasta')

annotation_dir = os.path.join(data_dir, 'annotation')
alignment_dir_auto = os.path.join(data_dir, 'alignment_auto')

loop_dir = os.path.join(data_dir, 'loops')

views_dir = os.path.join(data_dir, 'views')

neares_protein_data_dir = os.path.join(data_dir, 'nearest_protein_data')

# libraries
dssr_dir = os.path.join(lib_dir, 'DSSR')
motifscanx_dir = os.path.join(lib_dir, 'RNAMotifScanX-release')
pymol_py3_root = 'pymol-py3'
pymol_py3_dir = os.path.join(os.path.expanduser("~"), pymol_py3_root + '/lib/python')
fonts_dir = os.path.join(lib_dir, 'fonts')

#temporary directory
temp_dir = os.path.join(data_dir, 'temp')

cluster_source_dir = os.path.join(data_dir, 'supercluster_sources')

base_path = os.path.dirname(root_dir)
base_path_len = len(base_path) + 1

##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES ##### DIRECTORIES #####

pdbx_url = 'https://files.rcsb.org/download/'
# fasta_url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?compressionType=uncompressed&structureIdList='
fasta_url = 'https://www.rcsb.org/fasta/entry/'
fr3d_url = "http://rna.bgsu.edu/rna3dhub/pdb/XXXX/interactions/fr3d/all/csv"

