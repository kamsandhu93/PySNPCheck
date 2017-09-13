
""" Modules PySNP Check """
""" Main Script starts at 556 """

""" MODULE TO VALIDATE INPUT IS IN CORRECT FORMAT """
# This module is used to validate the primer pair line

def validate_input(line):
#    import re

    try:
        a,b,c,d = line.split(" ")


        # if not re.search(r"[^ATGC]", b) or not re.search(r"[^ATGC]", c):
        #     print("An ambiguos base has been found in one of the primers on line:\n"
        #           + line +"\nPlease check both primers and try again. This line " +
        #           "will be skipped for now.")
        #     return(False)

        # If either primer is longer than 50 bp return false
        if len(b) > 51 or len(c) > 51:
            (print("One or both of the primers found on the following line are too "+
                   "long:\n" + line +"\nPrimers should not be over 50 bases long." +
                   " The line will be skipped for now."))
            return(False)
            
        # If either primer is shorter than 15 bp return false
        if len(b) < 14 or len(c) < 14:
            (print("One or both of the primers found on the following line are too "+
                   "short:\n" + line +"\nPrimers should not be under 15 bases short." +
                   " The line will be skipped for now."))
            return(False)

        statement = ("The chromosome parameter on the following line may be " +
                         "inccorect:\n" + line + "\nChromsomes must be a number" +
                         " between 1-23, or either an X or Y (capitalised).")
        
        # Check the chromosme parameter is a number for 1:23, or X and Y
        # This section may need editing if data from species with different
        # chromosmes are used (more than 23 etc).
        try:
            if int(d) > 23 or int(d) < 1:
                print(statement)
                return(False)
        except(ValueError):
            if (d != "X" and d != "Y"):
                print(statement)
                return(False)

    # Valuse error will occur if 4 parameters are not specified 
    except(ValueError):
        print("The following input line is not in the correct format:\n" + line +
              "\nThe correct input format can be viewed in the manual using the" +
              " -h argument. This line will be skipped for now\n.")
        return(False)

    return(True)




""" LOCAL BLAST """
# Module for Local BLAST requiring command line BLAST  to be installed.
def local_blast(Primer_seq, chromosome, blastdb_path):

    from subprocess import Popen, PIPE
    from tempfile import NamedTemporaryFile
    from os import remove
    from os import path


    # Get current working directory, used for finding the path of the supplied 
    # Blast DB, Ireelevant when user uses own BLAST DB
    dir_path = path.dirname(path.realpath(__file__))


    # create a temp txt file for the primer input to be written to
    fp = NamedTemporaryFile(suffix=".txt", delete = False,  mode='w+t')
    fp.write('>PRIMER\n')
    fp.write(Primer_seq)
    fp.close()


    

 
    # Location of blastn, may need to be specified by the user 
    # This can be avoiding by adding BLAST to the environment variable
    #blast_loc = ('/usr/local/ncbi/blast/bin/blastn')


    import os

    # If user wants to use the supplied BLAST DB
    if blastdb_path == "use_hg38":
        blastDB_loc = (dir_path +'{1}blastdb{1}hg38{1}chr{0}{1}chr{0}.fa'.format(
                       chromosome, os.sep))

    else:
        # User will input path to own BLAST DB stored as the variable path from 
        # command line argument
        blastDB_loc = blastdb_path
 
        
    # Name of temp file created    
    input_primer = str(fp.name)


    # The command to be used by Popen and run on the command line
    # Here the parameters for BLAST are specified 
    command = ("blastn -query {0} -db {1} -max_hsps 1 -evalue 1000 -word_size \
               7").format(input_primer, blastDB_loc)

    # Executing the BLAST search as subprocess
    process = Popen(command, stdout=PIPE,stderr=PIPE, shell = True)
    # Get output and error
    # Should include error trapping here if to detect if there is an 
    # Error with the BLAST run. The error will be stored in the variable 
    # stderr
    stdout, stderr = process.communicate()


    # Makes script wait untill blast is finished
    # neccasary to prevent deletion of input temp file before
    # subprocess has finished
    process_wait = process.wait()


    # Parsing BLAST output
    # Split Top HSP result into a list
    blast_output_list = (stdout.decode("utf-8")).split("\n")

    # Prints error if error with BLAST
    if stderr.decode("utf-8"):
        print("There has been a problem with Local BLAST, please check the \
              error:/n")
        print(stderr.decode)
    # Removes temp file, neccassry for windows users
    remove(fp.name)


    # Initiate some variables to store BLAST output
    hits= None
    start = None
    end = None
    strand = None



    # Parse the list and get the resulst we need
    for line in blast_output_list:


        # Get bases matched
        if line.startswith(" Iden"):
            split_line = line.split(" ")
            hits = split_line[3]

        # Find which strand primer is on (will be opposite of this)
        if line.startswith(" Strand"):
            strand = line.strip(" ")

        # Get location
        if line.startswith("Sbjct"):
            split_line = line.split(" ")


            if int(split_line[2]) > int(split_line[6]):
                start = split_line[6]
                end = split_line[2]
            else:
                start = split_line[2]
                end = split_line[6]


    # Returns a list that wul be processed by the BLAST processing module
    return([hits, start, end, strand])



""" REMOTE BLAST """
# Module for onlien BLAST
def online_blast(primer, chromosome):
    from Bio.Blast import NCBIWWW 

    # Using the BioPython wraper to query the NCBI BLAST serevers
    # The BLAST parameters are specified here 
    result_handle = NCBIWWW.qblast("blastn", "refseq_genomic_human", primer,\
                                    entrez_query="txid9606[ORGN]",\
                                    expect="1000", word_size="7",
                                    format_type="Text", perc_ident="80",
                                    hitlist_size= "20")

    # Get results and split into list for each line
    results = result_handle.read()
    results_split = results.split('\n')


    extracted_hsp = []
    
    # Parse the results and extract the information required
    i = 0
    for line in results_split:
        if line.startswith(">") and "GRCh38" in line and "Primary Assembly" in line \
        and "chromosome {}".format(chromosome) in line:
            i = 100000
        if i > 99999:
            extracted_hsp.append(line)
        if i > 99999 and line.startswith("Sbjct"):
            break
        
    # In case where no BLAST resuts are returned 
    if not extracted_hsp:
        print("No BLAST resulst found for this primer on this chromosome")


        # Initiate some variables to store BLAST output
        hits= None
        start = None
        end = None
        strand = None
        # if no BLAST results, all variables stored as None, this will be 
        # detected in the main script and cause the program to carry on to 
        # the next query



    # Parse the list and get the resulst we need
    for line in extracted_hsp:


        # Get bases matched
        if line.startswith(" Iden"):
            split_line = line.split(" ")
            hits = split_line[3]

        # Find which strand primer is on (will be opposite of this)
        if line.startswith(" Strand"):
            strand = line.strip(" ")

        # Get location
        if line.startswith("Sbjct"):
            split_line = line.split(" ")


            if int(split_line[2]) > int(split_line[6]):
                start = split_line[6]
                end = split_line[2]
            else:
                start = split_line[2]
                end = split_line[6]


    return([hits, start, end, strand])



""" BLAST RESULT PROCCESSING """
# Takes in input from BLAST results for both primers and interprets them. 
def process_blast_results(primer1_out, primer2_out):

    reverse = {}
    forward = {}
    amplicon = {}

    reverse["type"] = "reverse"
    forward["type"]= "forward"
    amplicon["type"]= "amplicon"
    # Initialise some dictiobaries for storring of results for each region

    # Checks which strand the primer was found to deduce if forward or reverse
    # Then adds the results to the correct dictionary 
    # Basically- " If primer 1 was the reverse: 
    if "Minus" in primer1_out[3]:
        reverse["pos"] = primer2_out[1] + "-" + primer2_out[2]
        reverse["bscore"]=(primer2_out[0])


        forward["pos"] = primer1_out[1] + "-" + primer1_out[2]
        forward["bscore"]=(primer1_out[0])


        amplicon_start = int(primer2_out[2]) + 1
        amplicon_end = int(primer1_out[1]) - 1
        amplicon_length = amplicon_end - amplicon_start + 1

        amplicon["pos"] = (str(amplicon_start) + "-" + str(amplicon_end))
        amplicon["length"] = (amplicon_length)
    
    # Vice versa of previous block of code
    # Basically if primer 2 was the reverse:
    elif "Minus" in primer2_out[3]:
        reverse["pos"] = primer1_out[1] + "-" + primer1_out[2]
        reverse["bscore"]=(primer1_out[0])


        forward["pos"] = primer2_out[1] + "-" + primer2_out[2]
        forward["bscore"]=(primer2_out[0])


        amplicon_start = int(primer1_out[2]) + 1
        amplicon_end = int(primer2_out[1]) - 1
        amplicon_length = amplicon_end - amplicon_start + 1

        amplicon["pos"]=(str(amplicon_start) + "-" + str(amplicon_end))
        amplicon["length"] = (amplicon_length)

    # Reyurn output is 3 dictionaries for each region for variant seearching
    return(reverse, forward, amplicon)




""" ENTREZ SNP SEARCH """
# Using Entrez to find variants 
def entrez_snps(location, chromosome, email):
    output_list = []
    start, stop = location.split("-")
    # takes in parametrs from dicts in above module
    
    # Use Biopython Entrez wrapper to query NCBI servers
    from Bio import Entrez
    Entrez.email = email
    handle = Entrez.esearch(db = "SNP",
                           term="{0}[CHR] AND Homo sapiens[ORGN] AND {1}:{2}[CPOS]".format(
                                   chromosome, start, stop), retmax=200)
    records = Entrez.read(handle)
    
    
    # If no results return:
    if records.get("Count") == "0":
        output_list.append("No Variants found")
        return(output_list)
        
        
    # List of IDs returned from ENtrez search
    IdList = list(records['IdList'])
    str_IdList = ", ".join(str(e) for e in IdList)

    #print("test to see if script carries on when no variants found")
    handle.close()
    

    # Now use the ID list to fetch the actual records for each varianmt 
    # with the specified parameters 
    handle2 = Entrez.efetch(db="SNP", id=str_IdList, rettype="docset", retmode="text")
    records2 = handle2.read()
    handle2.close()

    hits_in_list = records2.split("\n\n")

    i = 1

    # Header line for output 
    pre_output_list = []
    pre_output_list.append('Hit')
    pre_output_list.append('RS_ID')
    pre_output_list.append('Position')
    pre_output_list.append('dbSNP_add')
    pre_output_list.append('dbSNP_update')
    pre_output_list.append('Valid')
    pre_output_list.append('Type')

    output_list.append(pre_output_list)
    
    
    # Itterate over hits and store them as a listed list 
    # outer list = each hit 
    # inner list = the results from each hit i.e hit no, rs id , posiition etc 
    # ecample: [[hit1, rsid, location],[hit2, rsid, location]]
    # This format is used to make lift over easier 
    
    for hit in hits_in_list[:-1]:
        if hit == "":
            continue

        hit_lines = hit.split("\n")

        for line in hit_lines:
            if line.startswith("rs"):
                try:
                     x = line.split(" ")
                     snp_id = x[0]
                except(IndexError):
                    line = snp_id;
            if line.startswith("CHROMOSOME BASE POSITION"):
                position = line.split("=")[1]


            if line.startswith("SNP_CLASS"):
                vartype = line.split("=")[1]
            if line.startswith("CREATE_BUILD_ID"):
                createbuild = line.split("=")[1]
            if line.startswith("MODIFIED_BUILD_ID"):
                modbuild = line.split("=")[1]
            if line.startswith("VALIDATED"):
                validated = line.split("=")[1]


        pre_output_list = []
        pre_output_list.append(i)
        pre_output_list.append(snp_id)
        pre_output_list.append(position)
        pre_output_list.append(createbuild)
        pre_output_list.append(modbuild)
        pre_output_list.append(validated)
        pre_output_list.append(vartype)

        output_list.append(pre_output_list)


        i+=1

    return(output_list)
   
    
    
""" SQL """
# SQL lite module
def sql_snps(db, location, chromosome):
    import sqlite3 as lite
    
    # get start and stop from input location 
    start, stop = location.split("-")
    
    output_list = []
    
    # connect to input databse
    con = lite.connect(db)
    
    
    # Search the database
    with con:
        cur = con.cursor()
      
        # SQL statment to find table for input chromosome and then look 
        # for results between start and stop
        
        cur.execute("SELECT * FROM chrom{0} WHERE coordinate BETWEEN {1}\
                    AND {2}".format(chromosome, str(start), str(stop)))
    
        # Result rows
        rows = cur.fetchall()
        
        i = 1
        
        if len(rows) == 0: 
            output_list.append("No Variants Found")
            return(output_list)
        
        # Header line
        pre_output_list = []
        pre_output_list.append("Hit")
        pre_output_list.append("RS_ID")
        pre_output_list.append("Position")
        pre_output_list.append("dbSNP_Add")
        pre_output_list.append("dbSNP_Update")
        pre_output_list.append("valid")   
        output_list.append(pre_output_list)
        
        
        # Itterate over row results and extract data in listed list formar 
        for row in rows:  
            rs = row[1]
            pos = row[0]
            dbAdd = row[4]
            dbUpdate=row[5]
            if row[3] == "Y":
                valid = "True"
            else:
                valid = "False"
                
            pre_output_list = []
            pre_output_list.append(i)
            pre_output_list.append("rs"+rs)
            pre_output_list.append(chromosome + ":" + str(pos))
            pre_output_list.append(dbAdd)
            pre_output_list.append(dbUpdate)
            pre_output_list.append(valid)
            
          
            output_list.append(pre_output_list)
            
            i+=1
            
            
    return(output_list)



""" LIFT OVER BLAST RESULTS """
# used to change coordinates of results hg19 


# Module that takes input of hg38 location, and gives back the hg19 location
def hg38_to_hg19(coordinate, chromosome, lo):
    from pyliftover import LiftOver


    converted = lo.convert_coordinate('chr'+chromosome, int(coordinate), '-')

    if len(converted) == 0:
        print("Coordinate at " + coordinate + "could not be converted from hg38" +
              "to hg19 as the locus may not exist in hg19 (see manual for more" +
              "detail. The hg38 coordinate will be displayed")
        return(coordinate)
    else:
        result = converted[0]

        new_coordinate = result[1]
        return(str(new_coordinate))


# Convers the locations in the BLAST dic returned for the BLAST results 
# processing module
        
def convertBlastDic(dic, chromosome, lo):
    hg38loc = dic["pos"]
    start38, stop38 = hg38loc.split("-")

    start19 = hg38_to_hg19(int(start38), chromosome, lo)
    stop19 = hg38_to_hg19(int(stop38), chromosome, lo)

    hg19loc = str(start19) + "-" + str(stop19)
    dic["pos"] = hg19loc



# Converts the cooridnates in the listed list results from the variant search 
# modules
    
def convertSNPsearch(SNP_list,chromosome,lo):
    if len(SNP_list) == 1: return
    else:
        for i in range(1, len(SNP_list)):
            working_list = SNP_list[i]
            pre_hg38loc = working_list[2]
            chromo, hg38loc = pre_hg38loc.split(":")

            pre_hg19loc = hg38_to_hg19(int(hg38loc), chromosome, lo)
            hg19loc = chromo + ":" + str(pre_hg19loc)

            working_list[2] = hg19loc
            
            









""" Main Script """
        

"""
Mainscript for PySNPCheck
"""
input_text = """The program takes two inputs. 
1) A single pair of primers with the four parameters: Primer 
pair name, primer one, primer two, and the letter or number of the chromosome 
the primers are found on. The following format must be used for this type of 
input: “(Primer pair name) (primer one) (primer two) (chromosome number/lette
r)” where each parameter is separated by a single whitespace, and typed in 
the order specified. 2) The path to a single “.txt” file containing several
 lines of primers in format specified in 1.""" 
 
rb_text = """Use remote BLAST via the NCBI servers for sequence searching."""
lb_text = """Use local BLAST against a custom BLAST database. The direct path 
to the custom BLAST database must be supplied after this argument."""
sql_text = """Use a custom SQlite databse for variant searching. The full path 
to the SQlite database must be specified after the argument. More information on 
the schema for a custom SQlite databse can be found in the README file in the 
root of the GitHub repository for PySNPCheck."""
hg19_text = """Outputs results with Hg19 cooridnates instead of the default 
Hg38. Note this argument should not be used in conjunction with –sql argument 
unless the specified database is in the Hg38 format."""
 

# ARGUMENTS 

import argparse

parser = argparse.ArgumentParser(description="LOCALSNPCHECK")
parser.add_argument('input', nargs='+', help=input_text)
parser.add_argument('-rb', '--RemoteBLAST', action='store_true', help=rb_text)
parser.add_argument('-lb', '--LocalBLAST', metavar='[Path to BLAST database]' ,
                    help=lb_text)
parser.add_argument('-sql', '--SQlite', metavar='[Path to SQlite database]' ,
                    help=sql_text)

parser.add_argument('-hg19', '--hg19', action='store_true', help=hg19_text)

args = parser.parse_args()

   

from datetime import datetime
startTime = datetime.now()
import modulesv2

""" Input parametrs """

# Decide which BLAST and Variant search module to use based on command line 
# input 

max_amp_size = 5000

if args.RemoteBLAST: seq_search = "blast_remote"
elif args.LocalBLAST: seq_search = "blast_local_user"
else: seq_search = "default_blast_local"



if args.SQlite: snp_search = "sql"
else: snp_search = "Entrez"


""" Input processing """
# Error trap each line of input using the input module

list_input = []

# If file supplied:
if len(args.input) == 1 and ".txt" in (args.input)[0]:
    f = open(args.input[0])
    for line in f:

        if modulesv2.validate_input(line[:-1]):# module to test eror trap input
            list_input.append(line[:-1])
    f.close()

# No file stanadard input
elif len(args.input) == 4:
    str_input = " ".join(str(e) for e in args.input)
    if modulesv2.validate_input(str_input):
        list_input.append(str_input)

else: # There has been an error in input
    print("Check input for error")
    #TERMINATE
query_num = 1



""" Loop over input list and perfom main code """
for line in list_input:
    
    # Split input line into 4 variables
    primer_name, primer1, primer2, chromosome = line.split(" ")


    print("###### Query {}:{} ######".format(query_num, primer_name))

    # Do the BLAST based on command line parameters
    """BLAST"""
    if seq_search == "default_blast_local":
        one_primer_blasted = modulesv2.local_blast(primer1, chromosome, "use_hg38")
        two_primer_blasted = modulesv2.local_blast(primer2, chromosome, "use_hg38")
    elif seq_search == "blast_remote":
        one_primer_blasted = modulesv2.online_blast(primer1, chromosome)
        two_primer_blasted = modulesv2.online_blast(primer2, chromosome)
    else:
        one_primer_blasted = modulesv2.local_blast(primer1, chromosome, args.LocalBLAST)
        two_primer_blasted = modulesv2.local_blast(primer2, chromosome, args.LocalBLAST)
      
    
    # Checks if BLAST was successful for both primers and skips if not.
    if one_primer_blasted[0] == None or two_primer_blasted[0] == None:
        print("BLAST was unsuccessful, skipping line")
        continue
        
    """PROCESS BLAST SEARCH"""
    reverse, forward, amplicon = modulesv2.process_blast_results(one_primer_blasted,
                                                              two_primer_blasted)
    



    """MAXIMUM AMPLICON LENGTH """
    # Checks max amplicon length i svalid
    

    if int(amplicon["length"]) > max_amp_size:
        print("Amplicon for " + primer_name +	 " is bigger than the max amplicon" +
              " size: " + str(max_amp_size) + ". Skipping...")
        continue


    """Variant SEARCH"""

    # Decide which variant search module to use based on command line 
    # arguments 
    
    if snp_search == "sql":
       sql_db =  args.SQlite
       reverse_SNPs =  modulesv2.sql_snps(sql_db, reverse["pos"], chromosome)
       forward_SNPs = modulesv2.sql_snps(sql_db, forward["pos"], chromosome)
       amplicon_SNPs = modulesv2.sql_snps(sql_db, amplicon["pos"], chromosome)


#    elif snp_search == "mv":
#       reverse_SNPs =  modulesv2.mv_snps(reverse["pos"], chromosome)
#       forward_SNPs = modulesv2.mv_snps(forward["pos"], chromosome)
#       amplicon_SNPs = modulesv2.mv_snps(amplicon["pos"], chromosome)


    else:
        # Use Entrez module, needs email specifying 
        with open("email.txt", "r") as emailF:
            j = 1
            for line in emailF:
                if j != 1:
                    break
                pre_email = line
                j+=1

            email = pre_email.strip()

        reverse_SNPs =  modulesv2.entrez_snps(reverse["pos"], chromosome, email)
        forward_SNPs = modulesv2.entrez_snps(forward["pos"], chromosome, email)
        amplicon_SNPs = modulesv2.entrez_snps(amplicon["pos"], chromosome, email)




    """LIFT OVER FROM HG38 to HG19"""

    # If lift over required 
    
    if args.hg19:
        from pyliftover import LiftOver
        lo = LiftOver('hg38ToHg19.over.chain')
    
        modulesv2.convertBlastDic(reverse, chromosome, lo)
        modulesv2.convertBlastDic(forward, chromosome, lo)
        modulesv2.convertBlastDic(amplicon, chromosome, lo)
    
        modulesv2.convertSNPsearch(reverse_SNPs, chromosome, lo)
        modulesv2.convertSNPsearch(forward_SNPs, chromosome, lo)
        modulesv2.convertSNPsearch(amplicon_SNPs, chromosome, lo)
    



    """OUTPUT RESULTS"""
    # Printing results 

    print("Reverse:")
    print(reverse["bscore"])
    print(reverse["pos"])
    for line in reverse_SNPs:
        output_line = " ".join(str(e) for e in line)
        print(output_line)
    print("\n")


    print("Forward:")
    print(forward["bscore"])
    print(forward["pos"])
    for line in forward_SNPs:
        output_line = " ".join(str(e) for e in line)
        print(output_line)
    print("\n")

    print("Amplicon:")
    print(amplicon["pos"])
    for line in amplicon_SNPs:
        output_line = " ".join(str(e) for e in line)
        print(output_line)
    print("\n")

    query_num +=1





import time


#Print assembly used
if args.hg19: print("Human Genome Assembly: hg19")
else: print("Human Genome Assembly: hg38")
print("dbSNP Build: 150")
now = time.strftime("%c")
print ("\nResults published: " + time.strftime("%c"))
print("Time taken: " + str(datetime.now() - startTime))
