"""

Manual Curation Helper Version 1.3.5
Novelties:
  * Removed the steps to create files with removed sequences
  * Modified/added messages about number of TEs processed in each module
  * Fixed some mistakes with verbose messages
  * Corrected bug with semi-auto flow
  * Corrected bug with fully-manual flow

"""

import sys
import os
import argparse
import pandas as pd
import multiprocessing
from pdf2image import convert_from_path
from pdf2image.exceptions import PDFPageCountError
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import cv2
from matplotlib import pyplot as plt
import psutil
import glob
import shutil
import time

dicc_orders = {2: "LTR", 3: "COPIA", 4: "GYPSY", 5: "BELPAO",
               6: "TRIM", 7: "LARD", 8: "LINE", 9: "SINE", 10: "TIR", 11: "MITE", 12: "HELITRON", 13: "MAVERICK",
               14: "PLE", 15: "DIRS", 16: "R2", 17: "RTE", 18: "JOCKEY", 19: "L1", 20: "I", 21: "TC1MARINER",
               22: "HAT", 23: "MUTATOR", 24: "MERLIN", 25: "TRANSIB", 26: "P", 27: "PIGGYBAC", 28: "PIFHARBINGER",
               29: "CACTA", 30: "CRYPTON", 31: "UNCLASSIFIED", 31: "UNKNOWN", 32: "RETROTRANSPOSON", 33: "TRANSPOSON"}


def create_fasta_files(ref_tes, keep_seqs, file_name, outputdir, orders):
    record_curated = []
    create_output_folders(outputdir)
    try:
        for sequence in SeqIO.parse(ref_tes, "fasta"):
            real_name = str(sequence.id).split('#')[0]
            if real_name in keep_seqs:
                if len(orders) > 0:
                    pos = keep_seqs.index(str(real_name))
                    classification = dicc_orders[orders[pos]]
                    if orders[pos] >= 3 and orders[pos] <= 5:
                        classification = "LTR/" + classification
                    elif orders[pos] >= 16 and orders[pos] <= 20:
                        classification = "LINE/" + classification
                    elif orders[pos] >= 21 and orders[pos] <= 30:
                        classification = "TIR/" + classification
                    elif orders[pos] == 31:
                        classification = "UNCLASSIFIED"
                    sequence.id = real_name + "#" + str(classification)
                    sequence.description = ''
                    record_curated.append(sequence)
                else:
                    record_curated.append(sequence)

        SeqIO.write(record_curated, file_name, "fasta")
    except FileNotFoundError:
        print(
            "FATAL ERROR: I couldn't find the files, please check: '" + ref_tes + "' and '" + file_name + "'. Path not found")
        sys.exit(0)
    except PermissionError:
        print(
            "FATAL ERROR: I couldn't access the files, please check: ' " + ref_tes + "' and '" + file_name + "'. I don't have permissions.")
        sys.exit(0)
    except Exception as ex:
        print("FATAL ERROR: There is a unknown problem writing sequences in : '" + outputdir + "'.")
        print(ex.args)
        sys.exit(0)


def write_sequences_file(sequences, filename):
    try:
        SeqIO.write(sequences, filename, "fasta")
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't find the file, please check: '" + filename + "'. Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't access the files, please check: '" + filename + "'. I don't have permissions.")
        sys.exit(0)
    except Exception as exp :
        print("FATAL ERROR: There is a unknown problem writing sequences in : '" + filename + "'.")
        print(exp)
        sys.exit(0)


def create_output_folders(folder_path):
    try:
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
    except FileNotFoundError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". Path not found")
        sys.exit(0)
    except PermissionError:
        print("FATAL ERROR: I couldn't create the folder " + folder_path + ". I don't have permissions.")
        sys.exit(0)


def delete_files(file_path):
    try:
        os.remove(file_path)
    except FileNotFoundError:
        print("WARNING: The file " + file_path + " couldn't be removed. File not found.")
    except PermissionError:
        print("WARNING: The file " + file_path + " couldn't be removed because I don't have permissions.")


def filter_flf(ref_tes, flf_file, minFullLenCopies, outputdir):
    flf_lines = open(flf_file, 'r').readlines()[1:]
    seqsID_with_flf = [s.split("\t")[0] for s in flf_lines]
    seqs_with_flf = []
    numCopies = {}
    for te in SeqIO.parse(ref_tes, "fasta"):
        seq_name = str(te.id).split("#")[0]
        if seq_name in seqsID_with_flf:
            numfullCopies = int([s.split("\t")[6] for s in flf_lines if s.split("\t")[0] == seq_name][0])
            numfullFrags = int([s.split("\t")[4] for s in flf_lines if s.split("\t")[0] == seq_name][0])
            if numfullFrags >= minFullLenCopies:  # number of copies greater or equal than threshold?
                seqs_with_flf.append(te)
                numCopies[seq_name] = [numfullFrags, numfullCopies]

    write_sequences_file(seqs_with_flf, outputdir + "/cons_flf.fa")
    return outputdir + "/cons_flf.fa", numCopies


def processing_gff(seqname, gff_files, pre):
    posBlan = []
    postBlax = []
    posBlax = []
    posTR = []

    gff_blastn = open(gff_files + '/' + pre + '_TE_BLRn.gff3', 'r').readlines()
    for anno in gff_blastn:
        columns = anno.split('\t')
        # if columns[0] == seqname and columns[2] == 'TE_BLRn':
        if columns[0] == seqname and columns[2] == 'exon':
            posBlan.append((int(columns[3]), int(columns[4]) - int(columns[3])))

    gff_tblastx = open(gff_files + '/' + pre + '_TE_BLRtx.gff3', 'r').readlines()
    for anno in gff_tblastx:
        columns = anno.split('\t')
        # if columns[0] == seqname and columns[2] == 'TE_BLRtx':
        if columns[0] == seqname and columns[2] == 'exon':
            postBlax.append((int(columns[3]), int(columns[4]) - int(columns[3])))

    gff_blastx = open(gff_files + '/' + pre + '_TE_BLRx.gff3', 'r').readlines()
    for anno in gff_blastx:
        columns = anno.split('\t')
        # if columns[0] == seqname and columns[2] == 'TE_BLRx':
        if columns[0] == seqname and columns[2] == 'exon':
            posBlax.append((int(columns[3]), int(columns[4]) - int(columns[3])))

    gff_tr = open(gff_files + '/' + pre + '_TR.gff3', 'r').readlines()
    for anno in gff_tr:
        columns = anno.split('\t')
        if columns[0] == seqname and columns[2] == 'exon':
            posTR.append((int(columns[3]), int(columns[4]) - int(columns[3])))

    return posBlan, postBlax, posBlax, posTR


def count_flf_fasta(ref_tes, genome, cores, outputdir):
    create_output_folders(outputdir)
    create_output_folders(outputdir + "/classifiedModule")

    flf_df = pd.DataFrame(columns=['TE', 'length', 'covg', 'frags', 'fullLgthFrags', 'copies', 'fullLgthCopies'])

    if not os.path.exists(genome + ".nhr"):
        output = subprocess.run(
            ['makeblastdb', '-in', genome, '-dbtype', 'nucl'], stdout=subprocess.PIPE, text=True)

    output = subprocess.run(
        ['blastn', '-query', ref_tes, '-db', genome, '-out', outputdir + "/classifiedModule/TEs_vs_genome.blast", '-num_threads',
         str(cores),
         "-outfmt", "6", "-evalue", "10e-8"], stdout=subprocess.PIPE, text=True)

    blastresult = pd.read_table(outputdir + "/classifiedModule/TEs_vs_genome.blast", sep='\t',
                                names=['qseqid', 'sseqid', 'pident',
                                       'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                                       'bitscore'])

    for te in SeqIO.parse(ref_tes, "fasta"):
        seq_name = str(te.id).split("#")[0]
        te_hits = blastresult[blastresult.qseqid == str(te.id)].reset_index()
        count_frag = 0
        count_flf = 0
        te_len = len(str(te.seq))
        for i in range(te_hits.shape[0]):
            len_hit = te_hits.loc[i, 'length']
            if len_hit >= (te_len * 0.94):  # it's a full-length fragment
                count_flf += 1
            else:
                count_frag += 1
        flf_df = flf_df.append(
            {'TE': seq_name, 'length': te_len, 'covg': 0, 'frags': count_frag, 'fullLgthFrags': count_flf,
             'copies': count_frag, 'fullLgthCopies': count_flf}, ignore_index=True)

    flf_df.to_csv(outputdir + "/classifiedModule/fullLengthFrag.txt", header=True, sep='\t', index=False)
    delete_files(outputdir + "/classifiedModule/TEs_vs_genome.blast")
    return outputdir + "/classifiedModule/fullLengthFrag.txt"


def check_classification_Userlibrary(user_library, outputdir):
    orders = ['LTR', 'TRIM', 'LARD', 'LINE', 'SINE', 'PLE', 'DIRS', 'TIR', 'MITE', 'HELITRON', 'MAVERICK', 'CRYPTON',
              'UNCLASSIFIED']
    superfamilies = ["COPIA", "GYPSY", "BELPAO", "PLE", "DIRS", "R2", "RTE", "JOCKEY", "L1", "I", 'SINE', "TC1MARINER",
                     "HAT", "MUTATOR", "MERLIN", "TRANSIB", "P", "PIGGYBAC", "PIFHARBINGER", "CACTA", "MITE",
                     "HELITRON", "MAVERICK", "CRYPTON", "UNCLASSIFIED", "UNKNOWN"]
    ids_no_contained = []
    reasons = []
    for te in SeqIO.parse(user_library, "fasta"):
        if "#" not in str(te.id):  # it doesn't have the correct format
            ids_no_contained.append(te.id)
            reasons.append(
                "The sequence ID doesn't have the character '#' needed to separate the ID to the classification")
        else:
            seq_name = str(te.id).split("#")[0]
            classification = str(te.id).split("#")[1].upper().replace("-", "")
            if "/" in classification:
                order_given = classification.split("/")[0]
                superfamily = classification.split("/")[1]
            else:
                order_given = classification
                superfamily = "NA"
            if superfamily not in superfamilies and order_given not in superfamilies:
                if order_given not in orders:
                    ids_no_contained.append(te.id)
                    reasons.append(
                        "The classification wasn't find in my classification system. Remember that I use the Wicker nomenclature.")

    if len(ids_no_contained) == 0:  # Congrats!! everything looks great!
        return 0
    else:
        output_file = open(outputdir + "/sequences_with_problems.txt", "w")
        output_file.write("Sequence ID\tProblem\n")
        for i in range(len(ids_no_contained)):
            output_file.write(ids_no_contained[i] + "\t" + reasons[i] + "\n")
        return -1


def check_repet_input_folder(repet_input_dir, proj_name):
    ref_tes = repet_input_dir + "/" + proj_name + "_refTEs.fa"
    flf_file = repet_input_dir + "/TEannot/" + proj_name + "_chr_allTEs_nr_noSSR_join_path.annotStatsPerTE_FullLengthFrag.txt"
    repet_table = repet_input_dir + "/" + proj_name + "_denovoLibTEs_PC.classif"
    plots_dir = repet_input_dir + "/plotCoverage"
    gff_files = repet_input_dir + "/gff_reversed"

    if not os.path.exists(repet_input_dir):
        valid = False
        reason = "FATAL ERROR: "+repet_input_dir+" does not exist. Please check the path and re-run the software"
    elif not os.path.exists(ref_tes):
        valid = False
        reason = "FATAL ERROR: "+ref_tes+" does not exist. Please check the path and re-run the software"
    elif os.path.getsize(ref_tes) == 0:
        valid = False
        reason = "FATAL ERROR: "+ref_tes+" is empty. Please check the file and re-run the software"
    elif not os.path.exists(flf_file):
        valid = False
        reason = "FATAL ERROR: "+flf_file+" does not exist. Please check the path and re-run the software"
    elif os.path.getsize(flf_file) == 0:
        valid = False
        reason = "FATAL ERROR: "+flf_file+" is empty. Please check the file and re-run the software"
    elif not os.path.exists(repet_table):
        valid = False
        reason = "FATAL ERROR: "+repet_table+" does not exist. Please check the path and re-run the software"
    elif os.path.getsize(repet_table) == 0:
        valid = False
        reason = "FATAL ERROR: "+repet_table+" is empty. Please check the file and re-run the software"
    elif not os.path.exists(plots_dir):
        valid = False
        reason = "FATAL ERROR: "+plots_dir+" does not exist. Please check the path and re-run the software"
    elif not os.path.exists(gff_files):
        valid = False
        reason = "FATAL ERROR: "+gff_files+" does not exist. Please check the path and re-run the software"
    else:
        # In this case, everything is ok to run the curation process !!
        valid = True
        reason = ""

    return valid, reason


def find_TRs(te, outputdir, polyA_path):
    lenTR = 0
    typeTR = ""
    seq_name = str(te.id).split("#")[0]
    middle = int(len(str(te.seq)) / 2)
    firstmid = SeqRecord(Seq(te.seq[0:middle]), id="fivePrime")
    secondmid = SeqRecord(Seq(te.seq[middle:len(str(te.seq))]), id="threeprime")
    write_sequences_file(te, outputdir + "/" + seq_name + ".fa")
    write_sequences_file(firstmid, outputdir + "/" + seq_name + "_fiveprime.fa")
    write_sequences_file(secondmid, outputdir + "/" + seq_name + "_threeprime.fa")

    output = subprocess.run(
        ['makeblastdb', '-in', outputdir + "/" + seq_name + "_threeprime.fa", '-dbtype', 'nucl'],
        stdout=subprocess.PIPE, text=True)

    output = subprocess.run(
        ['blastn', '-query', outputdir + "/" + seq_name + "_fiveprime.fa", '-db',
         outputdir + "/" + seq_name + "_threeprime.fa", '-out',
         outputdir + "/" + str(seq_name) + ".blast", '-num_threads', str(cores), "-outfmt", "6",
         "-evalue", "1e-20"], stdout=subprocess.PIPE, text=True)

    blastresult = pd.read_table(outputdir + "/" + str(seq_name) + ".blast", sep='\t',
                                names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                       'sstart', 'send', 'evalue', 'bitscore'])

    if blastresult.shape[0] > 0:
        lenTR = int(blastresult.loc[0, "qend"]) - int(blastresult.loc[0, "qstart"])
        if int(blastresult.loc[0, "sstart"]) > int(blastresult.loc[0, "send"]):  # it's a TIR
            typeTR = "termTIR"
        else:
            typeTR = "termLTR"
    else:
        # Here we need to put the PolyATail program
        output = subprocess.run([polyA_path+'/polyAtail', outputdir + "/" + seq_name + ".fa"],
                                stdout=subprocess.PIPE, text=True)
        if "termPolyA" in output.stdout:
            for line in output.stdout.split('\n'):
                if "termPolyA" in line and "non-termPolyA" not in line:
                    typeTR = "polyAtail"
                    lenTR = int(line.split("|")[1].split("=")[1])

    delete_files(outputdir + "/" + str(seq_name) + ".blast")
    delete_files(outputdir + "/" + seq_name + ".fa")
    if os.path.exists(outputdir + "/" + seq_name + ".fa.polyA.set"):
        delete_files(outputdir + "/" + seq_name + ".fa.polyA.set")
    delete_files(outputdir + "/" + seq_name + "_fiveprime.fa")
    delete_files(outputdir + "/" + seq_name + "_threeprime.fa")
    delete_files(outputdir + "/" + seq_name + "_threeprime.fa.nhr")
    delete_files(outputdir + "/" + seq_name + "_threeprime.fa.nin")
    delete_files(outputdir + "/" + seq_name + "_threeprime.fa.nsq")

    return lenTR, typeTR


def find_profiles(te, outputdir, ref_profiles):
    # domains:      LTR retrotransposons/DIRS/LINE
    domains_dict = {'_GAG_': 0, '_AP_': 0, '_INT_': 0, '_RT_': 0, '_RNaseH_': 0, '_ENV_': 0,
                    # PLE       TIRs         Helitron    Maverick
                    '_EN_': 0, '_Tase_': 0, '_HEL_': 0, '_Prp': 0, '_ATPase_': 0}
    seq_name = str(te.id).split("#")[0]
    write_sequences_file(te, outputdir + "/" + seq_name + "_putative_te.fa")

    output = subprocess.run(
        ['getorf', '-sequence', outputdir + "/" + seq_name + "_putative_te.fa", '-minsize', '300', '-outseq',
         outputdir + "/" + seq_name + "_putative_te_orf.fa"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    if not os.path.exists(ref_profiles+".h3m"):
        output = subprocess.run(['hmmpress', ref_profiles],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    output = subprocess.run(
        ['hmmscan', '--tblout', outputdir + "/" + seq_name + "_profiles_found.hmm", '-E', '10', '--noali', '--cpu',
         '1', ref_profiles, outputdir + "/" + seq_name + "_putative_te_orf.fa"],
         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    result = ""
    struc_domains = ""

    if os.path.exists(outputdir + "/" + seq_name + "_profiles_found.hmm"):
        command = 'grep -v "^#" ' + outputdir + '/' + seq_name + '_profiles_found.hmm | wc -l'
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        hits = int(process.stdout.read())
        if hits > 0:
            command = 'grep -v "^#" ' + outputdir + '/' + seq_name + '_profiles_found.hmm | sed "s/_RVT_/_RT_/g" | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}\' > ' + outputdir + "/" + seq_name + "_profiles_found.hmm_formatted"
            process = subprocess.Popen(command, shell=True)
            process.communicate()

            hmm_results = pd.read_table(outputdir + "/" + seq_name + "_profiles_found.hmm_formatted", sep='\t',
                                        names=['domain', 'acc', 'orf', 'acc_orf', 'evalue', 'score', 'A', 'B', 'C',
                                               'D'])
            hasDomains = False
            for domain in domains_dict.keys():
                best_hit = hmm_results[hmm_results.domain.str.contains(domain)].reset_index()
                if best_hit.shape[0] > 0:
                    hasDomains = True
                    domains_dict[domain] = best_hit.loc[0, 'evalue']
                    result += best_hit.loc[0, 'domain'] + ": " + str(best_hit.loc[0, 'evalue']) + ", "
                    struc_domains += domain.replace("_", "") + " "

            if hasDomains:
                result = "profiles: " + result[:-2]
            delete_files(outputdir + "/" + seq_name + "_profiles_found.hmm_formatted")
        delete_files(outputdir + "/" + seq_name + "_profiles_found.hmm")

    delete_files(outputdir + "/" + seq_name + "_putative_te_orf.fa")
    delete_files(outputdir + "/" + seq_name + "_putative_te.fa")

    return result, struc_domains


def find_blast_hits(te, outputdir, blastx_db, blastn_db):
    seq_name = str(te.id).split("#")[0]
    write_sequences_file(te, outputdir + "/" + seq_name + "_putative_te.fa")
    output = subprocess.run(
        ['blastx', '-query', outputdir + "/" + seq_name + "_putative_te.fa", '-db', blastx_db, '-out',
         outputdir + "/" + seq_name + "_putative_te.fa.blastx", '-num_threads', '1', "-outfmt",
         "6 sseqid pident evalue bitscore"],
        stdout=subprocess.PIPE, text=True)

    blastresult = pd.read_table(outputdir + "/" + seq_name + "_putative_te.fa.blastx", sep='\t',
                                names=['sseqid', 'pident', 'evalue', 'bitscore'])
    blastx_result = ""
    i = 0
    if blastresult.shape[0] > 0:
        blastx_result = "TE_BLRx: "
        while i < 3 and i < blastresult.shape[0]:
            blastx_result += blastresult.loc[0, "sseqid"] + ": " + str(blastresult.loc[0, "pident"]) + "%, "
            i += 1
        blastx_result = blastx_result[:-2]

    output = subprocess.run(
        ['tblastx', '-query', outputdir + "/" + seq_name + "_putative_te.fa", '-db', blastn_db, '-out',
         outputdir + "/" + seq_name + "_putative_te.fa.tblastx", '-num_threads', '1', "-outfmt",
         "6 sseqid pident evalue bitscore"],
        stdout=subprocess.PIPE, text=True)

    blastresult = pd.read_table(outputdir + "/" + seq_name + "_putative_te.fa.tblastx", sep='\t',
                                names=['sseqid', 'pident', 'evalue', 'bitscore'])
    blasttx_result = ""
    i = 0
    if blastresult.shape[0] > 0:
        blasttx_result = "TE_BLRtx: "
        while i < 3 and i < blastresult.shape[0]:
            blasttx_result += blastresult.loc[0, "sseqid"] + ": " + str(blastresult.loc[0, "pident"]) + "%, "
            i += 1
        blasttx_result = blasttx_result[:-2]

    delete_files(outputdir + "/" + seq_name + "_putative_te.fa")
    delete_files(outputdir + "/" + seq_name + "_putative_te.fa.tblastx")
    delete_files(outputdir + "/" + seq_name + "_putative_te.fa.blastx")
    return blastx_result, blasttx_result


def build_class_table_parallel(ref_tes, cores, outputdir, blastn_db, blastx_db, tools_path, ref_profiles, do_blast):
    tes = [te for te in SeqIO.parse(ref_tes, "fasta")]
    n = len(tes)
    seqs_per_procs = int(n / cores)
    remain = n % cores
    ini_per_thread = []
    end_per_thread = []
    for p in range(cores):
        if p < remain:
            init = p * (seqs_per_procs + 1)
            end = n if init + seqs_per_procs + 1 > n else init + seqs_per_procs + 1
        else:
            init = p * seqs_per_procs + remain
            end = n if init + seqs_per_procs > n else init + seqs_per_procs
        ini_per_thread.append(init)
        end_per_thread.append(end)

    # Run in parallel the checking
    pool = multiprocessing.Pool(processes=cores)
    localresults = [pool.apply_async(build_class_table,
                                     args=[tes[ini_per_thread[x]:end_per_thread[x]], ref_profiles, outputdir, blastn_db,
                                           blastx_db, tools_path, do_blast]) for x in range(cores)]

    local_dfs = [p.get() for p in localresults]
    class_df = pd.DataFrame(
        columns=['Seq_name', 'length', 'strand', 'confused', 'class', 'order', 'Wcode', 'sFamily', 'CI', 'coding',
                 'struct', 'other'])
    for i in range(len(local_dfs)):
        class_df = class_df.append(local_dfs[i], ignore_index=True)
    pool.close()

    class_df.to_csv(outputdir + "/denovoLibTEs_PC.classif", header=True, sep='\t', index=False)

    new_tes_user = []
    for te in SeqIO.parse(ref_tes, "fasta"):
        te.id = str(te.id).split("#")[0]
        te.description = ""
        new_tes_user.append(te)

    write_sequences_file(new_tes_user, outputdir + "/new_user_lib.fa")


def build_class_table(ref_tes, ref_profiles, outputdir, blastn_db, blastx_db, tools_path, do_blast):
    class_df = pd.DataFrame(
        columns=['Seq_name', 'length', 'strand', 'confused', 'class', 'order', 'Wcode', 'sFamily', 'CI', 'coding',
                 'struct', 'other'])
    orders = ['LTR', 'TRIM', 'LARD', 'LINE', 'SINE', 'PLE', 'DIRS', 'TIR', 'MITE', 'HELITRON', 'MAVERICK', 'CRYPTON']
    superfamilies = ["COPIA", "GYPSY", "BELPAO", "PLE", "DIRS", "R2", "RTE", "JOCKEY", "L1", "I", 'SINE', "TC1MARINER",
                     "HAT", "MUTATOR", "MERLIN", "TRANSIB", "P", "PIGGYBAC", "PIFHARBINGER", "CACTA", "MITE",
                     "HELITRON",
                     "MAVERICK", "CRYPTON"]
    for te in ref_tes:
        seq_name = str(te.id).split("#")[0]
        classification = str(te.id).split("#")[1].upper().replace("-", "")
        length = len(str(te.seq))
        if "/" in classification:
            order_given = classification.split("/")[0]
            superfamily = classification.split("/")[1]
        else:
            order_given = classification
            superfamily = "NA"

        if superfamily not in superfamilies:
            if order_given in superfamilies:
                superfamily = order_given

        classTE = ""
        if superfamily in superfamilies:
            # superfamily found !!!
            sFamily = superfamily
            order = ""
            sFamily_index = superfamilies.index(superfamily)
            if sFamily_index <= 2:
                order = "LTR"
            elif sFamily_index == 3:
                order = "PLE"
            elif sFamily_index == 4:
                order = "DIRS"
            elif sFamily_index >= 5 and sFamily_index <= 9:
                order = "LINE"
            elif sFamily_index == 10:
                order = "SINE"
            elif sFamily_index >= 11 and sFamily_index <= 19:
                order = "TIR"
            elif sFamily_index == 20:
                order = "MITE"
            elif sFamily_index == 21:
                order = "Helitron"
            elif sFamily_index == 22:
                order = "Maverick"
            elif sFamily_index == 23:
                order = "Crypton"

            if sFamily_index <= 10:
                classTE = "I"
            else:
                classTE = "II"
        elif order_given in orders:
            # Well, superfamily didn't find, but I found the order
            sFamily = "NA"
            order = order_given
            orden_index = orders.index(order)
            if orden_index <= 6:
                classTE = "I"
            else:
                classTE = "II"
        else:
            #print("WARNING: I didn't find the superfamily, neither the order: "+order_given+"/"+superfamily)
            classTE = "Unclassified"
            order = "Unclassified"
            sFamily = "Unclassified"

        lenTR, typeTR = find_TRs(te, outputdir, tools_path)
        profiles, struc_dom = find_profiles(te, outputdir, ref_profiles)
        if do_blast:
            blastx, blasttx = find_blast_hits(te, outputdir, blastx_db, blastn_db)
        else:
            blastx, blasttx = "", ""

        terminals = ''
        if lenTR > 0:
            terminals = 'TermRepeats: ' + typeTR + ': ' + str(lenTR) + ';'

        final_coding_string = ""
        if blasttx != "":
            final_coding_string = blasttx + "; "
        if blastx != '':
            final_coding_string += blastx + ";"
        if profiles != "":
            final_coding_string += profiles

        if blastx == "" and blasttx == "" and profiles == "":
            final_coding_string = "NA"

        class_df = class_df.append(
            {'Seq_name': seq_name, 'length': length, 'strand': '+', 'confused': 'False', 'class': classTE,
             'order': order, 'Wcode': 'NA', 'sFamily': sFamily, 'CI': 0,
             'coding': 'coding=(' + final_coding_string + ')',
             'struct': 'struct=(TElength: ' + str(length) + 'bps; ' + terminals + ' ' + struc_dom + ')',
             'other': 'other=(NA)'}, ignore_index=True)

    return class_df


def count_domains_by_order(profiles, order):
    right_doms = 0
    other_doms = 0
    if order == "LINE":
        right_doms = len(
            [x for x in profiles.split(",") if '_RT_' in x or '_EN_' in x or '_RNaseH_' in x or '_GAG_' in x])
        other_doms = len([x for x in profiles.split(",") if '_AP_' in x or '_INT_' in x or '_ENV_' in x or '_Tase_' in x
                          or '_HEL_' in x or '_Prp' in x or '_ATPase_' in x])
    elif order == "LTR":
        right_doms = len([x for x in profiles.split(",") if
                        '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x])
        other_doms = len([x for x in profiles.split(",") if '_EN_' in x or '_Tase_' in x or '_HEL_' in x or '_Prp' in x or '_ATPase_' in x])

    elif order == "TIR":
        right_doms = len([x for x in profiles.split(",") if '_Tase_' in x])
        other_doms = len([x for x in profiles.split(",") if '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x or '_EN_' in x
                          or '_HEL_' in x or '_Prp' in x or '_ATPase_' in x])

    elif order == "HELITRON":
        right_doms = len([x for x in profiles.split(",") if '_HEL_' in x])
        other_doms = len([x for x in profiles.split(",") if
                          '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x or '_EN_' in x
                          or '_Tase_' in x or '_Prp' in x or '_ATPase_' in x])

    elif order == "MAVERICK":
        right_doms = len([x for x in profiles.split(",") if '_Prp' in x or '_ATPase_' in x])
        other_doms = len([x for x in profiles.split(",") if
                          '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x or '_EN_' in x
                          or '_Tase_' in x or '_HEL_' in x])

    return right_doms, other_doms


def inferring_domains(input_profiles):
    inferred = False
    new_class = 0

    if len(input_profiles) > 0:
        profiles = input_profiles[0]
        # Class I or II ?
        class2_doms = len(
            [x for x in profiles.split(",") if '_Tase_' in x or '_HEL_' in x or '_Prp' in x or '_ATPase_' in x])
        class1_doms = len([x for x in profiles.split(",") if '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or
                          '_RNaseH_' in x or '_ENV_' in x or '_EN_' ])

        if class1_doms > 0 and class2_doms == 0:
            # Retrotransposon !
            inferred = True
            new_class = 32

            # Is it a LTR-RT?
            ltr_count = len([x for x in profiles.split(",") if '_ENV_' in x or '_INT_' in x])
            if ltr_count > 0:
                new_class = 2

        elif class1_doms == 0 and class2_doms > 0:
            # Transposon !
            inferred = True
            new_class = 33

            hel_count = len([x for x in profiles.split(",") if '_HEL_' in x])
            mav_count = len([x for x in profiles.split(",") if '_Prp' in x or '_ATPase_' in x])
            tir_count = len([x for x in profiles.split(",") if '_Tase_' in x])

            if hel_count > 0 and mav_count == 0 and tir_count == 0:
                new_class = 12
            elif hel_count == 0 and mav_count > 0 and tir_count == 0:
                new_class = 13
            elif hel_count == 0 and mav_count == 0 and tir_count > 0:
                new_class = 10

    return inferred, new_class


def decision_tree_rules(struc_table, profiles, i, keep_seqs, minDomLTR, num_copies, orders, minFLNA, kept_seqs_record, ref_tes, automatic):
    status = -1  # 0 = delete, 1 = keep, -1 = Manual Inspection classified module, -3 = unclassified module
    reason = ""
    superFamily = str(struc_table.at[i, "sFamily"])

    #### Class 1 Retrotransposons
    if str(struc_table.at[i, "order"]).upper() == 'LINE' or str(struc_table.at[i, "order"]).upper() == 'SINE':
        if len(profiles) > 0 and len(profiles[0].split(",")) >= 1:
            rigth_doms, other_doms = count_domains_by_order(profiles[0], "LINE")
            if rigth_doms > 0 and other_doms == 0:
                keep_seqs.append(struc_table.at[i, "Seq_name"])
                kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                         str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                if superFamily.upper() == 'R2':
                    orders.append(16)
                elif superFamily.upper() == 'RTE':
                    orders.append(17)
                elif superFamily.upper() == 'JOCKEY':
                    orders.append(18)
                elif superFamily.upper() == 'L1':
                    orders.append(19)
                elif superFamily.upper() == 'I':
                    orders.append(20)
                else:
                    orders.append(8)
                reason = "It's a LINE !!! Kept for having at least one domain !"
                status = 1
            else:
                status = -1
                if automatic != 'F':
                    reason = "Sent to manual inspection because having non-LINE domains !"
                else:
                    reason = "Marked as incomplete TE for having non-LINE domains"
        elif len(profiles) == 0 and "polyAtail" in struc_table.at[i, "struct"] and (
                num_copies[struc_table.at[i, "Seq_name"]][0] >= minFLNA or num_copies[struc_table.at[i, "Seq_name"]][1]
                >= minFLNA):
            keep_seqs.append(struc_table.at[i, "Seq_name"])
            orders.append(9)
            kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                     str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
            reason = "It's a SINE !!! Kept for having polyA tail, no domains and at least " + str(
                minFLNA) + " full-length fragments or copies !"
            status = 1
        else:
            if automatic != 'F':
                reason = "Sent to manual inspection for don't having any domains or PolyA"
            else:
                reason = "Marked as incomplete TE for don't having any domains or PolyA"
            status = -1
    elif str(struc_table.at[i, "order"]).upper() == 'LTR' or str(struc_table.at[i, "order"]).upper() == 'LARD' or \
            str(struc_table.at[i, "order"]).upper() == 'TRIM':
        status_ltr = False
        if "termLTR: " in struc_table.at[i, "struct"]:
            if len(profiles) > 0 and len(profiles[0].split(",")) >= minDomLTR:
                rigth_doms, other_doms = count_domains_by_order(profiles[0], "LTR")
                if rigth_doms > minDomLTR and other_doms == 0:
                    reason = "It's a LTR-RT !!! Kept for having at least " + str(minDomLTR) + " domains and LTRs !"
                    status_ltr = True
                    status = 1
                else:
                    status = -1
                    if automatic != 'F':
                        reason = "Sent to manual inspection for don't having non-LTR-RTs domains !"
                    else:
                        reason = "Marked as incomplete TE for don't having non-LTR-RTs domains !"
            elif len(profiles) > 0 and len(profiles[0].split(",")) > 0 and struc_table.at[i, "order"] == 'LTR':
                rigth_doms, other_doms = count_domains_by_order(profiles[0], "LTR")
                if rigth_doms > 0 and other_doms == 0:
                    reason = "It's a truncated LTR-RT !!! Kept for having " + str(
                        len(profiles[0].split(","))) + " domains and LTRs !"
                    status_ltr = True
                    status = 1
                else:
                    status = -1
                    if automatic != 'F':
                        reason = "Sent to manual inspection for don't having non-LTR-RTs domains !"
                    else:
                        reason = "Marked as incomplete TE for don't having non-LTR-RTs domains !"
            elif len(profiles) == 0 and (num_copies[struc_table.at[i, "Seq_name"]][0] >= minFLNA or
                                         num_copies[struc_table.at[i, "Seq_name"]][1] >= minFLNA):
                keep_seqs.append(struc_table.at[i, "Seq_name"])
                kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                         str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                if int(struc_table.at[i, "length"]) <= 2500:
                    order_nonaut = 6
                    name = "TRIM"
                else:
                    order_nonaut = 7
                    name = "LARD"
                orders.append(order_nonaut)
                reason = "It's a " + name + " !!! Kept for having LTRs, no domains and at least " + str(
                    minFLNA) + " full-length fragments or copies !"
                status = 1
        else:
            # no profiles and no LTRs
            if len(profiles) == 0:
                if automatic != 'F':
                    reason = "Sent to manual inspection for don't having any domains neither LTRs!"
                else:
                    reason = "Marked as incomplete TE for don't having any domains neither LTRs!"
            else:
                if automatic != 'F':
                    reason = "Sent to manual inspection for don't having at least " + str(minDomLTR) + " domains and LTRs!"
                else:
                    reason = "Marked as incomplete TE for don't having at least " + str(minDomLTR) + " domains and LTRs!"
            status = -1

        if status_ltr:
            if superFamily.upper() == 'GYPSY':
                orders.append(4)
            elif superFamily.upper() == 'COPIA':
                orders.append(3)
            elif superFamily.upper().replace("-", "") == 'BELPAO':
                orders.append(5)
            else:
                orders.append(2)
            keep_seqs.append(struc_table.at[i, "Seq_name"])
            kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                     str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
    elif str(struc_table.at[i, "order"]).upper() == 'PLE' or str(struc_table.at[i, "order"]).upper() == 'DIRS':
        # manual inspection in classified module
        if automatic != 'F':
            reason = "Sent to manual inspection"
        else:
            reason = "Marked as incomplete TE"
        status = -1

    #### Class 2 DNA Transposons
    elif str(struc_table.at[i, "order"]).upper() == 'TIR' or str(struc_table.at[i, "order"]).upper() == 'MITE':
        if "termTIR: " in struc_table.at[i, "struct"]:
            if len(profiles) > 0 and len(profiles[0].split(",")) > 0:
                rigth_doms, other_doms = count_domains_by_order(profiles[0], "TIR")
                if rigth_doms > 0 and other_doms == 0:
                    keep_seqs.append(struc_table.at[i, "Seq_name"])
                    kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                             str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                    if superFamily.upper().replace("-", "") == 'TC1MARINER':
                        orders.append(21)
                    elif superFamily.upper() == 'HAT':
                        orders.append(22)
                    elif superFamily.upper() == 'MUTATOR':
                        orders.append(23)
                    elif superFamily.upper() == 'MERLIN':
                        orders.append(24)
                    elif superFamily.upper() == 'TRANSIB':
                        orders.append(25)
                    elif superFamily.upper() == 'P':
                        orders.append(26)
                    elif superFamily.upper() == 'PIGGYBAC':
                        orders.append(27)
                    elif superFamily.upper().replace("-", "") == 'PIFHARBINGER':
                        orders.append(28)
                    elif superFamily.upper() == 'CACTA':
                        orders.append(29)
                    else:
                        orders.append(10)
                    reason = "It's a TIR !!! Kept for having at least one domain !"
                    status = 1
                else:
                    status = -1
                    if automatic != 'F':
                        reason = "Sent to manual inspection because having non-TIR domains !"
                    else:
                        reason = "Marked as incomplete TE because having non-TIR domains !"
            else:
                if num_copies[struc_table.at[i, "Seq_name"]][0] >= minFLNA or num_copies[struc_table.at[i, "Seq_name"]][
                1] >= minFLNA:
                    keep_seqs.append(struc_table.at[i, "Seq_name"])
                    kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                             str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                    orders.append(11)
                    reason = "It's a MITE !!! Kept for having TIRs, no domains and at least " + str(
                        minFLNA) + " full-length fragments or copies !"
                    status = 1
        else:
            # no profiles and no TIRs
            if len(profiles) == 0:
                status = -1
                if automatic != 'F':
                    reason = "Sent to manual inspection for don't having TIRs neither domains !"
                else:
                    reason = "Marked as incomplete TE for don't having TIRs neither domains !"
            else:
                if automatic != 'F':
                    reason = "Sent to manual inspection for don't having TIRs and for having domains !"
                else:
                    reason = "Marked as incomplete TE for don't having TIRs and for having domains !"
                status = -1
    elif str(struc_table.at[i, "order"]).upper() == 'CRYPTON':
        # manual inspection in classified module
        status = -1
        if automatic != 'F':
            reason = "Sent to manual inspection"
        else:
            reason = "Marked as incomplete TE"
    elif str(struc_table.at[i, "order"]).upper() == 'HELITRON':
        if len(profiles) > 0:
            rigth_doms, other_doms = count_domains_by_order(profiles[0], "HELITRON")
            if rigth_doms > 0 and other_doms == 0:
                keep_seqs.append(struc_table.at[i, "Seq_name"])
                kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                         str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                orders.append(12)
                reason = "It's a Helitron !!! Kept for having Helicase !"
                status = 1
            else:
                status = -1
                if automatic != 'F':
                    reason = "Sent to manual inspection for having non-Helitron domains"
                else:
                    reason = "Marked as incomplete TE for having non-Helitron domains"
        elif "helitronExtremities: " in struc_table.at[i, "struct"] :
            keep_seqs.append(struc_table.at[i, "Seq_name"])
            kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                     str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
            orders.append(12)
            reason = "It's a Helitron !!! Kept for having helitron Extremities !"
            status = 1
        else:
            status = -1
            if automatic != 'F':
                reason = "Sent to manual inspection for don't having Helitron Extremities neither helicase !"
            else:
                reason = "Marked as incomplete TE for don't having Helitron Extremities neither helicase !"
    elif str(struc_table.at[i, "order"]).upper() == 'MAVERICK':
        if len(profiles) > 0 and len(profiles[0].split(",")) >= 1:
            rigth_doms, other_doms = count_domains_by_order(profiles[0], "MAVERICK")
            if rigth_doms > 0 and other_doms == 0:
                keep_seqs.append(struc_table.at[i, "Seq_name"])
                kept_seqs_record.append([x for x in SeqIO.parse(ref_tes, "fasta") if
                                         str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])
                orders.append(13)
                reason = "It's a Maverick !!! Kept for having at least one domain !"
                status = 1
            else:
                status = -1
                if automatic != 'F':
                    reason = "Sent to manual inspection for having non-Maverick domains !"
                else:
                    reason = "Marked as incomplete TE for having non-Maverick domains !"
        else:
            status = -1
            if automatic != 'F':
                reason = "Sent to manual inspection for don't having domains !"
            else:
                reason = "Marked as incomplete TE for don't having domains !"
    return status, reason


def new_module1(plots_dir, ref_tes, gff_files, outputdir, pre, te_aid, automatic, minDomLTR, num_copies, minFLNA,
                verbose, repet, ext_nucl, max_nns, num_ite, blastn_db, blastx_db, tools_path, ref_profiles, library_path):
    kept_seqs_record = []
    seqs_to_module3 = []
    non_curated = []
    seqs_manu = 0
    seqs_module3 = 0
    keep_auto = 0
    ele_number = 0

    # obtain seqIDs of consensi with at least one full length fragment
    seqID_list = [str(te.id).split("#")[0] for te in SeqIO.parse(ref_tes, "fasta")]

    seqs_to_mi = []
    orders_seqs_mi = []
    ref_tes_bee = outputdir + "/classifiedModule/extended_cons.fa"
    if automatic != 'M':
        # Step 1 Extend all the consensus
        start_time = time.time()
        run_BEE_parallel(genome, ref_tes, ext_nucl, num_ite, outputdir + "/classifiedModule/", max_nns, cores)
        end_time = time.time()
        if verbose:
            print("MESSAGE: The sequences were extended successfully [" + str(end_time - start_time) + " seconds]")

        start_time = time.time()
        build_class_table_parallel(ref_tes_bee, cores, outputdir+'/classifiedModule/', blastn_db, blastx_db,
                                   tools_path, ref_profiles, False)
        struc_table = pd.read_csv(outputdir + "/classifiedModule/denovoLibTEs_PC.classif", sep='\t')
        end_time = time.time()
        if verbose:
            print("MESSAGE: TE Feature table was created [" + str(end_time - start_time) + " seconds]")

        # Step2 BLASTn with ref library
        start_time = time.time()
        keep_seqs, orders = run_blast(library_path, ref_tes_bee, cores, 80, 80)
        keep_auto += len(keep_seqs)
        for seq_id in keep_seqs:
            kept_seqs_record.append([x for x in SeqIO.parse(ref_tes_bee, "fasta") if x.id.split("#")[0] == seq_id][0])
        end_time = time.time()
        if verbose:
            print("MESSAGE: Total TEs kept because of hits with library: " + str(len(keep_seqs)) + " [" + str(end_time - start_time) + " seconds]")

        # Step 3 Structural Check
        totalTEs = struc_table.shape[0]
        for i in range(struc_table.shape[0]):
            status = -1
            if struc_table.at[i, "Seq_name"] not in keep_seqs and struc_table.at[i, "Seq_name"] in seqID_list:
                codings = struc_table.at[i, "coding"].replace("coding=(", " ").replace(")", "")
                codigs = codings.split(";")

                # Structural checks
                if str(struc_table.at[i, "class"]).upper() in ['UNCLASSIFIED', 'NA', 'NAN'] or str(struc_table.at[i, "order"]).upper() in ['UNCLASSIFIED', 'NA', 'NAN']:
                    seqs_to_module3.append(struc_table.at[i, "Seq_name"])
                    status = -3
                else:
                    profiles = [cod for cod in codigs if "profiles:" in cod]
                    status, reason = decision_tree_rules(struc_table, profiles, i, keep_seqs, minDomLTR, num_copies, orders, minFLNA, kept_seqs_record, ref_tes_bee, automatic)

                if status == 1:
                    keep_auto += 1
                elif status != -3:
                    all_class = [x.replace("–", "").upper() for x in dicc_orders.values()]
                    seqs_to_mi.append(struc_table.at[i, "Seq_name"])
                    if str(struc_table.at[i, "sFamily"]).upper() != "NAN":
                        class_mod2 = all_class.index(str(struc_table.at[i, "sFamily"]).upper()) + 2
                    else:
                        class_mod2 = all_class.index(str(struc_table.at[i, "order"]).upper()) + 2
                    orders_seqs_mi.append(class_mod2)

                if verbose:
                    print("[" + str(ele_number + 1) + "/" + str(totalTEs) + "] Seq Name: " + str(
                        struc_table.at[i, "Seq_name"]))
                    print("Length: " + str(struc_table.at[i, "length"]))
                    print("Class: CLASS " + str(struc_table.at[i, "class"]) + " / " + str(
                        struc_table.at[i, "order"]) + " / " + str(struc_table.at[i, "sFamily"]))
                    print("Copies: FLF = " + str(num_copies[struc_table.at[i, "Seq_name"]][0]) + " / FLC = " + str(
                        num_copies[struc_table.at[i, "Seq_name"]][1]))
                    print("Coding: ")
                    for cod in codigs:
                        print(cod)
                    print(struc_table.at[i, "struct"])
                    print(struc_table.at[i, "other"])
                    if status == 1:
                        print(reason)
                        print("---------------------------------------------------------------\n")
                    elif status == -1:
                        print(reason)
                        print("---------------------------------------------------------------\n")
                    elif status == -3:
                        print("Sent to unclassified module due to it's unclassified !")
                        print("---------------------------------------------------------------\n")
                        seqs_module3 += 1
            else:
                if verbose:
                    codings = struc_table.at[i, "coding"].replace("coding=(", " ").replace(")", "")
                    codigs = codings.split(";")
                    print("[" + str(ele_number + 1) + "/" + str(totalTEs) + "] Seq Name: " + str(
                        struc_table.at[i, "Seq_name"]))
                    print("Length: " + str(struc_table.at[i, "length"]))
                    print("Class: CLASS " + str(struc_table.at[i, "class"]) + " / " + str(
                        struc_table.at[i, "order"]) + " / " + str(struc_table.at[i, "sFamily"]))
                    print("Copies: FLF = " + str(num_copies[struc_table.at[i, "Seq_name"]][0]) + " / FLC = " + str(
                        num_copies[struc_table.at[i, "Seq_name"]][1]))
                    print("Coding: ")
                    for cod in codigs:
                        print(cod)
                    print(struc_table.at[i, "struct"])
                    print(struc_table.at[i, "other"])
                    print("TE Kept because it has homology with our databases !")
                    print("---------------------------------------------------------------\n")
            ele_number += 1
    else:
        ref_tes_bee = ref_tes
        start_time = time.time()
        build_class_table_parallel(ref_tes_bee, cores, outputdir + '/classifiedModule/', blastn_db, blastx_db,
                                   tools_path, ref_profiles, False)
        struc_table = pd.read_csv(outputdir + "/classifiedModule/denovoLibTEs_PC.classif", sep='\t')
        end_time = time.time()
        if verbose:
            print("MESSAGE: TE Feature table was created [" + str(end_time - start_time) + " seconds]")
        keep_seqs, orders = [], []
    # those seqs that after finished the iterations, don't have structural features yet. Manual Inspection if S mode
    # or marked as incomplete elements if A mode
    if te_aid == 'Y':
        run_te_aid_parallel(tools_path + "/TE-Aid-master/", genome, ref_tes_bee, outputdir + "/classifiedModule", cores)

    ele_number = 0
    for i in range(struc_table.shape[0]):
        if (automatic == 'M' or struc_table.at[i, "Seq_name"] in seqs_to_mi) and struc_table.at[i, "Seq_name"] in seqID_list:
            if automatic != 'F':
                codings = struc_table.at[i, "coding"].replace("coding=(", " ").replace(")", "")
                codigs = codings.split(";")
                print("[" + str(ele_number + 1) + "/" + str(struc_table.shape[0]) + "] Seq Name: " + str(struc_table.at[i, "Seq_name"]))
                print("Length: " + str(struc_table.at[i, "length"]))
                print("Class: CLASS " + str(struc_table.at[i, "class"]) + " / " + str(
                    struc_table.at[i, "order"]) + " / " + str(struc_table.at[i, "sFamily"]))
                print("Copies: FLF = " + str(num_copies[struc_table.at[i, "Seq_name"]][0]) + " / FLC = " + str(
                    num_copies[struc_table.at[i, "Seq_name"]][1]))
                print("Coding: ")
                for cod in codigs:
                    print(cod)
                print(struc_table.at[i, "struct"])
                print(struc_table.at[i, "other"])

                if te_aid == 'N':

                    fig = plt.figure(figsize=(20, 10))
                    # setting values to rows and column variables
                    rows = 1
                    columns = 3

                    # reading images
                    Image1 = cv2.imread(plots_dir + '/' + pre + '_' + struc_table.at[i, "Seq_name"] + '.png')

                    # Adds a subplot at the 1st position
                    fig.add_subplot(rows, columns, 1)

                    # showing image
                    plt.imshow(Image1)
                    plt.axis('off')
                    plt.title("REPET")

                    # obtaining the positions of Blastn, tBlastn y blastx from GFF:
                    posBltn, postBlasx, posBlastx, posTR = processing_gff(struc_table.at[i, "Seq_name"], gff_files, pre)

                    # ploting gff files
                    ax = fig.add_subplot(rows, columns, 2)
                    ax.broken_barh([(0, int(struc_table.at[i, "length"]))], (0, 2), facecolors='black')
                    ax.broken_barh(posBltn, (3, 2), facecolors='blue', alpha=0.2)
                    ax.broken_barh(postBlasx, (6, 2), facecolors='orange', alpha=0.2)
                    ax.broken_barh(posBlastx, (9, 2), facecolors='red', alpha=0.2)
                    ax.broken_barh(posTR, (12, 2), facecolors='green', alpha=0.2)
                    ax.set_ylim(0, 14)
                    ax.set_xlim(0, int(struc_table.at[i, "length"]))
                    ax.set_yticks([1, 4, 7, 10, 13])
                    ax.set_yticklabels(['TE seq', 'BLASTn', 'tBLASTx', 'BLASTx', 'TRs'])
                    ax.set_xlabel('TE consensus')
                    ax.set_title('REPET GFF files')

                    if os.path.exists(outputdir + '/classifiedModule/MSA_plots/' + struc_table.at[i, "Seq_name"] + '.copies.cialign_output.png'):
                        Image2 = cv2.imread(outputdir + '/classifiedModule/MSA_plots/' + struc_table.at[i, "Seq_name"] + '.copies.cialign_output.png')
                        fig.add_subplot(rows, columns, 3)
                        # showing image
                        plt.imshow(Image2)
                        plt.axis('off')
                        plt.title("MSA")
                    else:
                        ax = fig.add_subplot(rows, columns, 3)
                        ax.set_title('MSA')
                        ax.axis([0, 10, 0, 10])
                        ax.text(4, 5, 'Not enough TE copies', style='italic', fontsize=14, fontweight='bold')

                    plt.tight_layout()
                    plt.draw()
                    plt.pause(1)

                else:
                    # create figure
                    fig = plt.figure(figsize=(20, 10))

                    # setting values to rows and column variables
                    rows = 1
                    if repet:
                        columns = 3
                    else:
                        columns = 1
                    if automatic != 'M':
                        columns += 1
                    try:
                        pos_te_aid = 1
                        if repet:
                            # reading images
                            Image1 = cv2.imread(plots_dir + '/' + pre + '_' + struc_table.at[i, "Seq_name"] + '.png')

                            # Adds a subplot at the 1st position
                            fig.add_subplot(rows, columns, 1)

                            # showing image
                            plt.imshow(Image1)

                            plt.axis('off')
                            plt.title("REPET")

                            # obtaining the positions of Blastn, tBlastn y blastx from GFF:
                            posBltn, postBlasx, posBlastx, posTR = processing_gff(struc_table.at[i, "Seq_name"],
                                                                                  gff_files, pre)
                            # ploting gff files
                            ax = fig.add_subplot(rows, columns, 2)
                            ax.broken_barh([(0, int(struc_table.at[i, "length"]))], (0, 2), facecolors='black')
                            ax.broken_barh(posBltn, (3, 2), facecolors='blue', alpha=0.2)
                            ax.broken_barh(postBlasx, (6, 2), facecolors='orange', alpha=0.2)
                            ax.broken_barh(posBlastx, (9, 2), facecolors='red', alpha=0.2)
                            ax.broken_barh(posTR, (12, 2), facecolors='green', alpha=0.2)
                            ax.set_ylim(0, 14)
                            ax.set_xlim(0, int(struc_table.at[i, "length"]))
                            ax.set_yticks([1, 4, 7, 10, 13])
                            ax.set_yticklabels(['TE seq', 'BLASTn', 'tBLASTx', 'BLASTx', 'TRs'])
                            ax.set_xlabel('TE consensus')
                            ax.set_title('REPET GFF files')
                            pos_te_aid = 3

                        try:
                            pages = convert_from_path(
                                outputdir + '/classifiedModule/te_aid/' + struc_table.at[i, "Seq_name"] + '.fa.c2g.pdf')

                            # Saving pages in jpeg format
                            for page in pages:
                                page.save(
                                    outputdir + '/classifiedModule/te_aid/' + struc_table.at[i, "Seq_name"] + '.fa.c2g.jpeg',
                                    'JPEG')

                            Image2 = cv2.imread(
                                outputdir + '/classifiedModule/te_aid/' + struc_table.at[i, "Seq_name"] + '.fa.c2g.jpeg')

                            # Adds a subplot at the 2nd position

                            fig.add_subplot(rows, columns, pos_te_aid)

                            # showing image
                            plt.imshow(Image2)
                            plt.axis('off')
                            plt.title("TE-aid")
                            if automatic != 'M':
                                if os.path.exists(outputdir + '/classifiedModule/MSA_plots/' + struc_table.at[
                                    i, "Seq_name"] + '.copies.cialign_output.png'):
                                    Image2 = cv2.imread(
                                        outputdir + '/classifiedModule/MSA_plots/' + struc_table.at[i, "Seq_name"] + '.copies.cialign_output.png')
                                    fig.add_subplot(rows, columns, pos_te_aid + 1)
                                    # showing image
                                    plt.imshow(Image2)
                                    plt.axis('off')
                                    plt.title("MSA")
                                else:
                                    ax = fig.add_subplot(rows, columns, pos_te_aid + 1)
                                    ax.set_title('MSA')
                                    plt.axis('off')
                                    ax.axis([0, 10, 0, 10])
                                    ax.text(4, 5, 'Not enough TE copies', style='italic', fontsize=14, fontweight='bold')

                            plt.tight_layout()
                            plt.draw()
                            plt.pause(1)
                        except PDFPageCountError:
                            print("FATAL ERROR: the " + outputdir + "/classifiedModule/te_aid/" + struc_table.at[
                                i, "Seq_name"] + ".fa.c2g.pdf is empty. Please check the file.")

                        delete_files(outputdir + '/classifiedModule/te_aid/' + struc_table.at[i, "Seq_name"] + '.fa.c2g.jpeg')
                    except FileNotFoundError:
                        print("WARNING: Element " + struc_table.at[i, "Seq_name"] + "couldn't be processed")

                try:
                    keep = int(
                        input(
                            "Keep the sequence?[1: remove, 2: LTR, 3: LTR-Copia, 4: LTR-Gypsy, 5: LTR-Bel_pao,"
                            " 6: TRIM, 7: LARD, 8: LINE, 9: SINE, 10: TIR, 11: MITE, 12: Helitron, 13: Maverick"
                            ", 14: PLE, 15: DIRS,  16: R2, 17: RTE, 18: Jockey, 19: L1, 20: I, 21: Tc1–Mariner,"
                            " 22: hAT, 23: Mutator, 24: Merlin, 25: Transib, 26: P, 27: PiggyBac, 28: PIF–Harbinger, "
                            "29: CACTA, 30: CRYPTON, 31: Unclassified, 32: Retrotransposon, 33: Transposon] "))
                except ValueError:
                    keep = -1

                while keep < 1 or keep > 34:
                    try:
                        print('You must indicate a number between 1 and 33.')
                        keep = int(
                            input(
                                "Keep the sequence?[1: remove, 2: LTR, 3: LTR-Copia, 4: LTR-Gypsy, 5: LTR-Bel_pao,"
                                " 6: TRIM, 7: LARD, 8: LINE, 9: SINE, 10: TIR, 11: MITE, 12: Helitron, 13: Maverick"
                                ", 14: PLE, 15: DIRS,  16: R2, 17: RTE, 18: Jockey, 19: L1, 20: I, 21: Tc1–Mariner,"
                                " 22: hAT, 23: Mutator, 24: Merlin, 25: Transib, 26: P, 27: PiggyBac, 28: PIF–Harbinger,"
                                " 29: CACTA, 30: CRYPTON, 31: Unclassified, 32: Retrotransposon, 33: Transposon] "))
                    except ValueError:
                        keep = -1
                if keep == 31:
                    seqs_to_module3.append(struc_table.at[i, "Seq_name"])
                elif keep > 1:
                    keep_seqs.append(struc_table.at[i, "Seq_name"])
                    orders.append(keep)
                    kept_seqs_record.append([x for x in SeqIO.parse(ref_tes_bee, "fasta") if
                                             str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]][0])

                # increase the number of sequences manually analyzed
                seqs_manu += 1
                print("---------------------------------------------------------------\n")
            else:
                te = [x for x in SeqIO.parse(ref_tes_bee, "fasta") if
                      x.id.split("#")[0] == struc_table.at[i, "Seq_name"]][0]
                superFamilyLib = ""
                if str(struc_table.at[i, "sFamily"]).upper() in [x for x in dicc_orders.values()]:
                    superFamilyLib = "/" + str(struc_table.at[i, "sFamily"]).upper()
                te.id = str(struc_table.at[i, "Seq_name"]) + "_inc#" + str(struc_table.at[i, "order"]) + superFamilyLib
                te.description = ""
                non_curated.append(te)
        ele_number += 1

    print("")
    print("Total stats :")
    print("Sequences kept: " + str(len(keep_seqs) + len(non_curated)))
    print("-------------------------------------------")
    print("TEs detected as complete: " + str(len(kept_seqs_record)))
    print("TEs detected as incomplete: " + str(len(non_curated)))
    print("TEs manually analyzed : " + str(seqs_manu))
    print("TEs sent to Unclassified module : " + str(len(seqs_to_module3)))
    print("-------------------------------------------")
    print("")

    for index in range(len(kept_seqs_record)):
        classification = dicc_orders[orders[index]]
        if orders[index] >= 3 and orders[index] <= 5:
            classification = "LTR/" + classification
        elif orders[index] >= 16 and orders[index] <= 20:
            classification = "LINE/" + classification
        elif orders[index] >= 21 and orders[index] <= 29:
            classification = "TIR/" + classification
        new_name = keep_seqs[index] + "#" + classification
        kept_seqs_record[index].id = new_name
        kept_seqs_record[index].description = ""

    write_sequences_file(kept_seqs_record, outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa")
    write_sequences_file(non_curated, outputdir + "/classifiedModule/kept_seqs_classified_module_non_curated.fa")
    create_fasta_files(ref_tes_bee, seqs_to_module3, outputdir + "/classifiedModule/input_to_unclassified_module_seqs.fa", outputdir, [])

    delete_files(outputdir + "/classifiedModule/extended_cons.fa")
    delete_files(outputdir + "/classifiedModule/putative_TEs.fa")
    delete_files(outputdir + "/classifiedModule/new_user_lib.fa")


def run_BEE_parallel(genome, ref_tes, exe_nucl, num_ite, outputdir, max_nns, cores):
    create_output_folders(outputdir)

    print('MESSAGE: Starting with BEE step ...')

    # Check if the genome is already formatted to be a BLAST db
    if not os.path.exists(genome + ".n*"):  # Run makeblastdb
        output = subprocess.run(
            ['makeblastdb', '-in', genome, '-dbtype', 'nucl'], stdout=subprocess.PIPE, text=True)

    if os.path.getsize(genome) == 0:
        print("FATAL ERROR: "+genome+" is empty. Please check the file and re-run the software")
        sys.exit(0)
    elif os.path.getsize(ref_tes) == 0:
        print("FATAL ERROR: "+ref_tes+" is empty. Please check the file and re-run the software")
        sys.exit(0)

    tes = [te for te in SeqIO.parse(ref_tes, "fasta")]
    n = len(tes)
    seqs_per_procs = int(n / cores)
    remain = n % cores
    ini_per_thread = []
    end_per_thread = []
    for p in range(cores):
        if p < remain:
            init = p * (seqs_per_procs + 1)
            end = n if init + seqs_per_procs + 1 > n else init + seqs_per_procs + 1
        else:
            init = p * seqs_per_procs + remain
            end = n if init + seqs_per_procs > n else init + seqs_per_procs
        ini_per_thread.append(init)
        end_per_thread.append(end)

    # Run in parallel the checking
    pool = multiprocessing.Pool(processes=cores)
    localresults = [pool.apply_async(BEE,
                                     args=[genome, tes[ini_per_thread[x]:end_per_thread[x]], exe_nucl, num_ite,
                                           outputdir, max_nns, False, x]) for x in range(cores)]

    localChecks = [p.get() for p in localresults]
    final_seqs = []
    for i in range(len(localChecks)):
        for te in SeqIO.parse(outputdir + "/extended_cons.fa_" + str(i), "fasta"):
            final_seqs.append(te)
        delete_files(outputdir + "/extended_cons.fa_" + str(i))
    pool.close()

    write_sequences_file(final_seqs, outputdir + "/extended_cons.fa")

    create_output_folders(outputdir + "/MSA_plots/")
    pattern = "*.copies.cialign_output.png"
    files = glob.glob(outputdir + "/" + pattern)
    for file in files:
        file_name = os.path.basename(file)
        shutil.move(file, outputdir + "/MSA_plots/" + file_name)

    create_output_folders(outputdir + "/MSA_seeds/")
    pattern = "*.copies.cialign_with_consensus.fasta"
    files = glob.glob(outputdir + "/" + pattern)
    for file in files:
        file_name = os.path.basename(file)
        shutil.move(file, outputdir + "/MSA_seeds/" + file_name)


def BEE(genome, ref_tes, exe_nucl, num_ite, outputdir, max_nns, verbose, id_thread):
    final_results = []
    num_seq = 0

    for te in ref_tes:
        seq_name = str(te.id).split("#")[0]
        # Write the TE seq for the first iteration
        write_sequences_file([te], outputdir + "/" + str(seq_name) + ".copies.cialign_consensus.fasta")
        if verbose:
            print("   [" + str(num_seq) + "] Seq ID = " + str(te.id))

        for it in range(num_ite):
            # First step: BLASTn
            output = subprocess.run(
                ['blastn', '-query', outputdir + "/" + str(seq_name) + ".copies.cialign_consensus.fasta", '-db', genome, '-out',
                 outputdir + "/" + str(seq_name) + ".blast", '-num_threads', "1", "-outfmt", "6",
                 "-qcov_hsp_perc", "50", "-perc_identity", "80"], stdout=subprocess.PIPE, text=True)

            # Second step EXTRACT and EXTEND
            blastresult = pd.read_table(outputdir + "/" + str(seq_name) + ".blast", sep='\t',
                                        names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                                               'qend', 'sstart', 'send', 'evalue', 'bitscore'])
            blastresult = blastresult.sort_values(by=['bitscore', 'pident', 'length'], ascending=[False, False, False])

            result_file = open(outputdir + "/" + str(seq_name) + ".copies.fa", "w")
            original_cons = [consen for consen in SeqIO.parse(outputdir + "/" + str(seq_name) + ".copies.cialign_consensus.fasta", "fasta")][0]

            if verbose:
                print("        [" + str(it) + "] Length = " + str(len(original_cons.seq)))

            hit = 0
            while hit < blastresult.shape[0] and hit < 20:
                reverse = False
                subject_seq = blastresult.at[hit, 'sseqid']
                ini_hit = int(blastresult.at[hit, 'sstart'])
                end_hit = int(blastresult.at[hit, 'send'])
                scaff_seq = [x.seq for x in SeqIO.parse(genome, 'fasta') if str(x.id).split(" ")[0] == subject_seq]
                if len(scaff_seq) > 0:
                    scaff_seq = str(scaff_seq[0])
                    # for reversed hits
                    if end_hit < ini_hit:
                        tem = ini_hit
                        ini_hit = end_hit
                        end_hit = tem
                        reverse = True

                    if ini_hit - exe_nucl < 0:
                        ini_hit = 0
                    else:
                        ini_hit -= exe_nucl
                    if end_hit + exe_nucl > len(scaff_seq):
                        end_hit = len(scaff_seq)
                    else:
                        end_hit += exe_nucl

                    if reverse is False:
                        result_file.write(
                            ">copy_" + str(subject_seq) + "_" + str(ini_hit) + "_" + str(end_hit) + "\n" + scaff_seq[
                                                                                                           ini_hit:end_hit] + "\n")
                    else:
                        result_file.write(
                            ">copy_" + str(subject_seq) + "_" + str(ini_hit) + "_" + str(end_hit) + "\n" + str(
                                Seq(scaff_seq[ini_hit:end_hit]).reverse_complement()) + "\n")
                else:
                    print('WARNING: A subject id (sseqid in blast output file) was not found in the genome, please check:')
                    print(subject_seq)
                hit += 1
            result_file.close()

            if verbose:
                print("        [" + str(it) + "] BLAST hits = " + str(blastresult.shape[0]))

            if hit > 2:  # if there is at least two hit of the consensus and the genome
                # Third step: Multiple alignment, CIAlign, and Consensus creation
                output = subprocess.run(['mafft', '--auto', '--quiet', '--thread', '1', '--adjustdirection',
                                         outputdir + "/" + str(seq_name) + '.copies.fa'],
                                        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)

                mafft_comm_output = str(output.stdout)
                mafft_output = open(outputdir + "/" + str(seq_name) + ".copies.mafft", "w")
                mafft_output.write(mafft_comm_output)
                mafft_output.close()

                if len(mafft_comm_output) == 0:
                    if verbose:
                        print("         No mafft align there!!!")
                else:
                    output = subprocess.run(
                        ['CIAlign', '--infile', outputdir + "/" + str(seq_name) + ".copies.mafft", '--out', outputdir + "/" + str(seq_name) + ".copies.cialign",
                         '--plot_output', '--remove_short', '--make_consensus', '--remove_divergent', '--remove_divergent_minperc', '0.3', '--crop_ends',
                         '--remove_insertions', '--insertion_max_size', '200', '--silent'],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                    delete_files(outputdir + "/" + str(seq_name) + ".copies.cialign_log.txt")
                    if os.path.exists(outputdir + "/" + str(seq_name) + ".copies.cialign_removed.txt"):
                        delete_files(outputdir + "/" + str(seq_name) + ".copies.cialign_removed.txt")
                        delete_files(outputdir + "/" + str(seq_name) + ".copies.cialign_cleaned.fasta")

                delete_files(outputdir + "/" + str(seq_name) + ".copies.mafft")

        final_seq = [con for con in SeqIO.parse(outputdir + "/" + str(seq_name)+".copies.cialign_consensus.fasta", 'fasta')][0]

        # only keep those consensuses with maximum max_nns number of Ns in the final result
        perc_nns = (str(final_seq.seq).lower().count('n') * 100) / len(final_seq.seq)

        if verbose:
            print("      Final Length = " + str(len(final_seq.seq)))
            print("      Number of Ns = " + str(str(final_seq.seq).lower().count('n')) + " (" + str(perc_nns) + " %)")

        if perc_nns <= max_nns:
            final_seq.id = te.id
            final_seq.description = ""
            final_results.append(final_seq)
        else:
            final_results.append(te)

        num_seq += 1

        delete_files(outputdir + "/" + str(seq_name) + ".copies.cialign_consensus.fasta")
        if num_ite > 0:
            delete_files(outputdir + "/" + str(seq_name) + ".copies.fa")
            delete_files(outputdir + "/" + str(seq_name) + ".blast")

    write_sequences_file(final_results, outputdir + "/extended_cons.fa_" + str(id_thread))
    return num_seq


def filter_bad_candidates(new_ref_tes, perc_ssr, outputdir, tools_path, busco_library, RNAs_library, cores):

    if os.path.getsize(new_ref_tes) == 0:
        print("FATAL ERROR: " + new_ref_tes + " is empty. Please check the file and re-run the software")
        sys.exit(0)
    elif os.path.getsize(busco_library) == 0:
        print("FATAL ERROR: "+busco_library+" is empty. Please check the file and re-run the software")
        sys.exit(0)

    deleted_seqs = []
    kept_seqs = []

    current_path = os.getcwd()
    os.chdir(outputdir)
    try:
        output = subprocess.run(
            [tools_path + '/trf409.linux64', new_ref_tes, '2', '3', '5', '80', '10', '20', '15', '-h', '-d'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except Exception as exp:
        print("FATAL ERROR: I couldn't execute properly the TRF program. Please check the error: "+exp.args)

    os.chdir(current_path)
    # Search for single repeats in TEs
    dicc_sr = {}
    with open(outputdir + "/"+os.path.basename(new_ref_tes)+".2.3.5.80.10.20.15.dat", 'r') as inFile:
        nbInSeq = 0
        for line in inFile:
            row = line.split(" ")
            if len(row) > 1 and "Sequence:" in row[0]:
                nbInSeq += 1
                seqName = row[1][:-1]
            if len(row) >= 14 and not "Sequence:" in row[0]:
                start = row[0]
                end = row[1]
                if seqName in dicc_sr.keys():
                    dicc_sr[seqName] += int(end) - int(start)
                else:
                    dicc_sr[seqName] = int(end) - int(start)

    # Search for matches with References/BUSCO genes
    if not os.path.exists(busco_library+".h3m"):
        output = subprocess.run(['hmmpress', busco_library],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    output = subprocess.run(
        ['hmmscan', '--tblout', outputdir + "tes_vs_genes.hmm", '-E', '10', '--noali', '--cpu',
         str(cores), busco_library, new_ref_tes],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    command = 'grep -v "^#" ' + outputdir + 'tes_vs_genes.hmm | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}\' > ' + outputdir + "tes_vs_genes.hmm_formatted"
    process = subprocess.Popen(command, shell=True)
    process.communicate()
    hmm_results = pd.read_table(outputdir + "tes_vs_genes.hmm_formatted", sep='\t',
                                names=['target', 'acc', 'query', 'acc_query', 'evalue', 'score', 'A', 'B', 'C',
                                       'D'])

    tes_with_matches = []
    for x in range(hmm_results.shape[0]):
        if hmm_results.loc[x, "query"] not in tes_with_matches:
            tes_with_matches.append(hmm_results.loc[x, "query"])

    # Search for matches with tRNAs or rRNAs
    if not os.path.exists(RNAs_library + ".h3m"):
        output = subprocess.run(['hmmpress', RNAs_library],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    output = subprocess.run(
        ['hmmscan', '--tblout', outputdir + "/tes_vs_rnas.hmm", '-E', '10', '--noali', '--cpu',
         str(cores), RNAs_library, new_ref_tes],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    command = 'grep -v "^#" ' + outputdir + '/tes_vs_rnas.hmm | awk \'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}\' > ' + outputdir + '/tes_vs_rnas.hmm_formatted'
    process = subprocess.Popen(command, shell=True)
    process.communicate()
    hmm_results = pd.read_table(outputdir + '/tes_vs_rnas.hmm_formatted', sep='\t',
                                names=['target', 'acc', 'query', 'acc_query', 'evalue', 'score', 'A', 'B', 'C', 'D'])

    tes_with_rnas = []
    for x in range(hmm_results.shape[0]):
        if hmm_results.loc[x, "query"] not in tes_with_rnas:
            tes_with_rnas.append(hmm_results.loc[x, "query"])

    # Remove those TEs with SSR > 60% or those that match with genes
    for te in SeqIO.parse(new_ref_tes, "fasta"):
        if te.id not in tes_with_matches and te.id not in tes_with_rnas:  # It hasn't a match with Reference/BUSCO genes neither with RNAs
            if te.id in dicc_sr.keys():
                te_len = len(str(te.seq))
                lenSR = dicc_sr[te.id]
                if (lenSR * 100) / te_len < perc_ssr:  # It hasn't less than 60% of SSR in its sequence
                    kept_seqs.append(te)
                else:
                    deleted_seqs.append(te)
            else:  # It hasn't any SSR
                kept_seqs.append(te)
        else:
            deleted_seqs.append(te)

    write_sequences_file(kept_seqs, outputdir + "/putative_TEs.fa")
    delete_files(outputdir + "/" + os.path.basename(new_ref_tes) + ".2.3.5.80.10.20.15.dat")
    delete_files(outputdir + "/tes_vs_genes.hmm")
    delete_files(outputdir + "/tes_vs_genes.hmm_formatted")
    delete_files(outputdir + "/tes_vs_rnas.hmm")
    delete_files(outputdir + "/tes_vs_rnas.hmm_formatted")
    delete_files(busco_library + ".h3f")
    delete_files(busco_library + ".h3i")
    delete_files(busco_library + ".h3m")
    delete_files(busco_library + ".h3p")

    return outputdir + "/putative_TEs.fa", deleted_seqs


def module3(ref_tes, library_path, cores, outputdir, perc_ident, perc_cover, internal_library):
    # BLASTn with already curated TEs
    start_time = time.time()
    keep_seqs_records = []
    keep_seqs, orders = run_blast(library_path, ref_tes, cores, perc_ident, perc_cover)

    # BLASTn with internal reference TE database
    keep_seqs_internal, orders_internal = run_blast(internal_library, ref_tes, cores, perc_ident, perc_cover)

    for te_index in range(len(keep_seqs_internal)):
        if keep_seqs_internal[te_index] not in keep_seqs:
            keep_seqs.append(keep_seqs_internal[te_index])
            orders.append(orders_internal[te_index])

    for te_index in range(len(keep_seqs)):
        te = [x for x in SeqIO.parse(ref_tes, "fasta") if x.id.split("#")[0] == keep_seqs[te_index]][0]
        te.id = te.id.split("#")[0] + "#" + dicc_orders[orders[te_index]]
        te.description = ""
        keep_seqs_records.append(te)

    end_time = time.time()
    print("MESSAGE: BLASTn successfully run [" + str(end_time - start_time) + " seconds]")

    start_time = time.time()
    struc_table = pd.read_csv(outputdir + "/unclassifiedModule/denovoLibTEs_PC.classif", sep='\t')
    num_infered_tes = 0
    num_structural_tes = 0
    for i in range(struc_table.shape[0]):
        if struc_table.at[i, "Seq_name"] not in keep_seqs:
            codings = struc_table.at[i, "coding"].replace("coding=(", " ").replace(")", "")
            codigs = codings.split(";")
            profiles = [cod for cod in codigs if "profiles:" in cod]
            inferred, new_class = inferring_domains(profiles)
            if inferred:
                te_selected = [x for x in SeqIO.parse(ref_tes, "fasta") if x.id.split("#")[0] == struc_table.at[i, "Seq_name"]]
                te_selected[0].id = te_selected[0].id.split("#")[0] + "#" + dicc_orders[new_class]
                te_selected[0].description = ""
                keep_seqs_records.append(te_selected[0])
                num_infered_tes += 1
            else:
                if "termLTR: " in struc_table.at[i, "struct"]:
                    te_selected = [x for x in SeqIO.parse(ref_tes, "fasta") if x.id.split("#")[0] == struc_table.at[i, "Seq_name"]]
                    te_selected[0].id = te_selected[0].id.split("#")[0] + "_unconfirmed#LTR"
                    te_selected[0].description = ""
                    keep_seqs_records.append(te_selected[0])
                    num_structural_tes += 1
                elif "termTIR: " in struc_table.at[i, "struct"]:
                    te_selected = [x for x in SeqIO.parse(ref_tes, "fasta") if x.id.split("#")[0] == struc_table.at[i, "Seq_name"]]
                    te_selected[0].id = te_selected[0].id.split("#")[0] + "_unconfirmed#TIR"
                    te_selected[0].description = ""
                    keep_seqs_records.append(te_selected[0])
                    num_structural_tes += 1
    end_time = time.time()
    print("MESSAGE: Order inferring from structural features done [" + str(end_time - start_time) + " seconds]")

    if len(keep_seqs_records) > 0:
        write_sequences_file(keep_seqs_records, outputdir + "/unclassifiedModule/kept_seqs_unclassified_module.fa")
        print("")
        print("Total stats :")
        print("Unclassified elements recovered: " + str(len(keep_seqs_records)))
        print("-------------------------------------------")
        print("TEs recovered by homology: " + str(len(keep_seqs)))
        print("TEs recovered by domain inferring: " + str(num_infered_tes))
        print("TEs recovered by structural features: " + str(num_structural_tes))
        print("-------------------------------------------")
        print("")
    else:
        print("WARNING: unclassified module couldn't find any TE !")


def run_te_aid_parallel(te_aid_path, genome, ref_tes, outputdir, cores):
    if not os.path.exists(outputdir + "/te_aid"):

        if not os.path.exists(genome + ".nhr"):
            output = subprocess.run(
                ['makeblastdb', '-in', genome, '-dbtype', 'nucl'], stdout=subprocess.PIPE, text=True)

        tes = [te for te in SeqIO.parse(ref_tes, "fasta")]
        n = len(tes)
        seqs_per_procs = int(n / cores)
        remain = n % cores
        ini_per_thread = []
        end_per_thread = []
        for p in range(cores):
            if p < remain:
                init = p * (seqs_per_procs + 1)
                end = n if init + seqs_per_procs + 1 > n else init + seqs_per_procs + 1
            else:
                init = p * seqs_per_procs + remain
                end = n if init + seqs_per_procs > n else init + seqs_per_procs
            ini_per_thread.append(init)
            end_per_thread.append(end)

        # Run in parallel the checking
        pool = multiprocessing.Pool(processes=cores)
        localresults = [pool.apply_async(run_te_aid,
                                         args=[te_aid_path, genome, outputdir + "/te_aid_" + str(x),
                                               tes[ini_per_thread[x]:end_per_thread[x]], ref_tes]) for x in
                        range(cores)]

        create_output_folders(outputdir + "/te_aid")

        localChecks = [p.get() for p in localresults]
        for i in range(len(localChecks)):
            if localChecks[i] == 0:
                print("FATAL ERROR: TE+Aid in Parallel didn't execute well. Please re run.")
                sys.exit(0)
            else:
                pattern = "*.pdf"
                files = glob.glob(outputdir + "/te_aid_" + str(i) + "/" + pattern)

                for file in files:
                    # extract file name form file path
                    file_name = os.path.basename(file)
                    shutil.move(file, outputdir + "/te_aid/" + file_name)
                shutil.rmtree(outputdir + "/te_aid_" + str(i))
        pool.close()
    else:
        print("TE+aid already run!")


def run_te_aid(te_aid_path, genome, outputdir, tes, ref_tes_path):
    status = -1
    create_output_folders(outputdir)
    for sequence in tes:
        seq_name = str(sequence.id).split("#")[0]
        write_sequences_file([sequence], outputdir + "/" + str(seq_name) + ".fa")

        try:
            output = subprocess.run(
                [te_aid_path + '/TE-Aid', '-q', outputdir + "/" + str(seq_name) + ".fa", '-g', genome, '-o', outputdir],
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            delete_files(outputdir + "/" + str(seq_name) + ".fa")
            status = 1
        except:
            status = 0
    return status


def run_blast(library_path, ref_tes, cores, perc_identity, perc_cov):
    keep_seqs = []
    orders = []

    if not os.path.exists(library_path + ".nhr"):
        output = subprocess.run(
            ['makeblastdb', '-in', library_path, '-dbtype', 'nucl'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)

    if not os.path.exists(ref_tes + ".blast"):
        output = subprocess.run(
            ['blastn', '-query', ref_tes, '-db', library_path, '-out', ref_tes + ".blast", '-num_threads', str(cores),
             "-outfmt", "6", "-qcov_hsp_perc", str(perc_cov), "-perc_identity", str(perc_identity), "-max_hsps", "1"],
            stdout=subprocess.PIPE, text=True)
    else:
        print("WARNING: Blast output already exists, skipping BLASTn....")

    orders_superfamilies = [x for x in dicc_orders.values()]
    blastresult = open(ref_tes + ".blast", "r").readlines()

    for hit in blastresult:
        if hit.split("\t")[0].split("#")[0] not in keep_seqs:
            blasted_order = hit.replace("?", "").split("\t")[1].split("#")[1]
            if "/" in blasted_order:
                order_given = blasted_order.split("/")[0].upper()
                superfamily = blasted_order.split("/")[1].upper()
            else:
                order_given = blasted_order.upper()
                superfamily = "NA"

            if superfamily.replace('-', '').upper() in orders_superfamilies:
                # superfamily found !!!
                orders.append(orders_superfamilies.index(superfamily.replace('-', '').upper()) + 2)
                keep_seqs.append(hit.split("\t")[0].split("#")[0])
            elif order_given.replace('-', '').upper() in orders_superfamilies:
                # Well, superfamily didn't find, but I found the order
                orders.append(orders_superfamilies.index(order_given.replace('-', '').upper()) + 2)
                keep_seqs.append(hit.split("\t")[0].split("#")[0])
            else:
                print("WARNING: Neither order and superfamily found: " + blasted_order)
                orders.append(31)
                keep_seqs.append(hit.split("\t")[0].split("#")[0])

    delete_files(ref_tes + ".blast")
    delete_files(library_path + ".nhr")
    delete_files(library_path + ".nin")
    delete_files(library_path + ".nsq")
    return keep_seqs, orders


def run_cdhit(ref_tes, outputdir, cores):
    create_output_folders(outputdir)

    if not os.path.exists(outputdir + "/non_redundant_lib.fa"):
        output = subprocess.run(
            ['cd-hit-est', '-i', ref_tes, '-o', outputdir + '/non_redundant_lib.fa', '-c', '0.8', '-aS', '0.8', '-G', '0', '-g' , '1', '-b', '500', '-T', str(cores), '-M', '0'],
            stdout=subprocess.PIPE, text=True)
    else:
        print("WARNING: cd-hit-est output already exists, skipping cd-hit-est....")

    return outputdir + "/non_redundant_lib.fa"


if __name__ == '__main__':

    Installation_path = os.path.dirname(os.path.realpath(__file__))
    ####################################################################################################################
    # Needed by Classified module
    ####################################################################################################################
    tools_path = Installation_path + "/tools/"
    library_path = Installation_path + "/db/allDatabases.clustered_rename.fa"
    ref_profiles = Installation_path + "/db/Pfam35.0.hmm"
    blastn_db = Installation_path + "/db/repbase20.05_ntSeq_cleaned_TE.fa"
    blastx_db = Installation_path + "/db/repbase20.05_aaSeq_cleaned_TE.fa"
    rnas_db = Installation_path + "/db/rRNA_Eukaryota.hmm"
    minDomLTR = 3  # minimum number of domains to automatically kept an LTR retrotransposon in classified module
    minFLNA = 2  # minimum number of full length fragments of copies to automatically kept a non-autonomous element

    ####################################################################################################################
    # Needed by BEE method to extension
    ####################################################################################################################
    max_nns = 10  # threshold of the percentage of Ns allowed in the final extended consensus to keep it (1-100)

    ####################################################################################################################
    # Needed by Unclassified module
    ####################################################################################################################
    FLF_UNCLASS = 2  # minimum number of full length copies to consider a TE in unclassified module
    perc_ident = 70  # Percentage of identity to a known element to keep a TE in unclassified module
    perc_cover = 70  # Percentage of the hit coverage in the TE sequence to keep it in unclassified module
    ####################################################################################################################

    print("\n#########################################################################")
    print("#                                                                       #")
    print("#                             MCHelper                                  #")
    print("#   Instructions:  Execute the program, press q to close the graph and  #")
    print("#     press a number 2-33 if you want to keep the TE, or 1 otherwise.   #")
    print("#                                                                       #")
    print("#########################################################################\n")

    ### read parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--module', required=True, dest='module_user',
                        help='module of curation [A, C, U, E]. Required*')
    parser.add_argument('-i', '--input', required=False, dest='input_dir',
                        help='Directory with the files required to do the curation (REPET output directory). Required*')
    parser.add_argument('-g', '--genome', required=False, dest='genome',
                        help='Genome used to detect the TEs. Required*')
    parser.add_argument('-o', '--output', required=False, dest='outputdir',
                        help='Path to the output directory. Required*')
    parser.add_argument('--te_aid', required=False, dest='te_aid', default='Y',
                        help='Do you want to use TE-aid? [Y or N]. Default=Y')
    parser.add_argument('-a', required=False, dest='automatic', default='F',
                        help='Level of automation: F: fully automated, S: semi-automated, M: fully manual?. Default=F')
    parser.add_argument('-n', required=False, dest='proj_name',
                        help='REPET project name. Required*')
    parser.add_argument('-t', required=False, dest='cores', default=-1,
                        help='cores to execute some steps in parallel')
    parser.add_argument('-j', required=False, dest='module2_seqs_file',
                        help='Path to the sequences to be used in the extension module')
    parser.add_argument('-k', required=False, dest='module3_seqs_file',
                        help='Path to the sequences to be used in the unclassified module')
    parser.add_argument('-m', required=False, dest='ref_library_module3',
                        help='Path to the sequences to be used as references in the unclassified module')
    parser.add_argument('-v', required=False, dest='verbose', default='N',
                        help='Verbose? [Y or N]. Default=N')
    parser.add_argument('--input_type', required=False, dest='input_type',
                        help='Input type: fasta or REPET.')
    parser.add_argument('-l', required=False, dest='user_library',
                        help='User defined library to be used with input type fasta.')
    parser.add_argument('-b', required=False, dest='busco_library',
                        help='Reference/BUSCO genes to filter out TEs (HMMs expected). ')
    parser.add_argument('-c', required=False, dest='minFullLenCopies', default=1,
                        help='Minimum number of full-length copies to process an element.')
    parser.add_argument('-s', required=False, dest='perc_ssr',
                        help='Maximum length covered by single repetitions (in percentage between 0-100) allowed for a TE not to be removed')
    parser.add_argument('-e', required=False, dest='ext_nucl',
                        help='Number of nucleotides to extend each size of the element. Default=1000')
    parser.add_argument('-x', required=False, dest='num_ite',
                        help='Number of iterations to extend the elements')
    parser.add_argument('--version', action='version', version='%(prog)s v1.1.0')


    options = parser.parse_args()
    module_user = options.module_user
    input_dir = options.input_dir
    proj_name = options.proj_name
    genome = options.genome
    outputdir = options.outputdir
    te_aid = options.te_aid
    automatic = options.automatic
    cores = options.cores
    module2_seqs_file = options.module2_seqs_file
    module3_seqs_file = options.module3_seqs_file
    ref_library_module3 = options.ref_library_module3
    verbose = options.verbose
    input_type = options.input_type
    user_library = options.user_library
    busco_library = options.busco_library
    minFullLenCopies = int(options.minFullLenCopies)
    perc_ssr = options.perc_ssr
    ext_nucl = options.ext_nucl
    num_ite = options.num_ite

    ####################################################################################################################
    # Parameter validation
    ####################################################################################################################
    module = 0
    if module_user.upper() not in ['A', 'C', 'U', 'E', '3333']:  # 3333 for debugging only
        print('FATAL ERROR: module (-r parameter) must be A (all steps), U (unclassified module), or E (BEE extension)')
        sys.exit(0)
    else:
        if module_user.upper() == 'A':
            module = 123
        elif module_user.upper() == 'C':
            module = 1
        elif module_user.upper() == 'U':
            module = 3
        elif module_user.upper() == 'E':
            module = 2
        elif module_user.upper() == '3333':
            module = 3333

    if outputdir is None:
        outputdir = os.getcwd()
        print('MESSAGE: Output directory will be '+outputdir)
    elif not os.path.exists(outputdir):
        print('FATAL ERROR: Output folder ' + outputdir + " doesn't exist.")
        sys.exit(0)
    else:
        outputdir = os.path.abspath(outputdir)

    if verbose is None:
        verbose = False
        print('MESSAGE: Verbose is not activated')
    elif verbose.upper() not in ['Y', 'N']:
        print('FATAL ERROR: verbose (-v) must be Y or N')
        sys.exit(0)
    else:
        if verbose == 'Y':
            verbose = True
        else:
            verbose = False

    if cores is None or cores == -1:
        cores = int(psutil.cpu_count())
        print("MESSAGE: Missing threads parameter, using by default: " + str(cores))
    else:
        cores = int(cores)

    if perc_ssr is None:
        perc_ssr = 60
        print("MESSAGE: Missing -s parameter, using by default: " + str(perc_ssr))
    elif not perc_ssr.isnumeric():
        print('FATAL ERROR: -s must be numeric, but I received: '+str(perc_ssr))
        sys.exit(0)
    elif int(perc_ssr) < 0 or int(perc_ssr) > 100:
        print('FATAL ERROR: -s must be between 0 and 100')
        sys.exit(0)
    else:
        perc_ssr = int(perc_ssr)

    if ext_nucl is None:
        ext_nucl = 1000
        print("MESSAGE: Missing -e parameter, using by default: " + str(ext_nucl))
    elif not ext_nucl.isnumeric():
        print('FATAL ERROR: -e must be numeric, but I received: '+str(ext_nucl))
        sys.exit(0)
    elif int(ext_nucl) < 1 or int(perc_ssr) > 1500:
        print('FATAL ERROR: -e must be between 1 and 1500')
        sys.exit(0)
    else:
        ext_nucl = int(ext_nucl)

    if num_ite is None:
        num_ite = 1
        print("MESSAGE: Missing -x parameter, using by default: " + str(num_ite))
    elif not num_ite.isnumeric():
        print('FATAL ERROR: -x must be numeric, but I received: '+str(num_ite))
        sys.exit(0)
    elif int(num_ite) < 0 or int(num_ite) > 100:
        print('FATAL ERROR: -x must be between 0 and 100')
        sys.exit(0)
    else:
        num_ite = int(num_ite)

    ####################################################################################################################
    # Classified module
    ####################################################################################################################
    if module in [1, 123]:
        if te_aid is None:
            te_aid = 'Y'
            print('MESSAGE: Using by default te_aid = Y')
        elif te_aid.upper() not in ['Y', 'N']:
            print('FATAL ERROR: unknown value of --te_aid parameter: ' + te_aid + '. This parameter must be Y or N')
            sys.exit(0)
        else:
            te_aid = te_aid.upper()
        if automatic is not None:
            if automatic.upper() not in ['F', 'S', 'M']:
                print('FATAL ERROR: automation level (-a) must be F, S or M (see help.)')
                sys.exit(0)
            else:
                automatic = automatic.upper()
        if genome is None:
            print('FATAL ERROR: -g parameter must be specified in for the classified module')
            sys.exit(0)
        if not os.path.exists(genome):
            print("FATAL ERROR: Genome file " + genome + " doesn't exist.")
            sys.exit(0)
        if input_type == None:
            print(
                'FATAL ERROR: You need to specify the input type parameter (--input_type). Values can be fasta or repet')
            sys.exit(0)
        elif input_type.upper() not in ['FASTA', 'REPET']:
            print(
                'FATAL ERROR: Unknown value in the input type parameter (--input_type). Values must be fasta or repet')
            sys.exit(0)
        else:
            input_type = input_type.lower()
        if (module == 1 or module == 123) and busco_library is None:
            print('FATAL ERROR: -b parameter must be specified for the unclassified module')
            sys.exit(0)
        elif (module == 1 or module == 123) and not os.path.exists(busco_library):
            print("FATAL ERROR: Reference/BUSCO genes file " + busco_library + " doesn't exist.")
            sys.exit(0)

        ################################################################################################################
        # classified module's Pre-validations :
        ################################################################################################################
        print("MESSAGE: Starting with Classified module...")
        use_repet = True
        if input_type == 'repet':
            if input_dir is None:
                print('FATAL ERROR: -i parameter must be specified in for the classified module and input type REPET')
                sys.exit(0)
            if proj_name is None:
                print('FATAL ERROR: -n parameter must be specified in for the classified module and input type REPET')
                sys.exit(0)

            start_time = time.time()
            input_valid, reason_valid = check_repet_input_folder(input_dir, proj_name)
            end_time = time.time()
            if verbose:
                print("MESSAGE: REPET Input checking done: [" + str(end_time - start_time) + " seconds]")

            if input_valid:
                ref_tes = input_dir + "/" + proj_name + "_refTEs.fa"
                flf_file = input_dir + "/TEannot/" + proj_name + "_chr_allTEs_nr_noSSR_join_path.annotStatsPerTE_FullLengthFrag.txt"
                repet_table = input_dir + "/" + proj_name + "_denovoLibTEs_PC.classif"
                plots_dir = input_dir + "/plotCoverage"
                gff_files = input_dir + "/gff_reversed"
            else:
                print(reason_valid)
                sys.exit(0)
        elif input_type == 'fasta':
            use_repet = False
            if user_library is None:
                print('FATAL ERROR: -l parameter must be specified in for the classified module using input type fasta')
                sys.exit(0)
            if not os.path.exists(user_library):
                print("FATAL ERROR: TE library file " + user_library + " doesn't exist.")
                sys.exit(0)
            ref_tes = user_library

            start_time = time.time()
            if check_classification_Userlibrary(ref_tes, outputdir) == 0:
                flf_file = count_flf_fasta(ref_tes, genome, cores, outputdir)
                if not os.path.exists(outputdir + "/classifiedModule/denovoLibTEs_PC.classif") and not os.path.exists(
                        outputdir + "/classifiedModule/new_user_lib.fa"):
                    if automatic == 'M':
                        do_blast = True
                    else:
                        do_blast = False
                repet_table = outputdir + "/classifiedModule/denovoLibTEs_PC.classif"

                gff_files = ""
                plots_dir = ""
            else:
                print(
                    'FATAL ERROR: There are some sequences with problems in your library. Please check the file: ' + outputdir + '/sequences_with_problems.txt')
                sys.exit(0)

            end_time = time.time()
            if verbose:
                print("MESSAGE: Fasta pre-processing done: [" + str(end_time - start_time) + " seconds]")

        if not os.path.exists(outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa"):
            ############################################################################################################
            # First step: reduce redundancy with cd-hit
            ############################################################################################################
            start_time = time.time()
            ref_tes = run_cdhit(ref_tes, outputdir + "/classifiedModule/", cores)
            end_time = time.time()
            if verbose:
                print("MESSAGE: The library was reduced to " + str(len(list(SeqIO.parse(ref_tes, 'fasta')))) + " after clustering by cd-hit [" + str(end_time - start_time) + " seconds]")

            ############################################################################################################
            # Second step: extract seqs with at least one full length fragment and minimum minFullLenCopies
            ############################################################################################################
            start_time = time.time()
            new_ref_tes, num_copies = filter_flf(ref_tes, flf_file, minFullLenCopies, outputdir+"/classifiedModule/")
            end_time = time.time()
            if verbose:
                print("MESSAGE: The library was reduced to " + str(len(list(SeqIO.parse(new_ref_tes, 'fasta')))) + " after FLF filtering [" + str(end_time - start_time) + " seconds]")

            ############################################################################################################
            # Third step: SSR, BUSCO genes, tRNA and rRNA filters
            ############################################################################################################
            start_time = time.time()

            new_ref_tes, deleted_seqs = filter_bad_candidates(new_ref_tes, perc_ssr, outputdir+"/classifiedModule/", tools_path, busco_library, rnas_db, cores)
            end_time = time.time()
            if verbose:
                print("MESSAGE: The library was reduced to " + str(len(list(SeqIO.parse(new_ref_tes, 'fasta')))) + " after SSR, genes and RNA filtering [" + str(end_time - start_time) + " seconds]")

            ############################################################################################################
            # Fourth step: Structural checks, BEE method, and Visualize plots
            ############################################################################################################
            start_time = time.time()
            # If repet, we need to include the classification info into the library
            if use_repet:
                rename_tes = []
                struc_table = pd.read_csv(repet_table, sep='\t')
                for i in range(struc_table.shape[0]):
                    is_in_lib = [x for x in SeqIO.parse(new_ref_tes, "fasta") if str(x.id).split("#")[0] == struc_table.at[i, "Seq_name"]]
                    if len(is_in_lib) > 0:
                        bad_named_te = is_in_lib[0]
                        classif = str(struc_table.at[i, "order"]) + "/" + str(struc_table.at[i, "sFamily"])
                        bad_named_te.id = str(bad_named_te.id).split("#")[0] + "#" + classif
                        bad_named_te.description = ""
                        rename_tes.append(bad_named_te)

                write_sequences_file(rename_tes, new_ref_tes+"_tmp")
                delete_files(new_ref_tes)
                shutil.move(new_ref_tes+"_tmp", new_ref_tes)

            new_module1(plots_dir, new_ref_tes, gff_files, outputdir, proj_name, te_aid, automatic, minDomLTR,
                        num_copies, minFLNA, verbose, use_repet, ext_nucl, max_nns, num_ite, blastn_db, blastx_db, tools_path,
                        ref_profiles, library_path)

        else:
            print("MESSAGE: classified module was already run, please verify or change the output directory")

        if module == 123:
            new_ref_library = []
            if os.path.exists(outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa"):
                for te in SeqIO.parse(outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa", "fasta"):
                    new_ref_library.append(te)

            if os.path.exists(outputdir + "/classifiedModule/kept_seqs_classified_module_non_curated.fa"):
                for te in SeqIO.parse(outputdir + "/classifiedModule/kept_seqs_classified_module_non_curated.fa", "fasta"):
                    new_ref_library.append(te)

            if len(new_ref_library) > 0:
                write_sequences_file(new_ref_library, outputdir + "/classifiedModule/kept_seqs_classified_module.fa")
                ref_library_module3 = outputdir + "/classifiedModule/kept_seqs_classified_module.fa"
            else:
                print("WARNING: Classified module didn't find any TE")
                ref_library_module3 = ""
            module3_seqs_file = outputdir + "/classifiedModule/input_to_unclassified_module_seqs.fa"

    ####################################################################################################################
    # Unclassified module
    ####################################################################################################################
    if module in [3, 123]:
        if module == 3 and module3_seqs_file is None:
            print('FATAL ERROR: -k parameter must be specified for the unclassified module')
            sys.exit(0)
        if module == 3 and ref_library_module3 is None:
            print('FATAL ERROR: -m parameter must be specified for the unclassified module')
            sys.exit(0)
        if module == 3 and input_type is None:
            print('FATAL ERROR: --input_type parameter must be specified for the unclassified module')
            sys.exit(0)
        if module == 3 and genome is None:
            print('FATAL ERROR: -g parameter must be specified for the unclassified module')
            sys.exit(0)
        if module == 3 and busco_library is None:
            print('FATAL ERROR: -b parameter must be specified for the unclassified module')
            sys.exit(0)
        elif module == 3 and not os.path.exists(busco_library):
            print("FATAL ERROR: Reference/BUSCO genes file " + busco_library + " doesn't exist.")
            sys.exit(0)

        if module == 123 and (os.path.getsize(module3_seqs_file) == 0 or module3_seqs_file == ""):
            print("MESSAGE: There is no sequences to unclassified Module, skipping ...")
        else:
            if module == 3 and input_type == "fasta":
                flf_file = count_flf_fasta(module3_seqs_file, genome, cores, outputdir)
            elif module == 3 and input_type == "repet":
                if module == 3 and input_dir is None:
                    print('FATAL ERROR: -i parameter must be specified for the unclassified module')
                    sys.exit(0)
                if module == 3 and proj_name is None:
                    print('FATAL ERROR: -n parameter must be specified for the unclassified module')
                    sys.exit(0)
                flf_file = input_dir + "/TEannot/" + proj_name + "_chr_allTEs_nr_noSSR_join_path.annotStatsPerTE_FullLengthFrag.txt"

            if ref_library_module3 != "":
                create_output_folders(outputdir + "/unclassifiedModule")
                print("MESSAGE: Starting with Unclassified module...")
                ########################################################################################################
                # First step: Filter elements with not enough FLF copies in the genome
                ########################################################################################################
                start_time = time.time()
                module3_seqs_file, num_copies = filter_flf(module3_seqs_file, flf_file, FLF_UNCLASS, outputdir+"/unclassifiedModule/")
                end_time = time.time()
                if verbose:
                    print("MESSAGE: The library was reduced to " + str(len(list(SeqIO.parse(module3_seqs_file, 'fasta')))) + " after FLF filtering [" + str(end_time - start_time) + " seconds]")

                ########################################################################################################
                # Second step: run BEE in parallel
                ########################################################################################################
                start_time = time.time()
                run_BEE_parallel(genome, module3_seqs_file, ext_nucl, num_ite, outputdir + "/unclassifiedModule/", max_nns, cores)
                module3_seqs_file, num_copies = filter_flf(module3_seqs_file, flf_file, FLF_UNCLASS, outputdir + "/unclassifiedModule/")
                end_time = time.time()
                if verbose:
                    print("MESSAGE: The sequences were extended successfully [" + str(end_time - start_time) + " seconds]")
                ########################################################################################################
                # Third step: find structural features
                ########################################################################################################
                start_time = time.time()
                build_class_table_parallel(outputdir + "/unclassifiedModule/extended_cons.fa", cores, outputdir + '/unclassifiedModule/',
                                           blastn_db, blastx_db, tools_path, ref_profiles, False)
                end_time = time.time()
                if verbose:
                    print("MESSAGE: TE Feature table was created [" + str(end_time - start_time) + " seconds]")

                ########################################################################################################
                # Fourth step: process the sequences, BLASTn againts a given curated library, inferred classifications,
                # and create outputs
                ########################################################################################################
                module3(outputdir + "/unclassifiedModule/extended_cons.fa", ref_library_module3, cores, outputdir,
                   perc_ident, perc_cover, library_path)
            else:
                if module == 123:
                    print("WARNING: Unclassified module didn't run because classified module didn't find any TE")
                else:
                    print("WARNING: "+ref_library_module3+" is empty")

    ####################################################################################################################
    # BEE module
    ####################################################################################################################
    if module == 2:
        if user_library is None:
            print('FATAL ERROR: -l parameter must be specified for the BEE extension module using input type fasta')
            sys.exit(0)
        if not os.path.exists(user_library):
            print("FATAL ERROR: TE library file " + user_library + " doesn't exist.")
            sys.exit(0)
        if genome is None:
            print('FATAL ERROR: -g parameter must be specified in for the BEE extension module')
            sys.exit(0)
        if not os.path.exists(genome):
            print("FATAL ERROR: Genome file " + genome + " doesn't exist.")
            sys.exit(0)

        ########################################################################################################
        # First step: run BEE in parallel
        ########################################################################################################
        create_output_folders(outputdir + "/bee_module/")
        run_BEE_parallel(genome, user_library, ext_nucl, num_ite, outputdir + "/bee_module/", max_nns, cores)

    ####################################################################################################################
    # Debugging only !!!!
    if module in [3333]:
        print("Debugging....")

    ####################################################################################################################
    # Writing the final results
    final_seqs = []
    if module in [1, 12, 123]:
        try:
            for te in SeqIO.parse(outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa", "fasta"):
                final_seqs.append(te)
        except FileNotFoundError:
            print("WARNING: the file " + outputdir + "/classifiedModule/kept_seqs_classified_module_curated.fa was not created. Please check.")
        try:
            for te in SeqIO.parse(outputdir + "/classifiedModule/kept_seqs_classified_module_non_curated.fa", "fasta"):
                final_seqs.append(te)
        except FileNotFoundError:
            print("WARNING: the file " + outputdir + "/classifiedModule/kept_seqs_classified_module_non_curated.fa was not created. Please check.")
        delete_files(outputdir + "/classifiedModule/non_redundant_lib.fa")
        delete_files(outputdir + "/classifiedModule/non_redundant_lib.fa.clstr")
    if module in [3, 123]:
        try:
            for te in SeqIO.parse(outputdir + "/unclassifiedModule/kept_seqs_unclassified_module.fa", "fasta"):
                final_seqs.append(te)
        except FileNotFoundError:
            print("WARNING: the file " + outputdir + "/unclassifiedModule/kept_seqs_unclassified_module.fa was not created. Please check.")

    # Reduce the final redundancy
    if len(final_seqs) > 0:
        write_sequences_file(final_seqs, outputdir + "/curated_sequences_R.fa")
        non_redundat = run_cdhit(outputdir + "/curated_sequences_R.fa", outputdir, cores)
        shutil.move(non_redundat, outputdir + "/curated_sequences_NR.fa")

