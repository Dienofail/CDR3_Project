def analysis(Sequence_Reads, output_path, results, chain, with_statistics=True, with_reverse_complement_search=True, omitN=True):
    import numpy as np
    import decimal as dec
    import string
    import operator as op
    import collections as coll
    import os
    from Bio import SeqIO
    from acora import AcoraBuilder
    from time import time, clock
    from string import Template
    from operator import itemgetter, attrgetter
    import Levenshtein as lev
    print Sequence_Reads
    results_file = open(str(results) + '.txt', "w")
    # output_path = str(results) + '.txt'
    # J_results_file = open( str(Jresults) + '.txt', "w")
    Nseqs = 0
    # print("Path at terminal when executing this file")
    # print(os.getcwd() + "\n")
    chain = 'beta'
    print output_path
    

    v_half_split, j_half_split = [10, 6]  # Do not change - V tags are split at position 10, J at position 6, to look for half tags if no full tag is found.

    ################

    print 'Commencing analysis on a total of', Sequence_Reads, 'file(s)'

    ################
    print ('Importing known V, D and J gene segments and tags...')

    handle = open("human_TRBV_region.fasta", "rU")
    v_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRBJ_region.fasta", "rU")
    j_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    v_regions = []
    for j in range(0, len(v_genes)):
        v_regions.append(string.upper(v_genes[j].seq))

    j_regions = []
    for j in range(0, len(j_genes)):
        j_regions.append(string.upper(j_genes[j].seq))

    ##############
    # # Build keyword tries of V and J tags for fast assignment
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_trbv.txt", "rU"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("new_tags_trbj.txt", "rU"), j_half_split)   

    v_builder = AcoraBuilder()
    for i in range(0, len(v_seqs)):
        v_builder.add(str(v_seqs[i]))  # Add all V tags to keyword trie

    v_key = v_builder.build()

    j_builder = AcoraBuilder()
    for i in range(0, len(j_seqs)):
        j_builder.add(str(j_seqs[i]))  # Add all J tags to keyword trie

    j_key = j_builder.build()

    ##############
    # Build keyword tries for first and second halves of both V and J tags
    v_half1_builder = AcoraBuilder()
    for i in range(0, len(half1_v_seqs)):
        v_half1_builder.add(str(half1_v_seqs[i]))
    half1_v_key = v_half1_builder.build()

    v_half2_builder = AcoraBuilder()
    for i in range(0, len(half2_v_seqs)):
        v_half2_builder.add(str(half2_v_seqs[i]))
    half2_v_key = v_half2_builder.build()

    # j_half1_builder = AcoraBuilder()
    # for i in range(0,len(half1_j_seqs)):
        # j_half1_builder.add(str(half1_j_seqs[i]))
    # half1_j_key = j_half1_builder.build()

    # j_half2_builder = AcoraBuilder()
    # for i in range(0,len(half2_j_seqs)):
        # j_half2_builder.add(str(half2_j_seqs[i]))
    # half2_j_key = j_half2_builder.build()

    ###############
    # # Initialise variables
    assigned_count = 0  # this will just increase by one every time we correctly assign a seq read with all desired variables
    seq_count = 0  # this will simply track the number of sequences analysed in file
    Nseqs = 0  # Number of raw reads containing N nucleotides
    t0 = time()  # Begin timer
    assigned_v_count = 0
    assigned_j_count = 0
    total_count = 0
    total_quality = 0
    twenty_five = 0
    thirty = 0
    thirty_five = 0
    fourty = 0
    ###############
    # # Open .txt file created at the start of analysis
    stemplate = Template("$v\t$j\t$del_v\t$del_j\t$nt_insert\t$seqid\t$v_end\t$quality")  # Creates stemplate, a holder, for f. Each line will have the 5 variables separated by a space
    # output_template = Template('$total\t$all_found    $v_found    $j_found    $time    $v_end')
    ###############
    # # Begin analysing sequences

        
    print 'Importing sequences from', Sequence_Reads, ' and assigning V and J regions...'
    handle = open(Sequence_Reads, "rU")
    
    for record in SeqIO.parse(handle, "fastq"):
        total_count += 1
        if 'N' in str(record.seq) and omitN == True:
            Nseqs += 1
            #print str(record.seq)
            
        else:
            forward_phred = record.letter_annotations["phred_quality"]
            # print to_print_str
           # list_of_str = to_print_str.split("\s")
            # for idx, val in enumerate(list_of_str):
                # print 'Currently at ' + str(idx) + ' With value ' + str(val) + "\n"
            # print "The fifth element of to_print_str is " + str(to_print_str[4]) + "\n"
            found_seq_match = 0
            found_v_match = 0
            found_j_match = 0
            seq_count += 1
            average_quality = 0
            my_iterator = 0
            my_end_iterator = 0
            hold_v = v_key.findall(str(record.seq))
            hold_j = j_key.findall(str(record.seq))
 
            if hold_v:
                v_match = v_seqs.index(hold_v[0][0])  # Assigns V
                assigned_v_count += 1
                temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1  # Finds where the end of a full V would be
                my_iterator = hold_v[0][1] 
                my_end_iterator = my_iterator + 21
                
                if get_v_deletions(record.seq, v_match, temp_end_v, v_regions):  # If the number of deletions has been found
                    [ end_v, deletions_v] = get_v_deletions(record.seq, v_match, temp_end_v, v_regions)
            else:
                found_v_match = 0
                hold_v1 = half1_v_key.findall(str(record.seq))
                hold_v2 = half2_v_key.findall(str(record.seq))
                for i in range(len(hold_v1)):
                    indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v1[i][1]:hold_v1[i][1] + len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                            if lev.hamming(v_seqs[k], str(record.seq)[hold_v1[i][1]:hold_v1[i][1] + len(v_seqs[k])]) <= 1:
                                v_match = k
                                temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1  # Finds where the end of a full V would be
                                found_v_match += 1
                                # assigned_v_count += 1
                                my_iterator = hold_v1[i][1] 
                                my_end_iterator = my_iterator + 11
                                 
                for i in range(len(hold_v2)):
                    indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1] - v_half_split:hold_v2[i][1] + len(v_seqs[half2_v_seqs.index(hold_v2[i][0])]) - v_half_split]):
                            if lev.hamming(v_seqs[k], str(record.seq)[hold_v2[i][1] - v_half_split:hold_v2[i][1] + len(v_seqs[k]) - v_half_split]) <= 1:
                                v_match = k
                                temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1  # Finds where the end of a full V would be
                                found_v_match += 1
                                # assigned_v_count += 1
                                my_iterator = hold_v2[i][1]
                                my_end_iterator = my_iterator + 11
                if found_v_match == 1:
                    assigned_v_count += 1
                        
            if hold_j:
                j_match = j_seqs.index(hold_j[0][0])  # Assigns J
                assigned_j_count += 1
                temp_start_j = hold_j[0][1] - jump_to_start_j[j_match]  # Finds where the start of a full J would be
                if get_j_deletions(record.seq, j_match, temp_start_j, j_regions):
                    [ start_j, deletions_j] = get_j_deletions(record.seq, j_match, temp_start_j, j_regions)
            # else:
                # found_j_match = 0
                # hold_j1 = half1_j_key.findall(str(record.seq))
                # hold_j2 = half2_j_key.findall(str(record.seq))
                # for i in range(len(hold_j1)):
                    # indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                    # for k in indices:
                        # if len(j_seqs[k]) == len(str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                            # if lev.hamming( j_seqs[k], str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                # j_match = k
                                # temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                # found_j_match += 1
                # for i in range(len(hold_j2)):
                    # indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                    # for k in indices:
                        # if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])-j_half_split]):
                            # if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                # j_match = k
                                # temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                # found_j_match += 1

            if hold_v and hold_j:
                if get_v_deletions(record.seq, v_match, temp_end_v, v_regions) and get_j_deletions(record.seq, j_match, temp_start_j, j_regions):
                    quality_array = forward_phred[my_iterator:my_end_iterator]
                    counter = 0
                    for idx, val in enumerate(quality_array):
                        average_quality += val
                        counter += 1
                    true_quality = average_quality / counter
                    total_quality += true_quality
                    if true_quality < 25:
                        twenty_five += 1
                    elif true_quality <30:
                        thirty += 1
                    elif true_quality <35:
                        thirty_five += 1
                    elif true_quality =<40:
                        fourty += 1
                    f_seq = stemplate.substitute(v=str(v_match), j=str(j_match), del_v=str(deletions_v), del_j=str(deletions_j) , nt_insert=str(record.seq[end_v + 1:start_j]), seqid=str(record.seq), v_end=str(temp_end_v), quality=str(true_quality))
                    print >> results_file, f_seq  # Write to results_file (text file) the classification of the sequence
                    assigned_count += 1
                    found_seq_match = 1
            # elif hold_v and found_j_match == 1:
                # if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    # [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                    # f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id) )
                    # print >> results_file, f_seq
                    # assigned_count += 1
                    # found_seq_match = 1
            elif found_v_match == 1 and hold_j:
                if get_v_deletions(record.seq, v_match, temp_end_v, v_regions) and get_j_deletions(record.seq, j_match, temp_start_j, j_regions):
                    quality_array = forward_phred[my_iterator:my_end_iterator]
                    counter = 0
                    for idx, val in enumerate(quality_array):
                        average_quality += val
                        counter += 1
                    true_quality = average_quality / counter
                    if true_quality < 25:
                        twenty_five += 1
                    elif true_quality <30:
                        thirty += 1
                    elif true_quality <35:
                        thirty_five += 1
                    elif true_quality =<40:
                        fourty += 1
                    total_quality += true_quality
                    [ end_v, deletions_v] = get_v_deletions(record.seq, v_match, temp_end_v, v_regions)
                    f_seq = stemplate.substitute(v=str(v_match), j=str(j_match), del_v=str(deletions_v), del_j=str(deletions_j), nt_insert=str(record.seq[end_v + 1:start_j]), seqid=str(record.seq), v_end=str(temp_end_v), quality=str(true_quality))
                    print >> results_file, f_seq
                    assigned_count += 1
                    found_seq_match = 1
            # elif found_v_match == 1 and found_j_match == 1:
                # if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    # [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
                    # [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
                    # f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record.seq[end_v+1:start_j])+str(','), seqid = str(record.id) )
                    # print >> results_file, f_seq
                    # assigned_count += 1
                    # found_seq_match = 1
            if found_seq_match == 0 and with_reverse_complement_search == True:
                
                #####################
                # REVERSE COMPLEMENT
                #####################

                record_reverse = record.reverse_complement()
                reverse_phred = record_reverse.letter_annotations["phred_quality"]
                hold_v = v_key.findall(str(record_reverse.seq))
                hold_j = j_key.findall(str(record_reverse.seq))
                
                if hold_v:                
                    v_match = v_seqs.index(hold_v[0][0])  # Assigns V
                    assigned_v_count += 1
                    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1  # Finds where the end of a full V would be
                    my_iterator = hold_v[0][1] 
                    my_end_iterator = my_iterator + 21
                    if get_v_deletions(record_reverse.seq, v_match, temp_end_v, v_regions):  # If the number of deletions has been found
                        [ end_v, deletions_v] = get_v_deletions(record_reverse.seq, v_match, temp_end_v, v_regions)

                else:
                    found_v_match = 0
                    hold_v1 = half1_v_key.findall(str(record_reverse.seq))
                    hold_v2 = half2_v_key.findall(str(record_reverse.seq))
                    for i in range(len(hold_v1)):
                        indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1] + len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                if lev.hamming(v_seqs[k], str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1] + len(v_seqs[k])]) <= 1:
                                    v_match = k
                                    temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1  # Finds where the end of a full V would be
                                    found_v_match += 1
                                    my_iterator = hold_v1[i][1] 
                                    my_end_iterator = my_iterator + 11
                                    # assigned_v_count += 1
                    for i in range(len(hold_v2)):
                        indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1] - v_half_split:hold_v2[i][1] + len(v_seqs[half2_v_seqs.index(hold_v2[i][0])]) - v_half_split]):
                                if lev.hamming(v_seqs[k], str(record_reverse.seq)[hold_v2[i][1] - v_half_split:hold_v2[i][1] + len(v_seqs[k]) - v_half_split]) <= 1:
                                    v_match = k
                                    temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - v_half_split - 1  # Finds where the end of a full V would be
                                    found_v_match += 1
                                    my_iterator = hold_v2[i][1] 
                                    my_end_iterator = my_iterator + 11
                                    # assigned_v_count += 1
                    if found_v_match == 1:
                        assigned_v_count += 1
                if hold_j:
                    j_match = j_seqs.index(hold_j[0][0])  # Assigns J
                    assigned_j_count += 1
                    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match]  # Finds where the start of a full J would be
                    if get_j_deletions(record_reverse.seq, j_match, temp_start_j, j_regions):
                        [ start_j, deletions_j] = get_j_deletions(record_reverse.seq, j_match, temp_start_j, j_regions)
                # else:
                    # found_j_match = 0
                    # hold_j1 = half1_j_key.findall(str(record_reverse.seq))
                    # hold_j2 = half2_j_key.findall(str(record_reverse.seq))
                    # for i in range(len(hold_j1)):
                        # indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                        # for k in indices:
                            # if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                # if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                    # j_match = k
                                    # temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                    # found_j_match += 1
                    # for i in range(len(hold_j2)):
                        # indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                        # for k in indices:
                            # if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])-j_half_split]):
                                # if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]-j_half_split:hold_j2[i][1]+len(j_seqs[k])-j_half_split] ) <= 1:
                                    # j_match = k
                                    # temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - j_half_split # Finds where the start of a full J would be
                                    # found_j_match += 1
                
                if hold_v and hold_j:
                    if get_v_deletions(record_reverse.seq, v_match, temp_end_v, v_regions) and get_j_deletions(record_reverse.seq, j_match, temp_start_j, j_regions):
                        quality_array = reverse_phred[my_iterator:my_end_iterator]
                        counter = 0
                        for idx, val in enumerate(quality_array):
                            average_quality += val
                            counter += 1
                        true_quality = average_quality / counter
                        if true_quality < 25:
                            twenty_five += 1
                        elif true_quality <30:
                            thirty += 1
                        elif true_quality <35:
                            thirty_five += 1
                        elif true_quality =<40:
                            fourty += 1
                        total_quality += true_quality
                        f_seq = stemplate.substitute(v=str(v_match), j=str(j_match), del_v=str(deletions_v), del_j=str(deletions_j), nt_insert=str(record_reverse.seq[end_v + 1:start_j]), seqid=str(record.seq), v_end=str(temp_end_v), quality=str(true_quality))
                        print >> results_file, f_seq  # Write to results_file (text file) the classification of the sequence
                        assigned_count += 1
                        found_seq_match = 1
                # elif hold_v and found_j_match == 1:
                    # if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        # [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                        # f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id) )
                        # print >> results_file, f_seq
                        # assigned_count += 1
                        # found_seq_match = 1
                elif found_v_match == 1 and hold_j:
                    if get_v_deletions(record_reverse.seq, v_match, temp_end_v, v_regions) and get_j_deletions(record_reverse.seq, j_match, temp_start_j, j_regions):
                        quality_array = reverse_phred[my_iterator:my_end_iterator]
                        counter = 0
                        for idx, val in enumerate(quality_array):
                            average_quality += val
                            counter += 1
                        true_quality = average_quality / counter
                        if true_quality < 25:
                            twenty_five += 1
                        elif true_quality <30:
                            thirty += 1
                        elif true_quality <35:
                            thirty_five += 1
                        elif true_quality =<40:
                            fourty += 1
                        total_quality += true_quality
                        [ end_v, deletions_v] = get_v_deletions(record_reverse.seq, v_match, temp_end_v, v_regions)
                        f_seq = stemplate.substitute(v=str(v_match), j=str(j_match), del_v=str(deletions_v), del_j=str(deletions_j), nt_insert=str(record_reverse.seq[end_v + 1:start_j]), seqid=str(record.seq), v_end=str(temp_end_v), quality=str(true_quality))
                        print >> results_file, f_seq
                        assigned_count += 1
                        found_seq_match = 1
                # elif found_v_match == 1 and found_j_match == 1:
                    # if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        # [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                        # [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                        # f_seq = stemplate.substitute( v = str(v_match)+str(','), j = str(j_match)+str(','), del_v = str(deletions_v)+str(','), del_j = str(deletions_j)+str(','), nt_insert = str(record_reverse.seq[end_v+1:start_j])+str(','), seqid = str(record.id) )
                        # print >> results_file, f_seq
                        # assigned_count += 1
                        # found_seq_match = 1
                
    handle.close()
    results_file.close()
    average_total_quality = total_quality / assigned_count
    print 'Number of V genes assigned', str(assigned_v_count), 'Number of J genes assigned', str(assigned_j_count), 'Total count', str(total_count), "Average Total Quality", str(average_total_quality)
    # J_results_file.close()
    timed = time() - t0
    print str(seq_count), 'sequences were analysed'
    print str(assigned_count), 'sequences were successfully assigned'
    if omitN == True:
        print Nseqs, 'sequences contained ambiguous N nucleotides'
    print 'Time taken =', str(timed), 'seconds'
    # output_template = Template('$total    $all_found    $v_found    $j_found    $time')
    # output_template.substitute(total=str(total_count), all_found=str(assigned_count), v_found=str(assigned_v_count), j_found=str(assigned_j_count), time=str(timed))
    stringtoprint = str(Sequence_Reads) + "\t" + str(total_count) + "\t" + str(assigned_count) + "\t" + str(assigned_v_count) + "\t" + str(assigned_j_count) + "\t" + str(timed) + "\t" + str(Nseqs) + "\t" + str(average_total_quality) + "\t" + str(twenty_five) + "\t" + str(thirty) + "\t" + str(thirty_five) + "\t" + str(fourty) + "\n"
#         output_file = open(str(output_path), "a")
#         output.write(stringtoprint)
#         output_file.close()
    print str(stringtoprint)
        
    with open(output_path, "a") as myfile:
        myfile.write(stringtoprint)
def get_v_deletions(rc, v_match, temp_end_v, v_regions_cut):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] and str(v_regions_cut[v_match])[pos - 1] == str(rc)[function_temp_end_v - 1] and str(v_regions_cut[v_match])[pos - 2] == str(rc)[function_temp_end_v - 2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        return []

def get_j_deletions(rc, j_match, temp_start_j, j_regions_cut):
    # This function determines the number of J deletions in sequence rc
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in rc.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j + 2 < len(str(rc)):
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and str(j_regions_cut[j_match])[pos + 1] == str(rc)[function_temp_start_j + 1] and str(j_regions_cut[j_match])[pos + 2] == str(rc)[function_temp_start_j + 2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        return []

def get_v_tags(file_v, half_split):
    import string
    
    v_seqs = []  # Holds all V tags
    jump_to_end_v = []  # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n")  # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0])  # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1]))  # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = []  # Holds all J tags
    jump_to_start_j = []  # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]

def get_distinct_clones(handle, handle_results, with_count=False):

    # # LOOKS THROUGH TEXT FILE OF CLASSIFIERS AND WRITES NEW FILE CONTAINING ALL DISTINCT CLASSIFIERS, OPTIONALLY WITH COUNT OF ALL DISTINCT CLASSIFIERS
    # # with_count=True writes file with counts of all distinct classifiers
    
    from string import Template
    import collections as coll
    from operator import itemgetter, attrgetter

    write_to = open(str(handle_results) + '.txt', "w")

    if with_count == True:
        stemplate = Template('$count $element')
        d = coll.defaultdict(int)
        for line in handle:
            elements = line.rstrip("\n")
            d[elements] += 1
        d_sorted = sorted(d.items(), key=itemgetter(1), reverse=True)
        for k in d_sorted:
            kcount = k[1]
            f_seq = stemplate.substitute(count=str(kcount) + str(','), element=k[0])
            print >> write_to, f_seq
    else:
        stemplate = Template('$element')
        d = coll.defaultdict(int)
        for line in handle:
            elements = line.rstrip("\n")
            if elements not in d:
                d[elements] = 1
                f_seq = stemplate.substitute(element=elements)
                print >> write_to, f_seq
    
    handle.close()
    write_to.close()

def get_translated_sequences(handle, handle_results, chain, with_outframe=False, fullaaseq=False):

    # # TRANSLATES CLASSIFIERS TO AA SEQUENCES VIA THEIR NT SEQUENCE
    # # Default settings are -
    # # chain = "beta" or chain = "alpha"
    # # with_outframe=True or False: writes all aa seqeunces to file, including those that are out-of-frame (with stop codon symbol *)
    # # fullaaseq=True or False: True writes the whole V(D)J aa sequence to file, False, writes only the CDR3 region.

    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna
    import string
    import re

    handle_vb = open("human_TRBV_region.fasta", "rU")
    handle_jb = open("human_TRBJ_region.fasta", "rU")
    handle_va = open("human_TRAV_region.fasta", "rU")
    handle_ja = open("human_TRAJ_region.fasta", "rU")
    handle_vg = open("human_TRGV_region.fasta", "rU")
    handle_jg = open("human_TRGJ_region.fasta", "rU")
    handle_vd = open("human_TRDV_region.fasta", "rU")
    handle_jd = open("human_TRDJ_region.fasta", "rU")
    
    vb_raw = list(SeqIO.parse(handle_vb, "fasta"))
    handle_vb.close()
    jb_raw = list(SeqIO.parse(handle_jb, "fasta"))
    handle_jb.close()
    va_raw = list(SeqIO.parse(handle_va, "fasta"))
    handle_va.close()
    ja_raw = list(SeqIO.parse(handle_ja, "fasta"))
    handle_ja.close()
    vg_raw = list(SeqIO.parse(handle_vg, "fasta"))
    handle_vg.close()
    jg_raw = list(SeqIO.parse(handle_jg, "fasta"))
    handle_jg.close()
    vd_raw = list(SeqIO.parse(handle_vd, "fasta"))
    handle_vd.close()
    jd_raw = list(SeqIO.parse(handle_jd, "fasta"))
    handle_jd.close()

    vb_regions = []
    for i in range(0, len(vb_raw)):
        vb_regions.append(string.upper(vb_raw[i].seq))

    jb_regions = []
    for i in range(0, len(jb_raw)):
        jb_regions.append(string.upper(jb_raw[i].seq))

    va_regions = []
    for i in range(0, len(va_raw)):
        va_regions.append(string.upper(va_raw[i].seq))

    ja_regions = []
    for i in range(0, len(ja_raw)):
        ja_regions.append(string.upper(ja_raw[i].seq))

    vg_regions = []
    for i in range(0, len(vg_raw)):
        vg_regions.append(string.upper(vg_raw[i].seq))

    jg_regions = []
    for i in range(0, len(jg_raw)):
        jg_regions.append(string.upper(jg_raw[i].seq))

    vd_regions = []
    for i in range(0, len(vd_raw)):
        vd_regions.append(string.upper(vd_raw[i].seq))

    jd_regions = []
    for i in range(0, len(jd_raw)):
        jd_regions.append(string.upper(jd_raw[i].seq))

    write_to = open(str(handle_results) + '.txt', "w")

    if chain == "alpha":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 6:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ', ''))
            elif len(classifier) == 5:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''
                
            if delv != 0:
                used_v = va_regions[v][:-delv]
            elif delv == 0:
                used_v = va_regions[v]

            if delj != 0:
                used_j = ja_regions[j][delj:]
            elif delj == 0:
                used_j = ja_regions[j]

            seq = str(used_v + ins + used_j)
            start = (len(seq) - 1) % 3
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G', str(aaseq)) and re.findall('C', str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G', str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper + 4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "beta":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 6:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ', ''))
            elif len(classifier) == 5:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vb_regions[v][:-delv]
            elif delv == 0:
                used_v = vb_regions[v]

            if delj != 0:
                used_j = jb_regions[j][delj:]
            elif delj == 0:
                used_j = jb_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq) % 3 + 2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G', str(aaseq)) and re.findall('C', str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G', str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper + 4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "gamma":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 6:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ', ''))
            elif len(classifier) == 5:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vg_regions[v][:-delv]
            elif delv == 0:
                used_v = vg_regions[v]

            if delj != 0:
                used_j = jg_regions[j][delj:]
            elif delj == 0:
                used_j = jg_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq) % 3 + 2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G', str(aaseq)) and re.findall('C', str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G', str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper + 4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3

    if chain == "delta":
        for line in handle:
            elements = line.rstrip("\n")
            classifier = elements.split(',')

            if len(classifier) == 6:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = str(classifier[4].replace(' ', ''))
            elif len(classifier) == 5:
                v = int(classifier[0])
                j = int(classifier[1])
                delv = int(classifier[2])
                delj = int(classifier[3])
                ins = ''

            if delv != 0:
                used_v = vd_regions[v][:-delv]
            elif delv == 0:
                used_v = vd_regions[v]

            if delj != 0:
                used_j = jd_regions[j][delj:]
            elif delj == 0:
                used_j = jd_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq) % 3 + 2
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, elements + ', ' + str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, elements + ', ' + str(aaseq)
            else:     
                if re.findall('FG.G', str(aaseq)) and re.findall('C', str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G', str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper + 4]
                    if with_outframe == True:
                        print >> write_to, elements + ', ' + cdr3
                    elif '*' not in aaseq:
                        print >> write_to, elements + ', ' + cdr3
            
    handle.close()
    write_to.close()


def plot_v_usage(handle, chain, savefilename="Vusage", order="frequency"):

    # # PLOTS V GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain == "alpha":
        tags = open("tags_trav.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order == "frequency":        
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_v)
            percent_usage_v = [0] * num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1-1', 'V1-2', 'V10', 'V12-1', 'V12-2', 'V12-3', 'V13-1', 'V13-2', 'V14/D4', 'V16', 'V17', 'V18', 'V19', 'V2', 'V20', 'V21', 'V22', 'V23/D6', 'V24', 'V25', 'V26-1', 'V26-2', 'V27', 'V29/DV5', 'V3', 'V30', 'V34', 'V35', 'V36/DV7', 'V38-1', 'V38-2/DV8', 'V39', 'V4', 'V40', 'V41', 'V5', 'V6', 'V7', 'V8-1', 'V8-2/8-4', 'V8-3', 'V8-6', 'V9-1', 'V9-2', 'DV1', 'DV2', 'DV3')
            v_linked = [0] * len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0] * len(sorted_v)
            v_percents = [0] * len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.yticks(pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)
            
        elif order == "chromosome":
            total = sum(freq_vector_v)
            fv = [0] * num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1-1', 'V1-2', 'V10', 'V12-1', 'V12-2', 'V12-3', 'V13-1', 'V13-2', 'V14/D4', 'V16', 'V17', 'V18', 'V19', 'V2', 'V20', 'V21', 'V22', 'V23/D6', 'V24', 'V25', 'V26-1', 'V26-2', 'V27', 'V29/DV5', 'V3', 'V30', 'V34', 'V35', 'V36/DV7', 'V38-1', 'V38-2/DV8', 'V39', 'V4', 'V40', 'V41', 'V5', 'V6', 'V7', 'V8-1', 'V8-2/8-4', 'V8-3', 'V8-6', 'V9-1', 'V9-2', 'DV1', 'DV2', 'DV3')
            chromosome_order = [0, 1, 13, 24, 32, 35, 36, 37, 38, 42, 2, 3, 39, 40, 6, 4, 7, 8, 43, 5, 41, 9, 10, 11, 12, 14, 15, 16, 17, 44, 18, 19, 20, 22, 23, 25, 21, 26, 27, 28, 29, 30, 31, 33, 34, 45, 46]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize=8)
            ax.set_xticks(ind + 1 * width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "beta":
        tags = open("tags_trbv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order == "frequency":        
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_v)
            percent_usage_v = [0] * num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
            v_linked = [0] * len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0] * len(sorted_v)
            v_percents = [0] * len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.yticks(pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)
            
        elif order == "chromosome":
            total = sum(freq_vector_v)
            fv = [0] * num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
            chromosome_order = [ 14, 21, 23, 26, 31, 24, 25, 37, 32, 38, 44, 0, 3, 1, 4, 33, 39, 27, 34, 28, 40, 29, 35, 41, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22 ]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + 1 * width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "gamma":
        tags = open("tags_trgv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order == "frequency":        
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_v)
            percent_usage_v = [0] * num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V2', 'V3', 'V4', 'V5', 'V8', 'V9')
            v_linked = [0] * len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0] * len(sorted_v)
            v_percents = [0] * len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.yticks(pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)
            
        elif order == "chromosome":
            total = sum(freq_vector_v)
            fv = [0] * num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V2', 'V3', 'V4', 'V5', 'V8', 'V9')
            chromosome_order = [0, 1, 2, 3, 4, 5]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + 1 * width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "delta":
        tags = open("tags_trdv.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_v = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_v[int(elements.split(',')[0])] += 1

        if order == "frequency":        
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_v)
            percent_usage_v = [0] * num_genes
            for i in range(num_genes):
                percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')
            v_linked = [0] * len(percent_usage_v)
            for i in range(len(percent_usage_v)):
                v_linked[i] = (gene_list_v[i], percent_usage_v[i])
            sorted_v = sorted(v_linked, key=itemgetter(1))
            v_labels = [0] * len(sorted_v)
            v_percents = [0] * len(sorted_v)
            for j in range(len(sorted_v)):
                v_labels[j] = sorted_v[j][0]
                v_percents[j] = sorted_v[j][1]
            pos_v = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.yticks(pos_v, v_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_v, v_percents, align='center', color='yellow', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)

        elif order == "chromosome":
            total = sum(freq_vector_v)
            fv = [0] * num_genes
            for i in range(num_genes):
                fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
            gene_list_v = ('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')
            chromosome_order = [ 3, 5, 0, 4, 6, 7, 1, 2 ]
            gene_list_v = [ gene_list_v[i] for i in chromosome_order]
            fv = [fv[i] for i in chromosome_order]

            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fv, width, color='yellow')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + 1 * width)
            ax.set_xticklabels(gene_list_v)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()
    tags.close()

    
def plot_j_usage(handle, chain="beta", savefilename="Jusage", order="frequency"):

    # # PLOTS J GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    if chain == "alpha":
        
        tags = open("tags_traj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order == "frequency":
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_j)
            percent_usage_j = [0] * num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J20', 'J21', 'J22', 'J23', 'J24', 'J26', 'J27', 'J28', 'J29', 'J3', 'J30', 'J31', 'J32', 'J33', 'J34', 'J36', 'J37', 'J38', 'J39', 'J4', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J5', 'J50', 'J52', 'J53', 'J54', 'J56', 'J57', 'J6', 'J7', 'J8', 'J9')
            j_linked = [0] * len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0] * len(sorted_j)
            j_percents = [0] * len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.yticks(pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)

        elif order == "chromosome":
            total = sum(freq_vector_j)
            fj = [0] * num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J20', 'J21', 'J22', 'J23', 'J24', 'J26', 'J27', 'J28', 'J29', 'J3', 'J30', 'J31', 'J32', 'J33', 'J34', 'J36', 'J37', 'J38', 'J39', 'J4', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J5', 'J50', 'J52', 'J53', 'J54', 'J56', 'J57', 'J6', 'J7', 'J8', 'J9')
            chromosome_order = [45, 44, 43, 42, 41, 40, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 27, 26, 25, 24, 23, 22, 21, 20, 19, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 49, 48, 47, 46, 39, 28, 18]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "beta":
        
        tags = open("tags_trbj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order == "frequency":
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_j)
            percent_usage_j = [0] * num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
            j_linked = [0] * len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0] * len(sorted_j)
            j_percents = [0] * len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.yticks(pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)

        elif order == "chromosome":
            total = sum(freq_vector_j)
            fj = [0] * num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "gamma":
         
        tags = open("tags_trgj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order == "frequency":
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_j)
            percent_usage_j = [0] * num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1', 'J2', 'JP', 'JP1', 'JP2')
            j_linked = [0] * len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0] * len(sorted_j)
            j_percents = [0] * len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.yticks(pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)

        elif order == "chromosome":
            total = sum(freq_vector_j)
            fj = [0] * num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1', 'J2', 'JP', 'JP1', 'JP2')
            chromosome_order = [3, 2, 0, 4, 1]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "delta":
         
        tags = open("tags_trdj.txt", "rU")
        num_genes = 0
        for line in tags:
            num_genes += 1

        freq_vector_j = [0] * num_genes
        for line in handle:
            elements = line.rstrip("\n")
            freq_vector_j[int(elements.split(',')[1])] += 1
            
        if order == "frequency":
            plt.rcParams['figure.figsize'] = 10, 10
            total = sum(freq_vector_j)
            percent_usage_j = [0] * num_genes
            for i in range(num_genes):
                percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1', 'J2', 'J3', 'J4')
            j_linked = [0] * len(percent_usage_j)
            for i in range(len(percent_usage_j)):
                j_linked[i] = (gene_list_j[i], percent_usage_j[i])
            sorted_j = sorted(j_linked, key=itemgetter(1))
            j_labels = [0] * len(sorted_j)
            j_percents = [0] * len(sorted_j)
            for j in range(len(sorted_j)):
                j_labels[j] = sorted_j[j][0]
                j_percents[j] = sorted_j[j][1]
            pos_j = np.arange(num_genes) + 1
            plt.figure()
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.yticks(pos_j, j_labels)
            plt.xlabel('Frequency Usage')
            plt.barh(pos_j, j_percents, align='center', color='red', height=0.2)
            plt.grid(True)
            plt.savefig(str(savefilename) + '.png', dpi=300)

        elif order == "chromosome":
            total = sum(freq_vector_j)
            fj = [0] * num_genes
            for i in range(num_genes):
                fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
            gene_list_j = ('J1', 'J2', 'J3', 'J4')
            chromosome_order = [0, 3, 1, 2]
            gene_list_j = [ gene_list_j[i] for i in chromosome_order]
            fj = [fj[i] for i in chromosome_order]
            
            ind = np.arange(num_genes)
            width = 0.25

            fig = plt.figure()
            ax = fig.add_subplot(111)
            rects1 = ax.bar(ind, fj, width, color='red')

            ax.set_ylabel('Frequency', fontsize=10)
            ax.set_xticks(ind + width)
            ax.set_xticklabels(gene_list_j)
            plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=6)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=6)
            plt.grid(True)

            plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()
    tags.close()

def plot_del_v(handle, savefilename="Vdels"):

    # # PLOTS V GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_v = [0] * 50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_v[int(elements.split(',')[2])] += 1

    total = sum(deletions_v)
    for i in range(len(deletions_v)):
        deletions_v[i] = deletions_v[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_v, width, color='yellow')

    ax.set_ylabel('Frequency', fontsize=16)
    ax.set_xlabel('Number of V germline deletions', fontsize=16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=16)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=16)
    plt.ylim((0, 0.2))
    plt.xlim((0, 20))

    plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()

def plot_del_j(handle, savefilename="Jdels"):

    # # PLOTS J GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_j = [0] * 50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_j[int(elements.split(',')[3])] += 1

    total = sum(deletions_j)
    for i in range(len(deletions_j)):
        deletions_j[i] = deletions_j[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_j, width, color='red')

    ax.set_ylabel('Frequency', fontsize=16)
    ax.set_xlabel('Number of J germline deletions', fontsize=16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=16)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=16)
    plt.ylim((0, 0.2))
    plt.xlim((0, 20))

    plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()

def plot_vj_joint_dist(handle, chain="beta", savefilename="VJusage"):

    # # PLOTS VJ JOINT GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
        
    if chain == "alpha":
        
        tags_v = open("tags_trav.txt", "rU")
        tags_j = open("tags_traj.txt", "rU")

        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v, num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v, j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V1-1', 'V1-2', 'V10', 'V12-1', 'V12-2', 'V12-3', 'V13-1', 'V13-2', 'V14/D4', 'V16', 'V17', 'V18', 'V19', 'V2', 'V20', 'V21', 'V22', 'V23/D6', 'V24', 'V25', 'V26-1', 'V26-2', 'V27', 'V29/DV5', 'V3', 'V30', 'V34', 'V35', 'V36/DV7', 'V38-1', 'V38-2/DV8', 'V39', 'V4', 'V40', 'V41', 'V5', 'V6', 'V7', 'V8-1', 'V8-2/8-4', 'V8-3', 'V8-6', 'V9-1', 'V9-2', 'DV1', 'DV2', 'DV3')
        gene_list_j = ('J10', 'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J20', 'J21', 'J22', 'J23', 'J24', 'J26', 'J27', 'J28', 'J29', 'J3', 'J30', 'J31', 'J32', 'J33', 'J34', 'J36', 'J37', 'J38', 'J39', 'J4', 'J40', 'J41', 'J42', 'J43', 'J44', 'J45', 'J46', 'J47', 'J48', 'J49', 'J5', 'J50', 'J52', 'J53', 'J54', 'J56', 'J57', 'J6', 'J7', 'J8', 'J9')
            
        pos_v = np.arange(num_v) + 1
        pos_j = np.arange(num_j) + 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v - 0.5
        pos_ticks_j = pos_j - 0.5
        plt.yticks(pos_ticks_v, gene_list_v)
        plt.xticks(pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, rotation='vertical', fontsize='8')
        plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "beta":
        tags_v = open("tags_trbv.txt", "rU")
        tags_j = open("tags_trbj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v, num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v, j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
        gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')

        pos_v = np.arange(num_v) + 1
        pos_j = np.arange(num_j) + 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v - 0.5
        pos_ticks_j = pos_j - 0.5
        plt.yticks(pos_ticks_v, gene_list_v)
        plt.xticks(pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "gamma":
        tags_v = open("tags_trgv.txt", "rU")
        tags_j = open("tags_trgj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v, num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v, j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V2', 'V3', 'V4', 'V5', 'V8', 'V9')
        gene_list_j = ('J1', 'J2', 'JP', 'JP1', 'JP2')
        
        pos_v = np.arange(num_v) + 1
        pos_j = np.arange(num_j) + 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v - 0.5
        pos_ticks_j = pos_j - 0.5
        plt.yticks(pos_ticks_v, gene_list_v)
        plt.xticks(pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename) + '.png', dpi=300)

    if chain == "delta":
        tags_v = open("tags_trdv.txt", "rU")
        tags_j = open("tags_trdj.txt", "rU")
        
        num_v = 0
        for line in tags_v:
            num_v += 1

        num_j = 0
        for line in tags_j:
            num_j += 1
        
        joint_distribution = np.zeros((num_v, num_j))
        for line in handle:
            elements = line.rstrip("\n")

            v = int(elements.split(',')[0])
            j = int(elements.split(',')[1])

            joint_distribution[v, j] += 1

        joint_distribution = joint_distribution / sum(sum(joint_distribution))
        gene_list_v = ('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8')
        gene_list_j = ('J1', 'J2', 'J3', 'J4')
        
        pos_v = np.arange(num_v) + 1
        pos_j = np.arange(num_j) + 1
        
        plt.figure()
        plt.pcolor(joint_distribution)
        pos_ticks_v = pos_v - 0.5
        pos_ticks_j = pos_j - 0.5
        plt.yticks(pos_ticks_v, gene_list_v)
        plt.xticks(pos_ticks_j, gene_list_j)
        plt.colorbar()
        plt.pcolor(joint_distribution)
        yticklabels = plt.getp(plt.gca(), 'yticklabels')
        plt.setp(yticklabels, fontsize='8')
        xticklabels = plt.getp(plt.gca(), 'xticklabels')
        plt.setp(xticklabels, fontsize='8')
        plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()
    tags_v.close()
    tags_j.close()

def plot_insert_lengths(handle, savefilename="InsertLengths"):

    # # PLOTS DISTRIBUTION OF INSERT LENGTH BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    maxi = 500
    insert_lengths = [0] * maxi

    for line in handle:
        elements = line.rstrip("\n")

        classifier = elements.split(',')
        if len(classifier) == 6:
            insert_lengths[len(classifier[4].replace(' ', ''))] += 1
        else:
            insert_lengths[0] += 1

    total = sum(insert_lengths)
    for i in range(len(insert_lengths)):
        insert_lengths[i] = insert_lengths[i] / float(total)

    ind = np.arange(maxi)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, insert_lengths, width, color='blue')

    ax.set_ylabel('Frequency', fontsize=16)
    ax.set_xlabel('Number of nucleotides', fontsize=16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=16)
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=16)
    plt.xlim((0, 50))

    plt.savefig(str(savefilename) + '.png', dpi=300)

    handle.close()


# # Adaptation for running as an interactive command-line tool
    
print 'Welcome to DeCombinatoR_Edited, written by NICLAS THOMAS, edited by: ALVIN SHI'
filestoanalyze = []
resultslist = []
choose_file = raw_input('Enter path to the dir you wish to analyse. For example, enter C:/Users/Name/Desktop/')
chain = 'beta'
name = 'Full_Run_With_Qual_Scores'
overall_output = raw_input('Please enter the path to the outputfile, whose name will be PATH/total_output')
output = overall_output + 'Decombinator_output_run2.txt'
print output
rc = True

import os
import re
import platform
print platform.system()

os.chdir(choose_file)
for files in os.listdir("."):
    if files.endswith(".fastq"):
        filestoanalyze.append(files)
print(os.getcwd())
currentpath = os.getcwd()
print currentpath

for myfiles in filestoanalyze:
    toappend = re.sub('.fastq', '_output', myfiles)
    resultslist.append(toappend)

print str(filestoanalyze)
print str(resultslist)
    
if platform.system() == 'Windows':
    newpath = currentpath + '\\Results_' + str(name) + '\\'  # # Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)
elif platform.system() == 'Linux':
    newpath = currentpath + '/Results_' + str(choose_file) + '/'  # # Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)
elif platform.system() == 'Darwin':
    newpath = currentpath + '/Results_' + str(choose_file) + '/'  # # Ensure correct for specified platform
    if not os.path.exists(newpath):
        os.makedirs(newpath)
elif platform.system() == 'CYGWIN_NT-6.1': 
    newpath = currentpath + '/Results_' + str(name) + '/'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

# for i in range(len(ints)):
#    print i, ints[i]
#    

towritepath = newpath + str(output)
with open(towritepath, "a") as myfile:
    currentstringtoprint = "Filename\tTotalseqs\tAssignedseqs\tAssignedV\tAssignedJ\tTime\tNseqs\tAverage_Quality\tNum_Under_25\tNum_under_30\tNum_under_35\tNum_under_40\n"
    myfile.write(currentstringtoprint)

for idx, val in enumerate(filestoanalyze):
    print idx, val
    name_results = resultslist[idx]
    print ("Currently using this name on name_results:\n")
    print name_results
    analysis(val, towritepath, newpath + str(name_results), str(chain), with_statistics=True, with_reverse_complement_search=rc, omitN=True)

# seqs_found = 0
# seqcheck = open(newpath + str(name_results) + '.txt', "rU")
# for line in seqcheck:
#     seqs_found += 1
# seqcheck.close()
# 
# if seqs_found == 0:
#     print ' could not find any TcR', chain, 'sequences in the specified file.'
# else:
#     choose_plots = raw_input('Would you like to plot the results of your analysis? Enter (y/n): ')
# 
#     if choose_plots == 'n':
#         try:
#             choose_style = raw_input('Would you like to plot gene usage ordered by frequency or chromosomal position? Enter (frequency/chromosome): ')
#             print 'Plotting results of the analysis...'
#             plot_v_usage(open(newpath + str(name_results) + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + '_Vusage', order=str(choose_style))
#             plot_j_usage(open(newpath + str(name_results) + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + '_Jusage', order=str(choose_style))
#             plot_del_v(open(newpath + str(name_results) + '.txt', "rU"), savefilename=newpath + str(name_results) + '_Vdels')
#             plot_del_j(open(newpath + str(name_results) + '.txt', "rU"), savefilename=newpath + str(name_results) + '_Jdels')
#             plot_vj_joint_dist(open(newpath + str(name_results) + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + '_VJusage')
#             plot_insert_lengths(open(newpath + str(name_results) + '.txt', "rU"), savefilename=newpath + str(name_results) + '_InsertLengths')
#             print 'All plots successfully saved to: -'
#             print ''
#             print newpath
#             print ''
#         except:
#             print 'DeCombinatoR encountered an unexpected error while plotting your results.'
#             print 'If the problem persists, please contact niclas.thomas@gmail.com'
# 
#     choose_extras = raw_input('Would you like to use the extra functionality of DeCombinatoR to create additional files containing distinct clones and translated TcR sequences. Enter (y/n): ')
# 
#     if choose_extras == 'y':
#         try:   
#             choose_outframe = raw_input('Would you like to include out-of-frame TcR sequences in the file of translated TcR sequences? Enter (y/n): ')
#             choose_fullaaseq = raw_input('Would you like to translate the full TcR sequence or just the CDR3 region? Enter (full/cdr3): ')
#             choose_count = raw_input('Would you like to include the frequency of each distinct clone in the file of distinct clones? Note, to plot the features of distinct clones, you must select (n). Enter (y/n): ')
# 
#             if str(choose_outframe) == 'y':
#                 of = True
#             else:
#                 of = False
# 
#             if str(choose_fullaaseq) == 'full':
#                 fa = True
#             else:
#                 fa = False
# 
#             if str(choose_count) == 'y':
#                 c = True
#             else:
#                 c = False
#                 
#             get_distinct_clones(open(newpath + str(name_results) + '.txt', "rU"), handle_results=newpath + str('distinct_clones'), with_count=c)
#             get_translated_sequences(open(newpath + str(name_results) + '.txt', "rU"), handle_results=newpath + str('translated_sequences'), chain=str(chain), with_outframe=of, fullaaseq=fa)
#             print 'The requested file(s) have been successfully saved to: -'
#             print ''
#             print newpath
#             print ''
# 
#             if choose_count == 'n':
#                 
#                 print 'Plotting features of the distinct clones...'
#                 plot_v_usage(open(newpath + str('distinct_clones') + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + 'distinctclones_Vusage', order=str(choose_style))
#                 plot_j_usage(open(newpath + str('distinct_clones') + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + 'distinctclones_Jusage', order=str(choose_style))
#                 plot_del_v(open(newpath + str('distinct_clones') + '.txt', "rU"), savefilename=newpath + str(name_results) + 'distinctclones_Vdels')
#                 plot_del_j(open(newpath + str('distinct_clones') + '.txt', "rU"), savefilename=newpath + str(name_results) + 'distinctclones_Jdels')
#                 plot_vj_joint_dist(open(newpath + str('distinct_clones') + '.txt', "rU"), chain=str(chain), savefilename=newpath + str(name_results) + 'distinctclones_VJusage')
#                 plot_insert_lengths(open(newpath + str('distinct_clones') + '.txt', "rU"), savefilename=newpath + str(name_results) + 'distinctclones_InsertLengths')
#                 print 'All plots successfully saved to: -'
#                 print ''
#                 print newpath
#                 print ''
# 
#         except:
#             print 'DeCombinatoR encountered an unexpected error while using its extra functionality.'
#             print 'If the problem persists, please contact niclas.thomas@gmail.com'
        
raw_input('Press Enter to exit...')
