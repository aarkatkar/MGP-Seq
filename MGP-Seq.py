#########################################################################################
# Filename: 	MGP-Seq.py
# Written by: 	Peter Dornbos
# 				Graduate Assistant (LaPres Lab)
# 				Integrative Toxicology Program
# 				Department of Biochemistry & Molecular Biology
# 				Michigan State University
#               Anooj Arkatkar
#                               Undergraduate Assistant (LaPres Lab)
#                               Department of Biochemistry & Molecular Biology
#                               Michigan State University
# Prepared: 	July of 2019
#
# Purpose: Automated Gene Sequence Imputation Software
#########################################################################################

# The name of the Tkinter library is changed from Python 2 to Python 3
try:
    import Tkinter
except ImportError:
    # tkinter is used for GUI elements
    import tkinter as Tkinter
    from tkinter import messagebox as tkMessageBox
    from tkinter import ttk
    from tkinter import filedialog as tkFileDialog
    # urllib is used to connect to Sanger to download csv files
    from urllib.request import urlopen
    from urllib.parse import urlencode
else:
    import tkMessageBox
    import tkFileDialog
    import ttk
    from urllib import urlopen
    from urllib import urlencode
finally:
    # Used to load SNP and indel data
    import csv
    # Used to create hyperlinks
    import webbrowser
    # Used to find files in directory
    import os
    # Used to set delays when connecting to NCBI gene
    import time
    # Used to interpret files downloaded directly from Sanger
    import codecs
    # Used to parse strings more easily
    import re
    # Used to parse NCBI database searches
    import xml.etree.ElementTree as ET

# All outputs and downloads will go here (global variable)
output_directory = '.'
email = 'arkatkar@msu.edu'
tool = 'MGPSeq'

def sangerDownload(tag, gene, form, version='rel-1505'):
    # Form is 'SNP' or 'Indel'

    # sv, sn, and st are binary numbers that encode search criteria
    # sv = 31 (1111 in binary) searches for all 4 available structural variants
    query = {'sv': 31,  # Structural Variants
             'sn': 8589410303,  # SNP/Indel Type
             'st': 68719476735, # Strains
             'loc': tag,    # Locus
             'format': 'csv'}   # Output File
    
    url = 'https://www.sanger.ac.uk/sanger/Mouse_SnpViewer_Export/' +\
          '{}s/{}?'.format(form,version) + urlencode(query)
    
    f = urlopen(url)
    table = csv.reader(codecs.iterdecode(f, 'utf-8'))
    filename = os.path.join(output_directory,
                            'MouseGP.{}_{}s.csv'.format(gene, form))
    with open(filename, 'w') as f:
        outfile = csv.writer(f, lineterminator='\n')
        for row in table:
            outfile.writerow(row)
    f.close()
    tkMessageBox.showinfo("Download Successful",
                         "Downloaded {}s as {}".format(form, filename))
    return filename

gene_term = 'Mus Musculus[Orgn] AND {}[Gene]'


def esearch(gene, db='gene', wait=True):
    query = {'db': db,
             'term': gene_term.format(gene),
             'email': email,
             'tool': tool}
    if db == 'nucleotide':
        query['idtype'] = 'acc'

    # wait 0.37s before querying to avoid going over usage restrictions
    if wait:
        time.sleep(0.37)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' +\
          urlencode(query)

    with urlopen(url) as f:
        root = ET.fromstring(f.read())
        id_list = root[3]
        ids = ''
        for child in id_list:
            ids += child.text + ','
    return ids


def esummary(db, id_):    
    query = {'db': db,
             'id': id_,
             'email': email,
             'tool': tool}
    
    time.sleep(0.37)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' +\
    urlencode(query)
    table = []
    
    with urlopen(url) as f:
        root = ET.fromstring(f.read())
    return root


def efetch(filename, **query):
    time.sleep(0.37)
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' +\
    urlencode(query)

    with urlopen(url) as f:
        record = f.read().decode('utf-8')
        
    with open(filename, 'w') as outfile:
        outfile.write(record.rstrip('\n'))
        

def DNA_parse(root):
    table = []
    for child in root[0]:
        data = {}
        try:
            x = child.find('GenomicInfo')[0]
        except (IndexError, TypeError):
            continue
        for item in x:
            data[item.tag] = item.text
        name = child.find('Name').text
        description = child.find('Description').text
        start = int(data['ChrStart'])+1
        end = int(data['ChrStop'])+1
        acc = data['ChrAccVer']
        chromosome = data['ChrLoc']
        tag = chromosome + ':' + str(start-1) + '-' + str(end+1)
        results = name + ':\n' + description + '\n' + tag
        table.append((acc, tag, results, name, start, end))
    return table


def RNA_parse(summary, gene):
    d = {}
    for transcript in summary:
        key = ''
        acc = ''
        for child in transcript:
            if child.attrib.get('Name') == 'Title':
                title = child.text
                transcript_variant = title.split(',')[-2].strip()
                try:
                    gene = re.search('(\((.*?)\))', title).group(0)
                except AttributeError:
                    gene = ''
                if '(' in transcript_variant:
                    gene = ''
                key = transcript_variant + ' ' + gene
            elif child.attrib.get('Name') == 'Caption':
                acc = child.text
        d[key.strip()] = acc
    return d


def searchDNA(gene):
    ids = esearch(gene)
    summary = esummary('gene', ids)
    return DNA_parse(summary)


def searchRNA(gene):
    all_ids = esearch(gene, db='nucleotide')
    ids = ''
    for id_ in all_ids.split(','):
        if id_[:2] in ('NM', 'XM'):
            ids += id_ + ','
    summary = esummary('nucleotide', id_=ids)
    return RNA_parse(summary, gene)


def downloadDNA(acc, gene, start, end):
    filename = os.path.join(output_directory,
                            "NCBI_Gene.{}.fasta".format(gene))

    efetch(filename, db='nucleotide', id=acc, seq_start=start,
           seq_stop=end, rettype='fasta', retmode='text', strand=2,
           email=email, tool=tool)
    
    return filename


def downloadRNA(acc, gene):
    filename = os.path.join(output_directory,
                            "NCBI_RNA.{}.fasta".format(gene))
    efetch(filename, db='nucleotide', id=acc, rettype='fasta',
           retmode='text', email=email, tool=tool)
    return filename
        
base_dict = {"G": "C",
             "C":"G",
             "A":"T",
             "T":"A"}

bases = ['T', 'C', 'A', 'G']

def reverse_complement(seq):
    return [base_dict[b] for b in reversed(seq)]

class Sequence:
    def __init__(self, original_dna='', loc=0,
                 file=''):
        self.moving_seq = []
        self.original_dna = original_dna
        self.loc = loc
        self.alignment = Alignment(None, None)
        self.orientation = []
        self.imp_flip = False
        
        if file:
            self.load(file)
            
    def load(self, file):
        with open(file, 'r') as f:
            temp_seq = []
            string = f.read()
            original_dna = ''
            # First, we check for a FASTA header
            if string[0] == '#' or '|' in string:
                raise TypeError
            if string[0] != '>':
                for c in string:
                    if c in bases:
                        temp_seq.append(c)
                self.moving_seq = list(temp_seq)
                return (len(self.moving_seq), self.loc)

            # Next, we store the dna to self.dna
            sep = string.index('\n')
            header = string[:sep]
            for c in string[sep + 1:]:
                if c in bases:
                    temp_seq.append(c)
            self.moving_seq = list(temp_seq)

            # We try to parse the header for the starting locus
            indices = []
            for i, char in enumerate(header):
                if char in [":", "-", " "]:
                    indices.append(i)
            try:
                self.loc = int(header[indices[0]+1 : indices[1]])
                
            except ValueError:
                try:
                    self.loc = int(header[indices[1]+1 : indices[2]])
                except ValueError:
                    pass
                
            except IndexError:
                pass
            
            finally:
                return (len(self.moving_seq), self.loc)
        
    def transcribe(self):
        dna = self.moving_seq
        exons = self.alignment.exons_sorted
        spliced = []
        length = len(dna)
        if self.alignment.flip:
            for i in reversed(exons):
                if i < length:
                    temp = dna[i]
                    for bp in reversed(temp):
                        spliced.append(base_dict[bp])
        else:
            for i in exons:
                if i < length:
                    spliced.append(dna[i])
        self.moving_seq = spliced
    
    @property
    def peptide(self):
        seq = self.sequence
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))

        start_codon = seq.index('ATG')
        peptide = ''
        for t in range(start_codon, len(seq), 3):
            codon = seq[t:t + 3]
            amino_acid = codon_table.get(codon, '*')
            if amino_acid == '*':
                break
            peptide += amino_acid
        return peptide

    @property
    def peptide_longest(self):
        seq = self.sequence
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))

        peptides = []
        for frame in range(3):
            peptide = ''
            for t in range(frame, len(seq), 3):
                codon = seq[t:t + 3]
                if codon != 'ATG' and len(peptide) == 0:
                    continue
                amino_acid = codon_table.get(codon, '*')
                if amino_acid == '*':
                    peptides.append(peptide)
                    peptide = ''
                else:
                    peptide += amino_acid
                
        peptides.sort(key=lambda p: len(p), reverse=True)
        return peptides[0]

    def impute(self, imputations, snp_col, indel_col):
        snps = imputations.snps[snp_col]
        indels = imputations.indels[indel_col]
        temp_seq = list(self.moving_seq)
        exons = set(self.alignment.exons)
        for position in snps:
            if position < self.loc:
                continue
            try:
                temp_seq[position - self.loc] = snps[position]
            except IndexError:
                pass
            
        for position in indels:
            indel = indels[position][0]
            if position < self.loc:
                continue
            b_pos = position - self.loc
            b_range = indels[position][1]
            if bool(b_pos in exons) != bool(b_pos + b_range in exons):
                sign = bool(b_pos in exons)
                for mv in range(1, b_range):
                    if bool(mv+b_pos in exons) != sign:
                        if sign:
                            indel = indel[:mv]
                        else:
                            indel = indel[mv:]
                            b_pos += mv
                            b_range -= mv
                        break
            try:
                if b_pos >= 1:
                    temp_seq[b_pos] = indel
                    for extra in range(b_range):
                        temp_seq[b_pos + extra + 1] = ''
            except IndexError:
                pass
        self.moving_seq = list(temp_seq)

    @property
    def sequence(self):
        return ''.join(self.moving_seq)

    def deorient(self):
        if not self.imp_flip:
            return None
        
        self.moving_seq = reverse_complement(self.sequence)
    
    def orient_imputations(self, imputations):
        reference = imputations.reference
        retain = 0
        flip = 0
        for position in reference:
            gene_loc = position - self.loc
            if gene_loc < 0 or gene_loc >= len(self.moving_seq):
                continue
            predicted_base = reference[position]
            if self.moving_seq[gene_loc] == predicted_base:
                retain += 1
            if self.moving_seq[-1*gene_loc-1] == base_dict[predicted_base]:
                flip += 1

        if flip > retain:
            self.moving_seq = reverse_complement(self.moving_seq)
            self.imp_flip = True
            
        self.orientation = list(self.moving_seq)

    def reorient(self):
        self.moving_seq = list(self.orientation)
        
    def align(self, mRNA, k=12):
        self.alignment = Alignment(mRNA.sequence, self.sequence)
        self.alignment.align(k=k)
        self.imp_flip = False
            

class Alignment:
    def __init__(self, transcript=None, gene=None, emboss='', threshold=30):
        self.flip = False
        self.exons = set()
        self.exons_sorted = []
        self.transcript = transcript
        self.gene = gene
        if emboss:
            self.extractEMBOSS(emboss, threshold)
        
    def align(self, k=12):
        a = DB_align(self.transcript, self.gene, k=k)
        self.exons = a[0]
        self.flip = a[1]
        self.exons_sorted = sorted(self.exons)

class Imputations:
    def __init__(self):
        self.snps = {False:{}}
        self.indels = {False:{}}
        self.reference = {}
    
    def load_snps(self, source, snp_cols, position_col=1, reference_col=4,
                 confidence=False):
        with open(source, 'rU') as snp_source:
            data = csv.reader(snp_source, delimiter=',')
            data = list(data)
            for n, snp_col in enumerate(snp_cols):
                snps = {}
                if snp_col == False:
                    continue
                for i, row in enumerate(data):
                    if i == 0:
                        name = row[snp_col].replace("/", "_")
                        continue
                    try:
                        try:
                            snp = row[snp_col][0]
                            if n == 0:
                                reference = row[reference_col][0].upper()
                                self.reference[int(row[position_col])] = reference
                        except IndexError:
                            snp=''
                        if confidence and not snp.isupper():
                            continue
                        if snp.upper() in bases:
                            snps[int(row[position_col])] = snp.upper()
                    except ValueError:
                        pass
                self.snps[snp_col] = snps
            
    def load_indels(self, source, indel_cols, position_col=1, reference_col=4,
                 confidence=False):
        with open(source, 'rU') as snp_source:
            data = csv.reader(snp_source, delimiter=',')
            data = list(data)
            for n, indel_col in enumerate(indel_cols):
                indels = {}
                if indel_col == False:
                    continue
                for i, row in enumerate(data):
                    if i == 0:
                        name = row[indel_col].replace("/", "_")
                        continue
                            
                    indel = ''
                    for base in row[indel_col]:
                        if not base.upper() in bases:
                            break
                        indel += base
                    if confidence and not indel.isupper():
                        continue
                    if indel != '':
                        if n == 0:
                            reference = row[reference_col][0].upper()
                            self.reference[int(row[position_col])] = reference
                        indel_data = (indel.upper(), len(row[reference_col])-1)
                        indels[int(row[position_col])] = indel_data
                self.indels[indel_col] = indels


# A sorted list class
class priority_queue:
    def __init__(self):
        self.PQ = []
        self.meta = {}
        self.finished = set()
    def add_multiple(self, array):
        for a in array:
            self.PQ.append(a[0])
            self.meta[a[0]] = (a[1], a[2])
    def dequeue(self):
        temp = self.PQ.pop(-1)
        self.finished.add(temp)
        return temp
    def update(self, p, v):
        self.meta[p] = v
    def sort(self):
        self.PQ.sort(key=lambda x: self.meta[x][0], reverse=True)
    def info(self, node):
        if node in self.finished:
            return False
        else:
            return self.meta[node]


# Dijkstra's Shortest Path 
def path(start, end, connections):
    # Organize the connections into vertices(nodes) and edges(dictionary)
    nodes = set()
    for a, b,c in connections:
        nodes.add(a)
        nodes.add(b)
        
    d = dict(zip(nodes, ([] for _ in nodes)))
    for connection in connections:
        d[connection[0]].append((connection[1], connection[2]))
        
    dist = zip(nodes, (float('inf') for _ in nodes), (False for _ in nodes))
    PQ = priority_queue()
    PQ.add_multiple(dist)
    PQ.update(start, (0, False))
    PQ.sort()

    outpath = []
    while True:
        t = PQ.dequeue()
        current_distance = PQ.meta[t][0]
        if t == end:
            total = current_distance
            break
        
        for node in d[t]:
            old_distance = PQ.info(node[0])
            if old_distance is False:
                continue
            new_distance = current_distance + node[1]
            if new_distance < old_distance[0]:
                PQ.update(node[0], (new_distance, t))
        PQ.sort()

    visited = set()
    check = end
    while True:
        check = PQ.meta[check][1]
        if check in visited:
            print('Could not find alignment')
            return (None, None)
        visited.add(check)
        if check == start:
            break
        outpath.append(check)
    return (list(reversed(outpath)), total)

def needleman_wunsch(s1, s2):
    '''implementation of NW with affine gap penalty'''
    match = 5 #from DNAfull 
    mismatch = -4 #from DNAfull
    g = -4 #gap extension cost
    h = -16 #gap open cost
    extension_bonus = 2 #preference to alignments at beginning or end
    
    x = [c for c in s1]
    y = [c for c in s2]
    s = {True: match,
         False: mismatch}

    m = len(x) + 1
    n = len(y) + 1
    
    # Create three NW arrays (for affine gap penalty)
    M = [[(float('-inf'), ()) for _ in range(n)] for _ in range(m)]
    Ix = [[(float('-inf'), ()) for _ in range(n)] for _ in range(m)]
    Iy = [[(float('-inf'), ()) for _ in range(n)] for _ in range(m)]
    
    # hint: A -> array[column][row] (= [x][y])
    
    M[0][0] = (0, 0)
    for i in range(m):
        Iy[i][0] = (h + g*i, 2)
    for j in range(n):
        Ix[0][j] = (h + g*j, 1)

    # Fill in the NW Array with indel/mismatch = -1, match=+1
    for j in range(1, n):
        for i in range(1, m):
            # adjust M matrix
            match_score = s[x[i-1] == y[j-1]]
            v = M[i-1][j-1][0] + match_score
            d = 0
            if i == j or n-j == m-i:
                v += extension_bonus
            insertion_x = Ix[i-1][j-1][0] + match_score
            insertion_y = Iy[i-1][j-1][0] + match_score
            if insertion_y > v:
                v = insertion_y
                d = 2
            if insertion_x > v:
                v = insertion_x
                d = 1
            M[i][j] = (v,d)

           # adjust Ix matrix
            v = M[i][j-1][0] + h + g
            d = 0
            insertion_x = Ix[i][j-1][0] + g
            if insertion_x > v:
                v = insertion_x
                d = 1
            Ix[i][j] = (v,d)
            
            # adjust Iy matrix
            v = M[i-1][j][0] + h + g
            d = 0
            insertion_y = Iy[i-1][j][0] + g
            if insertion_y > v:
                v = insertion_y
                d = 2
            Iy[i][j] = (v,d)
            
    # Retrace back to beginning
    m -= 1
    n -= 1
    exons = []
    insertions = []
    
    c_val = max(M[m][n], Ix[m][n], Iy[m][n])
    if c_val == M[m][n]:
        cell = (0, m, n)
    elif c_val == Iy[m][n]:
        cell = (2, m, n)
    else:
        cell = (1, m, n)
    while cell[1:3] != (0, 0):
        pointer = cell[0]
        i = cell[1]
        j = cell[2]
        if pointer == 0:
            cell = (M[i][j][1], i-1, j-1)
            exons.append(j-1)
        elif pointer == 1:
            cell = (Ix[i][j][1], i, j-1)
        else:
            cell = (Iy[i][j][1], i-1, j)
            insertions.append(j-1)
    return exons, insertions


def DB_align(s1, s2, k=12, intron_cost=1, insertion_cost=3):
    '''Predicts each index of gDNA that will remain in cDNA'''
    
    # Begin with a simplified BLAT search for matching sequences
    # The query contains overlapping kmers
    query = []
    for i in range(len(s1)-k+1):
        query.append((s1[i:i+k], i))

    # The search contains non-overlapping kmers
    search = {}
    for i in range(0, len(s2)-k, k):
        temp = s2[i:i+k]

        if temp in search:
            search[temp].append(i)
        else:
            search[temp] = [i]
            
    # matches will be compared against {covered} to see if they
    # are already encompassed by a bigger match
    covered = set()
    
    for q in query:
        f = search.get(q[0])
        if f is None:
            continue

        # here we combine two BLAT steps
        # 1) we extend each match in both directions as far as they will go
        # 2) we pick the largest match if multiple are found
        best_match = (0, 0, 0)
        for match in f:
            c_start = q[1]
            c_end = q[1] + k

            g_start = match
            g_end = match + k
            for c in covered:
                # this checks if the match is already covered
                if c_start < c[1] and g_start - c_start == c[2]:
                    if c[1]-c[0] > best_match[1]-best_match[0]:
                        best_match = c
                    break
            else:
                # extend each match backwards
                while True:
                    if c_start <= 0 or g_start <= 0:
                        break
                    if s1[c_start-1] != s2[g_start-1]:
                        break
                    c_start -= 1
                    g_start -= 1

                # extend each match forwards
                while True:
                    try:
                        if s1[c_end + 1] != s2[g_end + 1]:
                            break
                    except IndexError:
                        break
                    c_end += 1
                    g_end += 1

                # replace "best_match" if this match is longer
                if c_end - c_start > best_match[1] - best_match[0]:
                    best_match = (c_start, c_end, g_start - c_start)
                    
        covered.add(best_match)

    # This block tells you the percentage of cDNA base pairs that match somewhere
    # Ideally, this number should be 1.0, but could be smaller
    # in the case of a PolyA tail, sequencing errors, etc.
    accounted = set()
    seeds = []
    for c in covered:
        seeds.append(c)
        for base in range(c[0], c[1]+1):
            accounted.add(base)
    ratio = len(accounted)/len(s1)

    convert = False
    # If less than 90% of the bases have been mapped, the sequence might be
    # anti-sense rather than sense. Try a flipped alignment
    if ratio < 0.9:
        check = check_complementary(s1, s2, ratio, search, k=k)
        if check[0]:
            print('flipping sequence')
            convert = True
            # I have gotten memory overload errors if I don't clear these
            # before starting up the complementary alignment
            seeds.clear()
            search.clear()
            query.clear()
            accounted.clear()

            # might as well take the exon seeds from check_complementary
            # rather than starting over again
            # the same is done with s1's reverse complement
            seeds = check[1]
            s1 = check[2]

    # Some helpful sorting here can narrow our alignment search space
    # since we know we are not going backwards on the cDNA strand
    seeds.sort(reverse=True)

    # connections <- [(exon1, exon2, linking_cost)...]
    connections = []
    final_index = len(s1)-1
    for j in range(len(seeds)):
        cDNA_start = seeds[j][0]
        cDNA_end = seeds[j][1]
        b = seeds[j][2]
        
        # Add connection to start
        connections.append(('e', j, final_index - cDNA_end)) 
        # Add connection to end
        connections.append((j, 's', cDNA_start))

        # See how well the exon connects to all other exon seeds
        for i in range(j+1, len(seeds)):
            i_cDNA_end = seeds[i][1]
            i_b = seeds[i][2]
            # i_b > b suggests an insertion into mRNA (probably polyadenylated)
            if i_b > b:
                gDNA_start = cDNA_start + b
                i_gDNA_start = seeds[i][0] + i_b
                if gDNA_start < i_gDNA_start:
                    continue
                else:
                    cost = cDNA_start - i_cDNA_end
                    if cost < insertion_cost:
                        cost = insertion_cost
                    connections.append((j, i, cost))
                    
            # this is a standard intron-exon boundary
            else:
                cost = cDNA_start - i_cDNA_end
                if cost < intron_cost:
                    cost = intron_cost
                connections.append((j, i, cost))

    # Run Dijkstra's Shortest Path to minimize alignment cost
    p = path('e', 's', connections)
    
    # Now, find optimal alignment given the observed path
    # Possible connections:
    # 1) Intron, 2) substitution, 3) insertion, 4) large gap
    # Intron - prefer splice signal
    # substitution - treat like one long exon
    # insertion - ignore
    # large gap - NW
    end_mimic = (len(s1), len(s1), len(s2)-len(s1))
    start_mimic = (0, 0, 0)
    exon_path = [end_mimic] + [seeds[k] for k in p[0]] + [start_mimic]
   
    exons = []
    end_cutoff = len(s1)-1

    # (index just before insertion, insertion bases)
    insertions = []
    for i in range(len(exon_path)-1):
        pre = exon_path[i]
        post = exon_path[i+1]
        if post[2] == pre[2]:
            # simple substitution - make into one long exon
            for k in range(post[1], end_cutoff+1):
                exons.append(k + pre[2])
            end_cutoff = post[1]-1
        elif post[1] >= pre[0]-1 and post[2] < pre[2]:
            # Intron connection
            if pre[0] == len(s1)-1:
                continue
            if post[1] == 0:
                for k in range(end_cutoff+1):
                    exons.append(k+pre[2])
                break
            
            for c in range(pre[0]-1, post[1]+1):
                intron_start = s2[c+1+post[2]:c+3+post[2]]
                intron_end = s2[c-1+pre[2]:c +1+ pre[2]]
                # attempt to set intron boundaries to match consensus
                if intron_start + intron_end in ("GTAG", "CTAC"):
                    break
            for k in range(c+1, end_cutoff+1):
                exons.append(k+pre[2])

            end_cutoff = c
        
        elif post[1] + post[2] >= pre[0] + pre[2]-1:
            # insertion connection - very rare; usually from experimental error
            insertions.append((pre[0]+pre[2], post[2]-pre[2]))
            for k in range(pre[0]-1, end_cutoff+1):
                exons.append(k + pre[2])
            end_cutoff = pre[0]-2
            
        else:
            # large gap - use NW to make alignment
            gap_start = post[1]+1
            if gap_start == 1:
                gap_start = 0
            gap_end = pre[0]
            align1 = s1[gap_start:gap_end]
            align2 = s2[gap_start+post[2]:gap_end+pre[2]]
            greedy_extension = needleman_wunsch(align1, align2)
            for g in greedy_extension[0]:
                exons.append(g + gap_start + post[2])
            if pre[1] != len(s1)-1: 
                for k in range(pre[0],end_cutoff + 1):
                    exons.append(k+pre[2])
            end_cutoff = post[1]
            
            if len(greedy_extension[1]) == 0:
                continue
            prev_i = greedy_extension[1][0]+gap_start
            prev_s = 1
            
            for k in greedy_extension[1][1:]:
                i_index = gap_start + k
                if i_index != prev_i - 1:
                    insertions.append((prev_i, prev_s))
                prev_i = i_index
                prev_s += 1
                
            insertions.append((prev_i, prev_s))

    return (exons, convert)

# check_complementary is near-identical to the first half of DB_align
def check_complementary(s1, s2, ratio, search, k=12):
    s1 = ''.join(reverese_complement(s1))
            
    query = []
    for i in range(len(s1)-k+1):
        query.append((s1[i:i+k], i))

    covered = set()
    for q in query:
        f = search.get(q[0])
        if f is None:
            continue
        best_match = (0, 0, 0)
        for match in f:
            c_start = q[1]
            c_end = q[1] + k

            g_start = match
            g_end = match + k
            for c in covered:
                if c_start < c[1] and g_start - c_start == c[2]:
                    if c[1]-c[0] > best_match[1]-best_match[0]:
                        best_match = c
                    break
            else:
                while True:
                    if c_start <= 0 or g_start <= 0:
                        break
                    if s1[c_start-1] != s2[g_start-1]:
                        break
                    c_start -= 1
                    g_start -= 1
                while True:
                    try:
                        if s1[c_end + 1] != s2[g_end + 1]:
                            break
                    except IndexError:
                        break
                    c_end += 1
                    g_end += 1
                if c_end - c_start > best_match[1] - best_match[0]:
                    best_match = (c_start, c_end, g_start - c_start)
        covered.add(best_match)
        
    accounted = set()
    seeds = []
    for c in covered:
        seeds.append(c)
        for base in range(c[0], c[1]+1):
            accounted.add(base)
            
    r = len(accounted)/len(s1)
    accounted.clear()
    covered.clear()
    search.clear()
    query.clear()
    if r > ratio:
        return (True, seeds, s1)
    else:
        return (False, [], '')

"""Error messages of various types"""
# Some programmers may cringe at this unnecessary OOP, so here it remains
class References:
    def __init__(self):
        pass
error = References()
# MISSING FIELD
error.missing_field = "Missing Field"
error.missing_alignment = "Please select an alignment file"
error.missing_loc = "Please enter a numeric gene location"
error.missing_ref = "Please select a reference gene"
error.no_ref = "No Reference Gene"

# redText warnings
warning = References()
warning.no_header = """No header found.
Please manually enter the gene location"""

# ALIGNMENT ERRORS
error.alignment = "Exon Alignment Error"
error.no_path = "Failed to align mRNA to specified gDNA\n" +\
              "Try submitting an EMBOSS alignment file instead"

# A definition is used so that the GUI can call the script at any time
def impute(s_data, i_data, ref_data, loc, separator,
           strain_set, confidence, output_key, threshold, alignment, gene_name,
           p_type):
    start_time = time.time()
    print('starting imputation')
    # Setup for all of the file-saving
    if not output_key:
            output_type = "gene"
    elif output_key == 1:
        output_type = "mRNA"
    else:
        output_type = "protein"

    if gene_name != "":
        gene_name = "_" + gene_name
    
    filename = output_directory
    
    gene = Sequence(file=ref_data, loc=loc)

    if output_key:
        mRNA = Sequence(file=alignment)

        if len(strain_set)==0:
            # For testing, allows the code to work even without SNPs or Indels
            gene.align(mRNA, k=threshold)
            
    # Load in the Imputations
    imputations = Imputations()
    snp_cols = [column[1] for column in strain_set]
    indel_cols = [column[2] for column in strain_set]
    if s_data:
        imputations.load_snps(s_data, snp_cols=snp_cols, confidence=confidence)
    if i_data:
        imputations.load_indels(i_data, indel_cols=indel_cols, confidence=confidence)
    
    for i, j in enumerate(strain_set):
        strain = j[0].replace('/', '_')
        # These (slow) steps only need to happen the first time
        if i == 0:
            if s_data or i_data:
                # Flip the sequence if it doesn't match the imputation files
                gene.orient_imputations(imputations)
            # Add an alignment to our gene for transcription
            if output_key:
                gene.align(mRNA, k=threshold)

##        print(strain, gene.sequence[:10])
        gene.impute(imputations, j[1], j[2])

        if output_key:
            gene.transcribe()
        else:
            gene.deorient()

        if output_key == 2:
            if p_type == 0:
                output = gene.peptide
            else:
                output = gene.peptide_longest
        else:
            output = gene.sequence
        # Give each strain a unique file name
        text_name = strain + "_" + output_type + gene_name
        output_dir = os.path.join(filename, text_name + ".txt")

        with open(output_dir, "w") as outfile:
            outfile.write(">{}\n".format(strain) + output)
            
        gene.reorient()
##        print(strain, gene.sequence[:10])
    # Repeat the steps for the C57BL_6J reference output
    if output_key:
        gene.transcribe()
    else:
        gene.deorient()
    if output_key == 2:
        if p_type == 0:   
            output = gene.peptide
        else:
            output = gene.peptide_longest
    else:
        output = gene.sequence
    
    text_name = "C57BL_6J_" + output_type + gene_name
    ref_output = os.path.join(filename, text_name + ".txt")
    with open(ref_output, "w") as r_out:
        r_out.write(">C57BL_6J\n" + output)
    
    print("Process Completed")
    print("Output directory: " + filename)
    file_no = len(strain_set)
    print("Completed in {:.2f}s".format(time.time()-start_time))
    if file_no == 1:
        tkMessageBox.showinfo("Process Completed",
                              "{} file was imputed".format(file_no))
    else:
        tkMessageBox.showinfo("Process Completed",
                              "{} files were imputed".format(file_no))


def call_sanger(event):
    webbrowser.open_new(r"https://www.sanger.ac.uk/sanger/Mouse_SnpViewer/")
    

def call_reference(event):
    webbrowser.open_new(r"https://www.ncbi.nlm.nih.gov/gene")

def call_EMBOSS(event):
    webbrowser.open_new(r"https://www.ebi.ac.uk/Tools/psa/emboss_stretcher/")

# Lists files in directory for the drop-down menus
# Bump up files that are more likely to the top of the list based on name
fasta_files = []
txt_files = []
snp_preferred = []
snp_regular = []
indel_preferred = []
indel_regular = []

directories = ["*Current Directory*", "*Open File Picker*"]
csv_files = []
for x in os.listdir('.'):
    extension = os.path.splitext(x)[1]
    if extension == ".csv":
        csv_files.append(x)
        if "SNP" in x.upper():
            snp_preferred.append(x)
        else:
            snp_regular.append(x)
        if "INDEL" in x.upper():
            indel_preferred.append(x)
        else:
            indel_regular.append(x)
    elif extension == ".txt":
        txt_files.append(x)
    elif extension == ".fasta":
        fasta_files.append(x)
    elif os.path.isdir(x):
        directories.append(x)
snp_sources = ["None", "*Open File Picker*"] + snp_preferred + snp_regular
indel_sources = ["None", "*Open File Picker*"] + indel_preferred + indel_regular
ref_files = ["*Open File Picker*"] + fasta_files + txt_files

b_color = "white"
a_color = "#e7ebf6"
    
    
# Everything from this point down is used for the GUI
class GUI:
    def __init__(self, master):
        # Organize the window into a bunch of frames
        self.master = master

        self.search_frame = Tkinter.Frame(master, background=b_color)
        self.search_frame.pack(anchor='n')

        # Frame containing parameters on the left
        self.m = Tkinter.Frame(master, padx=20, background=b_color)
        self.m.pack(side="left", anchor="n")

        # Separator line between right frame and left frame
        self.sep = ttk.Separator(master, orient="vertical")
        self.sep.pack(side="left", expand=True, fill="y", pady=15)

        # Frame containing checkboxes on the right
        self.c = Tkinter.Frame(master, padx=20, background=b_color)
        self.c.pack(pady=10)
        self.snp_open = {}
        self.indel_open = {}
        self.indel_cols = {}
        self.snp_cols = {}
        self.top_frame = Tkinter.Frame(self.m, pady=10, background=b_color)
        self.top_frame.pack(anchor="w")
        self.mid_frame = Tkinter.Frame(self.m, background=b_color)
        self.mid_frame.pack()
        self.bottom_frame = Tkinter.Frame(self.m, background=b_color)
        self.bottom_frame.pack()
        self.checks_s = []
        self.buttons_s = False
        self.checks_i = []
        self.buttons_i = False
        self.p_active = False
        self.p_type_active = False
        self.master = master
        self.left1_frame = Tkinter.Frame(self.c, background=b_color)
        self.left1_frame.grid(row=1, column=0, sticky="n", padx=5)
        self.right1_frame = Tkinter.Frame(self.c, background=b_color)
        self.right1_frame.grid(row=1, column=1, sticky="n", padx=5)
        self.button1_frame = Tkinter.Frame(self.c, background=b_color)
        self.button1_frame.grid(row=0, columnspan=2)
        self.left2_frame = Tkinter.Frame(self.c, background=b_color)
        self.left2_frame.grid(row=1, column=2, sticky="n", padx=5)
        self.right2_frame = Tkinter.Frame(self.c, background=b_color)
        self.right2_frame.grid(row=1, column=3, sticky="n", padx=5)
        self.button2_frame = Tkinter.Frame(self.c, background=b_color)
        self.button2_frame.grid(row=0, column=2, columnspan=2)
        self.rest = Tkinter.Frame(self.c, background=b_color)
        self.rest.grid(row=1, column=4)

        # Make the initial labels
        self.s_label = Tkinter.Label(self.button1_frame, bg=b_color,
                                     text="No SNP Data Selected")
        self.s_label.grid(row=0, columnspan=2, padx=(0,10))

        self.indel_label = Tkinter.Label(self.button2_frame, bg=b_color,
                                         text="No Indel Data Selected")
        self.indel_label.grid(row=0, columnspan=2)

        # Select SNP File
        self.s_file = Tkinter.StringVar(master)
        self.disp_s = Tkinter.StringVar(master)
        self.snp_menu = ttk.OptionMenu(self.top_frame, self.disp_s,
                                       "Select SNP Data", *snp_sources,
                                       command=self.reset_s_list, style="My.TMenubutton")
        self.snp_menu.grid(row=2, column=1, sticky="W")
        self.snp_label = Tkinter.Label(self.top_frame, text="SNPs:   ", bg=b_color)
        self.snp_label.grid(row=2, column=0, sticky="W")

        # Select Indel File
        self.i_file = Tkinter.StringVar(master)
        self.disp_i = Tkinter.StringVar(master)
        self.indel_menu = ttk.OptionMenu(self.top_frame, self.disp_i,
                                         "Select Indel Data", *indel_sources,
                                         command=self.reset_i_list, style="My.TMenubutton")
        self.indel_menu.grid(row=3, column=1, sticky="W", pady=2)
        self.i_label = Tkinter.Label(self.top_frame, text="Indels:   ", bg=b_color)
        self.i_label.grid(row=3, column=0, sticky="W")

        
        # Confidence Level
        self.conf = Tkinter.IntVar(master)
        self.conf_button = ttk.Checkbutton(self.top_frame, text="High Confidence Only",
                                           variable=self.conf, style="My.TCheckbutton")
        self.conf.set(0)
        self.conf_button.grid(row=4, columnspan=2, sticky="w", pady=(0,10))

        # Select Reference Gene
        self.disp_r = Tkinter.StringVar(master)
        self.r_file = Tkinter.StringVar(master)
        self.ref_menu = ttk.OptionMenu(self.top_frame, self.disp_r,
                                       "Select Reference Gene", style="My.TMenubutton",
                                       command=self.parse, *ref_files)
        self.ref_menu.grid(row=1, column=1, sticky="W", pady=2)
        self.r_label = Tkinter.Label(self.top_frame, text="Reference Gene:    ",
                                     bg=b_color)
        self.r_label.grid(row=1, column=0, sticky="W")

        # Enter Filename to Save 
        self.name = Tkinter.StringVar(master)
        self.name.set("*Current Directory*")
        self.entry = ttk.OptionMenu(self.top_frame, self.name,
                                    "*Current Directory*", *directories, style="My.TMenubutton",
                                    command=self.get_filename)
        self.entry.grid(row=0, column=1, sticky="W", pady=(0,10), padx=(0,0))
        self.e1_label = Tkinter.Label(self.top_frame, text="Output Directory:    ",
                                      bg=b_color)
        self.e1_label.grid(row=0, column=0, sticky="W", pady=(0,10))

        # Protein, Gene, or RNA
        self.disp_a = Tkinter.StringVar(master)
        self.a_file = Tkinter.StringVar(master)
        self.threshold = Tkinter.IntVar(master)
        self.threshold.set(12)
        self.p_or_g = Tkinter.IntVar(master)
        self.p_or_g.set(0)
        self.radio_frame = Tkinter.Frame(self.bottom_frame, background=b_color)
        self.radio_frame.grid(row=1, columnspan=2)
        self.gene_button = ttk.Radiobutton(self.radio_frame, text="Gene",
                                           variable=self.p_or_g, value=0,
                                           style="My.TRadiobutton",
                                           command=self.p_del)
        self.gene_button.grid(row=1, column=0)
        self.protein_button = ttk.Radiobutton(self.radio_frame, text="Protein",
                                              variable=self.p_or_g, value=2,
                                              style="My.TRadiobutton",
                                              command=self.p_info)
        self.protein_button.grid(row=1, column=2)
        self.protein_button = ttk.Radiobutton(self.radio_frame, text="mRNA",
                                              variable=self.p_or_g, value=1,
                                              style="My.TRadiobutton",
                                              command=self.p_info)
        self.protein_button.grid(row=1, column=1, padx=10)

        # Peptide-Type
        self.p_type = Tkinter.IntVar(master)
        self.p_type.set(0)

        # Enter Gene Location
        self.locus = Tkinter.IntVar(master)
        self.locus.set("")
        self.entry2 = ttk.Entry(self.mid_frame, textvariable=self.locus,
                                width=12)
        self.entry2.grid(row=4, column=1, sticky="W")
        self.e2_label = Tkinter.Label(self.mid_frame, text="Gene Location:",
                                      bg=b_color)
        self.e2_label.grid(row=4, column=0, sticky="W", padx=(0,25))

        # Query Genes
        self.db_query = Tkinter.StringVar(master)
        self.testdb = ttk.Entry(self.search_frame, textvariable=self.db_query)
        self.testdb.bind("<Return>", lambda event: self.db_fill(self.db_query.get()))
        self.testdb.grid(row=0, column = 1, sticky="W")
        self.db_title = Tkinter.Label(self.search_frame, text="Search Gene Symbol:",
                                      bg=b_color)
        self.db_title.grid(row=0, column=0, sticky="E")

        # Warning Label
        self.label = Tkinter.Label(self.mid_frame, text="",
                                   fg="red3", justify="left", bg=b_color)
        self.label.grid(row=6, columnspan=2)

        # Button to Start Imputation     
        self.button = ttk.Button(self.m, text="IMPUTE", command=self.start)
        self.button.pack(pady=(10,0))

        # Hyperlinks
        self.links = Tkinter.Frame(self.m, background=b_color)
        self.links.pack(pady=(40, 10))

        self.sanger_link = Tkinter.Label(self.links, text="Sanger SNP/Indel Database",
                                        fg="blue", cursor="hand2", bg=b_color)
        self.sanger_link.pack(pady=5)
        self.sanger_link.bind("<Button-1>", call_sanger)

        self.ref_link = Tkinter.Label(self.links, text="NCBI Reference Gene Database",
                                        fg="blue", cursor="hand2", bg=b_color)
        self.ref_link.pack(pady=5)
        self.ref_link.bind("<Button-1>", call_reference)
        
        self.EMBOSS_link = Tkinter.Label(self.links, text="EMBOSS Stretcher Alignment",
                                        fg="blue", cursor="hand2", bg=b_color)
        self.EMBOSS_link.pack(pady=5)
        self.EMBOSS_link.bind("<Button-1>", call_EMBOSS)
        
    def p_info(self):
        if not self.p_active:
            self.p_frame = Tkinter.Frame(self.bottom_frame, background=b_color)
            self.p_frame.grid(row=2, columnspan=2)
            self.a_menu = ttk.OptionMenu(self.p_frame, self.disp_a,
                                         "Select Alignment File", *ref_files,
                                         style="My.TMenubutton", command=self.a_load)
            self.a_menu.grid(row=2, column=1, sticky="E", pady=(10,5))
            self.a_label = Tkinter.Label(self.p_frame, text="Alignment:   ", bg=b_color)
            self.a_label.grid(row=2, column=0, sticky="W")
            self.t_entry = ttk.Entry(self.p_frame, textvariable=self.threshold,
                                     width=5)
            self.t_entry.grid(row=3, column=1, pady=5)
            self.t_label = Tkinter.Label(self.p_frame, text="Exon Seed Length:",
                                         bg=b_color)
            self.t_label.grid(row=3, column=0, sticky="W", padx = (0,10))
            
        if self.p_or_g.get() == 2 and not self.p_type_active:
            self.utr = ttk.Radiobutton(self.p_frame,
                                       text="Minimize 5' UTR",
                                       variable=self.p_type, value=0,
                                       style="My.TRadiobutton")
            self.orf = ttk.Radiobutton(self.p_frame,
                                       text="Maximize ORF",
                                       variable=self.p_type, value=1,
                                       style="My.TRadiobutton")
            
            self.utr.grid(row=4, column=0, sticky='W', padx=(0,10))
            self.orf.grid(row=4, column=1, sticky='E', padx=(0,10))
            
            self.p_type_active = True

        elif self.p_or_g.get() == 1 and self.p_type_active:
            self.utr.destroy()
            self.orf.destroy()
            self.p_type_active = False
            
        self.p_active = True
            
    def p_del(self):
        if self.p_active:
            self.p_frame.destroy()
        self.p_active = False
        self.p_type_active = False

        
    # Function for the "All" button for the SNPs
    def sel_all_s(self):
        for g in self.snp_open:
            self.snp_open[g].set(1)
            
    # Function for the "None" button for the SNPs
    def sel_none_s(self):
        for g in self.snp_open:
            self.snp_open[g].set(0)

    # Function for the "All" button for the Indels
    def sel_all_i(self):
        for g in self.indel_open:
            self.indel_open[g].set(1)

    # Function for the "None" button for the Indels
    def sel_none_i(self):
        for g in self.indel_open:
            self.indel_open[g].set(0)
    
    # Resets the list each time a file is chosen
    def reset_s_list(self, selection):
        # For some reason the button is sending a second argument to the function
        # Problem is fixed by allowing the function to take extra *args
        
        for widget in self.checks_s:
            widget.destroy()
        if selection == "None":
            # Clears everything if "None" is chosen as the file
            try:
                self.all_s.destroy()
                self.none_s.destroy()
                self.buttons_s = False
                self.snp_open.clear()
                self.snp_cols.clear()
                self.s_label.configure(text="No SNP Data Selected")
            except AttributeError:
                pass
            self.s_file.set("None")
        elif selection.split()[0] == "Download":
            self.disp_s.set("Downloading...")
            self.snp_menu.update()
            download = sangerDownload(self.db[1], self.db[3], "SNP")
            self.reset_s_list(download)
            return None
        elif selection == "*Open File Picker*":
            snp_file = (tkFileDialog.askopenfilename(title="Select SNP File",
                                                     filetypes=(("CSV files", "*.csv"),
                                                     ("all files", "*.*"))))
            if snp_file == '':
                self.reset_s_list('None')
                return None
            self.s_file.set(snp_file)
            self.disp_s.set(os.path.split(snp_file)[1])
        else:
            self.s_file.set(selection)
            self.disp_s.set(os.path.split(selection)[1])
            self.snp_menu.update()
        
        # Place the "All" and "None" buttons if a file is chosen
        if not self.buttons_s and self.s_file.get() not in ["None", "Select SNP Data"]:
            self.all_s = ttk.Button(self.button1_frame, text="All",
                                    command=self.sel_all_s, style="My.TButton")
            self.all_s.grid(row=1, column=0)
            self.none_s = ttk.Button(self.button1_frame, text="None",
                                     command=self.sel_none_s, style="My.TButton")
            self.none_s.grid(row=1, column=1, padx=20)
            self.buttons_s = True
            self.s_label.configure(text="SNP Strains")
        
        if self.s_file.get() not in ["None", "Select SNP Data"]:
            self.snp_open.clear()
            self.snp_cols.clear()
            # Collect the strain names from header
            with open(self.s_file.get(), "rt") as snp_source:
                data = csv.reader(snp_source)
                for new_datarow in data:
                    headers = new_datarow
                    break
                for i, s in enumerate(headers):
                    if i > 4 and i % 2 == 1 and s:
                        s = s.replace("_", "/")
                        self.snp_open[s] = 1
                        self.snp_cols[s] = i
                        
            # Store the checkbuttons themselves in a list to track their properties
            self.checks_s = []
            position_counter = 0
            # Cut the list in half for even length of column 1 and column 2 checkboxes
            length = int(len(self.snp_open) / 2) - 1
            for strain in sorted(self.snp_open):
                self.snp_open[strain] = Tkinter.IntVar()
                self.snp_open[strain].set(1)
                if position_counter > length:
                    t_frame = self.right1_frame
                else:
                    t_frame = self.left1_frame

                self.l = ttk.Checkbutton(t_frame, text=strain, style="My.TCheckbutton",
                                         variable=self.snp_open[strain], state="normal")
                self.l.pack(anchor="w")
                self.checks_s.append(self.l)
                position_counter += 1

    def get_filename(self, selection):
        if selection == "*Current Directory*":
            val = '.'
        elif selection == "*Open File Picker*":
            val = tkFileDialog.askdirectory()
        else:
            val = selection
        global output_directory
        output_directory = val
        if val != '.':
            self.name.set(os.path.split(val)[-1])
                
    # Same as reset_s_list but for indels
    def reset_i_list(self, selection):
        for widget in self.checks_i:
            widget.destroy()
        if selection == "None":
            try:
                self.all_i.destroy()
                self.none_i.destroy()
                self.buttons_i = False
                self.indel_open.clear()
                self.indel_cols.clear()
                self.indel_label.configure(text="No Indel Data Selected")
            except AttributeError:
                pass
            self.i_file.set("None")
        elif selection.split()[0] == "Download":
            self.disp_i.set("Downloading...")
            self.indel_menu.update()
            download = sangerDownload(self.db[1], self.db[3], "Indel")
            self.reset_i_list(download)
            return None
        elif selection == "*Open File Picker*":
            indel_file = (tkFileDialog.askopenfilename(title="Select Indel File",
                                                       filetypes=(("CSV files", "*.csv"),
                                                       ("all files", "*.*"))))
            if indel_file == '':
                self.reset_i_list('None')
                return None
            self.i_file.set(indel_file)
            self.disp_i.set(os.path.split(indel_file)[1])
        else:
            self.i_file.set(selection)
            self.disp_i.set(os.path.split(selection)[1])
            self.indel_menu.update()
        if not self.buttons_i and self.i_file.get() not in ["None", "Select Indel Data"]:
            self.all_i = ttk.Button(self.button2_frame, text="All",
                                    command=self.sel_all_i, style="My.TButton")
            self.all_i.grid(row=1, column=0)
            self.none_i = ttk.Button(self.button2_frame, text="None",
                                     command=self.sel_none_i, style="My.TButton")
            self.none_i.grid(row=1, column=1, padx=20)
            self.buttons_i = True
            self.indel_label.configure(text="Indel Strains")
        if self.i_file.get() not in ["None", "Select SNP Data"]:
            self.indel_open.clear()
            self.indel_cols.clear()
            with open(self.i_file.get(), "rt") as i_source:
                data = csv.reader(i_source)
                for new_datarow in data:
                    headers = new_datarow
                    break
                for i, s in enumerate(headers):
                    if i > 4 and i % 2 == 1 and s:
                        s = s.replace("_", "/")
                        self.indel_open[s] = 1
                        self.indel_cols[s] = i
            for w in self.checks_i:
                w.pack_forget()
            self.checks_i = []
            position_counter = 0
            length = int(len(self.indel_open) / 2) - 1
            for strain in sorted(self.indel_open):
                self.indel_open[strain] = Tkinter.IntVar()
                self.indel_open[strain].set(1)
                if position_counter > length:
                    t_frame = self.right2_frame

                else:
                    t_frame = self.left2_frame

                self.l = ttk.Checkbutton(t_frame, text=strain, variable=self.indel_open[strain],
                                         state="normal", style="My.TCheckbutton")
                self.l.pack(anchor="w")
                self.checks_i.append(self.l)
                position_counter += 1

    def a_load(self, selection):
        if selection == "*Open File Picker*":
            a_file = (tkFileDialog.askopenfilename(title="Select mRNA File",
                                                   filetypes=(("Text files", "*.txt"),
                                                   ("FASTA files", "*.fasta"),
                                                   ("all files", "*.*"))))
            if a_file == '':
                return None
            self.a_file.set(a_file)
            self.disp_a.set(os.path.split(a_file)[1])
        elif selection.split()[0] == "Download:":
            acc = self.transcripts[selection[10:]]
            self.disp_a.set("Downloading...")
            self.a_menu.update()
            download = downloadRNA(acc, self.db[3])
            self.a_file.set(download)
            self.disp_a.set(os.path.split(download)[1])
            self.a_menu.update()
        else:
            self.a_file.set(selection)
            self.disp_a.set(selection)

    def db_fill(self, gene):

        try:
            search = searchDNA(gene)
            index = self.selectDNA(search)
            if type(index) is not int:
                return None
            self.db = search[index]
            self.transcripts = searchRNA(self.db[3])
        except RuntimeError:
            tkMessageBox.showerror("Database Error",
                                   "Could not find gene symbol in database")
            return None
        except:
            tkMessageBox.showerror("Database Error",
                                   "Could not connect to NCBI Gene")
            return None
        
        official_gene = self.db[3]
        new_snp = "Download {} SNPS".format(official_gene)
        new_indel = "Download {} Indels".format(official_gene)
        snp_options = ["None", "*Open File Picker*", new_snp]
        self.update_menu(self.snp_menu, self.reset_s_list, snp_options)
        indel_options = ["None", "*Open File Picker*", new_indel]
        self.update_menu(self.indel_menu, self.reset_i_list, indel_options)
        r_gene = "Download {} Sequence".format(official_gene)
        r_options = ["None", "*Open File Picker*", r_gene]
        self.update_menu(self.ref_menu, self.parse, r_options)
        transcript_text = ["Download: " + x for x in self.transcripts]
        a_options = ["*Open File Picker*"]+transcript_text
        try:
            self.update_menu(self.a_menu, self.a_load, a_options)
        except AttributeError:
            global ref_files
            ref_files = a_options
    def update_menu(self, menu, c_func, options):
        m = menu["menu"]
        m.delete(0, "end")
        for string in options:
            m.add_command(label=string,
                          command=lambda value=string: c_func(value))

    def selectDNA(self, summary_info):
        if len(summary_info) == 1:
            tkMessageBox.showinfo("Gene Information",
                                  "We found the following gene:\n" +\
                                  summary_info[0][2])
            return 0
        if len(summary_info) == 0:
            tkMessageBox.showerror("Database Error",
                                  "No sequences available for queried gene")
            return 0
        inputDialog = DNA_SelectionWindow(root, summary_info)
        root.wait_window(inputDialog.top)
        return inputDialog.index
        
        
    # Parse FASTA Header for gene location check if complementary
    def parse(self, selection):
        if selection.split()[0] == "Download":
            self.disp_r.set("Downloading...")
            self.ref_menu.update()
            download = downloadDNA(self.db[0], self.db[3], self.db[4], self.db[5])
            self.r_file.set(download)
            self.disp_r.set(os.path.split(download)[1])
            self.ref_menu.update()
            self.parse(download)
            return None
        elif selection == "*Open File Picker*":
            r_file = (tkFileDialog.askopenfilename(title="Select mRNA File",
                                                   filetypes=(("Text files", "*.txt"),
                                                   ("FASTA files", "*.fasta"),
                                                   ("all files", "*.*"))))
            if r_file == '':
                return None
            self.r_file.set(r_file)
            self.disp_r.set(os.path.split(r_file)[1])
            self.parse(r_file)
            return None
        else:
            self.r_file.set(selection)
        try:
            self.label.configure(text="")
            with open(self.r_file.get(), "r") as f:
                string = f.read()
                # First, the FASTA header is separated and parsed for starting location
                sep = string.index("\n")
                header = string[:sep]
                indices = []
                for i, char in enumerate(header):
                    if char in [":", "-", " "]:
                        indices.append(i)
                # If the first number has a 'c' before it, complementary is made True
                try:
                    self.locus.set(int(header[indices[0] + 1:indices[1]]))
                except ValueError:
                    self.locus.set(int(header[indices[1] + 1:indices[2]]))
                    # Starting location of gene identified
                self.separator = sep
        except (ValueError, IndexError) as header_not_found:
            # Happens when there is no header
            self.label.configure(text=warning.no_header)
            self.locus.set("")
            self.separator = -1

    def start(self):
        proceed = True
        ref_file = self.r_file.get()
        snp_file = self.s_file.get()
        indel_file = self.i_file.get()
        try:
            locus = self.locus.get()
        except ValueError:
            locus = ""
        if ref_file == "Select Reference Gene":
            tkMessageBox.showerror(error.no_ref,
                                   error.missing_ref)
        else:
            # File becomes empty string if "None" is selected
            if snp_file == "None" or snp_file == "Select SNP Data":
                snp_file = ""
            indel_file = self.i_file.get()
            if indel_file == "None" or indel_file == "Select Indel Data":
                indel_file = ""
            if (snp_file or indel_file) and locus == "":
                proceed = False
                tkMessageBox.showerror(error.missing_field,
                                       error.missing_missing_loc)
                    
            if self.p_or_g.get() != 0 and self.a_file.get() == "Select Alignment File":
                if proceed:
                    proceed = False
                    tkMessageBox.showerror(error.missing_field,
                                           error.missing_alignment)
    
            # The impute function below triggers the actual script
            if proceed:
                confidence = bool(self.conf.get())
                t_strain = []
                # Collect combined list of selected strains without duplicates
                for g in set(self.snp_cols) | set(self.indel_cols):
                    s_col, i_col = False, False
                    if g in self.snp_open:
                        if self.snp_open[g].get():
                            s_col = self.snp_cols[g]
                    if g in self.indel_open:
                        if self.indel_open[g].get():
                            i_col = self.indel_cols[g]
                    # Tup stores the strain name and which column to look at in
                    # the SNP or Indel csv file
                    tup = g, s_col, i_col
                    if tup[1] or tup[2]:
                        t_strain.append(tup)
                # Send the user input to the imputer function                    
                impute(snp_file, indel_file, ref_file, locus, self.separator,
                       t_strain, confidence, self.p_or_g.get(),
                       self.threshold.get(), self.a_file.get(), self.db_query.get(),
                       self.p_type.get())

class DNA_SelectionWindow:
    def __init__(self, parent, summary_info):
        top = self.top = Tkinter.Toplevel(parent, bg=b_color)
        self.myLabel = Tkinter.Label(top, text='Select the Desired Gene',
                                     bg=b_color)
        self.myLabel.grid(row=0, columnspan=2, pady=(5, 10))
        self.chosenIndex = Tkinter.IntVar(parent)
        for i, s in enumerate(summary_info):
            self.gene_button = ttk.Radiobutton(top, text=s[2],
                                               variable=self.chosenIndex, value=i,
                                               style="My.TRadiobutton")
            self.gene_button.grid(row=i+1, columnspan=2, pady=5, padx=5,
                                  sticky="W")
        self.mySubmitButton = ttk.Button(top, text='Okay', command=self.send)
        self.mySubmitButton.grid(row=i+2, sticky='E', column=0, pady=10)
        self.myCancelButton = ttk.Button(top, text='Cancel', command=self.cancel)
        self.myCancelButton.grid(row=i+2, sticky='W', column=1)

    def send(self):
        self.index = self.chosenIndex.get()
        self.top.destroy()
    def cancel(self):
        self.username = None
        self.top.destroy()

# Activate GUI
root = Tkinter.Tk()
root.configure(bg=b_color)
root.style = ttk.Style()
root.style.configure("My.TCheckbutton", background=b_color)
root.style.configure("My.TMenubutton", background=a_color)
root.style.configure("My.TRadiobutton", background=b_color)
root.title("Gene Coordinates Imputation")
main_ui = GUI(root)
root.mainloop()

#########################################################################################
# .				   -Peter Dornbos-
#               /hs.            .sh/               
#             /mN-                -Nm/             
#           `yMM-                  -MMy`           
#           dMMy                    yMMd           
#          oMMMs       .-::-.       sMMMo          
#          NMMMm  .+ymMMMMMMMMmy+.  mMMMN          
#          MMMMMs-MMMMmhsssshmMMMM-sMMMMM          
#       `+hMMMMMMh/y/`        `/y/hMMMMMMh+`       
#     -hMMMMMMMMMMMh/`        `/hMMMMMMMMMMMh-     
#   `yMMMMMMNNNMMMMMMMmhs:/shmMMMMMMMNNNMMMMMMy`   
#  `mMMmo:. :yyo-/smMMMd:``:hMMMms/-syy: .:omMMm`  
#  dMh-     oMMM`   :dM.    .Mm:   `MMMo     -hMd  
# /M+       /MMM:     :o:``:o:     :MMM/       +M/ 
# ys        `mMMm`     yMMMMh     `mMMm         sy 
#  `         -NMMm:    -MMMM-    :mMMN-         `y 
#            .dMMMd/` -MMMM- `/dMMMd.          `/ 
#               :dMMMM+oMMMMo+MMMMd:            
#                 .+yy/MMMMMM/yy+.                 
#                   `sMMMMMMMMs`               
#                `:yNMMMMMMMMMMNy:`            
#        `+syyyhdmMMMMMMmy::smMMMMMMmdhyyys+`       
#           `:/+ooo+/-`      `-/+ooo+/:.                                                                 
#########################################################################################
