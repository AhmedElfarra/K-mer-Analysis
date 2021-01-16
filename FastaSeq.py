import os
import csv
import sys
from itertools import *
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import SeqIO
from sklearn import cluster
import numpy as np
import shutil
import scipy.stats

#helper function to return folder path (outside class)
#works for a filepath or folderpath
def getFolder(filepath):
    if sys.platform[:3] == "win":
        filename_list = os.path.splitext(filepath)[0].split('\\')
        folder_path = ''
        if '.' in filepath:
            for x in range(len(filename_list)-1):
                folder_path = folder_path + filename_list[x] + '\\'
            folder_path = folder_path[:-1]
        else:
            for x in range(len(filename_list)):
                folder_path = folder_path + filename_list[x] + '\\'
            folder_path = folder_path[:-1]
        return folder_path
    else:
        filename_list = os.path.splitext(filepath)[0].split('/')
        folder_path = ''
        if '.' in filepath:
            for x in range(len(filename_list)-1):
                folder_path = folder_path + filename_list[x] + '/'
            folder_path = folder_path[:-1]
        else:
            for x in range(len(filename_list)):
                folder_path = folder_path + filename_list[x] + '/'
            folder_path = folder_path[:-1]
        return folder_path
    
class FastaSeq:
    #initializer: takes a fasta file and converts it to a FastaSeq object (string)
    #default format is 'dsDNA': analyzes both the input sequence and reverse complement
    #can specify format to be 'ssDNA' (only looks at input sequence), no other valid formats
    def __init__(self, seqfile, format = 'dsDNA'):
        try:
            s1 = SeqIO.read(seqfile, 'fasta')
            if sys.platform[:3] == "win":
                filename_list = os.path.splitext(seqfile)[0].split('\\')
            else:
                filename_list = os.path.splitext(seqfile)[0].split('/')            
            self.seq_name = filename_list[len(filename_list)-1]
            self.folder_path = ''
            #check if command window has been set with os.chdir or if the whole filepath is specified, create folder_path attribute accordingly
            if 'Users' in seqfile:
                self.folder_path = getFolder(seqfile)
            else:   
                if sys.platform[:3] == "win":
                    self.folder_path = getFolder(os.getcwd() + '\\' + seqfile)
                else:
                    self.folder_path = getFolder(os.getcwd() + '/' + seqfile)
            self.seq_name = filename_list[len(filename_list)-1]
            #origin is required for FASTA file I/O, as we need our header information when writing files
            self.origin = s1
            self.sequence = s1.seq
            #if len(self.sequence) > 500000:
                #self.sequence = self.sequence[:500000]
            if format == 'dsDNA':
                self.format = 'dsDNA'
                self.sequence = self.sequence + 'N' + self.sequence.reverse_complement()
            elif format == 'ssDNA':
                self.format = 'ssDNA'
            else:
                print('Not a valid format - please choose from 1 (dsDNA - default) or 2 (ssDNA)')
        except:
            #try:
                if sys.platform[:3] == "win":
                    filename_list = os.path.splitext(seqfile)[0].split('\\')
                else:
                    filename_list = os.path.splitext(seqfile)[0].split('/')            
                self.seq_name = filename_list[len(filename_list)-1]
                self.folder_path = ''
                #check if command window has been set with os.chdir or if the whole filepath is specified, create folder_path attribute accordingly
                if 'Users' in seqfile:
                    self.folder_path = getFolder(seqfile)
                else:   
                    if sys.platform[:3] == "win":
                        self.folder_path = getFolder(os.getcwd() + '\\' + seqfile)
                    else:
                        self.folder_path = getFolder(os.getcwd() + '/' + seqfile)
                self.seq_name = filename_list[len(filename_list)-1]
                #origin is required for FASTA file I/O, as we need our header information when writing files
                record_iter = 0
                running_seq = ''
                for record in SeqIO.parse(seqfile, 'fasta') :
                    if record_iter == 0:
                        self.origin = record
                        running_seq = record.seq
                    else :
                        running_seq = running_seq + "N" + record.seq
                    record_iter += 1
                self.sequence = running_seq
                if format == 'dsDNA':
                    self.format = 'dsDNA'
                    self.sequence = self.sequence + 'N' + self.sequence.reverse_complement()
                elif format == 'ssDNA':
                    self.format = 'ssDNA'
                else:
                    print('Not a valid format - please choose from 1 (dsDNA - default) or 2 (ssDNA)')   
            #except:
                #print(seqfile, "not a valid file")
    
    #getKmers: for a kmer length of your choosing, return the frequency of each kmer found in the sequence as a dictionary
    #must install BioPython (via pip) to run this
    def getKmers(self, kmer):
        self.kmer_dict = {}
        kmer_tuples = list(product(['A', 'C', 'G', 'T'], repeat = kmer))
        kmer_list = []
        for x in range(len(kmer_tuples)):
            new_kmer = ''
            for y in range(kmer):
                new_kmer += kmer_tuples[x][y]
            kmer_list.append(new_kmer)
        for x in range(len(kmer_list)):
            self.kmer_dict[kmer_list[x]] = 0
        for substr in kmer_list:
            self.kmer_dict[substr] = self.sequence.count_overlap(substr)
        return self.kmer_dict
        
    #exportKmers: uses getKmers to export kmer frequencies to a .csv or .xlsx named after the seq file and kmer size, located in the folder specified
    #default parameter results in export to current directory folder, can specift output folder
    def exportKmers(self, kmer, exportfolder = None):
        if exportfolder == None:
            exportfolder = self.folder_path
        export_dict = self.getKmers(kmer)
        if sys.platform[:3] == "win":
            with open(exportfolder + '\\' + self.seq_name + str(kmer) + 'mers.csv', mode = 'w+') as csvfile:
                writer = csv.writer(csvfile, lineterminator = '\n')
                for key,item in export_dict.items():
                    writer.writerow([key, item])
        else:
            with open(exportfolder + '/' + self.seq_name + str(kmer) + 'mers.csv', mode = 'w+') as csvfile:
                writer = csv.writer(csvfile, lineterminator = '\n')
                for key,item in export_dict.items():
                    writer.writerow([key, item])
    
    #drawCGRPlot: draws a CGR plot, saves as a PNG in a specified folder
    #uncomment plt.show and comment plt.close if you prefer to view the images as they are created (will save upon closing the window)
    #default parameter results in export to current directory folder, can specift output folder
    #default point size is 0.1 (ideal for larger sequeces), if images are too sparse this can be increased to ~0.3
    def drawCGRPlot(self, exportfolder = None, pointsize = 0.1):
        if exportfolder == None:
            exportfolder = self.folder_path
        seq = self.sequence
        x = 0
        y = 0
        pointlist = [(x,y)]
        cnum = 0
        gnum = 0
        anum = 0
        tnum = 0
        totnum = 0
        cper = 0
        gper = 0
        aper = 0
        tper = 0
        for c in range(len(seq)):
            if seq[c] == 'A':
                anum += 1
                x = (x-1)/2
                y = (y-1)/2
            elif seq[c] == 'C':
                cnum += 1
                x = (x-1)/2
                y = (y+1)/2
            elif seq[c] == 'G':
                gnum += 1
                x = (x+ 1)/2
                y = (y+1)/2
            elif seq[c] == 'T': 
                tnum += 1
                x = (x+1)/2
                y = (y-1)/2
            pointlist.append((x,y))
        xs, ys = zip(*pointlist)
        plt.figure(figsize=(7,7))
        plt.plot(xs, ys, marker = 'o', color = 'black', markersize = pointsize, linewidth = 0)
        plt.title(self.seq_name.capitalize())
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.text(-1.1, -1.1, "A", fontsize = 16)
        plt.text(-1.1, 1.05, "C", fontsize = 16)
        plt.text(1.05, 1.05, "G", fontsize = 16)
        plt.text(1.05, -1.1, "T", fontsize = 16)
        
        totnum = anum + cnum + gnum + tnum
        aper = anum/totnum
        cper = cnum/totnum
        gper = gnum/totnum
        tper = tnum/totnum
        
        plt.text(1.02, 0.15, 'A: ' + str(aper*100)[:5] + r'%', fontsize = 8)
        plt.text(1.02, 0.05, 'C: ' + str(cper*100)[:5] + r'%', fontsize = 8)
        plt.text(1.02, -0.05, 'G: ' + str(gper*100)[:5] + r'%', fontsize = 8)
        plt.text(1.02, -0.15, 'T: ' + str(tper*100)[:5] + r'%', fontsize = 8)
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.gca().axes.get_xaxis().set_visible(False)
        if sys.platform[:3] == "win":
            plt.savefig(exportfolder + '\\' + self.seq_name + 'CGRPlot.png')
        else:
            plt.savefig(exportfolder + '/' + self.seq_name + 'CGRPlot.png')
        #uncomment .show() and comment .close() if you do not prefer autosaving without viewing the image
        #comment .savefig() if you do not want to save a file and simply want to view the plot
        #plt.show()
        plt.close()
        
    #drawFCGRPlot: draws an FCGR plot, saves as a .PNG in a specified folder
    #can change to view the image before saving as in CGR plot
    #can take a kmer input of 1-4, larger sizes not supported
    def drawFCGRPlot(self, kmer, exportfolder = None):
        import seaborn as sns
        if exportfolder == None:
            exportfolder = self.folder_path
        kmer_dict = self.getKmers(kmer)
        kmer_list = []
        label_list = []
        font_size = 0
        if kmer == 1:
            font_size = 15 
            label_list = ['C', 'G', 'A', 'T']
        elif kmer == 2:
            font_size = 12
            label_list = ['CC', 'GC', 'CG', 'GG', 'AC', 'TC', 'AG', 'TG', 'CA', 'GA', 'CT', 'GT', 'AA', 'TA', 'AT', 'TT']
        elif kmer == 3:
            font_size = 9
            label_list = ['CCC', 'GCC', 'CGC', 'GGC', 'CCG', 'GCG', 'CGG', 'GGG', \
                          'ACC', 'TCC', 'AGC', 'TGC', 'ACG', 'TCG', 'AGG', 'TGG', \
                          'CAC', 'GAC', 'ATC', 'CTC', 'CAG', 'GAG', 'CTG', 'GTG', \
                          'AAC', 'TAC', 'GTC', 'TTC', 'AAG', 'TAG', 'ATG', 'TTG', \
                          'CCA', 'GCA', 'CGA', 'GGA', 'CCT', 'GCT', 'CGT', 'GGT', \
                          'ACA', 'TCA', 'AGA', 'TGA', 'ACT', 'TCT', 'AGT', 'TGT', \
                          'CAA', 'GAA', 'CTA', 'GTA', 'CAT', 'GAT', 'CTT', 'GTT', \
                          'AAA', 'TAA', 'ATA', 'TTA', 'AAT', 'TAT', 'ATT', 'TTT']
        elif kmer == 4:
            label_list =  ['CCCC', 'GCCC', 'CGCC', 'GGCC', 'CCGC', 'GCGC', 'CGGC', 'GGGC', 'CCCG', 'GCCG', 'CGCG', 'GGCG', 'CCGG', 'GCGG', 'CGGG', 'GGGG']
            label_list += ['ACCC', 'TCCC', 'AGCC', 'TGCC', 'ACGC', 'TCGC', 'AGGC', 'TGGC', 'ACCG', 'TCCG', 'AGCG', 'TGCG', 'ACGG', 'TCGG', 'AGGG', 'TGGG']
            label_list += ['CACC', 'GACC', 'CTCC', 'GTCC', 'CAGC', 'GAGC', 'CTGC', 'GTGC', 'CACG', 'GACG', 'CTCG', 'GTCG', 'CAGG', 'GAGG', 'CTGG', 'GTGG']
            label_list += ['AACC', 'TACC', 'ATCC', 'TTCC', 'AAGC', 'TAGC', 'ATGC', 'TTGC', 'AACG', 'TACG', 'ATCG', 'TTCG', 'AAGG', 'TAGG', 'ATGG', 'TTGG']
            label_list += ['CCAC', 'GCAC', 'CGAC', 'GGAC', 'CATC', 'GATC', 'CCTC', 'GCTC', 'CCAG', 'GCAG', 'CGAG', 'GGAG', 'CCTG', 'GCTG', 'CGTG', 'GGTG']
            label_list += ['ACAC', 'TCAC', 'AGAC', 'TGAC', 'AATC', 'TATC', 'ACTC', 'TCTC', 'ACAG', 'TCAG', 'AGAG', 'TGAG', 'ACTG', 'TCTG', 'AGTG', 'TGTG']                       
            label_list += ['CAAC', 'GAAC', 'CTAC', 'GTAC', 'CGTC', 'GGTC', 'CTTC', 'GTTC', 'CAAG', 'GAAG', 'CTAG', 'GTAG', 'CATG', 'GATG', 'CTTG', 'GTTG']                                                
            label_list += ['AAAC', 'TAAC', 'ATAC', 'TTAC', 'AGTC', 'TGTC', 'ATTC', 'TTTC', 'AAAG', 'TAAG', 'ATAG', 'TTAG', 'AATG', 'TATG', 'ATTG', 'TTTG']                                                                        
            label_list += ['CCCA', 'GCCA', 'CGCA', 'GGCA', 'CCGA', 'GCGA', 'CGGA', 'GGGA', 'CCCT', 'GCCT', 'CGCT', 'GGCT', 'CCGT', 'GCGT', 'CGGT', 'GGGT']                                                                                                
            label_list += ['ACCA', 'TCCA', 'AGCA', 'TGCA', 'ACGA', 'TCGA', 'AGGA', 'TGGA', 'ACCT', 'TCCT', 'AGCT', 'TGCT', 'ACGT', 'TCGT', 'AGGT', 'TGGT']                                                                                                                        
            label_list += ['CACA', 'GACA', 'CTCA', 'GTCA', 'CAGA', 'GAGA', 'CTGA', 'GTGA', 'CACT', 'GACT', 'CTCT', 'GTCT', 'CAGT', 'GAGT', 'CTGT', 'GTGT']
            label_list += ['AACA', 'TACA', 'ATCA', 'TTCA', 'AAGA', 'TAGA', 'ATGA', 'TTGA', 'AACT', 'TACT', 'ATCT', 'TTCT', 'AAGT', 'TAGT', 'ATGT', 'TTGT']                        
            label_list += ['CCAA', 'GCAA', 'CGAA', 'GGAA', 'CCTA', 'GCTA', 'CGTA', 'GGTA', 'CCAT', 'GCAT', 'CGAT', 'GGAT', 'CCTT', 'GCTT', 'CGTT', 'GGTT']
            label_list += ['ACAA', 'TCAA', 'AGAA', 'TGAA', 'ACTA', 'TCTA', 'AGTA', 'TGTA', 'ACAT', 'TCAT', 'AGAT', 'TGAT', 'ACTT', 'TCTT', 'AGTT', 'TGTT']
            label_list += ['CAAA', 'GAAA', 'CTAA', 'GTAA', 'CATA', 'GATA', 'CTTA', 'GTTA', 'CAAT', 'GAAT', 'CTAT', 'GTAT', 'CATT', 'GATT', 'CTTT', 'GTTT']
            label_list += ['AAAA', 'TAAA', 'ATAA', 'TTAA', 'AATA', 'TATA', 'ATTA', 'TTTA', 'AAAT', 'TAAT', 'ATAT', 'TTAT', 'AATT', 'TATT', 'ATTT', 'TTTT']                        
            font_size = 5
        for x in range(len(label_list)):
            for key, item in kmer_dict.items():
                if key == label_list[x]:
                    kmer_list.append(item)
        r_and_c = 2 ** kmer
        kmer_values = np.array(kmer_list).reshape(r_and_c,r_and_c)
        labels = np.array(label_list).reshape(r_and_c,r_and_c)
        if kmer < 5:
            ax = sns.heatmap(kmer_values, square = True, cbar = True, yticklabels = False, xticklabels = False, annot = labels, fmt = '', annot_kws={"size": font_size}, cmap = ['#FFFFFF', '#F5F5F5', '#DCDCDC', '#D3D3D3', '#C0C0C0', '#A9A9A9', '#808080', '#696969', '#000000'], linewidths = 0.2, linecolor = 'black')
            ax.set_title(self.seq_name)
            if sys.platform[:3] == "win":
                plt.savefig(exportfolder + '\\' + self.seq_name + str(kmer) + 'merFCGRPlot.png')
            else:
                plt.savefig(exportfolder + '/' + self.seq_name + str(kmer) + 'merFCGRPlot.png')
            plt.close()
        else:
             print('This program does not support FCGR generation with kmers greater than 4')