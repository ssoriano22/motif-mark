#!/usr/bin/env python

import cairo
import math
import re
import bioinfo
import argparse

class Exon:

    def __init__(self,exon_seq,x1,y1,x2,y2):
        self.eseq = exon_seq
        self.start_x = x1
        self.start_y = y1
        self.end_x = x2
        self.end_y = y2
    
    def draw_exon(self,ctx):
        '''Draw a exon (thicker line) with pycairo. Enter context object to draw on.'''
        #ctx.rectangle(self.start_x,self.start_y,self.end_x,self.end_y) #(x1,y1,x2,y2)
        ctx.set_line_width(100)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to(self.start_x,self.start_y)  #(x1,y1)
        ctx.line_to(self.end_x,self.end_y) #(x2,y2)
        ctx.stroke()
        print("Exon successfully drawn!")

class Intron:
    
    def __init__(self,intron_seq,line_x1,line_y1,line_x2,line_y2):
        self.iseq = intron_seq
        self.Lx1 = line_x1
        self.Ly1 = line_y1
        self.Lx2 = line_x2
        self.Ly2 = line_y2

    def draw_intron(self,ctx2):
        '''Draw a intron (line) with pycairo. Enter context object to draw on.'''
        ctx2.set_line_width(10)
        ctx2.set_source_rgb(0, 0, 0)
        ctx2.move_to(self.Lx1,self.Ly1)  #(x1,y1)
        ctx2.line_to(self.Lx2,self.Ly2) #(x2,y2)
        ctx2.stroke()
        print("Intron successfully drawn!")

class Gene:

    def __init__(self,gene_Name,gene_Chr,gene_Start,gene_End,gene_Seq):
        self.name = gene_Name
        self.chr = gene_Chr
        self.startLoc = gene_Start
        self.endLoc = gene_End
        self.seq = gene_Seq
        self.motifmap = {}

    def get_Xons(self):
        '''Parse input fasta sequence into exon (uppercase seq), and upstream/downstream introns (lowercase seq). Returns a tuple(intron_upstream,exon,intron_downstream).'''
        intron_up = re.split("[A-Z]+",self.seq)[0]
        intron_down = re.split("[A-Z]+",self.seq)[1]
        exon = re.split("[a-z]+",self.seq)[1]
        Xons = (intron_up, exon, intron_down)
        return Xons
    
    def draw_motifmap(self,ctx3,c_height,margin,the_UpIntron,the_Exon,the_DownIntron):
        '''Draw motifmap for gene with pycairo. Enter context object to draw on.'''
        print("Full Gene Motif Map")
        print(self.motifmap)
        print("Intron/Exon/Intron Motif Map")
        mcolor_r = (0.1,0.2,0.4,0.5,0.6,0.8)
        mcolor_g = (0.5,0.2,0.8,0.1,0.4,0.6)
        mcolor_b = (0.4,0.5,0.6,0.8,0.1,0.2)
        for gkey,gseq in self.motifmap.items():
            #0 = up_intron, 1 = exon, 2 = down_intron
            print("Gene Level")
            print(gkey)
            print(gseq)
            i = 0 # Motif Color index - use diff color for each motif (max 5)
            for mkey,motseq in gseq.items():
                #Motif (key) and dict (value) of each coords motif was found
                print("Motif Level")
                print(mkey) #Motif
                print(motseq) #Dict: {(coords):motif}
                for ckey,mstr in motseq.items():
                    #Coords:Motif Dict
                    print("Coords Level")
                    print(ckey)
                    print(mstr)
                    print(i)
                    ctx3.set_line_width(100)
                    ctx3.set_source_rgb(mcolor_r[i],mcolor_g[i],mcolor_b[i])
                    if gkey == 0: #Up-Intron
                        ctx3.move_to(margin+ckey[0],c_height)  #(x1,y1) #margin = n (100 in this case), c_height altered per gene in main code
                        ctx3.line_to(margin+ckey[1],c_height) #(x2,y2)
                        ctx3.stroke()
                    elif gkey == 1: #Exon
                        ctx3.move_to(margin+len(the_UpIntron.iseq)+ckey[0],c_height)  #(x1,y1) #100 = n (margin), c_height altered per gene in main code
                        ctx3.line_to(margin+len(the_UpIntron.iseq)+ckey[1],c_height) #(x2,y2)
                        ctx3.stroke()
                    else: #Down-Intron
                        ctx3.move_to(margin+len(the_UpIntron.iseq)+len(the_Exon.eseq)+ckey[0],c_height)  #(x1,y1) #100 = n (margin), c_height altered per gene in main code
                        ctx3.line_to(margin+len(the_UpIntron.iseq)+len(the_Exon.eseq)+ckey[1],c_height) #(x2,y2)
                        ctx3.stroke()
                i += 1

        print("Motif Map successfully drawn!")

class Motif:

    def __init__(self, motif_seq):
        self.mseq = motif_seq
    
    def searchMotif(self,cseq):
        '''Check input sequence for specified motifs.'''
        found_motif_loc_dict = {}
        rx_mseq: str = self.mseq
        if "Y" in self.mseq:
            rx_mseq = re.sub("Y","[CTU]",rx_mseq)
        if "R" in self.mseq:
            rx_mseq = re.sub("R","[AG]",rx_mseq)
        if "U" in self.mseq:
            rx_mseq = re.sub("U","[TU]",rx_mseq)
        if "T" in self.mseq:
            rx_mseq = re.sub("T","[TU]",rx_mseq)
        if "W" in self.mseq:
            rx_mseq = re.sub("W","[ATU]",rx_mseq)
        if "S" in self.mseq:
            rx_mseq = re.sub("S","[CG]",rx_mseq)
        if "M" in self.mseq:
            rx_mseq = re.sub("M","[AC]",rx_mseq)
        if "K" in self.mseq:
            rx_mseq = re.sub("K","[GTU]",rx_mseq)
        if "B" in self.mseq:
            rx_mseq = re.sub("B","[CGTU]",rx_mseq)
        if "D" in self.mseq:
            rx_mseq = re.sub("D","[AGTU]",rx_mseq)
        if "H" in self.mseq:
            rx_mseq = re.sub("H","[ACTU]",rx_mseq)
        if "V" in self.mseq:
            rx_mseq = re.sub("V","[ACG]",rx_mseq)
        if "N" in self.mseq:
            rx_mseq = re.sub("N","[ACGTU]",rx_mseq)
        # print(cseq)
        # print(rx_mseq)

        for match in re.finditer(rx_mseq,cseq.upper()):
            # print(match.group())
            found_motif_loc_dict[match.span()] = self.mseq

        return found_motif_loc_dict

#Argparse
def get_args():
    '''Argparse: Retrieves user-input arguments from command line.'''
    parser = argparse.ArgumentParser(description="A program to create a genetic diagram (PNG) of motif locations in a given FASTA file. Inputs: -f INPUT_FILENAME.fasta -m MOTIFS_FILENAME.txt")

    parser.add_argument("-f","--file",help="Input FASTA file",type=str)
    parser.add_argument("-m","--motifs",help="Input file with motifs",type=str)
    return parser.parse_args()

#Get argparse variables
args = get_args()
input_FASTA = args.file
input_motifs = args.motifs

#Read in motifs.txt - store motif in known_motifs_set
known_motif_set = set()
with open(input_motifs,"r") as fh_motif:
    for motif in fh_motif:
        current_motif = Motif(motif.strip().upper())
        known_motif_set.add(current_motif)
#print(known_motif_set)

#Read in fasta and transform to onelinefasta with same name prefix
input_filename = input_FASTA.split(".fasta")[0]
output_filename = input_filename + "_temp.fasta"
bioinfo.oneline_fasta(input_FASTA,output_filename)

#Initialize pycairo canvas coordinates for display - png format
canvas_width, canvas_height = 1200,2000
surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, canvas_width, canvas_height)
#Initialize coordinates for drawing + add background
context = cairo.Context(surface)
context.set_source_rgb(1, 1, 1)
context.paint()

#Read in oneline fasta
i = 0 #Initialize counter - for testing
n = 100 #X-axis counter - to adjust x coords for drawing per gene
g = canvas_height/10 #Gene counter - adjusts y coords for drawing each gene
with open(output_filename,"r") as fh_fasta:
    while True:
        i += 1 #Increment record counter (starts w/ 1)
        # if i>1: #FOR TESTING - one gene
        #     break #FOR TESTING 
        #Get current gene header and split into name, chr, start/end loci
        current_header = fh_fasta.readline().strip()
        #EOF case:
        if current_header == "":
            break
        current_genename = re.split(">",current_header)[1].split(" ")[0]
        current_geneend = re.split(">",current_header)[1].split(" ")[1].split("-")[1]
        current_genechr = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[0]
        current_genestart = re.split(">",current_header)[1].split(" ")[1].split("-")[0].split(":")[1]
        #Get current gene sequence + create Gene object
        current_geneseq = fh_fasta.readline().strip()
        currentGene = Gene(current_genename,current_genechr,current_genestart,current_geneend,current_geneseq)
        #Break current sequence into introns and exon. Format: tuple(up_intron,exon,down_intron)
        current_Xons = currentGene.get_Xons()
        #Create Intron/Exon objects for each of current gene (2 Intron, 1 Exon)
        currentUpIntron = Intron(current_Xons[0],n,g,(n+len(current_Xons[0])),g)
        currentExon = Exon(current_Xons[1],(n+len(current_Xons[0])),g,(n+len(current_Xons[0]+current_Xons[1])),g)
        currentDownIntron = Intron(current_Xons[2],(n+len(current_Xons[0]+current_Xons[1])),g,(n+len(current_Xons[0]+current_Xons[1]+current_Xons[2])),g)
        #Draw introns and exon objects to context/pycairo surface
        currentUpIntron.draw_intron(context)
        currentExon.draw_exon(context)
        currentDownIntron.draw_intron(context)

        #Check each sequence for each motif object, compiling all motifs per seq, and then all seq per gene
        mcoords_gene_dict = {}
        s = 0 #Sequence counter: 0=up_intron, 1=exon, 2=down_intron
        for sequence in current_Xons:
            mcoords_seq_dict = {}
            for mot in known_motif_set:
                mcoords_seq_dict[mot.mseq] = mot.searchMotif(sequence)
            mcoords_gene_dict[s] = mcoords_seq_dict
            s += 1
        #print(mcoords_gene_dict)
        #Add gene motif map to current Gene object
        currentGene.motifmap = mcoords_gene_dict
        currentGene.draw_motifmap(context,g,n,currentUpIntron,currentExon,currentDownIntron)
        #Increment n according to total x length drawn for this gene
        g += (canvas_height/10)

#Write surface to png - name using same prefix as input FASTA file
surface.write_to_png(input_filename+".png")
print("Motif Mark PNG generated successfully! Image Name:",input_filename+".png")