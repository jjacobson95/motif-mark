#!/usr/bin/env python
from tracemalloc import start
import cairo
import random

#from itertools import count
import argparse
from ctypes.wintypes import HACCEL
from operator import truediv
from pickle import FALSE
import sys

motifs_file = open("./Fig_1_motifs.txt", "r")
gene_file_1 = open("./Figure_1.fasta", "r")


#add argparse

parser = argparse.ArgumentParser(description = "Argparse")
parser.add_argument('-f', '--files', help="Input fasta files in corresponding order to motif files", type=argparse.FileType('r'), nargs='+', required =True)
parser.add_argument('-m', '--motifs', help="Input motif files in corresponding order to fasta files", type=argparse.FileType('r'), nargs='+', required =True)
#parser.add_argument('-h','--help', help="Display long help message", required = False)

args = parser.parse_args()
fasta_files = args.files
motif_files = args.motifs

if len(motif_files) != len(fasta_files):
    # exits the program
    sys.exit("Number of Motif files must be equal to the number of Fasta files.")    
else:
    print("Number of Motif files is equal to the number of Fasta Files. Proceeding.")

class Gene:
    def __init__(self, the_sequence, gene_name, figure_number):
        """Class Gene: This contains the sequence, gene name, and figure number.
            Methods: seq_count: defermine length of sequence and add as attribute. Remove this. just use len in attrib def.
                     IE_pos: Locate and store the points where introns and exons exist.
        """
        self.sequence = the_sequence
        self.start_as = ""
        self.swap_points = []
        self.seq_length = 0
        self.name = gene_name
        self.figure_ID = figure_number
    def seq_count(self):
        self.seq_length = len(self.sequence)
        #print(self.sequence)
        return self.seq_length
    def IE_pos(self):
        counter = 0
        the_query = list(self.sequence)
        counter = -1
        for item in the_query:
            counter += 1
            if item.islower():
                the_query[counter] = "LOW"
                #print(item)
            if item.isupper():
                the_query[counter] = "UP"
                #print(item)
        #print(the_query)
        if the_query[0] == "LOW":
            self.start_as = "Intron"
        if the_query[0] == "UP":
            self.start_as = "Exon"
        holder_val = the_query[0]
        for index, item in enumerate(the_query):
            if (item != holder_val):
                self.swap_points.append(index)
                holder_val = (the_query[index])
        self.swap_points.append(self.seq_length)


class Motif:
    def __init__(self, m_sequence, gene_seq, gene_name, figure_number):
        """
        Class Motif: This contains the motif sequence, motif length/"thickness", gene sequence, gene name, 
        and figure number. 
        Methods: scanner: Find and store all locations of motifs in the gene sequence. 
                This "scanner" function uses a sliding window approach to determine motif presence & location.
        """
        self.type = m_sequence.upper()
        self.gene_sequence = gene_seq.upper()
        self.position = []
        self.thickness = len(m_sequence)
        self.associated_gene = gene_name
        self.associated_figure = figure_number
        self.colors = []
    # def scanner(self):
    #     kmer_length = len(self.type)
    #     for index, chunk in enumerate(self.gene_sequence):
    #         if self.gene_sequence[index:(index+kmer_length)] == self.type:
    #             self.position.append(index)
    def scanner(self):
        kmer_length = len(self.type)
        for index, chunk in enumerate(self.gene_sequence):
            kmer = self.gene_sequence[index:(index+kmer_length)]
            count = 0
            for item in kmer:
                if item in nucleotide_dict[self.type[count]]:
                    count +=1
                    if count == kmer_length:
                        self.position.append(index)
                else:
                    break
    def random_color_generator(self):
        self.colors.append(random.random())
        self.colors.append(random.random())
        self.colors.append(random.random())
                

nucleotide_dict = { "A" : ["A", "W", "M", "R", "D", "H", "V", "N"],
                    "G" : ["G", "R", "S", "K", "B", "D", "V", "N"],
                    "T" : ["T", "U", "W", "K", "Y", "B", "D", "H", "N"],
                    "C" : ["C", "S", "M", "Y", "B", "H", "V", "N"],
                    "U" : ["T"],
                    "W" : ["A","T"],
                    "S" : ["C","G"],
                    "M" : ["A","C"],
                    "K" : ["G","T"],
                    "R" : ["A","G"],
                    "Y" : ["T","C"],
                    "B" : ["G","T","C"],
                    "D" : ["A","G","T"],
                    "H" : ["A","C","T"],
                    "V" : ["A","C","G"],
                    "N" : ["A","C","T","G"] }



# Create motifs dictionary. Key = motif_file_X, value = list(file pointer)
counter = 1
number_of_motif_files = []
for x in motif_files:
    number_of_motif_files.append(counter)
    counter +=1
    
motif_files_dict = {}
counter = 0
for i in number_of_motif_files:
    motif_files_dict['motif_file_%s' % i] = [motif_files[counter]]
    counter +=1

#Update motifs dictionary. Key = motif_file_X, value = list(file pointer, [motifs list])
for value in motif_files_dict.values():
    motif_list = []
    for line in value[0]:
        line = line.strip()
        motif_list.append(line)
    value.append(motif_list)


#print(motif_files_dict)


# Create Genes dictionary. Key = fasta_file_X, value = list(file pointer)
counter = 1
number_of_genes_files = []
for x in fasta_files:
    number_of_genes_files.append(counter)
    counter +=1
    
fasta_files_dict = {}
counter = 0
for i in number_of_genes_files:
    fasta_files_dict['fasta_file_%s' % i] = [fasta_files[counter]]
    counter +=1


#Update Genes dictionary. Key = fasta_file_X, value = list(file pointer, [gene objects])
fig_num = 0
for value in fasta_files_dict.values():
    current_sequence = ""
    fig_num +=1
    fig_x_genes = [] 
    counter = -1   
    for line in value[0]:
        list_of_lines_in_file = line.split()
        #print(list_of_lines_in_file)
        for i in range(len(list_of_lines_in_file)):
            if i ==0:
                #print(list_of_lines_in_file[i])
                if line.startswith(">"):
                    #print(line)
                    gene_name = (list_of_lines_in_file[i][1:])
                    
                    #print(current_sequence)
                    #print(gene)
                    fig_x_genes.append(Gene(current_sequence,gene_name,fig_num))
                    current_sequence = ""
                    counter +=1      
                if line.startswith(">") == False:
                    current_sequence += (list_of_lines_in_file[i])
                    #print(current_sequence)
                    fig_x_genes[counter].sequence = current_sequence
    value.append(fig_x_genes)
                

# print(fasta_files_dict["fasta_file_2"][1][3].name)

#Find attribute values for all gene instances.
for value in fasta_files_dict.values():
    for gene in value[1]:
         gene.seq_count()
         gene.IE_pos()

# print(fasta_files_dict["fasta_file_2"][1][3].name)

motif_holder = [] 
counter = 0
for value in motif_files_dict.values():
    for motif in value[1]:
        #print(motif)
        for fasta in fasta_files_dict.values():
            for gene in fasta[1]:
                motif_holder.append((Motif(motif,gene.sequence,gene.name,gene.figure_ID)))



# print(motif_holder[1].associated_figure)

for motif in motif_holder:
    motif.scanner()
    motif.random_color_generator()


for value in fasta_files_dict.values():
    #print(value[0])
    for the_gene in value[1]:
        #print(the_gene.name)
        for motif in motif_holder:
            if motif.associated_gene == the_gene.name and motif.associated_figure == the_gene.figure_ID:
                if motif.position != []:
                    print(motif.type)
                    print(motif.associated_gene)
                    print(motif.position)




#print("\n")

# print(motif_files_dict["motif_file_2"][1].type)
# print(motif_files_dict["motif_file_2"][1].position)
# print(motif_files_dict["motif_file_2"][1].thickness)
# print(motif_files_dict["motif_file_2"][1].associated_gene)
# print(motif_files_dict["motif_file_2"][1].associated_figure)




#print(fasta_files_dict["fasta_file_1"][1][0])
# print(fasta_files_dict["fasta_file_2"][1][3].name)

class Figure:
    def __init__(self, figure_number):
        self.fig_num = figure_number
        self.gene_list = []
        self.motif_list = []
    def find_gene_and_motif_lists(self):
        for value in fasta_files_dict.values():
            for the_gene in value[1]:
                if the_gene.figure_ID == self.fig_num:
                    self.gene_list.append(the_gene)
                for motif in motif_holder:
                    if motif.associated_gene == the_gene.name and motif.associated_figure == the_gene.figure_ID == self.fig_num:
                        if motif.position != []:
                            self.motif_list.append(motif)
    def beautiful_creation(self):
        spacer = 0
        height = (200 * len(self.gene_list))
        surface = cairo.PDFSurface("test_plot.pdf",1500, height)
        context = cairo.Context(surface)
        context.set_line_width(1)
        start_x = 50
        start_y = 150
        for gene in self.gene_list:
            context.set_font_size(40)
            context.move_to(40 , 50)
            figure_name = ("Figure " + str(gene.figure_ID))
            context.show_text(figure_name)
            #diplay gene name here
            current_x = start_x 
            current_y = start_y + spacer
            number_of_segments = range(len(gene.swap_points))
            context.move_to(current_x+50 , current_y)
            context.set_font_size(20)
            context.show_text(gene.name)
            current_y = start_y + spacer + 40
            if gene.start_as == "Intron":
                intron_at_current_step= True
            else:
                intron_at_current_step = False
            for segment in number_of_segments:
                #stretch_factor = 500/gene.seq_length
                for motif in self.motif_list:
                    if motif.associated_gene == gene.name and motif.associated_figure == gene.figure_ID:
                        for position in motif.position:
                            context.set_source_rgb(motif.colors[0],motif.colors[1],motif.colors[2])
                            context.rectangle(position+start_x,current_y-15,motif.thickness,30)
                            context.fill()
                context.move_to(current_x,current_y)
                if intron_at_current_step == True:
                    context.set_source_rgb(.5, .5, .1)
                    current_x = start_x + gene.swap_points[segment]
                    context.line_to(current_x,current_y)
                    context.stroke()
                    intron_at_current_step = False
                    continue
                if intron_at_current_step == False:
                    context.set_source_rgb(0,0,0)
                    context.rectangle(current_x,current_y-15,gene.swap_points[segment]-current_x+start_x,30)
                    context.fill()
                    current_x = start_x + gene.swap_points[segment]
                    intron_at_current_step = True
                    continue
                
            spacer += 120
                #Add motifs in previous sections

        surface.write_to_png("test_plot.png")
        surface.finish()

            


        # context.move_to(450,25)
        # context.line_to(450,325)
        # context.rectangle(50,50,300,350)
        # context.stroke()
        # surface.write_to_png("plot.png")
        # surface.finish()


            


# self.sequence = the_sequence
#         self.start_as = ""
#         self.swap_points = []
#         self.seq_length = 0
#         self.name = gene_name
#         self.figure_ID = figure_number




        #create a practice sideways boxplot.
        #Then figure out how to create these in a loop.
        #Then figure out how to correspond the genes and motifs to each.
        #Color should maybe be in the motif class. 
        #Last, create a legend

figures = []
for stuff in range(len(fasta_files_dict)):
    stuff +=1
    figures.append(Figure(stuff))

for figure in figures:
    figure.find_gene_and_motif_lists()


# print(figures[1].gene_list)
# print(figures[1].motif_list)

for figure in figures:
    figure.beautiful_creation()



# surface = cairo.PDFSurface("plotxx.png", 1000, 1000)
# context = cairo.Context(surface)
# context.set_line_width(1)
# context.move_to(450,25)
# context.line_to(450,325)
# context.rectangle(50,50,300,350)
# context.stroke()
# surface.write_to_png("plotxx.png")
# surface.finish()


    # Motif Attributes
        # self.type = m_sequence
        # self.gene_sequence = gene_seq
        # self.position = []
        # self.thickness = len(self.type)
        # self.associated_gene = gene_name
        # self.associated_figure = figure_number

    # Gene Attributes:
        # self.sequence = the_sequence
        # self.start_as = ""
        # self.swap_points = []
        # self.seq_length = 0
        # self.name = gene_name
        # self.figure_ID = figure_number


