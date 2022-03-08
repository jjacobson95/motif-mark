#!/usr/bin/env python
import cairo
import random
import argparse
import sys
import os


#I don't know where half of the imported modules came from. But they are there and it works so I'm keeping them.


#Hello, the code begins here. 
#First import fasta and motif files.
parser = argparse.ArgumentParser(description = "Argparse")
parser.add_argument('-f', '--files', help="Input fasta files in corresponding order to motif files", type=argparse.FileType('r'), nargs='+', required =True)
parser.add_argument('-m', '--motifs', help="Input motif files in corresponding order to fasta files", type=argparse.FileType('r'), nargs='+', required =True)

args = parser.parse_args()
fasta_files = args.files
motif_files = args.motifs


print("""
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░███╗   ███╗ █████╗ ████████╗██╗███████╗  ███╗   ███╗ █████╗ ██████╗ ██╗  ██╗░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░████╗ ████║██╔══██╗╚══██╔══╝██║██╔════╝  ████╗ ████║██╔══██╗██╔══██╗██║ ██╔╝░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░██╔████╔██║██║  ██║   ██║   ██║█████╗    ██╔████╔██║███████║██████╔╝█████═╝░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░██║╚██╔╝██║██║  ██║   ██║   ██║██╔══╝ ░  ██║╚██╔╝██║██╔══██║██╔══██╗██╔═██╗░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░██║ ╚═╝ ██║╚█████╔╝ ░ ██║ ░ ██║██║   ░░  ██║ ╚═╝ ██║██║  ██║██║  ██║██║ ╚██╗░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░╚═╝ ░ ░ ╚═╝ ╚════╝  ░ ╚═╝ ░ ╚═╝╚═╝  ░░░  ╚═╝ ░ ░ ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
""")


#progress / exit statements
print("-------------\n----------\n-------\n----\n")
print("Motif Mark can process any number of Fasta and Motif files with any number of Genes and up to 18 Motifs. \n")
print("Unused Motifs will not appear in figure legends.")
print("Output Figures will be named in order of matching files read in.\n")
if len(motif_files) != len(fasta_files):
    print("Number of Motif files must be equal to the number of Fasta files. \n Number of Fasta files found: ", len(fasta_files), "\n Number of Motif files found: ", len(motif_files))
    sys.exit("Please Try Again.")
else:
    print("Number of Motif files is equal to the number of Fasta Files. \n Proceeding. \n Number of Fasta files found: ", len(fasta_files), "\n Number of Motif files found: ", len(motif_files))


#
# Here are reference lists / dictionaries.

#Dictionary is two directional for nucleotides. So weird nucleotides can be in motif or sequence.
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

#A list of colors, a tad bit sloppy at the moment
X = 255
my_color_list = [[220/X,20/X,60/X], [255/X,215/X,0/X], [34/X,139/X,34/X], [0/X,206/X,209/X],
                 [100/X,149/X,237/X], [138/X,43/X,226/X], [255/X,20/X,147/X], [139/X,69/X,19/X], [119/X,136/X,153/X], [240/X,255/X,240/X],
                 [random.random(),random.random(),random.random()], [random.random(),random.random(),random.random()],
                 [random.random(),random.random(),random.random()], [random.random(),random.random(),random.random()],
                 [random.random(),random.random(),random.random()], [random.random(),random.random(),random.random()],
                 [random.random(),random.random(),random.random()], [random.random(),random.random(),random.random()]]


#
# My Classes are here.
# Three Classes are below: Gene, Motif, Figure

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
    def color_generator(self, col_1,col_2,col_3):
        self.colors.append(col_1)
        self.colors.append(col_2)
        self.colors.append(col_3)


class Figure:
    def __init__(self, figure_number, figure_name):
        self.fig_num = figure_number
        self.fig_name = figure_name
        self.gene_list = []
        self.motif_list = []
        self.longest_seq = 0
        """
        Class Figure: This contains the gene list, and motif list. These lists hold instances from the gene and motif class.
        Methods: find_gene_and_motif_lists: Generates gene and motif subsets for this figure from preexisting comprehensive lists.
                 find_longest_sequence: This is used to generate the figure width.
                 beautiful_creation: Create a figure based on all of the previously gathered information.
        """
    def find_gene_and_motif_lists(self):
        for gene in fasta_files_list[self.fig_num-1]:
            self.gene_list.append(gene)
            for motif_list in motif_holder:
                for motif in motif_list:
                    if motif.associated_gene == gene.name and motif.associated_figure == gene.figure_ID == self.fig_num:
                        if motif.position != []:
                            self.motif_list.append(motif)
    def find_longest_sequence(self):
        for gene in fasta_files_list[self.fig_num-1]:
            if gene.seq_length > self.longest_seq:
                self.longest_seq = gene.seq_length 
    def beautiful_creation(self):
        spacer = 0
        height = (175 * len(self.gene_list) + 100)
        width = self.longest_seq + 400
        figure_file_name = self.fig_name + ".png"
        surface = cairo.ImageSurface(self.fig_num,width, height)
        context = cairo.Context(surface)
        context.set_line_width(1)
        start_x = 50
        start_y = 150
            #Genes
        for gene in self.gene_list:
            context.set_font_size(40)
            context.move_to(40 , 50)
            context.set_source_rgb(1,1,1)
            context.show_text(self.fig_name)
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
                reset = False
                #stretch_factor = 500/gene.seq_length
                context.move_to(current_x,current_y)
                context.set_source_rgb(1,1,1)
                    #Create Introns
                if intron_at_current_step == True and reset == False:
                    context.set_line_width(1)
                    current_x = start_x + gene.swap_points[segment]
                    context.line_to(current_x,current_y)
                    context.stroke()
                    intron_at_current_step = False
                    reset = True
                    #Create Exons
                if intron_at_current_step == False and reset == False:
                    context.rectangle(current_x,current_y-15,gene.swap_points[segment]-current_x+start_x,30)
                    context.set_source_rgb(1,1,1)
                    context.fill_preserve()
                    context.set_source_rgb(.721,.525,.0431)
                    context.set_line_width(1.75)
                    context.stroke()
                    current_x = start_x + gene.swap_points[segment]
                    intron_at_current_step = True
                    reset = True
                    #Create Motifs
                for motif in self.motif_list:
                    if motif.associated_gene == gene.name and motif.associated_figure == gene.figure_ID ==self.fig_num:
                        context.set_source_rgb(motif.colors[0],motif.colors[1],motif.colors[2])
                        for position in motif.position:
                            context.rectangle(position+start_x,current_y-10,motif.thickness,20)
                            context.fill()
            spacer += 140 
           #Create Legend 
        legend_region_x = width - 300
        legend_region_y = 100       
        legend_x_position =  legend_region_x  
        legend_y_position = legend_region_y
        used_motifs = []
        for motif in self.motif_list:
            if motif.type not in used_motifs:
                context.set_source_rgb(motif.colors[0],motif.colors[1],motif.colors[2])
                context.move_to(legend_x_position + 20, legend_y_position +20)
                context.show_text(motif.type)
                context.rectangle(legend_x_position,legend_y_position,12,20)
                context.fill()
                legend_y_position += 45
                used_motifs.append(motif.type)
        surface.write_to_png(figure_file_name)
        surface.finish()
        context.set_source_rgb(0.0, 0.0, 0.0)
        context.set_operator(cairo.OPERATOR_CLEAR)



#
# Code to get all of the objects to do what they are supposed to do begins here.
# I'll first create gene objects
# Then I'll create motif objects
# Then, I'll join them together in figure objects
# Last I'll run my figure.beautiful creation function to generate PNG images.


# GENES:
# Create Genes dictionary. Key = fasta_file_X, value = list(file pointer)

number_of_genes_files = []
for x in fasta_files:
    number_of_genes_files.append(x.name)
    
fasta_files_dict = {}
counter = 0
for i in number_of_genes_files:
    fasta_files_dict[os.path.splitext(i)[0]] = [fasta_files[counter]]
    counter +=1

# Update Genes dictionary. Key = fasta_file_X, value = list(file pointer, [gene objects])
#I also ended up creating a genes list here because the dictionary was hard to reference in Figure instances. 
fig_num = 0
fasta_files_list = []
for value in fasta_files_dict.values():
    current_sequence = ""
    fig_num +=1
    fig_x_genes = [] 
    counter = -1   
    for line in value[0]:
        list_of_lines_in_file = line.split()
        for i in range(len(list_of_lines_in_file)):
            if i ==0:
                if line.startswith(">"):
                    gene_name = (list_of_lines_in_file[i][1:])
                    fig_x_genes.append(Gene(current_sequence,gene_name,fig_num))
                    current_sequence = ""
                    counter +=1      
                if line.startswith(">") == False:
                    current_sequence += (list_of_lines_in_file[i])
                    fig_x_genes[counter].sequence = current_sequence
    value.append(fig_x_genes)
    fasta_files_list.append(fig_x_genes)
                
#Find attribute values for all gene instances.
for value in fasta_files_dict.values():
    for gene in value[1]:
         gene.seq_count()
         gene.IE_pos()


print(" Gene library generated.")

# MOTIFS:
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
#I also ended up creating a genes list here because the dictionary was hard to reference in Figure instances. 
motif_files_list = []
for value in motif_files_dict.values():
    motif_list = []
    for line in value[0]:
        line = line.strip()
        motif_list.append(line)
    value.append(motif_list)
    motif_files_list.append(motif_list)

#Generate motif instances
motif_holder = [] 
iterations_to_do = range(len(fasta_files_list))
for number in iterations_to_do:
    for gene in fasta_files_list[number]:
        temp_list = []
        for motif in motif_files_list[number]:
            temp_list.append((Motif(motif,gene.sequence,gene.name,gene.figure_ID)))
        motif_holder.append(temp_list)

#Generate motif attributes
reset_color = False
counter = 0
for motif_list in motif_holder:
    if reset_color == True:
        counter = 0
    for motif in motif_list:
        motif.scanner()
        motif.color_generator(my_color_list[counter][0],my_color_list[counter][1],my_color_list[counter][2])
        counter +=1
    reset_color = True

print(" Motif library generated.")
#FIGURES:
#Create list that holds figure instances
figures = []
counter = 1
for f_files in fasta_files_dict:
    figures.append(Figure(counter,f_files))
    counter +=1

#Creation of end product.
counter = 1
for figure in figures:
    figure.find_gene_and_motif_lists()
    figure.find_longest_sequence()
    figure.beautiful_creation()
    print(" ", figure.fig_name, " complete!")
    counter +=1
    
print("\nCurrent bug: Figures beyond the first may be in black and white.")
    #To DO
     #Order of things to work on:
        #2: make figure 2 have color.   


print("\n----\n-------\n----------\n-------------")
print("Author: Jeremy Jacobson\n") 