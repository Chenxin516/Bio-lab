# File:             geneMatch_original.py
# Author:           Akaash Venkat, Audi Liu

from selenium import webdriver
from os.path import expanduser
import glob
import math
import random
import os
import subprocess
import sys
import time

GENE_DATABASE_FILE = "info_files/gene_database.txt"
GENE_GROUP_FILE = "info_files/gene_group.txt"
UNIDENTIFIABLE_GENE_FILE = "info_files/unidentifiable_genes.txt"
CHANGED_NAME_GENE_FILE = "info_files/changed_name_genes.txt"

GENE_LIST = []
UNIDENTIFIABLE_LIST = []
CHANGED_NAME = {}
GROUP = {}
B_D_PAIR = {}




    
def readDatabase():
    with open(GENE_DATABASE_FILE) as database_file:
        for line_content in database_file:
            if line_content != "\n":
                gene_info = []
                line_content = line_content.replace(" ", "")
                line_content = line_content.replace(")", "")
                line_content = line_content.replace("\n", "")
                temp_list = line_content.split("-",1)
                main_gene = temp_list[0]
                connecting_genes_list = temp_list[1].split(",")
                neighbors = {}
                for connecting_gene in connecting_genes_list:
                    name = connecting_gene.split("(")[0] 
                    num = float(connecting_gene.split("(")[1])
                    neighbors[name] = num
                gene_info.append(main_gene)
                gene_info.append(neighbors)
                GENE_LIST.append(gene_info)
    database_file.close()


def readUnidentifiable():
    with open(UNIDENTIFIABLE_GENE_FILE) as unidentifiable_file:
        for line_content in unidentifiable_file:
            line_content = line_content.replace(" ", "")
            line_content = line_content.replace("\n", "")
            if line_content != "" and "followinggenescannotbefound" not in line_content:
                UNIDENTIFIABLE_LIST.append(line_content)
    unidentifiable_file.close()


def readChangedName():
    with open(CHANGED_NAME_GENE_FILE) as changed_name_file:
        for line_content in changed_name_file:
            line_content = line_content.replace(" ", "")
            line_content = line_content.replace("\n", "")
            if line_content != "" and "followinggeneshavebeenrenamed" not in line_content:
                orig_name = line_content.split("=>")[0]
                new_name = line_content.split("=>")[1]
                CHANGED_NAME[orig_name] = new_name
    changed_name_file.close()


def writeToDatabase():
    os.system('rm ' + GENE_DATABASE_FILE)
    os.system('touch ' + GENE_DATABASE_FILE)
    database_file = open(GENE_DATABASE_FILE, "w")
    for counter in range(0, len(GENE_LIST)):
        gene_info = GENE_LIST[counter]
        main_gene = gene_info[0]
        connecting_genes_list = gene_info[1]
        line_content = main_gene + " - "
        for key, value in sorted(connecting_genes_list.items() ):
            line_content = line_content + key + "(" + str(value) + "), " 
        line_content = line_content[:-2]
        database_file.write(line_content + "\n\n")
    database_file.close()


def writeGeneGroups():
    os.system('touch ' + GENE_GROUP_FILE)
    grouping_file = open(GENE_GROUP_FILE, "w")

    groups = ["A", "B", "C", "D"]
    descriptions = ["Input gene that has direct connection with another input gene", "Input gene that is indirectly connected to another input gene, via an intermediate gene", "Input gene that is not directly or indirectly connected to another input gene", "Intermediate gene that connects Group B genes with Group A or other Group B genes"]

    for counter in range(0, len(groups)):
        group_id = groups[counter]
        description = descriptions[counter]
        grouping_file.write("Group " + group_id + ": " + description + "\n")
        grouping_file.write("---\n")
        cluster = getListForGroup(group_id)
        for gene in cluster:
            grouping_file.write(gene + "\n")
        grouping_file.write("\n\n\n")

    grouping_file.close()


def writeUnidentifiable():
    os.system('touch ' + UNIDENTIFIABLE_GENE_FILE)
    cleaned_unidentifiable_list = []
    for gene in UNIDENTIFIABLE_LIST:
        if gene not in cleaned_unidentifiable_list:
            cleaned_unidentifiable_list.append(gene)
    if len(cleaned_unidentifiable_list) != 0:
        unidentifiable_file = open(UNIDENTIFIABLE_GENE_FILE, "w")
        unidentifiable_file.write("The following genes cannot be found on the online STRING database, and will not be used in this program:\n\n")
        for gene in cleaned_unidentifiable_list:
            unidentifiable_file.write(gene + "\n")
        unidentifiable_file.close()
    else:
        os.system('rm ' + UNIDENTIFIABLE_GENE_FILE)


def writeChangedName():
    os.system('touch ' + CHANGED_NAME_GENE_FILE)
    if not CHANGED_NAME:
        os.system('rm ' + CHANGED_NAME_GENE_FILE)
    else:
        changed_name_file = open(CHANGED_NAME_GENE_FILE, "w")
        changed_name_file.write("The following genes have been renamed, as per the online STRING database:\n\n")
        for key, value in CHANGED_NAME.iteritems():
            changed_name_file.write(key + " => " + value + "\n")
        changed_name_file.close()


def initialize_connections():
    for gene_info in GENE_LIST:
        gene = gene_info[0]
        GROUP[gene] = "C"


def identifyGroupA(gene_list):
    for i in range(0, len(gene_list)):
        gene = gene_list[i][0]
        gene_neighbors = gene_list[i][1]
        
        for j in range(0, len(gene_list)):
            if i == j:
                continue
            
            other_gene = gene_list[j][0]
            
            if other_gene in gene_neighbors.keys():
                GROUP[gene] = "A"


def identifyGroupB(gene_list):
    for i in range(0, len(gene_list)):
        
        content_list = []
        gene = gene_list[i][0]
        gene_neighbors = gene_list[i][1]
        
        
        if GROUP[gene] == "A":
            continue
        
        for j in range(0, len(gene_list)):
            if i == j:
                continue

            other_gene = gene_list[j][0]
            other_gene_neighbors = gene_list[j][1]

            for inter_gene in gene_neighbors.keys():
                if inter_gene in other_gene_neighbors.keys():
                    if gene_neighbors[inter_gene] > other_gene_neighbors[inter_gene]:
                        content_list.append([other_gene_neighbors[inter_gene], inter_gene, other_gene])
                    else:
                        content_list.append([gene_neighbors[inter_gene], inter_gene, other_gene])
                                             
        best_match = max(content_list)
        if GROUP[gene] == "C":
            GROUP[gene] = "B"
            GROUP[best_match[1]] = "D"
            B_D_PAIR[gene] = best_match[1]


def parseInput():
    os.system('clear')
    arg = ""
    content = []
    input_genes = []
    input_list = []
    if (len(sys.argv) == 1):
        print("Please try running program again, this time adding an argument.")
        sys.exit(0)
    else:
        arg = sys.argv[1]
    if (arg[-4:] == '.txt'): #checking the file extensions
        with open(arg) as arg_input:
            content = arg_input.readlines()
    else:
        print("Please try running program again, this time passing in an argument of type '.txt'")
        sys.exit(0)

    for gene in content:
        input_list.append(gene.replace(" ", "").replace("\n", ""))

    for gene in input_list:
        if gene not in input_genes:
            input_genes.append(gene)

    current_gene_list = GENE_LIST

    for gene in input_genes:
        already_present = False
        for iter in range(0, len(current_gene_list)):
            existing_gene = current_gene_list[iter][0]
            if existing_gene == gene:
                already_present = True
                break

        if already_present == False:
            
            if gene in CHANGED_NAME:
                break
            if gene in UNIDENTIFIABLE_LIST:
                break
            
            gene_info = []
            gene_neighbors = find_neighbor(gene)

            if gene_neighbors == -1:
                UNIDENTIFIABLE_LIST.append(gene)
            else:
                if isinstance(gene_neighbors, basestring):
                    correct_gene = gene_neighbors
                    CHANGED_NAME[gene] = correct_gene
                    gene = correct_gene
                    gene_neighbors = find_neighbor(gene)

                gene_info.append(gene)
                if "" in gene_neighbors:
                    time.sleep(1)
                    gene_info.append(find_neighbor(gene))
                else:
                    gene_info.append(gene_neighbors)
                GENE_LIST.append(gene_info)
        GENE_LIST.sort()
        writeToDatabase()

    initialize_connections()
    identifyGroupA(GENE_LIST)
    identifyGroupB(GENE_LIST)



def getListForGroup(group_id):
    cluster = []
    for gene in GROUP:
        if GROUP[gene] == group_id:
            cluster.append(gene)
    return cluster



def find_neighbor(input_gene):
    gene_connectors = {}
    driver = webdriver.Chrome()
    driver.get("http://string-db.org/")
    driver.find_element_by_id("search").click()
    driver.find_element_by_id("primary_input:single_identifier").send_keys(input_gene)
    driver.find_element_by_id("species_text_single_identifier").send_keys("Homo sapiens")
    driver.find_element_by_xpath("//*[@id='input_form_single_identifier']/div[4]/a").click()
    if "Sorry, STRING did not find a protein" in driver.page_source:
        return -1
    if "Please select one from the list below" in driver.page_source:
        driver.find_element_by_xpath("//*[@id='proceed_form']/div[1]/div/div[2]/a[2]").click()
    driver.find_element_by_id("bottom_page_selector_settings").click()
    driver.find_element_by_xpath("//*[@id='bottom_page_selector_legend']").click()
    page_data = driver.page_source
    split1 = page_data.split("<td class=\"td_name middle_row first_row last_row\" onclick=")
    split2 = split1[1].split("</td>")
    split3 = split2[0].split("\">")
    correct_gene_name = split3[1]
    if input_gene != correct_gene_name:
        return str(correct_gene_name)
    driver.find_element_by_xpath("//*[@id='bottom_page_selector_table']").click()
    driver.find_element_by_id("bottom_page_selector_settings").click()
    driver.find_element_by_xpath("//*[@id='standard_parameters']/div/div[1]/div[3]/div[2]/div[2]/div[1]/label").click()
    driver.find_element_by_xpath("//select[@name='limit']/option[text()='custom value']").click()
    driver.find_element_by_id("custom_limit_input").clear()
    driver.find_element_by_id("custom_limit_input").send_keys("500")
    time.sleep(5)
    driver.find_element_by_xpath("//*[@id='standard_parameters']/div/div[1]/div[5]/a").click()
    time.sleep(10)
    driver.find_element_by_id("bottom_page_selector_table").click()
    driver.find_element_by_xpath("//*[@id='bottom_page_selector_legend']").click()
    connectors = driver.find_elements_by_class_name("linked_item_row")
    for connector in connectors:
        neighbor = str(connector.text.split(' ')[0].split('\n')[0])
        confidence_value = str(connector.text.split(' ')[-1].split('\n')[-1])
        gene_connectors[neighbor] = float(confidence_value)
    return gene_connectors



def download_svg(gene_list):

    if len(gene_list) < 2:
        return -1
    
    SVG_STRING = ""
    for gene in gene_list:
        SVG_STRING = SVG_STRING + gene + "\n"
    
    driver = webdriver.Chrome()
    driver.get("http://string-db.org/")
    driver.find_element_by_id("search").click()
    driver.find_element_by_id("multiple_identifiers").click()
    driver.find_element_by_id("primary_input:multiple_identifiers").send_keys(SVG_STRING)
    driver.find_element_by_id("species_text_multiple_identifiers").send_keys("Homo sapiens")
    driver.find_element_by_xpath("//*[@id='input_form_multiple_identifiers']/div[5]/a").click()
    time.sleep(5)
    if "The following proteins in" in driver.page_source and "appear to match your input" in driver.page_source:
        driver.find_element_by_xpath("//*[@id='proceed_form']/div[1]/div/div[2]/a[2]").click()
    time.sleep(20)
    driver.find_element_by_id("bottom_page_selector_table").click()
    time.sleep(5)
    driver.find_element_by_id("bottom_page_selector_settings").click()
    time.sleep(15)
    driver.find_element_by_id("confidence").send_keys(" ")
    time.sleep(10)
    driver.find_element_by_id("block_structures").send_keys(" ")
    time.sleep(10)
    driver.find_element_by_xpath("//*[@id='standard_parameters']/div/div[1]/div[5]/a").click()
    time.sleep(15)
    driver.find_element_by_xpath("//*[@id='bottom_page_selector_legend']").click()
    time.sleep(10)
    driver.find_element_by_id("bottom_page_selector_table").click()
    time.sleep(25)
    driver.find_element_by_xpath("//*[@id='bottom_page_selector_table_container']/div/div[2]/div/div[3]/div[2]/a").click()
    time.sleep(10)

    terminal_output = subprocess.Popen(['pwd'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = terminal_output.communicate()
    pwd = stdout[:-1]
    downloads_dir = expanduser("~") + '/Downloads/*'
    downloaded_files = glob.glob(expanduser("~") + '/Downloads/*')
    recent_file = max(downloaded_files, key=os.path.getctime)
    move_command = 'mv "' + recent_file + '" "' + pwd + '/svg_files/_base.svg"'
    os.system(move_command)




def getGene(str):
    lst1 = str.split('>')
    lst2 = lst1[1].split('<')
    return lst2[0]



def writeToFile(content, file_name):
    os.system('touch ' + file_name)
    os.system('rm ' + file_name)
    os.system('touch ' + file_name)
    file = open(file_name, "w")
    for counter in range(0, len(content)):
        file.write(str(content[counter]))
    file.close()


"""
    New Idea for Tripod of A's:
    "Hardcode" positions for 3 main genes (ABCA4, RPE65, RHO)
    Split area into 3 zones and use the random generating in each zone
    Venn Diagram into 7 parts (A, B, C, AB, BC, AC, ABC)
        3 extra loser zones
        3 extra zones for the main ones (increase radius maybe, change color maybe)
"""

    
    


def modify_base_svg():
    grey = "(217,217,217)"
    yellow = "(255,255,0)"
    red = "(255,0,0)"
    green = "(0,153,0)"
    blue = "(52,152,219)"

    
    with open("svg_files/_base.svg") as base:
        content = base.readlines()

    content[0] = "<svg class=\"notselectable\" height=\"2250\" id=\"svg_network_image\" width=\"2400\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"

    height = 2250
    width = 2400   #

    old_pos_dict = {}

    length = len(content)
    for i in range(length): #change color
        gene_name = ""
        
        if "ellipse cx=" in content[i]:
            content[i] = ""
        
        if "circle class" in content[i]:
            
            gene_name = getGene(content[i+3])
            rgb_val = content[i].split("rgb")[1].split('"')[0]
            position_data = content[i].split("cx=")[1].split(" display")[0].replace("cy=", "").replace("\"", "")
            old_pos_dict[position_data] = gene_name
            
            if gene_name == "ABCA4" or gene_name == "RPE65" or gene_name == "RHO":
                content[i] = content[i].replace(rgb_val, blue)
            elif GROUP[gene_name] == "A":
                content[i] = content[i].replace(rgb_val, grey)
            elif GROUP[gene_name] == "B":
                content[i] = content[i].replace(rgb_val, yellow)
            elif GROUP[gene_name] == "C":
                content[i] = content[i].replace(rgb_val, red)
            elif GROUP[gene_name] == "D":
                content[i] = content[i].replace(rgb_val, green)

    a_list = getListForGroup("A")
    b_list = getListForGroup("B")
    c_list = getListForGroup("C")


    base_b_y = int(height / 6)
    lowest_d_y = 0
    new_pos_dict = {}
    

    b_length = len(b_list)
    for i in range(b_length):#put gene B and D new position into new_post_dict
        b_x = int((i + 1) * width / (b_length + 1))
        b_y = base_b_y + (80 / b_length) * i * math.pow(-1, i) + (b_length / 4.75) * i * math.pow((b_length / 8.3), i) #8.3
        d_y = b_y + int(height / 12)
        new_pos_dict[b_list[i]] = [str(b_x), str(b_y)]
        d_gene = B_D_PAIR[b_list[i]]
        if d_y > lowest_d_y:
            lowest_d_y = d_y
        new_pos_dict[d_gene] = [str(b_x), str(d_y)]

    c_length = len(c_list)
    for i in range(c_length): # c new positions
        c_x = int((i + 1) * width / (c_length + 1))
        c_y = int(height / 12)
        new_pos_dict[c_list[i]] = [str(c_x), str(c_y)]
        
    inv_old_pos_dict = {v: k for k, v in old_pos_dict.iteritems()}
    
    '''
    #The bounds to put the A's in
    new_x_list = []
    new_y_list = []

    x_left_bound = int(width/(b_length + 1))
    x_right_bound = int(width/(b_length + 1)*b_length) 
        
    y_up_bound =  int(height*11/12)
    y_low_bound = int(lowest_d_y + 100)
    
    '''
    '''
    a_list2 = a_list
    
    a_list2.remove("ABCA4")
    a_list2.remove("RPE65")
    a_list2.remove("RHO")
    '''
    
    new_pos_dict["ABCA4"] = [str(1200), str(775)] # 
    new_pos_dict["RPE65"] = [str(675), str(1500)] 
    new_pos_dict["RHO"] = [str(1725), str(1500)] 
    #center for :  1200,775 ((2250/3 - 50 + 2250/3 + 100)/2)
    #center for this rectangle: 675, 1500
    #center for this rectangle: 1725, 1500

   
    l = count_A_group() # list now has a b c ab bc ac abc loner list of group a genes
    height = 2250
    
    #a only : x  750 to 1650, y height/3 - 50 to height/3 + 100
    #b only : x 600 to 750, y height/3 + 100 to  height/3 + 1400
    #c only : x 1650 to 1800, y height/3 + 100 to  height/3 + 1400
    
    #ab: x 750 to 950, y height/3 + 100 to height/3 + 1300,
    #bc: x 950 to 1450, y height/3 + 1300 to height/3 +1400
    #ac: x 1450 to 1650, height/3 + 100 to height/3 + 1300,
    
    #abc: x 950 to 1450 y height/3 + 100 to height/3 + 1300
    
    #loners_left: 100 500, height/3 to  height/3 + 1500
    #loners_right:1900 2300  height/3 to  height/3 + 1500
    
    #Changing the location in new_pos_dict
    
    if l[0][0] in new_pos_dict.keys():
        print("not there")
    distribute_points(new_pos_dict, l[0],750, 1650 ,height/3 - 50 , height/3 +100 , True, 1200, 775)
    if l[0][0] in new_pos_dict.keys():
        print("now here")
        
    distribute_points(new_pos_dict, l[1],600, 750 ,height/3 + 100 , height/3 + 1400 , True, 675, 1500)
    distribute_points(new_pos_dict, l[2], 1650, 1800, height/3 + 100 , height/3 + 1400 , True, 1725, 1500)
    
    
    distribute_points(new_pos_dict, l[3], 750, 950 , height/3 + 100 , height/3 + 1300 , False, 0, 0)
    distribute_points(new_pos_dict, l[4], 950, 1450, height/3 + 1300, height/3 +1400, False, 0, 0)
    distribute_points(new_pos_dict, l[5], 1450, 1650, height/3 + 100, height/3 + 1300 , False, 0, 0)
    
    
    distribute_points(new_pos_dict, l[6],950 , 1450 ,height/3 + 100 ,height/3 + 1300 , False,0,0)
    distribute_points(new_pos_dict, l[7], 100, 500, height/3, height/3 + 1500, False, 0, 0)
    distribute_points(new_pos_dict, l[8], 1900, 2300, height/3, height/3 + 1500, False, 0, 0)

    
    

    #actually changing the content in the file 
    for i in range(length):

            
        if "line class=\"nw_edge\"" in content[i] and ".0\" stroke=" in content[i]:
            content[i] = ""
        
        
        if "line class=\"nw_edge\"" in content[i] and ".1\" stroke=" in content[i]: #updating the edge pos
            old_x1 = content[i].split("x1=\"")[1].split("\"")[0]
            old_y1 = content[i].split("y1=\"")[1].split("\"")[0]
            old_x2 = content[i].split("x2=\"")[1].split("\"")[0]
            old_y2 = content[i].split("y2=\"")[1].split("\"")[0]
            
            mod_old_x1 = str(float(old_x1) - 0.5)
            mod_old_x2 = str(float(old_x2) - 0.5)
            mod_old_y1 = str(float(old_y1) - 0.5)
            mod_old_y2 = str(float(old_y2) - 0.5)
            
            if ".0" in mod_old_x1:
                mod_old_x1 = mod_old_x1[:-2]
            if ".0" in mod_old_x2:
                mod_old_x2 = mod_old_x2[:-2]
            if ".0" in mod_old_y1:
                mod_old_y1 = mod_old_y1[:-2]
            if ".0" in mod_old_y2:
                mod_old_y2 = mod_old_y2[:-2]
            
            gene1_name = old_pos_dict[str(mod_old_x1) + " " + str(mod_old_y1)]
            gene2_name = old_pos_dict[str(mod_old_x2) + " " + str(mod_old_y2)]
            
            if gene1_name in new_pos_dict and gene2_name in new_pos_dict:
                new_pos1 = new_pos_dict[gene1_name]
                new_pos2 = new_pos_dict[gene2_name]
                updated_new_pos1 = [str(float(new_pos1[0]) + 0.5), str(float(new_pos1[1]) + 0.5)]
                updated_new_pos2 = [str(float(new_pos2[0]) + 0.5), str(float(new_pos2[1]) + 0.5)]
                first_half = content[i].split(" x1=")[0]
                second_half = content[i].split("/>")[1]
                content[i] = first_half + " x1=\"" + updated_new_pos1[0] + "\" y1=\"" + updated_new_pos1[1] + "\" x2=\"" + updated_new_pos2[0] + "\" y2=\"" + updated_new_pos2[1] + "\" />" + second_half
                
        
        if "<circle cx" in content[i]: #update the node pos
            old_x = content[i].split("cx=\"")[1].split("\"")[0]
            old_y = content[i].split("cy=\"")[1].split("\"")[0]
            gene_name = old_pos_dict[str(old_x) + " " + str(old_y)]
            if gene_name in new_pos_dict:
                new_pos = new_pos_dict[gene_name]
                first_half = content[i].split(" cx=")[0]
                second_half = content[i].split("fill=")[1]
                content[i] = first_half + " cx=\"" + new_pos[0] + "\" cy=\"" + new_pos[1] + "\" fill=" + second_half
        
        if "<circle class" in content[i]: #update the node pos too, need to change multiple places
            old_x = content[i].split("cx=\"")[1].split("\"")[0]
            old_y = content[i].split("cy=\"")[1].split("\"")[0]
            gene_name = old_pos_dict[str(old_x) + " " + str(old_y)]
            if gene_name in new_pos_dict:
                new_pos = new_pos_dict[gene_name]
                first_half = content[i].split(" cx=")[0]
                second_half = content[i].split("display=")[1]
                content[i] = first_half + " cx=\"" + new_pos[0] + "\" cy=\"" + new_pos[1] + "\" display=" + second_half
                    
        if "<text " in content[i]: #update the text pos
            old_text_x = content[i].split("x=\"")[1].split("\"")[0]
            old_text_y = content[i].split("y=\"")[1].split("\"")[0]
            old_x = str(float(old_text_x) - 18)
            old_y = str(float(old_text_y) + 18)
            
            if ".0" in old_x:
                old_x = old_x[:-2]
            if ".0" in old_y:
                old_y = old_y[:-2]
            
            gene_name = old_pos_dict[str(old_x) + " " + str(old_y)]
            if gene_name in new_pos_dict:
                new_pos = new_pos_dict[gene_name]
                new_text_pos = [str(float(new_pos[0]) + 18), str(float(new_pos[1]) - 18)]
                
                first_half = content[i].split("x=")[0]
                second_half = "x=\"" + new_text_pos[0] + "\" y=\"" + new_text_pos[1] + "\">" + gene_name + "</text>\n"
                
                content[i] = first_half + second_half
    
    
    #Visualize the boxes
    index = 0
    for i in range(length):
         if "g id=\"nodes\"" in content[i]:
             index = i - 1
             break
    
    #b only : x 600 to 750, y height/3 + 100 to  height/3 + 1400
    #c only : x 1650 to 1800, y height/3 + 100 to  height/3 + 1400
    #a only : x  750 to 1650, y height/3 - 50 to height/3 + 100
    
    #ab: x 750 to 950, y height/3 + 100 to height/3 + 1300,
    #ac: x 1450 to 1650, height/3 + 100 to height/3 + 1300,
    #bc: x 950 to 1450, y height/3 + 1300 to height/3 +1400
    
    #abc: x 950 to 1450 y height/3 + 100 to height/3 + 1300
    
    #left loners 
    insert_line(100, height/3, 500, height/3, content, index)
    insert_line(500, height/3, 500, height/3 + 1500, content, index)
    insert_line(100, height/3 + 1500, 500, height/3 + 1500, content, index)
    insert_line(100, height/3, 100, height/3 + 1500, content, index)
    
    
    #right loners
    insert_line(1900, height/3, 2300, height/3, content, index)
    insert_line(2300, height/3, 2300, height/3 + 1500, content, index)
    insert_line(2300, height/3 + 1500, 1900, height/3 + 1500, content, index)
    insert_line(1900, height/3 + 1500, 1900, height/3, content, index)
    
    
    #A middle
    insert_line(750, height/3 - 50, 1650, height/3 - 50, content, index)
    insert_line(1650, height/3 - 50, 1650, height/3 + 1300, content, index)
    insert_line(1650, height/3 + 1300, 750, height/3 + 1300, content, index)
    insert_line(750, height/3 + 1300, 750, height/3 - 50, content, index)
    
    
    
    #B left 
    insert_line(600, height/3 + 100, 1450, height/3 + 100, content, index)
    insert_line(1450, height/3 + 100, 1450, height/3 + 1400, content, index)
    insert_line(1450, height/3+1400, 600, height/3 + 1400, content, index)
    insert_line(600, height/3 +1400, 600, height/3 + 100, content, index)
    


    #C right
    insert_line(950, height/3 + 100, 1800, height/3 + 100, content, index)
    insert_line(1800, height/3 + 100, 1800, height/3 + 1400, content, index)
    insert_line(1800, height/3+1400, 950, height/3 + 1400, content, index)
    insert_line(950, height/3 +1400, 950, height/3 + 100, content, index)
    
    writeToFile(content, "svg_files/_modified_base.svg")


    #the function to distribute the points in a rectangle

def distribute_points(new_pos_dict, a_list, x_left_bound, x_right_bound, y_low_bound, y_up_bound, iso_zone, x_center, y_center):
    new_x_list = []
    new_y_list = []
    key_num = random.SystemRandom()

    for i in range(len(a_list)): 
        a_x_val = -1
        a_y_val = -1
        
        while True: #the statement in the while loop is just for one pick
            count = 0 #number of times we've checked the items ti make sure not overlapping with chosen points

            A_new_x = key_num.randint( x_left_bound, x_right_bound )
            A_new_y = key_num.randint( y_low_bound, y_up_bound )
            
            if iso_zone: 
                dist2 = math.sqrt( abs(x_center -  A_new_x) ** 2 + abs(y_center -  A_new_y) **2  )
                if dist2 < 42:
                    continue

            for j in range(len(new_y_list)):
                dist = math.sqrt(abs(new_x_list[j] -  A_new_x) * abs(new_x_list[j] -  A_new_x) + abs(new_y_list[j] -  A_new_y) * abs(new_y_list[j] -  A_new_y))
                if (dist < 57):
                    break #failed, generate a new number
                count = count + 1
                
            if count == len(new_y_list):
                a_x_val = A_new_x
                a_y_val = A_new_y
                break #when count in 
                
            
        new_x_list.append(a_x_val)
        new_y_list.append(a_y_val)
        new_pos_dict[a_list[i]] = [str(a_x_val), str(a_y_val)] 
        

   
def insert_line(x1, y1, x2, y2, content, i):
    first = "\t\t<g class=\"nwlinkwrapper\" id=\"edge.1111111.1111111\">\n"
    second =	"	\t\t<line class=\"nw_edge\" id=\"line.1111111.1111111.1\" stroke=\"#FF4500\" stroke-dasharray=\"none\" stroke-opacity=\"0.362\" stroke-width=\"7\" style=\"\" x1=\"" + str(x1) + "\" y1=\""+ str(y1)+ "\" x2=\"" + str(x2) + "\" y2=\"" + str(y2)+ "\" />\n"
    
    #index = n
    #second_2 =  "n" + "\" stroke=\"#0000FF\" stroke-dasharray=\"none\" stroke-opacity=\"0.362\" stroke-width=\"5\" style=\"\" x1=\"1609.5\" y1=\"1353.5\" x2=\"489.5\" y2=\"2051.5\" />"
    third = "		</g>\n"
    content.insert(i,third)
    content.insert(i,second)
    content.insert(i,first)
    
    
    

def create_svg(b_gene):

    with open("svg_files/_modified_base.svg") as modified_base:
        content = modified_base.readlines()
    
    naughty_b_list = []
    naughty_d_list = []
    b_list = getListForGroup("B")
    for b in b_list:
        if b != b_gene:
            naughty_b_list.append(b)
            naughty_d_list.append(B_D_PAIR[b])

    naughty_b_numbers = []
    naughty_d_numbers = []

    length = len(content)
    for i in range(length):
        if "<text fill=\"white\"" in content[i]:
            gene_line = content[i].split("</text>")[0].split(">")[1]
            
            if gene_line in naughty_b_list or gene_line in naughty_d_list:
                naughty_x = content[i-1].split("cx=\"")[1].split("\"")[0]
                naughty_y = content[i-1].split("cy=\"")[1].split("\"")[0]
                
                naughty_x = str(float(naughty_x) + 0.5)
                naughty_y = str(float(naughty_y) + 0.5)
                
                if gene_line in naughty_b_list:
                    naughty_b_numbers.append(naughty_x + " " + naughty_y)
                if gene_line in naughty_d_list:
                    naughty_d_numbers.append(naughty_x + " " + naughty_y)
                    for j in range(i - 4, i + 3):
                        content[j] = ""
    for i in range(length):
        if "<line class" in content[i]:
            x1 = content[i].split("x1=\"")[1].split("\"")[0]
            y1 = content[i].split("y1=\"")[1].split("\"")[0]
            pos1 = x1 + " " + y1
            x2 = content[i].split("x2=\"")[1].split("\"")[0]
            y2 = content[i].split("y2=\"")[1].split("\"")[0]
            pos2 = x2 + " " + y2
            if pos1 in naughty_b_numbers or pos2 in naughty_b_numbers or pos1 in naughty_d_numbers or pos2 in naughty_d_numbers:
                for j in range(i - 1, i + 2):
                    content[j] = ""

    writeToFile(content, "svg_files/" + b_gene + ".svg")


def count_A_group(): #ABCA4, RPE65, RHO 
    a = []
    b = []
    c = []
    ab = []
    bc = []
    ac = []
    abc = []
    loner = []
    
    A_genes = getListForGroup("A")
    print(len(A_genes))
    for gene in A_genes:
        for l in GENE_LIST:
            if gene == l[0]:
                #print(l[0])
                #print(l[1])

                genes_connected = l[1].keys()
                A =  "ABCA4" in genes_connected 
                B = "RPE65" in genes_connected 
                C = "RHO" in genes_connected 
                if gene != "ABCA4" and gene != "RPE65" and gene != "RHO":
                    if (A and B and C):
                        abc.append(gene)
                    elif (A and B and not C):
                        ab.append(gene)
                    elif (A and C and not B):
                        ac.append(gene)
                    elif (B and C and not A):
                        bc.append(gene)
                    elif (A and not B and not C):
                        a.append(gene)
                    elif (B and not A and not C):
                        b.append(gene)
                    elif (C and not A and not B):
                        c.append(gene)
                    else:
                        loner.append(gene)
    
    print("A:" +  str(len(a)))
    print("B:" + str(len(b)))
    print("C:" +  str(len(c)))
    print("AB:" + str(len(ab)))
    print("BC:" + str(len(bc)))
    print("AC:" + str(len(ac)))
    print("ABC:" + str(len(abc)))
    print("loner:" + str(len(loner)))
    
    print(len(a) + len(b) +len(c) + len(ab) + len(bc) + len(ac) + len(abc) + len(loner) )
    l = []
    l.append(a)
    l.append(b)
    l.append(c)
    l.append(ab)
    l.append(bc)
    l.append(ac)
    l.append(abc)
    
    half = len(loner)/2
    loners_left = loner[:half]
    loners_right = loner[half:]

    l.append(loners_left)
    l.append(loners_right)
    return l 
                    
    #print(GENE_LIST) # a list of lists (that contains the gene name and the dictionary)

def main():   
    os.system('mkdir info_files')
    os.system('mkdir svg_files')
    os.system('touch ' + GENE_DATABASE_FILE)
    os.system('touch ' + UNIDENTIFIABLE_GENE_FILE)
    os.system('touch ' + CHANGED_NAME_GENE_FILE)
    readDatabase()
    readUnidentifiable()
    readChangedName()
    
    
    parseInput()
    writeGeneGroups()
    writeUnidentifiable()
    writeChangedName()

    
    '''
    entire_list = []
    entire_list.extend(getListForGroup("A"))
    entire_list.extend(getListForGroup("B"))
    entire_list.extend(getListForGroup("C"))
    entire_list.extend(getListForGroup("D"))
    
    download_svg(entire_list)
    '''
    
    
    modify_base_svg()
    
    
    b_list = getListForGroup("B")
    for b_gene in b_list:
        create_svg(b_gene)
    
    
    
if __name__ == "__main__":
    main()
