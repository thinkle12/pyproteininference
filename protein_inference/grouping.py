#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:51:50 2017

@author: hinklet
"""

import collections
###Should I just import ProteinInference instead?
import datastore
from physical import ProteinGroup

class Grouper(object):
    """
    Parent Grouper class for all grouper subset classes
    """
    
    def __init__(self):
        None
        
class SimpleSubsetting(Grouper):
    """
    Class that performs simple subset grouping of Protein objects from protein_inference.physical.Protein()

    If we have two proteins A and B, and protein B's peptides are a COMPLETE subset of protein A's peptides.
    Then we create a protein group between the two proteins with protein A leading the group.

    Example: protein_inference.grouping.SimpleSubsetting(data_class = data)

    Where data is a DataStore object

    """
    def __init__(self,data_class):
        #scored_data is scored data from PI_picker class or from PI_scoring class
        if data_class.picked_proteins_scored:
            self.scored_data = data_class.picked_proteins_scored
        else:
            self.scored_data = data_class.scored_proteins
        self.data_class = data_class

        
    def execute(self):
        
        #Get the protein to peptide dictionary from datastore method
        ptp = datastore.ProteinToPeptideDictionary(self.data_class)
        ptp.execute()
        #Stores the dictionary as protein_peptide_dictionary in data_class
        prot_pep_dict = self.data_class.protein_peptide_dictionary
        
        #Get the higher or lower variable
        #try to see if high_low_better is already a variable...
        if not self.data_class.high_low_better:
            #If not run HigherOrLower method
            hl = datastore.HigherOrLower(self.data_class)
            hl.execute()
            higher_or_lower = self.data_class.high_low_better
        else:
            higher_or_lower = self.data_class.high_low_better

        
        #Generate the number of peptides for all protein objects (This can probably be moved to when we create the pre score data
        for protein_objects in self.scored_data:
            cur_identifier = protein_objects.identifier
            protein_objects.num_peptides = len(prot_pep_dict[cur_identifier])
            
        #Sort by the number of peptides for each protein with the most peptides at the top
        sprots = sorted(self.scored_data, key = lambda k: k.num_peptides, reverse=True)
        #Create two protein lists for grouping
        prot_obj_list = [x for x in sprots]
        prot_obj_list2 = [x for x in sprots]
        prot_name_list = [x.identifier for x in prot_obj_list2]

        #Reverse the list to be lowest number of peptides to highest number of peptides
        prot_obj_list.reverse()
        prot_obj_list2.reverse()
        prot_name_list.reverse()

        #We do two backwards for-loops so we can delete from the list as we go
        #Group from Top to bottom - doing bottom to top backwards (as we reverse the sorting above)
        #By doing this, subsets can only be below you in the list
        #Once something has been grouped, remove it from the list... no need to try it, its already in a group
        grouped = []
        #Loop backwards...
        for i in range(len(prot_obj_list)-1, -1, -1):
            #Here check to see if the current protein has been grouped already... if it has pass
            if prot_obj_list[i].identifier in prot_name_list:
                #Not grouped, so create a group for it
                grouped.append([prot_obj_list[i]])
                #Find its current peptides from the peptide dictionary
                cur_peps = prot_pep_dict[prot_obj_list[i].identifier]
                #Start the interior reverse loop
                for j in range(len(prot_obj_list2)-1, -1, -1):
                    #Make sure the protein we test is not identical to the current lead in the protein group
                    if prot_obj_list2[j].identifier!=prot_obj_list[i].identifier:
                        #If its not get the peptides from the inner protein
                        test_peps = prot_pep_dict[prot_obj_list2[j].identifier]
                        #Test to see if the inner proteins peptides are a subset of the outer proteins peptides
                        if test_peps.issubset(cur_peps):
                            #If it is, append the inner protein to the current group
                            grouped[-1].append(prot_obj_list2[j])
                            #Delete the inner protein from its list as it cannot be in another group, no reason to loop through it
                            del prot_name_list[j]
                            del prot_obj_list2[j]
            else:
                pass
            print len(prot_name_list)


        #Now get the list of scored proteins
        scored_proteins = list(self.scored_data)
        #Also get a list of the protein identifiers of all scored proteins
        #This is neccesary as protein picker deletes a target or decoy depending on which scored worse
        #So some proteins in the group no longer exist in 'scored_proteins'
        protein_finder = [x.identifier for x in scored_proteins]

        scores_grouped = []
        group_id = 0
        for groups in grouped:
            sub_groups = []
            group_id = group_id+1
            for prots in groups:
                try:
                    #For the proteins in this group see if it exists in the scored proteins
                    #If it does exist put its score and infoprmation in this sub group
                    pindex = protein_finder.index(prots.identifier)
                    cur_protein = scored_proteins[pindex]
                    cur_protein.group_id.append(group_id)
                    sub_groups.append(cur_protein)
                except ValueError:
                    #If the protein does no longer exist, pass (if it got removed by protein picker)
                    pass

            #Sort the sub groups by score....
            if len(sub_groups)>0:
                if higher_or_lower=='lower':
                    sub_groups = sorted(sub_groups, key = lambda k: k.score, reverse=False)
                if higher_or_lower=='higher':
                    sub_groups = sorted(sub_groups, key = lambda k: k.score, reverse=True)

                scores_grouped.append(sub_groups)

        #Sort all groups by the lead protein score of each group...
        if higher_or_lower=='lower':
            scores_grouped = sorted(scores_grouped, key = lambda k: k[0].score, reverse=False)
        if higher_or_lower=='higher':
            scores_grouped = sorted(scores_grouped, key = lambda k: k[0].score, reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped



class GlpkSetup(Grouper):
    """
    This class is used to setup the glpk file for analysis.

    Example: protein_inference.grouping.GlpkSetup(data_class = data ,glpkin_filename='glpkin_example.mod'))

    Class will use attributes from DataStore object data.
    Class will also write the glpkin filename to be used in protein_inference.grouping.GlpkRunner()

    The Bulk of the glpk input file looks as follows:
    s.t. c1: y[5658] >=1;
    s.t. c2: y[14145]+y[4857]+y[4858]+y[10143]+y[2966] >=1;
    s.t. c3: y[320]+y[4893]+y[4209]+y[911]+y[2767]+y[2296]+y[10678]+y[3545] >=1
    """
    def __init__(self,data_class,glpkin_filename='glpkin.mod'):
        self.data_class = data_class
        self.glpkin_filename = glpkin_filename

    def execute(self):

        #Here we get the peptide to protein dictionary
        peptoprot = datastore.PeptideToProteinDictionary(self.data_class)
        peptoprot.execute()
        pep_prot_dict = self.data_class.peptide_protein_dictionary


        #Here we get the list of all proteins
        plist = []
        for peps in pep_prot_dict.keys():
            for prots in list(pep_prot_dict[peps]):
                plist.append(prots)

        #Here we get the unique proteins
        unique_prots = list(set(plist).union())

        #Setup default dictionaries
        dd_num = collections.defaultdict(list)
        dd_prot_nums = collections.defaultdict(list)

        #For all the unique proteins from the search create a number to protein dictionary and a protein to number dictionary
        #Here we essentially assign a number to each protein
        #This is important as the glpk analysis just treats proteins as numbers...
        for p in range(len(unique_prots)):
            dd_num[unique_prots[p]].append(p)
            dd_prot_nums[p].append(unique_prots[p])

        #Store this data as glpk_protein_number_dictionary and glpk_number_protein_dictionary
        #The numbers are important as they are used in the GLPK input and we need to know what number in the GLPK output corresponds with which protein from the search
        self.data_class.glpk_protein_number_dictionary = dd_num
        self.data_class.glpk_number_protein_dictionary = dd_prot_nums

        #Create the GLPK input file
        fileout = open(self.glpkin_filename,'w')

        #Not sure if this header string is correct or if it needs to be here...
        fileout.write('/* sets */'+'\n'+'set PROTEINS;'+'\n'+'\n'+'\n')
        fileout.write('/* decision variables: yi, i in {1,..,5}. yi = 1 -> protein i is selected */'+'\n')
        fileout.write('var y {i in PROTEINS} binary >=0;'+'\n')
        fileout.write('/* objective function */'+'\n')
        fileout.write('minimize z: sum{i in PROTEINS} y[i];'+'\n'+'\n')
        fileout.write('/* Constraints */'+'\n')

        #Here we create the bulk of the input file which needs to look as follows:
        #s.t. c1: y[5658] >=1;
        #s.t. c2: y[14145]+y[4857]+y[4858]+y[10143]+y[2966] >=1;
        #s.t. c3: y[320]+y[4893]+y[4209]+y[911]+y[2767]+y[2296]+y[10678]+y[3545] >=1;
        #Each of the lines (constants, c1,c2,c3) is a peptide and each of the y[x] is a protein
        tot_peps = pep_prot_dict.keys()
        for j in range(len(pep_prot_dict)):
            combine = ['y['+str(dd_num[x][0])+']' for x in list(pep_prot_dict[tot_peps[j]])]
            fileout.write('s.t. c'+str(j+1)+': '+'+'.join(combine)+' >=1;'+'\n')

        #Finish writing the rest of the file and close it
        fileout.write('\n')
        fileout.write('data;'+'\n')
        numlist = [str(dd_num[x][0]) for x in unique_prots]
        strlist = ' '.join(numlist)
        #End the file with listing the entire set of proteins... (as its number identifier)
        fileout.write('set PROTEINS := '+strlist+' ;'+'\n'+'\n')

        fileout.write('end;')
        fileout.close()

class GlpkRunner(Grouper):
    """
    The GlpkRunner class takes a path to glpsol, the glpk input file from protein_inference.grouping.GlpkSetup(), a glpkout filename as well as a file_override option

    Example: protein_inference.grouping.GlpkRunner(path_to_glpsol = '/glpsol',glpkin='glpkin_example.mod',glpkout='glpkout_example.sol',file_override=False)

    path to glpsol on rescomp3 is: '/gne/research/apps/protchem/glpk/bin/glpsol'

    Typically set file_override to false unless you know what you are doing (IE you have a specific glpk solution file you want to use)

    Important output of this class is the glpk output solution file to be used in protein_inference.grouping.GlpkGrouper
    """

    def __init__(self,path_to_glpsol = '/glpsol',glpkin='glpkin.mod',glpkout='glpkout.sol',file_override=False):
        self.path_to_glpsol = path_to_glpsol
        self.file_override = file_override
        self.glpkin = glpkin
        self.glpkout = glpkout

    def execute(self):
        #If there is no file_override (mainly for offline testing)
        if not self.file_override:
            import subprocess
            #Run GLPK with the following command
            p = subprocess.Popen(str(self.path_to_glpsol)+' -m '+str(self.glpkin)+' -o '+str(self.glpkout), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

            output = p.communicate()

            print 'Start Command line Stdout'
            print output[0]
            print 'End Command line Stdout'
            print 'Start Command line Stderr'
            print output[1]
            print 'End Command line Stderr'

            if output[0]=='':
                raise ValueError('Glpk did not produce any output... See potential error output above')
        else:
            pass



class GlpkGrouper(Grouper):
    """
    This class takes a digest class object, a glpk solution file, as well as an option for swissprot override (protein naming convention).

    Example: protein_inference.grouping.GlpkGrouper(data_class = data,digest_class = digest,swissprot_override="soft",glpksolution_filename='glpkout_example.sol')

    Where data is a DataStore object, digest is a Digest.insilicodigest.InSilicoDigest() class object and glpksolution_filename is a glpk solution file.

    Finally, swissprot_override is a lead protein override for naming convention which can have 3 options: False, "soft", or "hard".
    selecting False skips the override.

    selecting "soft" will cycle through all protein leads and groups. If a protein lead is an unreviewed hit (trembl) then all proteins in the lead proteins group are inspected.
    If any of these intra group proteins are a reviewed hit (swissprot) and the lead trembl proteins peptides are a complete subset of the swissprot proteins peptides, then we select
    to swap the trembl and swissprot hits so that the swissprot is now the new lead protein.

    We opt to use this soft override scheme as default because given the redundancy of our databases. Out of GLPK we see a lot of unreviewed hits being called lead proteins.
    If proteins share the same peptides or share a very similar set of peptides we are unsure how GLPK selects which protein to supply as the lead.
    As such, we get many reviewed and unreviewed proteins as lead proteins.
    Performing this "soft" override switches many of these lead trembl hits to reviewed swissprot hits.
    This is important as group members are used to seeing swissprot identifiers as opposed to trembl identifiers

    selecting "hard" will perform the same swapping as the "soft" override. However,
    the "hard" override will swap a trembl lead with the swissprot protein with the highest number of peptides in its group (so long as its already not a lead protein itself)
    even if the unreviewed has a unique peptide relative to the peptides of all proteins in its group.
    This setting is not recommended given that you can potentially lose out on Unique peptides

    """
    ###Define indicies better and make this code more readable...
    ###Define indicies in init and show commented examples of how the data looks...

    def __init__(self,data_class,digest_class,swissprot_override="soft",glpksolution_filename='glpkout.sol'):
        self.data_class = data_class
        self.digest_class = digest_class
        self.swissprot_override = swissprot_override
        if data_class.picked_proteins_scored:
            self.scored_data = data_class.picked_proteins_scored
        else:
            self.scored_data = data_class.scored_proteins
        self.data_class = data_class
        self.glpksolution_filename = glpksolution_filename
        self.lead_protein_set = None

    def execute(self):
        print('this is working')
        import os
        glpk_out = open(self.glpksolution_filename,'r')

        #Get the number protein dictionary from glpk_setup
        dd_prot_nums = self.data_class.glpk_number_protein_dictionary

        glpk_out = glpk_out.read()
        glpk_out = glpk_out.split('\n')


        #Cant find a better way to do this... there are modules out there that work with glpk...
        start = glpk_out.index('   No. Column name       Activity     Lower bound   Upper bound')


        newlist = []
        #Fix this -13 and +2... not really sure how
        #Based on the output file we should start two lines after the start index and stop 13 before the end of the file...
        for lines in range(start+2,len(glpk_out)-13):
            new = [x.strip() for x in glpk_out[lines].split(' ')]
            res = []
            for stuff in new:
                if stuff!='':
                    res.append(stuff)
            newlist.append(res)

        #Use 1 here as 1 is the location in the line of the protein number
        #3 is the location of the binary (which indicates whether or not the protein number has a unique peptide making it a lead protein)
        numbers = [x[1].split(']')[0].split('[')[-1] for x in newlist]
        binary = [x[3] for x in newlist]

        self.numbers = numbers

        #Here we extract out the lead proteins...
        lead_proteins = []
        for k in range(len(numbers)):
            if binary[k]=='1':
                try:
                    passing_protein_number = int(numbers[k])
                    lead_proteins.append(dd_prot_nums[passing_protein_number][0])
                except IndexError:
                    print("No Protein for Protein Number"+str(passing_protein_number))

        lead_protein_set = set(lead_proteins)
        self.lead_protein_set = lead_protein_set

        print 'Number of lead proteins = '+str(len(lead_proteins))



        scored_proteins = list(self.scored_data)
        protein_finder = [x.identifier for x in scored_proteins]


        lead_protein_objects = []
        lead_protein_identifiers = []
        for proteins in lead_proteins:
            if proteins in protein_finder:
                p_ind = protein_finder.index(proteins)
                protein_object = scored_proteins[p_ind]
                lead_protein_objects.append(protein_object)
                lead_protein_identifiers.append(protein_object.identifier)

        self.lead_protein_objects = lead_protein_objects
        #Delete input and output files...
        #os.remove('glpkin.mod')
        #os.remove(self.glpksolution_filename)

        #Now we have the lead Proteins so we need to get the peptides for each lead protein
        #Then we find all proteins that share at least 1 peptide with each lead protein
        #If they share at least 1 peptide then assign that protein to the group...
        self.data_class.glpk_lead_proteins = lead_protein_objects

        prottopep = datastore.ProteinToPeptideDictionary(self.data_class)
        prottopep.execute()
        prot_pep_dict = self.data_class.protein_peptide_dictionary

        #Get the higher or lower variable
        if not self.data_class.high_low_better:
            hl = datastore.HigherOrLower(self.data_class)
            hl.execute()
            higher_or_lower = self.data_class.high_low_better
        else:
            higher_or_lower = self.data_class.high_low_better




        #####The Nnum_unique_peptide thing isnt working in this code block below I think...
        #Its saying that a bunch of protein groups dont have unique peptides which doesnt make sense...
        #Need to go over all proteins and all peptides in the protein and check if the peptide is unique
        #do this by taking in_silico_peptides_to_proteins (entire mapping of peptides to proteins)
        #And seeing if the length of the dictionary value (which are proteins) is equal to 1...
        #If its equal to 1 then that peptide (key) has only 1 protein it maps to potentially...

        list_of_prots_not_in_db = []
        list_of_peps_not_in_db = []
        print 'Grouping Proteins...'
        in_silico_peptides_to_proteins = self.digest_class.peptide_to_protein_dictionary
        grouped = []
        for protein_objects in lead_protein_objects:
            protein_objects.peptides = [x for x in list(prot_pep_dict[protein_objects.identifier]) if x in self.data_class.restricted_peptides]
            sub_group = [protein_objects]
            cur_peptides = list(prot_pep_dict[protein_objects.identifier])
            # Use a default dict set here for automatic union
            all_potential_proteins = collections.defaultdict(set)
            all_potential_proteins['current'] = set([])
            for peptides in cur_peptides:  #Probably put an if here... if peptides is in the list of peptides after being restricted by datastore.RestrictMainData
                if peptides in self.data_class.restricted_peptides:
                    # Get the proteins that map to the current peptide using in_silico_peptides_to_proteins
                    # First make sure our peptide is formatted properly...
                    if not peptides.isupper() or not peptides.isalpha():
                        # If the peptide is not all upper case or if its not all alphabetical...
                        rm = datastore.RemoveMods(peptides)
                        # Then run remove mods on it...
                        rm.execute()
                        # Redefine the peptide as the stripped version below...
                        peptides = rm.stripped_peptide
                    potential_protein_list = list(in_silico_peptides_to_proteins[peptides])
                    if len(potential_protein_list) == 0:
                        list_of_prots_not_in_db.append(peptides)
                        list_of_peps_not_in_db.append(protein_objects.identifier)
                        print 'Protein '+str(protein_objects.identifier)+' and Peptide '+str(peptides)+' is not in database...'
                    #If these proteins are not a lead add them to all potential protein default dict
                    for prots in potential_protein_list:
                        if prots not in lead_protein_set:
                            try:
                                #Try to find its object using protein_finder (list of identifiers) and scored_proteins (list of Protein Objects)
                                cur_index = protein_finder.index(prots)
                                current_protein_object = scored_proteins[cur_index]
                                other_peptides = list(prot_pep_dict[current_protein_object.identifier])
                                # Assign the peptides only if they are in restricted_peptides....
                                current_protein_object.peptides = [x for x in other_peptides if x in self.data_class.restricted_peptides]
                                all_potential_proteins['current'].add(current_protein_object)
                            except ValueError:
                                # Put a print statement here showing which proteins werent found... Meaning they dont have a protein object
                                # I think this could be due to the protein having only 1 peptide and it gets filtered out...
                                # Check to see though...
                                pass
            #Next append them to the current sub group... (Which is just the lead protein object)
            sub_group = sub_group+list(all_potential_proteins['current'])
            # sub_group at first is just the lead protein object...
            # We then try apply grouping by looking at all peptides from the lead...
            # For all of these peptides look at all other non lead proteins and try to assign them to the group...
            # We assign the entire protein object as well... in the above try/except
            #Then append this sub group to the main list
            #The variable grouped is now a list of lists which each element being a Protein object and each list of protein objects corresponding to a group
            grouped.append(sub_group)

        self.list_of_prots_not_in_db = list_of_prots_not_in_db
        self.list_of_peps_not_in_db = list_of_peps_not_in_db

        #Here are the proteins from the database that are reviewed (swissprot)
        sp_protein_set = set(self.digest_class.swiss_prot_protein_dictionary['swiss-prot'])

        #Here dont do swissprot override...
        if not self.swissprot_override:
            print 'Applying Group IDs...'
            # Here we create group ID's for all groups and do some sorting
            scores_grouped = []
            list_of_group_objects = []
            group_id = 0
            for groups in grouped:
                sub_groups = []
                group_id = group_id + 1
                #Create a ProteinGroup object...
                pg = ProteinGroup(group_id)
                print str(group_id)
                for prots in groups:
                    try:
                        # The following loop assigns group_id's, reviewed/unreviewed status, and number of peptides...
                        pindex = protein_finder.index(prots.identifier)
                        cur_protein = scored_proteins[pindex]
                        if group_id not in cur_protein.group_identification:
                            cur_protein.group_identification.append(group_id)
                        if prots.identifier in sp_protein_set:
                            cur_protein.reviewed = True
                        else:
                            cur_protein.unreviewed = True
                        cur_identifier = prots.identifier
                        cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                        # Here append the number of peptides... so we can use this as secondary sorting...
                        sub_groups.append(cur_protein)
                    except ValueError:
                        pass
                #Assign the protein objects to their group object...
                pg.proteins = sub_groups
                scores_grouped.append(sub_groups)
                list_of_group_objects.append(pg)


        if self.swissprot_override=='hard':

            #Hard swissprot_override will exchange any lead unreviewed hit with the best reviewed hit in its group... even if the unreviewed has a unique peptide relative to the peptides of all proteins in its group...
            print 'Applying Group IDs... and Executing Hard Swissprot Override...'
            #Here we create group ID's for all groups and do some sorting
            scores_grouped = []
            group_id = 0
            leads = set()
            lead_replaced_prot_pairs = []
            list_of_group_objects = []
            for groups in grouped:
                sub_groups = []
                group_id = group_id + 1
                pg = ProteinGroup(group_id)
                print str(group_id)
                for prots in groups:
                    try:
                        # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
                        pindex = protein_finder.index(prots.identifier)
                        cur_protein = scored_proteins[pindex]
                        if group_id not in cur_protein.group_identification:
                            cur_protein.group_identification.append(group_id)
                        if prots.identifier in sp_protein_set:
                            cur_protein.reviewed = True
                        else:
                            cur_protein.unreviewed = True
                        cur_identifier = prots.identifier
                        cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                        #Here append the number of unique peptides... so we can use this as secondary sorting...
                        sub_groups.append(cur_protein)
                        #Sorted groups then becomes a list of lists...


                    except ValueError:
                        #Here we pass if the protein does not have a score...
                        #Potentially it got 'picked' (removed) by protein picker...
                        pass


                #Sort the groups based on higher or lower indication, secondarily sort the groups based on number of unique peptides
                #We use the index [1:] as we do not wish to sort the lead protein... from GLPK
                if higher_or_lower=='lower':
                    sub_groups[1:] = sorted(sub_groups[1:], key = lambda k: (float(k.score),-float(k.num_peptides)), reverse=False)
                if higher_or_lower=='higher':
                    sub_groups[1:] = sorted(sub_groups[1:], key = lambda k: (float(k.score),float(k.num_peptides)), reverse=True)

                #scores_grouped is the MAIN list of lists with all protein objects
                scores_grouped.append(sub_groups)
                #If the lead is reviewed append it to leads and do nothing else...

                if sub_groups[0].reviewed:
                    if '-' in sub_groups[0].identifier:
                        pure_id = sub_groups[0].identifier.split('-')[0]
                        # Start to loop through sub_groups which is the current group...
                        for potential_replacement in sub_groups[1:]:
                            isoform_override = potential_replacement
                            if isoform_override.identifier==pure_id and isoform_override.identifier not in leads and set(sub_groups[0].peptides).issubset(set(isoform_override.peptides)):
                                isoform_override_index = scores_grouped[-1].index(isoform_override)
                                cur_iso_lead = scores_grouped[-1][0]
                                print cur_iso_lead.identifier
                                scores_grouped[-1][0], scores_grouped[-1][isoform_override_index] = scores_grouped[-1][isoform_override_index], scores_grouped[-1][0]
                                scores_grouped[-1][isoform_override_index], scores_grouped[-1][0]
                                new_iso_lead = scores_grouped[-1][0]
                                print new_iso_lead.identifier
                                lead_replaced_prot_pairs.append([cur_iso_lead, new_iso_lead])
                                leads.add(sub_groups[0].identifier)
                #If the lead is unreviewed then try to replace it with the best reviewed hit
                if not sub_groups[0].reviewed:
                    #If the lead is unreviewed attempt to replace it...
                    #Start to loop through sorted_groups which is the current sub group... (sub_groups)
                    for hits in sub_groups[1:]:
                        #Find the first reviewed hit... if its not a lead protein already then swap positions in scores_grouped and break...
                        if hits.reviewed:
                            best_swiss_prot_prot = hits
                            if best_swiss_prot_prot.identifier not in leads:
                                #We use -1 as the index of scores_grouped because the current 'sub_groups' is the last entry appended to scores_grouped
                                #Essentially scores_grouped[-1]==sub_groups
                                #We need this syntax so we can switch the location of the unreviewed lead identifier with the best reviewed identifier in scores_grouped
                                swiss_prot_override_index = scores_grouped[-1].index(best_swiss_prot_prot)
                                cur_tr_lead = scores_grouped[-1][0]
                                print cur_tr_lead.identifier
                                scores_grouped[-1][0], scores_grouped[-1][swiss_prot_override_index] = scores_grouped[-1][swiss_prot_override_index], scores_grouped[-1][0]
                                new_sp_lead =  scores_grouped[-1][0]
                                print new_sp_lead.identifier
                                lead_replaced_prot_pairs.append([cur_tr_lead,new_sp_lead])
                                #Append new_sp_lead protein to leads, to make sure we dont repeat leads
                                leads.add(new_sp_lead.identifier)
                                break
                            else:
                                #If no reviewed and none not in leads then pass...
                                pass
                        else:
                            pass

                pg.proteins = sub_groups
                list_of_group_objects.append(pg)
            self.data_class.lead_replaced_proteins = lead_replaced_prot_pairs



        if self.swissprot_override=='soft':
            print 'Applying Group IDs... and Executing Soft Swissprot Override...'
            # Here we create group ID's for all groups and do some sorting
            scores_grouped = []
            group_id = 0
            leads = set()
            lead_replaced_prot_pairs = []
            list_of_group_objects = []
            for groups in grouped:
                sub_groups = []
                group_id = group_id + 1
                pg = ProteinGroup(group_id)
                print str(group_id)
                for prots in groups:
                    try:
                        # The following loop assigns group_id's, reviewed/unreviewed status, and number of unique peptides...
                        pindex = protein_finder.index(prots.identifier)
                        cur_protein = scored_proteins[pindex]
                        if group_id not in cur_protein.group_identification:
                            cur_protein.group_identification.append(group_id)
                        if prots.identifier in sp_protein_set:
                            cur_protein.reviewed = True
                        else:
                            cur_protein.unreviewed = True
                        cur_identifier = prots.identifier
                        cur_protein.num_peptides = len(prot_pep_dict[cur_identifier])
                        # Here append the number of unique peptides... so we can use this as secondary sorting...
                        sub_groups.append(cur_protein)
                        # Sorted groups then becomes a list of lists... of protein objects


                    except ValueError:
                        # Here we pass if the protein does not have a score...
                        # Potentially it got 'picked' (removed) by protein picker...
                        pass


                # Sort the groups based on higher or lower indication, secondarily sort the groups based on number of unique peptides
                # We use the index [1:] as we do not wish to sort the lead protein... from GLPK
                if higher_or_lower == 'lower':
                    sub_groups[1:] = sorted(sub_groups[1:],
                                            key=lambda k: (float(k.score), -float(k.num_peptides)),reverse=False)
                if higher_or_lower == 'higher':
                    sub_groups[1:] = sorted(sub_groups[1:],
                                            key=lambda k: (float(k.score), float(k.num_peptides)), reverse=True)

                # scores_grouped is the MAIN list of lists with grouped protein objects
                scores_grouped.append(sub_groups)
                # If the lead is reviewed append it to leads and do nothing else...
                # If the lead is unreviewed then try to replace it with the best reviewed hit
                if not sub_groups[0].reviewed:
                    # If the lead is unreviewed attempt to replace it...
                    # Start to loop through sub_groups which is the current group...
                    for hits in sub_groups[1:]:
                        # Find the first reviewed it... if its not a lead protein already then do score swap and break...
                        if hits.reviewed:
                            best_swiss_prot_prot = hits
                            #If the lead proteins peptides are a subset of the best swissprot.... then swap the proteins... (meaning equal peptides or the swissprot completely covers the tremble reference)
                            if best_swiss_prot_prot.identifier not in leads and set(sub_groups[0].peptides).issubset(set(best_swiss_prot_prot.peptides)):
                                # We use -1 as the idex of scores_grouped because the current 'sub_groups' is the last entry appended to scores grouped
                                # Essentially scores_grouped[-1]==sub_groups
                                # We need this syntax so we can switch the location of the unreviewed lead identifier with the best reviewed identifier in scores_grouped
                                swiss_prot_override_index = scores_grouped[-1].index(best_swiss_prot_prot)
                                cur_tr_lead = scores_grouped[-1][0]
                                print cur_tr_lead.identifier
                                scores_grouped[-1][0], scores_grouped[-1][swiss_prot_override_index] = scores_grouped[-1][swiss_prot_override_index], scores_grouped[-1][0]
                                scores_grouped[-1][swiss_prot_override_index], scores_grouped[-1][0]
                                new_sp_lead = scores_grouped[-1][0]
                                print new_sp_lead.identifier
                                lead_replaced_prot_pairs.append([cur_tr_lead, new_sp_lead])
                                # Append new_sp_lead protein to leads, to make sure we dont repeat leads
                                leads.add(new_sp_lead.identifier)
                                break
                            else:
                                # If no reviewed and none not in leads then pass...
                                pass
                        else:
                            pass
                ###NEW ISOFORM OVERRIDE
                if sub_groups[0].reviewed:
                    if '-' in sub_groups[0].identifier:
                        pure_id = sub_groups[0].identifier.split('-')[0]
                        # Start to loop through sub_groups which is the current group...
                        for potential_replacement in sub_groups[1:]:
                            isoform_override = potential_replacement
                            if isoform_override.identifier==pure_id and isoform_override.identifier not in leads and set(sub_groups[0].peptides).issubset(set(isoform_override.peptides)):
                                isoform_override_index = scores_grouped[-1].index(isoform_override)
                                cur_iso_lead = scores_grouped[-1][0]
                                print cur_iso_lead.identifier
                                scores_grouped[-1][0], scores_grouped[-1][isoform_override_index] = scores_grouped[-1][isoform_override_index], scores_grouped[-1][0]
                                scores_grouped[-1][isoform_override_index], scores_grouped[-1][0]
                                new_iso_lead = scores_grouped[-1][0]
                                print new_iso_lead.identifier
                                lead_replaced_prot_pairs.append([cur_iso_lead, new_iso_lead])
                                leads.add(sub_groups[0].identifier)

                pg.proteins = sub_groups
                list_of_group_objects.append(pg)
            self.data_class.lead_replaced_proteins = lead_replaced_prot_pairs

        print 'Sorting Results based on lead Protein Score'
        if higher_or_lower=='lower':
            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key = lambda k: float(k.proteins[0].score), reverse=False)
        if higher_or_lower=='higher':
            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score), reverse=True)


        # Sometimes we have cases where:
        # protein a maps to peptides 1,2,3
        # protein b maps to peptides 1,2
        # protein c maps to a bunch of peptides and peptide 3
        # Therefore, in the model proteins a and b are equivalent in that they map to 2 peptides together - 1 and 2. peptide 3 maps to a but also to c...
        # Sometimes the model (glpk) will spit out protein b as the lead... we wish to swap protein b as the lead with protein a because it will likely have a better score...
        print('Potentially Reassigning leads...')
        lead_protein_set = set([x.proteins[0].identifier for x in list_of_group_objects])
        for i in range(len(list_of_group_objects)):
            for j in range(1, len(list_of_group_objects[i].proteins)):  # Loop over all sub proteins in the group...
                # if the lead proteins peptides are a subset of one of its proteins in the group, and the secondary protein is not a lead protein and its score is better than the leads... and it has more peptides...
                if higher_or_lower == 'higher':
                    if set(list_of_group_objects[i].proteins[0].peptides).issubset(set(list_of_group_objects[i].proteins[j].peptides)) and list_of_group_objects[i].proteins[j].identifier not in lead_protein_set and list_of_group_objects[i].proteins[0].score < list_of_group_objects[i].proteins[j].score and len(list_of_group_objects[i].proteins[0].peptides) < len(list_of_group_objects[i].proteins[j].peptides):
                        print('protein ' + str(list_of_group_objects[i].proteins[j].identifier) + ' will replace protein ' + str(list_of_group_objects[i].proteins[0].identifier) + ' as lead, with index '+str(j))
                        new_lead = list_of_group_objects[i].proteins[j]
                        old_lead = list_of_group_objects[i].proteins[0]
                        lead_protein_set.add(new_lead.identifier)
                        lead_protein_set.remove(old_lead.identifier)
                        # Swap their positions in the list
                        list_of_group_objects[i].proteins[0], list_of_group_objects[i].proteins[j] = new_lead, old_lead
                        print(j)
                        break


                if higher_or_lower == 'lower':
                    if set(list_of_group_objects[i].proteins[0].peptides).issubset(set(list_of_group_objects[i].proteins[j].peptides)) and list_of_group_objects[i].proteins[j].identifier not in lead_protein_set and list_of_group_objects[i].proteins[0].score > list_of_group_objects[i].proteins[j].score and len(list_of_group_objects[i].proteins[0].peptides) < len(list_of_group_objects[i].proteins[j].peptides):
                        print('protein ' + str(list_of_group_objects[i].proteins[j].identifier) + ' will replace protein ' + str(list_of_group_objects[i].proteins[0].identifier) + ' as lead, with index '+str(j))
                        new_lead = list_of_group_objects[i].proteins[j]
                        old_lead = list_of_group_objects[i].proteins[0]
                        lead_protein_set.add(new_lead.identifier)
                        lead_protein_set.remove(old_lead.identifier)
                        # Swap their positions in the list
                        list_of_group_objects[i].proteins[0], list_of_group_objects[i].proteins[j] = new_lead, old_lead
                        break


        print 'Re Sorting Results based on lead Protein Score'
        if higher_or_lower=='lower':
            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0].score), reverse=False)
            list_of_group_objects = sorted(list_of_group_objects, key = lambda k: float(k.proteins[0].score), reverse=False)
        if higher_or_lower=='higher':
            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0].score), reverse=True)
            list_of_group_objects = sorted(list_of_group_objects, key=lambda k: float(k.proteins[0].score), reverse=True)

        self.data_class.grouped_scored_proteins = scores_grouped
        self.data_class.protein_group_objects = list_of_group_objects


#class multi_subsetting(grouper):
#    ###Here we only do simple subset grouping
#    ###If we have two proteins A and B, and protein B's peptides are a complete subset of protein A's peptides
#    ###Then we create a protein group between the two proteins with protein A leading the group.
#    ###scored_data is sscored data from PI_picker class or from PI_scoring class
#    def __init__(self,data_class):
#        if data_class.picked_proteins_scored:
#            self.scored_data = data_class.picked_proteins_scored
#        else:
#            self.scored_data = data_class.scored_proteins
#
#        self.data_class = data_class
#
#    def execute(self):
#
#        #Get the protein to peptide dictionary from PI_data_store method
#        ptp = PI_data_store.protein_to_peptide_dictionary(self.data_class)
#        ptp.execute()
#        prot_pep_dict = self.data_class.protein_peptide_dictionary
#
#        #Get the higher or lower variable
#        if not self.data_class.high_low_better:
#            hl = PI_data_store.higher_or_lower(self.data_class)
#            hl.execute()
#            higher_or_lower = self.data_class.high_low_better
#        else:
#            higher_or_lower = self.data_class.high_low_better
#
#
#        #Create a list of lists of the number of unique peptides with the associated protein
#        pep_prot_list = []
#        for keys in prot_pep_dict.keys():
#            pep_prot_list.append([len(prot_pep_dict[keys]),keys])
#
#        #Sort by the number of peptides unique for each protein with the most peptides at the top
#        sprots = sorted(pep_prot_list,reverse=True)
#        #Create two protein lists for grouping
#        prot_list = [x[1] for x in sprots]
#        prot_list2 = [x[1] for x in sprots]
#
#        #Reverse the list to be lowest number of peptides to highest number of peptides
#        prot_list.reverse()
#        prot_list2.reverse()
#
#        #We do two backwards for loops so we can delete from the list as we go
#        #Group from Top to bottom - doing bottom to top backwards (as we reverse the sorting above)
#        #By doing this, subsets can only be below you in the list
#        #Once something has been grouped, remove it from the list... no need to try it, its already in a group
#        grouped = []
#        for i in range(len(prot_list)-1, -1, -1):
#            #Here check to see if the current protein has been grouped already... if it has pass
#            if prot_list[i] in prot_list2:
#                #Not grouped, so create a group for it
#                grouped.append([prot_list[i]])
#                #Find its current peptides from the peptide dictionary
#                cur_peps = prot_pep_dict[prot_list[i]]
#                #Start the interior reverse loop
#                for j in range(len(prot_list2)-1, -1, -1):
#                    #Make sure the protein we test is not identical to the current lead in the protein group
#                    if prot_list2[j]!=prot_list[i]:
#                        #If its not get the peptides from the inner protein
#                        test_peps = prot_pep_dict[prot_list2[j]]
#                        #Test to see any of the inner proteins peptides are a subset of the outer proteins peptides
#                        if len(test_peps-cur_peps)!=len(cur_peps):
#                            cur_peps=test_peps-cur_peps
#                            #If it is, append the inner protein to the current group
#                            grouped[-1].append(prot_list2[j])
#                            #Delete the inner protein from its list as it cannot be in another group, no reason to loop through it
#                            del prot_list2[j]
#
#            else:
#                pass
#            print len(prot_list2)
#
#        scored_proteins = list(self.scored_data)
#        protein_finder = [x[0] for x in scored_proteins]
#
#        leads = []
#        scores_grouped = []
#        group_id = 0
#        for groups in grouped:
#            sorted_groups = []
#            group_id = group_id+1
#            for prots in groups:
#                if prots not in leads:
#                    try:
#                        pindex = protein_finder.index(prots)
#                        cur_score = scored_proteins[pindex]
#                        cur_score.append(group_id)
#                        sorted_groups.append(cur_score)
#                    except ValueError:
#                        pass
#
#            if len(sorted_groups)>0:
#                if higher_or_lower=='lower':
#                    sorted_groups = sorted(sorted_groups, key = lambda k: float(k[1]), reverse=False)
#                if higher_or_lower=='higher':
#                    sorted_groups = sorted(sorted_groups, key = lambda k: float(k[1]), reverse=True)
#                leads.append(sorted_groups[0][0])
#                scores_grouped.append(sorted_groups)
#
#        if higher_or_lower=='lower':
#            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0][1]), reverse=False)
#        if higher_or_lower=='higher':
#            scores_grouped = sorted(scores_grouped, key = lambda k: float(k[0][1]), reverse=True)
#
#        print scores_grouped
#
#        self.data_class.grouped_scored_proteins = scores_grouped