#!/usr/bin/env python2
"""
Get data from BRENDA

Copyright (C) 2017  Martin Engqvist Lab
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
from os.path import join, dirname, basename, exists, isdir

### Load environmental variables from the project root directory ###
# find .env automagically by walking up directories until it's found

path = '.'
while '.env' not in os.listdir(path):
	path = '../' + path

dotenv_path = os.path.abspath(path) + '/'



# set up directory paths
CURRENT_DIR = os.getcwd()
PROJ = dirname(dotenv_path) # project root directory

DATA = join(PROJ, 'data') #data directory
RAW_EXTERNAL = join(DATA, 'raw_external') # external data raw directory
RAW_INTERNAL = join(DATA, 'raw_internal') # internal data raw directory
INTERMEDIATE = join(DATA, 'intermediate') # intermediate data directory
FINAL = join(DATA, 'final') # final data directory

RESULTS = join(PROJ, 'results') # output directory
FIGURES = join(RESULTS, 'figures') # figure output directory
PICTURES = join(RESULTS, 'pictures') # picture output directory


# make folders specific for certain data
folder_name = 'brenda_data_2019_1'
if folder_name != '':
	#make folders if they don't exist
	if not exists(join(RAW_EXTERNAL, folder_name)):
		os.makedirs(join(RAW_EXTERNAL, folder_name))

	if not exists(join(INTERMEDIATE, folder_name)):
		os.makedirs(join(INTERMEDIATE, folder_name))

	if not exists(join(FINAL, folder_name)):
		os.makedirs(join(FINAL, folder_name))


#### Your code here ####


# from SOAPpy import SOAPProxy ## for usage without WSDL file
import time
import string
import hashlib
import re
import json
#
# # BRENDA_USERNAME and BRENDA_PASS needs to be registered in the BRENDA website and then added to the .env file of this repository.
# parameters = BRENDA_USERNAME + ',' + hashlib.sha256(BRENDA_PASS).hexdigest()
#
#
# endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
# client = SOAPProxy(endpointURL)


def get_substrates(filepath):
	'''
	Download data regarding the substrates for each EC
	'''
	# Get the names of all organisms
	resultString = client.getEcNumbersFromSubstrate(parameters)

	ec_list = resultString.split('!')
	print('There are %s EC with substrate data' % len(ec_list))

	data = {}

	counter = 0
	for ec in ec_list:
		counter +=1
		if counter % 250 == 0:
			print(counter)

		data[ec] = client.getSubstrate(parameters + ',' + "ecNumber*%s" % ec)

		time.sleep(1)

	# save to file
	print('saving', filepath)
	with open(filepath, 'w') as outfile:
		json.dump(data, outfile)



def get_natural_substrates(filepath):
	'''
	Download data regarding the natural substrates for each EC
	'''
	# Get the names of all organisms
	resultString = client.getEcNumbersFromNaturalSubstrate(parameters)

	ec_list = resultString.split('!')
	print('There are %s EC with natural substrate data' % len(ec_list))

	data = {}

	counter = 0
	for ec in ec_list:
		counter +=1
		if counter % 250 == 0:
			print(counter)

		data[ec] = client.getNaturalSubstrate(parameters + ',' + "ecNumber*%s" % ec)

		time.sleep(1)

	# save to file
	with open(filepath, 'w') as outfile:
		json.dump(data, outfile)




def get_temperature_optima(filepath):

	# Get the names of all organisms with temperature data associated with them.
	resultString = client.getOrganismsFromTemperatureOptimum(parameters)

	org_list = resultString.split('!')
	print('There are %s organisms with enzyme temperature data' % len(org_list))


	#now get the temps for those organisms
	data = {}

	counter = 0
	for organism in org_list:
		counter +=1
		if counter % 250 == 0:
			print(counter)
		elif organism == 'Shiraia sp.':
			continue

		#get data for organism
		org_optima_result = client.getTemperatureOptimum(parameters + ',' + "organism*%s" % organism) #return data containing all the
		#print(org_optima_result)

        #process the data if there is an actual return result
		if org_optima_result != '':
			organism = organism.lower().replace(' ', '_').replace('.', '')
			organism = '_'.join(organism.split('_')[:2])

			if data.get(organism) is None: #this organism was not there before
				data[organism] = org_optima_result

			else: #organism is present, merge data
				data[organism] = data[organism] + '!' + org_optima_result

		time.sleep(1)

	# save to file
	with open(filepath, 'w') as outfile:
		json.dump(data, outfile)


def parse_temp(filepath):
	'''
	Parses the temperature records retrieved from each organism.
	The temperatures are stored on a per EC number basis.
	'''

	# get data
	with open(filepath, 'r') as infile:
		org_list_result = json.reads(infile.read())

	#process the output to get temperatures
	data={}
	for organism in org_list_result.keys():

		if len(organism.split('_')) != 2: #get rid of stuff that does not have two names
			continue

		elif organism.startswith('unidentified') or organism.startswith('uncultured') or organism.startswith('synthetic') or organism.startswith('unknown'):
			continue

		else:
			ec_nums = {} #for temporarily storing the EC numbers and their temp data

			#get record for a given organism
			record = org_list_result[organism]

			#split up the record to prepare it for temperature extraction
			org_optima_list = record.split('!')

			#for each of the entries for this organism
			for item in org_optima_list:

				#find ec number and add it to dictionary if it is not present
				m = re.search('ecNumber[*]([0-9]+.[0-9]+.[0-9]+.[B0-9]+)#', item)
				if m is None:
					print(item)
				else:
					ec = m.group(1)
					if ec_nums.get(ec) is None:
						ec_nums[ec] = {}
						ec_nums[ec]['temps'] = []

					#now find the temperature and add that to the ec number
					m = re.search('(temperatureOptimum[*])([-0-9.]+)#', item)
					if m is None:
						print(item)
					else:
						if m.group(2) != '-999': #-999 seems to be some sort of N/A value
							ec_nums[ec]['temps'].append(float(m.group(2)))


			#now go through the ec numbers for this organism and make the average
			for ec in ec_nums.keys():
				if len(ec_nums[ec]['temps']) != 0:
					ec_nums[ec]['average'] = int(round(sum(ec_nums[ec]['temps'])/float(len(ec_nums[ec]['temps']))))
					data[organism] = ec_nums


	return data




def parse_substrates(filepath):
    """
    Process the BRENDA data and retain the first molecule named as the main substrate.
    """

    with open(filepath, 'r') as f:
        in_data = json.loads(f.read())

    out_data = {}
    for ec in sorted(in_data.keys()):
        out_data[ec] = {'first_substrate':set([]), 'other_substrates':set([]), 'first_product':set([]), 'other_products':set([])}

        for entry in in_data[ec].split('!'):
            for item in entry.split('#'):
                #print(item)

                if item.startswith(u'reactionPartners'):
                    reaction = item.replace(u'reactionPartners*', u'')

                    # first attempt to get all substrates and all products separated
                    try:
                        reactants, products = reaction.split('=')
                    except:
                        reactants = reaction
                        products = ''

                    # split the reactants
                    try:
                        reactant_list = reactants.split(' + ')
                    except:
                        reactant_list = reactants

                    # take the first one
                    out_data[ec]['first_substrate'].add(re.sub('^[0-9]+?[ ]', '', reactant_list[0].strip()))

                    # get the others
                    out_data[ec]['other_substrates'].update([re.sub('^[0-9]+?[ ]', '', i.strip()) for i in reactant_list[1:]])


                    # split the products
                    try:
                        product_list = products.split(' + ')
                    except:
                        product_list = products

                    # take the first one
                    out_data[ec]['first_product'].add(re.sub('^[0-9]+?[ ]', '', product_list[0].strip()))

                    # get the others
                    out_data[ec]['other_products'].update([re.sub('^[0-9]+?[ ]', '', i.strip()) for i in product_list[1:]])


    # get rid of "more" as a substrate
    for ec in out_data.keys():
        for key in out_data[ec].keys():
            out_data[ec][key] = sorted(list(out_data[ec][key] - set(['more'])))


    return out_data


def parse_natural_substrates(filepath):
    """
    Process the BRENDA data and retain the first molecule named as the main substrate.
    """

    with open(filepath, 'r') as f:
        in_data = json.loads(f.read())

    out_data = {}
    for ec in sorted(in_data.keys()):
        out_data[ec] = {'first_substrate':set([]), 'other_substrates':set([]), 'first_product':set([]), 'other_products':set([])}

        for entry in in_data[ec].split('!'):
            for item in entry.split('#'):
                #print(item)

                if item.startswith(u'naturalReactionPartners'):
                    reaction = item.replace(u'naturalReactionPartners*', u'')

                    # first attempt to get all substrates and all products separated
                    try:
                        reactants, products = reaction.split('=')
                    except:
                        reactants = reaction
                        products = ''

                    # split the reactants
                    try:
                        reactant_list = reactants.split(' + ')
                    except:
                        reactant_list = reactants

                    # take the first one
                    out_data[ec]['first_substrate'].add(re.sub('^[0-9]+?[ ]', '', reactant_list[0].strip()))

                    # get the others
                    out_data[ec]['other_substrates'].update([re.sub('^[0-9]+?[ ]', '', i.strip()) for i in reactant_list[1:]])


                    # split the products
                    try:
                        product_list = products.split(' + ')
                    except:
                        product_list = products

                    # take the first one
                    out_data[ec]['first_product'].add(re.sub('^[0-9]+?[ ]', '', product_list[0].strip()))

                    # get the others
                    out_data[ec]['other_products'].update([re.sub('^[0-9]+?[ ]', '', i.strip()) for i in product_list[1:]])


    # get rid of "more" as a substrate
    for ec in out_data.keys():
        for key in out_data[ec].keys():
            out_data[ec][key] = sorted(list(out_data[ec][key] - set(['more'])))


    return out_data
#
# get_substrates(filepath=join(RAW_EXTERNAL, folder_name, 'substrates_data.txt'))
#
#
# get_natural_substrates(filepath=join(RAW_EXTERNAL, folder_name, 'natural_substrates_data.txt'))


subst = parse_substrates(join(RAW_EXTERNAL, folder_name, 'substrates_data.txt'))
with open(join(INTERMEDIATE, folder_name, 'substrates_filtered.json'), 'w') as f:
	f.write(json.dumps(subst))


nat_subst = parse_natural_substrates(join(RAW_EXTERNAL, folder_name, 'natural_substrates_data.txt'))
with open(join(INTERMEDIATE, folder_name, 'natural_substrates_filtered.json'), 'w') as f:
	f.write(json.dumps(nat_subst))

#get_temperature_optima('temp_optima.json')
#parse_temp('temp_optima.json')
