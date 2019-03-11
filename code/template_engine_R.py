#!/usr/bin/env python3
"""
template_engine_R.py is used to generate new R scripts with useful boilerplate code

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

import sys

# Checking if argparse is installed
try:
	import argparse
except ImportError:
	sys.stderr.write("[Error] The python module \"argparse\" is not installed\n")
	sys.stderr.write("[--] Would you like to install it now using 'sudo easy_install' [Y/N]? ")
	answer = sys.stdin.readline()
	if answer[0].lower() == "y":
		sys.stderr.write("[--] Running \"sudo easy_install argparse\"\n")
		from subprocess import call
		call(["sudo", "easy_install", "argparse"])
	else:
		sys.exit("[Error] Exiting due to missing dependency \"argparser\"")

parser = argparse.ArgumentParser(prog="template_engine_R.py")
parser.add_argument("name", nargs='?', help="Set the name for the script to create", default="new_script.R")
args = parser.parse_args()


def main():
	output_file = open(args.name, "w")
	output_file.write('''#!/usr/bin/env Rscript
"
Modify this line to briefly discribe the functionality of %s

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
"

#make sure the packages are installed
wants <- c('tidyverse')
has   <- wants %%in%% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has], repos="http://ftp.acc.umu.se/mirror/CRAN/")

#load all the packages
lapply(wants, require, character.only = TRUE)


#find the "code" directory, this should be one folder away from the project root
while (basename(getwd()) != 'code'){
  setwd("../")
  if (getwd() == '/'){
    quit(status = 1)
  }
}

#go one folder up and set as project root directory
setwd("../")
PROJ <- getwd()
setwd(paste0(PROJ, '/data/'))

#set folder paths relative to the project root
CURRENT <- getwd()
DATA <- getwd()
RAW_EXTERNAL <- paste0(DATA, '/raw_external/')
RAW_INTERNAL <- paste0(DATA, '/raw_internal/')
INTERMEDIATE = paste0(DATA, '/intermediate/')
FINAL = paste0(DATA, '/final/')

RESULTS <- paste0(PROJ, '/results/')
FIGURES <- paste0(RESULTS, 'figures/')
PICTURES <- paste0(RESULTS, 'pictures/')


#### Your code here ###


''' % args.name)



if __name__ == "__main__":
	main()
