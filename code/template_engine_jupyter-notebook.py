#!/usr/bin/env python3
"""
template_engine_python.py is used to generate new python scripts with useful boilerplate code

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

parser = argparse.ArgumentParser(prog="template_engine_jupyter-notebook.py")
parser.add_argument("name", nargs='?', help="Set the name for the script to create", default="new_script.ipynb")
args = parser.parse_args()


def main():
	output_file = open(args.name, "w")

	output_file.write('''
{
 "cells": [
   {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Modify this line to briefly discribe the functionality of %s<br/>",
    "<br/>",
    "Copyright (C) 2017  Martin Engqvist Lab<br/>",
    "This program is free software: you can redistribute it and/or modify<br/>",
    "it under the terms of the GNU General Public License as published by<br/>",
    "the Free Software Foundation, either version 3 of the License, or<br/>",
    "(at your option) any later version.<br/>",
    "This program is distributed in the hope that it will be useful,<br/>",
    "but WITHOUT ANY WARRANTY; without even the implied warranty of<br/>",
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the<br/>",
    "GNU General Public License for more details.<br/>",
    "You should have received a copy of the GNU General Public License<br/>",
    "along with this program.  If not, see <http://www.gnu.org/licenses/>."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\\n",
    "from dotenv import load_dotenv, find_dotenv\\n",
    "from os.path import join, dirname, basename, exists, isdir\\n",
    "\\n",
    "### Load environmental variables from the project root directory ###\\n",
    "# find .env automagically by walking up directories until it's found\\n",
    "dotenv_path = find_dotenv()\\n",
    "\\n",
    "# load up the entries as environment variables\\n",
    "load_dotenv(dotenv_path)\\n",
    "\\n",
    "# now you can get the variables using their names\\n",
    "\\n",
    "# Check whether a network drive has been specified\\n",
    "DATABASE = os.environ.get(\\"NETWORK_URL\\")\\n",
    "if DATABASE == 'None':\\n",
    "    pass\\n",
    "else:\\n",
    "    pass\\n",
    "    #mount network drive here\\n",
    "\\n",
    "# set up directory paths\\n",
    "CURRENT_DIR = os.getcwd()\\n",
    "PROJ = dirname(dotenv_path) # project root directory\\n",
    "\\n",
    "DATA = join(PROJ, 'data') #data directory\\n",
    "RAW_EXTERNAL = join(DATA, 'raw_external') # external data raw directory\\n",
    "RAW_INTERNAL = join(DATA, 'raw_internal') # internal data raw directory\\n",
    "INTERMEDIATE = join(DATA, 'intermediate') # intermediate data directory\\n",
    "FINAL = join(DATA, 'final') # final data directory\\n",
    "\\n",
    "RESULTS = join(PROJ, 'results') # output directory\\n",
    "FIGURES = join(RESULTS, 'figures') # figure output directory\\n",
    "PICTURES = join(RESULTS, 'pictures') # picture output directory\\n",
    "\\n",
    "\\n",
    "# make folders specific for certain data\\n",
    "folder_name = ''\\n",
    "if folder_name != '':\\n",
    "    #make folders if they don't exist\\n",
    "    if not exists(join(RAW_EXTERNAL, folder_name)):\\n",
    "        os.makedirs(join(RAW_EXTERNAL, folder_name))\\n",
    "\\n",
    "    if not exists(join(INTERMEDIATE, folder_name)):\\n",
    "        os.makedirs(join(INTERMEDIATE, folder_name))\\n",
    "\\n",
    "    if not exists(join(FINAL, folder_name)):\\n",
    "        os.makedirs(join(FINAL, folder_name))",
	"\\n\\n",
	"print('Standard variables loaded, you are good to go!')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "## Your code here ##"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.3"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}


''' % args.name)



if __name__ == "__main__":
	main()
