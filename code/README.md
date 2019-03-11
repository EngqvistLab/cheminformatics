Folder containing all programmatic code relating to the project.


The file "template_engine.py" can be used to generate python scripts with standard code snippets inside, whereas "template_engine_R.py" and "template_engine_Rmd.py" are used to make R and R-markdown scripts.


One of the main purposes of the boilerplate code inside the generated scripts is to locate where on the system a project is located, to construct folder paths (to data folders etc.) using this information, and to provide these paths in defined variables. This should make the project portable, enabling several people to collaborate on code. Kindly __make use of these files when creating new scripts__.


`python3 template_engine.py test_file.py`


`python3 template_engine_R.py test_file.R`


`python3 template_engine_Rmd.py test_file.Rmd`
