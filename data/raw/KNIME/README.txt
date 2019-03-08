
## From DOI: 10.1021/acs.jcim.8b00640


KNIME Annotated Workflow for Drug Discovery Maps application
------------------------------------------------------------

To import in KNIME: simply double-click the .knar archive and a workflow group will be created 
containing all needed data and the workflow in your local KNIME workspace.

----------

To test: run 'execute all' for an example application using Erlotinib. The output prediction can be found in the Interactive Table node.

----------

To use: draw or paste a molecule of your choosing in the Marvin Sketch node. Optionally tune the parameters using the dedicated nodes. Execute all nodes. The output prediction is in the Interactive Table node.

----------

Dependencies:
This workflow was created in the KNIME version 3.6 but is expected to work in older build of version 3.x. 
Several plugin dependencies are in place (RDKit, Erlwood and others) but these will be identified and installed by KNIME upon opening the workflow (after approval by user).
Python 3 needs to be installed and configured for use in KNIME including the relevant Python modules (pandas, numpy, scipy, scikit-learn)
