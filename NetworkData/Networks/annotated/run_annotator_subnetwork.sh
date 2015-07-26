# PF AD
python ../../../Code/network_annotater.py ../PF_AD_links_for_cytoscape.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/PF_AD_RAS_deg2/PF_AD_RAS_deg2.txt.sif.nodes
python ../../../Code/network_annotater.py ../PF_AD_links_for_cytoscape.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/PF_AD_RAS_deg3/PF_AD_RAS_deg3.txt.sif.nodes

# PF Normal
python ../../../Code/network_annotater.py ../PF_normal_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/PF_normal_RAS_deg2/PF_normal_RAS_deg2.sif.nodes
python ../../../Code/network_annotater.py ../PF_normal_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/PF_normal_RAS_deg3/PF_normal_RAS_deg3.sif.nodes

# Delimited Blood
python ../../../Code/network_annotater.py ../delimited_blood_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_blood_deg2/delimited_blood_deg2.sif.nodes
python ../../../Code/network_annotater.py ../delimited_blood_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_blood_deg3/delimited_blood_deg3.sif.nodes

# Delimited pan intestine
python ../../../Code/network_annotater.py ../delimited_pan_intestine_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_pan_intestine_deg2/delimited_pan_intestine_deg2.sif
python ../../../Code/network_annotater.py ../delimited_pan_intestine_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_pan_intestine_deg3/delimited_pan_intestine_deg3.sif

# Delimited risk
python ../../../Code/network_annotater.py ../delimited_risk_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_risk_deg2/delimited_risk_deg2.sif.nodes
python ../../../Code/network_annotater.py ../delimited_risk_network.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/delimited_risk_deg3/delimited_risk_deg3.sif.nodes

# Ileum links
python ../../../Code/network_annotater.py ../ileum_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/ileum_links_deg2/ileum_links_deg2.sif.nodes
python ../../../Code/network_annotater.py ../ileum_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/ileum_links_deg3/ileum_links_deg3.txt.sif.nodes

# Omental links
python ../../../Code/network_annotater.py ../omental_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/omental_links_RAS_deg2/omental_links_RAS_deg2.txt.sif.nodes
python ../../../Code/network_annotater.py ../omental_links_space_delimited.txt /Users/kwj/Projects/GelbRotation/NetworkData/Subnetworks/omental_links_deg3/omental_links_deg3.txt.sif.nodes

mv ../*anno*txt .
