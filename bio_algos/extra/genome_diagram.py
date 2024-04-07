from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SeqIO.read("NC_005816.gb", "genbank")
print(record)

# create empty diagram
gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
# add empty track to diagram
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
# add empty feature set to track
gd_feature_set = gd_track_for_features.new_set()

print("----------------------------------------------------------------------------")

for feature in record.features:
	if feature.type != "gene":
		continue
	if len(gd_feature_set) % 2 == 0:
		color = colors.blue
	else:
		color = colors.lightblue
	gd_feature_set.add_feature(feature, color=color, label=True, sigil="ARROW", label_size=8, label_angle=30)

	gd_diagram.draw(
		format="circular", # linear or circular
		orientation="landscape",
		pagesize="A4",
		fragments=1,
		start=0,
		end=len(record),
	)

	gd_feature_set.add_feature(feature, sigil="ARROW", color="orange", arrowhead_length=1)

	gd_diagram.write("./diagrams/plasmid_linear.pdf", "PDF")
	gd_diagram.write("./diagrams/plasmid_linear.eps", "EPS")
	gd_diagram.write("./diagrams/plasmid_linear.svg", "SVG")

