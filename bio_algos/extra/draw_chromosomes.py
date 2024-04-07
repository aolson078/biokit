from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome

entries = [
	("Chr I", 30432564),
	("Chr II", 19705359),
	("Chr III", 23470805),
	("Chr IV", 18585042),
	("Chr V", 26992728),
]

max_len = 30432563
telomere_length = 1000000

chr_diagram = BasicChromosome.Organism()
chr_diagram.page_size = (29.7 * cm, 21 * cm)  # A4 landscape

for name, length in entries:
	cur_chromosome = BasicChromosome.Chromosome(name)
	# set scale to max len plus 2 telemores
	# want same scale for all 5 so they can be compared
	cur_chromosome.scale_num = max_len + 2 * telomere_length