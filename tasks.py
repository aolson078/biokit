"""Celery tasks for the BioKit application."""

from flask_bio_app import create_app
from models import Record, Report, db
from bio_algos.phylo_tree import generate_tree
from bio_algos import dot_plot, heat_map, stacked_bar_chart
from bio_algos.gc_content_line import gc_line_graph
from bio_algos.gc_skew import gc_skew_plot
from bio_algos.nucleotide_pie import nucleotide_pie_chart

app = create_app()
celery = app.celery

@celery.task()
def compile_report_task(employee_id):
    records = Record.query.all()
    if len(records) <= 1:
        return None

    nucleotide_ids = []
    organisms = []
    nucleotides = []
    for record in records:
        nucleotide_ids.append(record.nucleotide_id)
        organisms.append(record.organism)
        nucleotides.append(record.nucleotides)

    report = Report(
        nucleotide_ids=nucleotide_ids,
        organisms=organisms,
        nucleotides=nucleotides,
        employee_id=employee_id
    )
    db.session.add(report)
    db.session.commit()

    report.associated_records.extend(records)
    db.session.commit()

    count = {}
    for i in range(len(organisms)):
        if organisms[i] in count:
            count[organisms[i]] += 1
            organisms[i] = f"{organisms[i]}{count[organisms[i]]}"
        else:
            count[organisms[i]] = 1

    phylo_tree_path = f"./static/graphs/phylo_tree/tree{report.id}.png"
    generate_tree(nucleotides, phylo_tree_path, organisms)
    report.phylo_tree = phylo_tree_path

    dot_paths = []
    heat_paths = []
    for i in range(len(organisms)):
        for j in range(i + 1, len(organisms)):
            sequences = [nucleotides[i], nucleotides[j]]
            dot_path = f"./static/graphs/dot_plot/dot{report.id}-{i}.png"
            dot_plot.dot_plot(sequences, [organisms[i], organisms[j]], dot_path)
            dot_paths.append(dot_path)
            heat_path = f"./static/graphs/heat_map/heat{report.id}-{i}.png"
            heat_map.heat_map(sequences, [organisms[i], organisms[j]], heat_path)
            heat_paths.append(heat_path)

    report.dot_line_graph = dot_paths
    report.heat_map = heat_paths

    bar_chart_path = f"./static/graphs/stacked_bar/bar{report.id}.png"
    stacked_bar_chart.stacked_bar_chart(nucleotides, organisms, bar_chart_path)
    report.bar_chart = bar_chart_path

    gc_paths = []
    skew_paths = []
    pie_paths = []
    for idx, seq in enumerate(nucleotides):
        line_path = f"./static/graphs/gc_line/line{report.id}-{idx}.png"
        gc_line_graph(seq, output=line_path)
        gc_paths.append(line_path)

        skew_path = f"./static/graphs/gc_skew/skew{report.id}-{idx}.png"
        gc_skew_plot(seq, output=skew_path)
        skew_paths.append(skew_path)

        pie_path = f"./static/graphs/nuc_pie/pie{report.id}-{idx}.png"
        nucleotide_pie_chart(seq, output=pie_path)
        pie_paths.append(pie_path)
    report.gc_line_graphs = gc_paths
    report.gc_skew_graphs = skew_paths
    report.nuc_pie_charts = pie_paths

    db.session.commit()
    return report.id
