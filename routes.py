import os
import re

from flask import (
	render_template,
	redirect,
	flash,
	url_for,
	request, jsonify
)

from sqlalchemy.exc import (
	InvalidRequestError)

from flask_bcrypt import check_password_hash

from flask_login import (
	login_user,
	logout_user,
	login_required, current_user,
)

import bio_algos.dot_plot
import models
from app import create_app, db, login_manager, bcrypt
from bio_algos.phylo_tree import generate_tree
from bio_algos.sequence_profile import amino_acid_composition, hydrophobicity, ss_propensity
from bio_algos.siRNA import *
from bio_algos.utilities import *
from bio_algos import stacked_bar_chart, heat_map
from forms import login_form, register_form
from models import fetch_records, User

app = create_app()


# Keeps track of current user object
# User loader function
@login_manager.user_loader
def load_user(user_id):
	return User.query.get(int(user_id))


@app.route('/profile')
@login_required
def profile():
	if current_user.is_authenticated:
		user_id = current_user.id
		# Use the user_id as needed
		return f'User ID: {user_id}'
	else:
		return 'User not authenticated'


# Home

# @app.route('/', methods=['GET'])
# def index():
#     # Your logic for handling the root URL
#     return render_template('index.html')

@app.route('/', methods=("GET", "POST"), strict_slashes=False)
def index():
	return render_template('index.html', title="Home", active_page='home')


@app.route("/index.html")
def home():
	return redirect(url_for("index"))


@app.route("/index")
def home1():
	return redirect(url_for("index"))


@app.route("/search", methods=["POST"], strict_slashes=False)
def search():
	# Get data from employee form
	query = request.form.get("search")
	results = fetch_records(query)
	return render_template('results.html', query=query, results=results)


@app.route('/compile_report', methods=["Get", "POST"], strict_slashes=False)
def compile_report():
	try:
		records = models.Record.query.all()
		nucleotide_ids = []
		organisms = []
		nucleotides = []
		for record in records:
			nucleotide_ids.append(record.nucleotide_id)
			organisms.append(record.organism)
			nucleotides.append(record.nucleotides)

		# create the Report object
		report = models.Report(
			nucleotide_ids=nucleotide_ids,
			organisms=organisms,
			nucleotides=nucleotides,
			employee_id=current_user.id
		)
		db.session.add(report)
		db.session.commit()

		# associate the records with the report
		report.associated_records.extend(records)
		db.session.commit()

		# create phylo tree
		organisms = report.organisms
		count = {}
		# Iterate through the list and update the count of each element
		for i in range(len(organisms)):
			if organisms[i] in count:
				count[organisms[i]] += 1
				organisms[i] = f"{organisms[i]}{count[organisms[i]]}"
			else:
				count[organisms[i]] = 1

		phylo_tree_path = f"./static/images/graphs/phylo_tree/tree{report.id}.png"
		generate_tree(nucleotides, phylo_tree_path, organisms)
		report.phylo_tree = phylo_tree_path

		dot_paths = []
		heat_paths = []
		for i in range(len(organisms)):
			for j in range(i + 1, len(organisms)):
				sequences = [nucleotides[i], nucleotides[j]]
				dot_paths.append(f"./static/images/graphs/dot_plot/dot{report.id}-{i}.png")
				# create dot graph
				bio_algos.dot_plot.dot_plot(sequences, [organisms[i], organisms[j]],
											f"./static/images/graphs/dot_plot/dot{report.id}-{i}.png")
				# create heat map
				heat_paths.append(f"./static/images/graphs/heat_map/heat{report.id}-{i}.png")
				heat_map.heat_map(sequences, [organisms[i], organisms[j]],
								  f"./static/images/graphs/heat_map/heat{report.id}-{i}.png")

		# input lists of paths into db
		report.dot_line_graph = dot_paths
		report.heat_map = heat_paths

		# create stacked bar chart
		bar_chart_path = f"./static/images/graphs/stacked_bar/bar{report.id}.png"
		stacked_bar_chart.stacked_bar_chart(nucleotides, organisms, bar_chart_path, )
		report.bar_chart = bar_chart_path

		db.session.commit()

		return redirect(url_for('employee'))  # Redirect to the employee route

	except Exception as e:
		return f"Error compiling report: {str(e)}"


# Display Report
@app.route('/display_report')
@app.route('/display_report/<int:report_id>', methods=("GET", "POST"), strict_slashes=False)
def display_report(report_id):
	report = models.Report.query.get(report_id)
	if report:
		graph_folders = os.listdir('./static/images/graphs')
		# Create a dictionary to store the image filenames for each folder
		graph_images = {}
		for folder in graph_folders:
			folder_path = os.path.join('./static/images/graphs', folder)
			graph_images[folder] = [image for image in os.listdir(folder_path) if str(report_id) in image]

		filtered_dict = {folder: images for folder, images in graph_images.items() if images}

		return render_template('display_report.html', report=report, graph_folders=filtered_dict.keys(),
							   graph_images=filtered_dict)
	else:
		return "Report not found", 404


# Employee
@app.route('/employee.html/', methods=("GET", "POST"), strict_slashes=False)
@app.route('/employee.html/<selected_result>', methods=("GET", "POST"), strict_slashes=False)
@login_required  # Make sure the user is logged in
def employee(selected_result=None):
	employee_id = current_user.id  # Get the employee ID from the current user

	if selected_result:
		# extract gene data from the selected_result from URL param, and extract relevant information
		lines = selected_result.split(" ")
		nucleotide_id = lines[0][1::]
		lines = lines[1::]

		nucleotides = re.sub(r'[^ACTG]', '', lines[-1])

		lines = lines[0:-1]
		organism = lines[0] + " " + lines[1]
		for i in range(len(lines)):
			if lines[i] == 'of':
				lines[i] = ''
		lines = lines[2:]

		gene_info = " ".join(lines)
		nucleotides = ''.join(x for x in nucleotides if not x.islower())

		# convert DNA nuc string to RNA
		rna_seq = dna_to_rna(nucleotides)

		# calculate siRNA target and gc content
		siRNA_target, gc_content = select_target_sequence(rna_seq)
		# compute sense and antisense strands https://en.wikipedia.org/wiki/Sense_(molecular_biology)
		sense_strand, antisense_strand = create_rna_strands(siRNA_target)
		# calculate similarity between sense and antisense strands (Uses Jaccard similarity coefficient)
		# minimizes off target effects and maximizes siRNA efficiency
		sense_similarity = calculate_similarity(sense_strand, antisense_strand)
		mole_weight = calculate_molecular_weight(rna_seq)
		melting_temp = calculate_melting_temp(rna_seq)

		# generates potential siRNA sequences from RNA
		siRNA_candidates = design_siRNA(rna_seq)

		efficiency_scores = {}

		for candidate in siRNA_candidates:
			efficiency_scores[candidate] = predict_efficiency(candidate)

		final_candidates = []

		max_score = max(efficiency_scores.values())

		for candidate, score in efficiency_scores.items():
			if score == max_score:
				final_candidates.append(candidate)

		siRNA_choice = final_candidates[0]

		amino_acids = transcribe_dna(nucleotides)

		amino_dict = amino_acid_composition(rna_seq)

		hydrophobicity_score = hydrophobicity(amino_dict)

		secondary_structure = ss_propensity(amino_dict)

		# create the Record object
		record = models.Record(
			nucleotide_id=nucleotide_id,
			organism=organism,
			gene_info=gene_info,
			nucleotides="".join(str(part) for part in nucleotides),  # Concatenate nucleotides

			# perform single record calculations
			gc_content=gc_content,
			amino_acids=amino_acids,
			sense_similarity=sense_similarity,
			molecular_weight=mole_weight,
			melting_temp=melting_temp,
			siRNA=siRNA_choice,
			hydrophobicity=hydrophobicity_score,
			secondary_structure_prediction=secondary_structure
		)

		# Add the objects to the session and commit to the database
		db.session.add(record)
		db.session.commit()

	reports = models.Report.query.filter_by(employee_id=employee_id).all()

	return render_template('employee.html', title="Employee", active_page='employee', reports=reports)


# Manager
@app.route('/manager.html/', methods=("GET", "POST"), strict_slashes=False)
def manager():
	reports = models.Report.query.all()
	return render_template('manager.html', title="Manager", active_page='manager', reports=reports)

@app.route('/get_employee_reports/<int:employee_id>', methods=['GET'])
def get_employee_reports(employee_id):
	employee_id = current_user.id
	if employee_id:
		reports = [
			{
				'id': report.id,
				'name': report.name,
			}
			for report in employee_id.reports
		]
		return jsonify({'reports': reports})
	else:
		return jsonify({'error': 'Employee not found'}), 404

# Admin
@app.route('/admin.html/', methods=("GET", "POST"), strict_slashes=False)
def admin():
	users = User.query.all()
	return render_template('admin.html', title="Admin", users=users, active_page='admin')


# Login route
@app.route("/login.html/", methods=("GET", "POST"), strict_slashes=False)
def login():
	form = login_form()

	if form.validate_on_submit():
		try:
			user = User.query.filter_by(email=form.email.data).first()
			if check_password_hash(user.password, form.password.data):
				login_user(user)
				return redirect(url_for('index'))
			else:
				flash("Invalid username or password!", "danger")
		except Exception as e:
			flash(e, "danger")

	return render_template("auth.html", form=form, active_page='login')


# Register route
@app.route("/register", methods=("GET", "POST"), strict_slashes=False)
def register():
	form = register_form()
	if form.validate_on_submit():
		try:
			email = form.email.data
			password = form.password.data
			username = form.username.data

			new_user = User(
				username=username,
				email=email,
				password=bcrypt.generate_password_hash(password),
			)

			db.session.add(new_user)
			db.session.commit()
			flash(f"Account created", "success")
			return redirect(url_for('admin'))

		except InvalidRequestError:
			db.session.rollback()
			flash(f"Something went wrong!", "danger")

	return render_template("auth.html", form=form)


@app.route("/logout")
@login_required
def logout():
	logout_user()
	return redirect(url_for('login'))


# Routes for the db ----------------------------------------------------------------------------------------------------

@app.route("/delete_user/<int:user_id>", methods=["DELETE"])
def deleteUser(user_id):
	try:
		user = User.query.get(user_id)
		if user:
			db.session.delete(user)
			db.session.commit()
			return f"User ID {user_id} deleted"
		else:
			return f"User ID {user_id} not found", 404
	except Exception as e:
		return f"Error deleting user: {str(e)}", 500


@app.route("/change_username/<int:user_id>", methods=["PUT"])
def changeUsername(user_id):
	try:
		data = request.get_json()  # Retrieve data from request body
		new_username = data.get('new_username')  # Get the new username

		user = User.query.get(user_id)
		if user:
			user.username = new_username  # Update the username
			db.session.commit()
			return f"User ID {user_id} username changed to {new_username} successfully"
		else:
			return f"User ID {user_id} not found", 404
	except Exception as e:
		return f"Error changing username: {str(e)}", 500


@app.route("/change_password/<int:user_id>", methods=["PUT"])
def changePassword(user_id):
	try:
		data = request.get_json()  # Retrieve data from request body
		new_password = data.get('new_password')

		user = User.query.get(user_id)
		if user:
			user.password = bcrypt.generate_password_hash(new_password),
			db.session.commit()
			return f"User ID {user_id} password changed successfully"
		else:
			return f"User ID {user_id} not found", 404
	except Exception as e:
		return f"Error changing password: {str(e)}", 500


if __name__ == "__main__":
	app.run(debug=True)
