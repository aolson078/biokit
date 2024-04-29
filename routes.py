import os
import re

from flask import (
	render_template,
	redirect,
	flash,
	url_for,
	request, jsonify, session,  send_file, abort, make_response
)


from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

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
from custom_exceptions import NoRecordsError
from forms import login_form, register_form
from models import *
from bio_algos.utilities import fetch_records

# create instance of flask app
app = create_app()


# keeps track of current user object
# user loader function
@login_manager.user_loader
def load_user(user_id):
	return User.query.get(int(user_id))

# home/index route
@app.route('/', methods=("GET", "POST"), strict_slashes=False)
def index():
	return render_template('index.html', title="Home", active_page='home')

# home/index route
@app.route("/index")
def home():
	return redirect(url_for("index"))

# route for when user lacks permission to access a route/resource
@app.route("/denied", methods=["GET"], strict_slashes=False)
def denied():
	return render_template("denied.html")





# route to handle the user query and form data
@app.route("/search", methods=["POST"], strict_slashes=False)
@login_required
def search():
	# Get data from employee form
	query = request.json.get('search')
	if not query:
		return jsonify({'error': 'Empty search query'}), 400

	results = fetch_records(query)
	return jsonify(results)



# route to handle processing the records into a report and generating/saving graphs
@app.route('/compile_report', methods=["GET", "POST"], strict_slashes=False)
@login_required
def compile_report():
	try:
		records = models.Record.query.all()
		if len(records) <= 1:
			raise NoRecordsError("Two or more records needed to generate report.")

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


	except NoRecordsError as e:
		print("No records to compile report")
		return redirect(url_for("employee"))


# route to handle building the visible report for the end user.
@app.route('/display_report')
@app.route('/display_report/<int:report_id>', methods=("GET", "POST"), strict_slashes=False)
@login_required
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


# route to handle the employee user interactions, compiling the report, displaying report
@app.route('/employee.html/', methods=("GET", "POST"), strict_slashes=False)
@app.route('/employee.html/<selected_result>', methods=("GET", "POST"), strict_slashes=False)
@login_required
def employee():
	# Get the employee ID from the current user
	employee_id = current_user.id

	if request.method == "POST":
		selected_result = request.json.get("selected_result")
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
				secondary_structure_prediction=secondary_structure,

				employee_id=employee_id
			)

			# Add the objects to the session and commit to the database
			db.session.add(record)
			db.session.commit()
			return redirect(url_for('employee'))

	# retrieve all reports made by the current logged in user.
	reports = models.Report.query.filter_by(employee_id=employee_id).all()

	# retrieve all records from current user by descending order of id
	all_records = models.Record.query.order_by(models.Record.id.desc()).all()

	return render_template('employee.html', title="Employee", active_page='employee', reports=reports,
	                       all_records=all_records)

# route to handle retrieving the logged in employees reports from the database
@app.route('/get_user_reports/<int:user_id>', methods=['GET'])
# @login_required
def get_user_reports(user_id):
    reports = Report.query.filter_by(employee_id=user_id).all()
    if reports:
        # Extract relevant data from each report
        report_data = [
            {
                'id': report.id,
                'nucleotide_ids': report.nucleotide_ids,
                'organisms': report.organisms,
                # Add other fields as needed
            }
            for report in reports
        ]
        return jsonify({'reports': report_data})
    else:
        return jsonify({'error': 'No reports found for this user'}), 404


''' This could work too for deleting reports (sola)
@app.route("/delete_report/<int:report_id>", methods=["DELETE"])
def delete_report(report_id):
    # Retrieve the report
    report = Report.query.filter_by(id=report_id).first()

    if report:
        try:
            # Remove all associated records
            report_record.report_id.clear()
            # report.record.clear()
            # Delete the report from the database
            db.session.delete(report)
            db.session.commit()

            return jsonify({'message': 'Report deleted successfully'}), 200
        except Exception as e:
            db.session.rollback()
            return jsonify({'error': f'Error deleting report: {str(e)}'}), 500
    else:
        return jsonify({'error': 'Report not found'}), 404
'''

# Route for downloading a PDF from the report table
@app.route('/download_report/<int:report_id>', methods=['GET'])
def download_report(report_id):
    # Retrieve the report data based on the report_id
    report = Report.query.get(report_id)

    if report:
        # Generate the PDF report
        buffer = BytesIO()
        pdf = canvas.Canvas(buffer, pagesize=letter)
        pdf.drawString(100, 750, f'Report ID: {report.id}')
        pdf.drawString(100, 730, f'Employee ID: {report.employee_id}')
        pdf.drawString(100, 710, f'Nucleotide IDs: {", ".join(report.nucleotide_ids)}')
        pdf.drawString(100, 690, f'Organisms: {", ".join(report.organisms)}')
        pdf.drawString(100, 670, f'Nucleotides: {", ".join(report.nucleotides)}')
        pdf.drawString(100, 650, f'Phylogenetic Tree: {report.phylo_tree}')
        pdf.drawString(100, 630, f'Dot Line Graph: {report.dot_line_graph}')
        pdf.drawString(100, 610, f'Heat Map: {report.heat_map}')
        pdf.drawString(100, 590, f'Bar Chart: {report.bar_chart}')

        # Save the PDF to the buffer
        pdf.save()
        buffer.seek(0)

        # Create a response with the PDF content
        response = make_response(buffer.getvalue())
        response.headers['Content-Type'] = 'application/pdf'

        # Set the Content-Disposition header to specify the filename as an attachment
        response.headers['Content-Disposition'] = f'attachment; filename=report_{report.id}.pdf'

        return response
    else:
        abort(404)


### MANAGER ###

# route to handle manager user interactions, changing report settings, save/print/deleting reports etc
@app.route('/manager.html/', methods=("GET", "POST"), strict_slashes=False)
# @login_required
# @models.is_manager
def manager():
	users = User.query.all()
	reports = Report.query.all()
	return render_template('manager.html', title="Manager", active_page='manager', reports=reports, users=users)


#Will (delete reports)
@app.route('/delete_report/<int:report_id>', methods=['DELETE'])
@login_required
def delete_report(report_id):
    report = models.Report.query.get(report_id)
    if report:
        db.session.delete(report)
        db.session.commit()
        return jsonify({'message': 'Report deleted successfully'}), 200
    else:
        return jsonify({'error': 'Report not found'}), 404

@app.route('/update_permissions', methods=['POST'])
def update_permissions():
	if request.method == 'POST':
		# Extract data from the form submission
		user_id = request.form.get('userId')
		view_reports = 'viewReports' in request.form
		delete_reports = 'deleteReports' in request.form
		print_reports = 'printReports' in request.form
		change_reports = 'changeReports' in request.form

		# Find the user by ID
		user = User.query.get(user_id)
		if user:
			# Update the user's permissions
			user.view_reports = view_reports
			user.delete_reports = delete_reports
			user.print_reports = print_reports
			user.change_reports = change_reports
			db.session.commit()
			flash('Permissions updated successfully', 'success')
			return redirect(url_for('admin'))
		else:
			flash('User not found', 'error')
			return redirect(url_for('admin'))


### ADMIN ###

# route to handle admin user interactions, creating new users, changing
@app.route('/admin.html/', methods=("GET", "POST"), strict_slashes=False)
# @login_required
# @models.is_admin
def admin():
	users = User.query.all()
	return render_template('admin.html', title="Admin", users=users, active_page='admin')


# route to handled user login and authentication
@app.route("/login.html/", methods=("GET", "POST"), strict_slashes=False)
def login():
	form = login_form()

	if form.validate_on_submit():
		try:
			user = User.query.filter_by(email=form.email.data).first()
			if check_password_hash(user.password, form.password.data):
				login_user(user)
				session['logged_in'] = True
				return redirect(url_for('index'))
			else:
				flash("Invalid username or password!", "danger")
		except Exception as e:
			flash(e, "danger")

	return render_template("auth.html", form=form, active_page='login')


# route to handle registering a new user (admin)
@app.route("/register", methods=("GET", "POST"), strict_slashes=False)
# @login_required
def register():
	form = register_form()
	if form.validate_on_submit():
		try:
			email = form.email.data
			password = form.password.data
			username = form.username.data
			role = form.role.data

			new_user = User(
				username=username,
				email=email,
				password=bcrypt.generate_password_hash(password),
				role=role
			)

			db.session.add(new_user)
			db.session.commit()
			flash(f"Account created", "success")
			return redirect(url_for('admin'))

		except InvalidRequestError:
			db.session.rollback()
			flash(f"Something went wrong!", "danger")

	return render_template("auth.html", form=form)


# route to handle logging current user out
@app.route("/logout")
@login_required
def logout():
	logout_user()
	session.pop('logged_in', None)
	return redirect(url_for('login'))


# Routes for the db ----------------------------------------------------------------------------------------------------
# route to handle deleting user from database (admin)
# @app.route("/delete_user/<int:user_id>", methods=["DELETE"])
# @login_required
# def deleteUser(user_id):
# 	try:
# 		user = User.query.get(user_id)
# 		if user:
# 			db.session.delete(user)
# 			db.session.commit()
# 			return f"User ID {user_id} deleted"
# 		else:
# 			return f"User ID {user_id} not found", 404
# 	except Exception as e:
# 		return f"Error deleting user: {str(e)}", 500

#added
@app.route("/delete_user/<int:user_id>", methods=["DELETE"])
@login_required
def delete_user(user_id):
    try:
        user = User.query.get(user_id)
        if user:
            db.session.delete(user)
            db.session.commit()
            return jsonify({'message': f'User ID {user_id} deleted successfully'}), 200
        else:
            return jsonify({'error': f'User ID {user_id} not found'}), 404
    except Exception as e:
        return jsonify({'error': f'Error deleting user: {str(e)}'}), 500



# route to handle changing username in database (admin)
@app.route("/change_username/<int:user_id>", methods=["PUT"])
@login_required
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

# route to handle changing password in database (admin)
@app.route("/change_password/<int:user_id>", methods=["PUT"])
@login_required
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

@app.route('/process_selected_user', methods=['POST'])
def process_selected_user():
    data = request.get_json()
    user_id = data.get('userId')
    username = data.get('username')
    return jsonify({'message': f'Selected user {username} with ID {user_id} processed successfully'})


if __name__ == "__main__":
	app.run(debug=True)
