"""
Upgraded Flask Bioinformatics Application
A modern, secure, and scalable bioinformatics web application.
"""

import os
import re
import logging
from typing import Dict, List, Optional, Tuple, Any
from datetime import datetime, timedelta
from functools import wraps
from dataclasses import dataclass
from pathlib import Path

from flask import (
    Flask, Blueprint, render_template, redirect, flash,
    url_for, request, jsonify, session, send_file,
    abort, make_response, current_app
)
from flask_login import (
    LoginManager, login_user, logout_user,
    login_required, current_user
)
from flask_wtf import FlaskForm
from flask_wtf.csrf import CSRFProtect
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from werkzeug.security import generate_password_hash, check_password_hash
from werkzeug.utils import secure_filename
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy.orm import joinedload
from marshmallow import Schema, fields, validate, ValidationError
from celery import Celery

from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Import your bio algorithms (assuming they're refactored)
from bio_algos import (
    dot_plot, phylo_tree, sequence_profile,
    siRNA, utilities, stacked_bar_chart, heat_map
)
from bio_algos.gc_content_line import gc_line_graph
from bio_algos.gc_skew import gc_skew_plot
from bio_algos.nucleotide_pie import nucleotide_pie_chart

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Initialize extensions early so decorators can use them
db = SQLAlchemy()
login_manager = LoginManager()
csrf = CSRFProtect()
cache = Cache()
limiter = Limiter(key_func=get_remote_address)
celery = Celery(__name__)


# Configuration classes
class Config:
    """Base configuration."""
    SECRET_KEY = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL', 'sqlite:///bioinformatics.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    
    # Security
    SESSION_COOKIE_SECURE = True
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'
    WTF_CSRF_TIME_LIMIT = None
    
    # File upload
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB max file size
    UPLOAD_FOLDER = 'uploads'
    ALLOWED_EXTENSIONS = {'txt', 'fasta', 'fa', 'fastq'}
    
    # Celery
    CELERY_BROKER_URL = os.environ.get('REDIS_URL', 'redis://localhost:6379/0')
    CELERY_RESULT_BACKEND = os.environ.get('REDIS_URL', 'redis://localhost:6379/0')
    
    # Cache
    CACHE_TYPE = 'redis'
    CACHE_REDIS_URL = os.environ.get('REDIS_URL', 'redis://localhost:6379/1')
    
    # NCBI
    NCBI_EMAIL = os.environ.get('NCBI_EMAIL', 'your-email@example.com')
    NCBI_API_KEY = os.environ.get('NCBI_API_KEY')
    
    # Rate limiting
    RATELIMIT_STORAGE_URL = os.environ.get('REDIS_URL', 'redis://localhost:6379/2')


# Validation Schemas
class SearchSchema(Schema):
    """Validation schema for search requests."""
    search = fields.Str(required=True, validate=validate.Length(min=1, max=200))


class ReportSettingsSchema(Schema):
    """Validation schema for report settings."""
    hide_nucleotide_id = fields.Bool(load_default=False)
    hide_organism = fields.Bool(load_default=False)
    hide_nucleotide = fields.Bool(load_default=False)
    hide_phylogenetic = fields.Bool(load_default=False)


class UserPermissionsSchema(Schema):
    """Validation schema for user permissions."""
    user_id = fields.Int(required=True)
    view_reports = fields.Bool(load_default=False)
    delete_reports = fields.Bool(load_default=False)
    print_reports = fields.Bool(load_default=False)
    change_reports = fields.Bool(load_default=False)


# Custom decorators
def role_required(role):
    """Decorator to check if user has required role."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if not current_user.is_authenticated:
                return redirect(url_for('auth.login'))
            if current_user.role != role:
                abort(403)
            return f(*args, **kwargs)
        return decorated_function
    return decorator


def validate_json(schema_class):
    """Decorator to validate JSON input."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            schema = schema_class()
            try:
                data = schema.load(request.get_json())
                request.validated_data = data
            except ValidationError as err:
                return jsonify({'errors': err.messages}), 400
            return f(*args, **kwargs)
        return decorated_function
    return decorator


# Service classes
class BioinformaticsService:
    """Service class for bioinformatics operations."""
    
    @staticmethod
    def process_sequence_data(selected_result: str) -> Dict[str, Any]:
        """Process selected sequence data."""
        try:
            lines = selected_result.split(" ")
            nucleotide_id = lines[0][1::]
            
            # Extract organism name
            organism_parts = []
            gene_parts = []
            in_gene_section = False
            
            for i, part in enumerate(lines[1:]):
                if part == 'of' or in_gene_section:
                    in_gene_section = True
                    gene_parts.append(part)
                elif len(organism_parts) < 2:
                    organism_parts.append(part)
                else:
                    gene_parts.append(part)
            
            organism = " ".join(organism_parts)
            gene_info = " ".join(gene_parts)
            
            # Extract nucleotides
            nucleotides = re.sub(r'[^ACTG]', '', lines[-1])
            nucleotides = ''.join(x for x in nucleotides if not x.islower())
            
            # Perform bioinformatics calculations
            rna_seq = utilities.dna_to_rna(nucleotides)
            siRNA_target, gc_content = siRNA.select_target_sequence(rna_seq)
            sense_strand, antisense_strand = siRNA.create_rna_strands(siRNA_target)
            sense_similarity = siRNA.calculate_similarity(sense_strand, antisense_strand)
            mole_weight = siRNA.calculate_molecular_weight(rna_seq)
            melting_temp = siRNA.calculate_melting_temp(rna_seq)
            
            # Generate siRNA candidates
            siRNA_candidates = siRNA.design_siRNA(rna_seq)
            efficiency_scores = {
                candidate: siRNA.predict_efficiency(candidate)
                for candidate in siRNA_candidates
            }
            
            # Select best candidate
            max_score = max(efficiency_scores.values())
            best_candidates = [
                candidate for candidate, score in efficiency_scores.items()
                if score == max_score
            ]
            siRNA_choice = best_candidates[0] if best_candidates else None
            
            # Additional analyses
            amino_acids = utilities.transcribe_dna(nucleotides)
            amino_dict = sequence_profile.amino_acid_composition(rna_seq)
            hydrophobicity_score = sequence_profile.hydrophobicity(amino_dict)
            secondary_structure = sequence_profile.ss_propensity(amino_dict)
            
            return {
                'nucleotide_id': nucleotide_id,
                'organism': organism,
                'gene_info': gene_info,
                'nucleotides': nucleotides,
                'gc_content': gc_content,
                'amino_acids': amino_acids,
                'sense_similarity': sense_similarity,
                'molecular_weight': mole_weight,
                'melting_temp': melting_temp,
                'siRNA': siRNA_choice,
                'hydrophobicity': hydrophobicity_score,
                'secondary_structure_prediction': secondary_structure
            }
            
        except Exception as e:
            logger.error(f"Error processing sequence data: {e}")
            raise ValueError(f"Failed to process sequence data: {str(e)}")


class ReportGenerator:
    """Service class for report generation."""
    
    @staticmethod
    def generate_graphs(report_id: int, nucleotides: List[str], organisms: List[str]) -> Dict[str, Any]:
        """Generate all graphs for a report."""
        graph_paths = {
            'phylo_tree': None,
            'dot_plots': [],
            'heat_maps': [],
            'bar_chart': None,
            'gc_lines': [],
            'gc_skews': [],
            'nuc_pies': []
        }
        
        try:
            # Ensure graph directories exist
            graph_dirs = [
                'static/graphs/phylo_tree',
                'static/graphs/dot_plot',
                'static/graphs/heat_map',
                'static/graphs/stacked_bar',
                'static/graphs/gc_line',
                'static/graphs/gc_skew',
                'static/graphs/nuc_pie'
            ]
            for dir_path in graph_dirs:
                Path(dir_path).mkdir(parents=True, exist_ok=True)
            
            # Generate phylogenetic tree
            phylo_path = f"static/graphs/phylo_tree/tree{report_id}.png"
            phylo_tree.generate_tree(nucleotides, phylo_path, organisms)
            graph_paths['phylo_tree'] = phylo_path
            
            # Generate pairwise comparisons
            for i in range(len(organisms)):
                for j in range(i + 1, len(organisms)):
                    sequences = [nucleotides[i], nucleotides[j]]
                    org_names = [organisms[i], organisms[j]]
                    
                    # Dot plot
                    dot_path = f"static/graphs/dot_plot/dot{report_id}-{i}-{j}.png"
                    dot_plot.dot_plot(sequences, org_names, dot_path)
                    graph_paths['dot_plots'].append(dot_path)
                    
                    # Heat map
                    heat_path = f"static/graphs/heat_map/heat{report_id}-{i}-{j}.png"
                    heat_map.heat_map(sequences, org_names, heat_path)
                    graph_paths['heat_maps'].append(heat_path)
            
            # Generate stacked bar chart
            bar_path = f"static/graphs/stacked_bar/bar{report_id}.png"
            stacked_bar_chart.stacked_bar_chart(nucleotides, organisms, bar_path)
            graph_paths['bar_chart'] = bar_path

            # Generate GC content line graphs, GC skew plots, and nucleotide pie charts
            for idx, seq in enumerate(nucleotides):
                line_path = f"static/graphs/gc_line/line{report_id}-{idx}.png"
                gc_line_graph(seq, output=line_path)
                graph_paths['gc_lines'].append(line_path)

                skew_path = f"static/graphs/gc_skew/skew{report_id}-{idx}.png"
                gc_skew_plot(seq, output=skew_path)
                graph_paths['gc_skews'].append(skew_path)

                pie_path = f"static/graphs/nuc_pie/pie{report_id}-{idx}.png"
                nucleotide_pie_chart(seq, output=pie_path)
                graph_paths['nuc_pies'].append(pie_path)

            return graph_paths
            
        except Exception as e:
            logger.error(f"Error generating graphs for report {report_id}: {e}")
            raise


class PDFReportGenerator:
    """Service class for PDF report generation."""
    
    @staticmethod
    def generate_pdf(report) -> BytesIO:
        """Generate a professional PDF report."""
        buffer = BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter)
        story = []
        styles = getSampleStyleSheet()
        
        # Title
        title = Paragraph(f"Bioinformatics Report #{report.id}", styles['Title'])
        story.append(title)
        story.append(Spacer(1, 12))
        
        # Report metadata
        metadata = [
            ['Report ID:', str(report.id)],
            ['Generated:', datetime.now().strftime('%Y-%m-%d %H:%M')],
            ['Employee:', report.employee.username],
            ['Number of Sequences:', str(len(report.nucleotide_ids))]
        ]
        
        metadata_table = Table(metadata, colWidths=[120, 400])
        metadata_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.grey),
            ('TEXTCOLOR', (0, 0), (0, -1), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
        ]))
        story.append(metadata_table)
        story.append(Spacer(1, 20))
        
        # Sequences section
        story.append(Paragraph("Analyzed Sequences", styles['Heading2']))
        story.append(Spacer(1, 12))
        
        for i, (nuc_id, org) in enumerate(zip(report.nucleotide_ids, report.organisms)):
            seq_data = [
                ['Nucleotide ID:', nuc_id],
                ['Organism:', org],
                ['Length:', f"{len(report.nucleotides[i])} bp"],
                ['GC Content:', f"{utilities.gc_content(report.nucleotides[i]):.2f}%"]
            ]
            
            seq_table = Table(seq_data, colWidths=[120, 400])
            seq_table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ]))
            story.append(seq_table)
            story.append(Spacer(1, 10))
        
        # Analysis results
        story.append(Paragraph("Analysis Results", styles['Heading2']))
        story.append(Spacer(1, 12))
        
        analysis_text = f"""
        This report contains the following analyses:
        • Phylogenetic tree showing evolutionary relationships
        • Dot plots for sequence similarity visualization
        • Heat maps for detailed sequence comparison
        • Nucleotide composition bar chart
        • GC content trends and GC skew plots
        Generated graphs are available in the web interface.
        """
        story.append(Paragraph(analysis_text, styles['Normal']))
        
        # Build PDF
        doc.build(story)
        buffer.seek(0)
        return buffer


# Celery tasks
def make_celery(app):
    """Configure the global Celery instance."""
    celery.conf.update(
        broker_url=app.config['CELERY_BROKER_URL'],
        result_backend=app.config['CELERY_RESULT_BACKEND'],
        **app.config
    )
    
    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)
    
    celery.Task = ContextTask
    return celery


# Blueprint for authentication routes
auth_bp = Blueprint('auth', __name__)

@auth_bp.route('/login', methods=['GET', 'POST'])
@limiter.limit("5 per minute")
def login():
    """Handle user login."""
    from forms import login_form
    if current_user.is_authenticated:
        return redirect(url_for('main.index'))

    form = login_form()
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        
        if user and check_password_hash(user.password, form.password.data):
            login_user(user, remember=form.remember_me.data)
            next_page = request.args.get('next')
            
            # Security: validate next URL to prevent open redirect
            if next_page and is_safe_url(next_page):
                return redirect(next_page)
            
            return redirect(url_for('main.index'))
        
        flash('Invalid email or password', 'danger')
        logger.warning(f"Failed login attempt for email: {form.email.data}")
    
    return render_template('auth/login.html', form=form)


@auth_bp.route('/logout')
@login_required
def logout():
    """Handle user logout."""
    logout_user()
    flash('You have been logged out successfully', 'info')
    return redirect(url_for('auth.login'))


# Blueprint for main routes
main_bp = Blueprint('main', __name__)

@main_bp.route('/')
def index():
    """Home page."""
    return render_template('index.html', title="Home")


@main_bp.route('/search', methods=['POST'])
@login_required
@validate_json(SearchSchema)
@limiter.limit("10 per minute")
def search():
    """Handle NCBI search requests."""
    try:
        query = request.validated_data['search']
        
        # Use caching for repeated searches
        cache_key = f"search:{query}"
        cached_results = cache.get(cache_key)
        
        if cached_results:
            return jsonify(cached_results)
        
        # Perform search
        results = utilities.fetch_records(query)
        
        # Cache results for 1 hour
        cache.set(cache_key, results, timeout=3600)
        
        return jsonify(results)
        
    except Exception as e:
        logger.error(f"Search error: {e}")
        return jsonify({'error': 'Search failed'}), 500


@main_bp.route('/process_sequence', methods=['POST'])
@login_required
def process_sequence():
    """Process selected sequence and create record."""
    try:
        selected_result = request.json.get('selected_result')
        if not selected_result:
            return jsonify({'error': 'No sequence selected'}), 400
        
        # Process sequence data
        sequence_data = BioinformaticsService.process_sequence_data(selected_result)
        
        # Create record
        record = Record(
            **sequence_data,
            employee_id=current_user.id
        )
        
        db.session.add(record)
        db.session.commit()
        
        return jsonify({
            'message': 'Record created successfully',
            'record_id': record.id
        })
        
    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        logger.error(f"Error processing sequence: {e}")
        db.session.rollback()
        return jsonify({'error': 'Failed to process sequence'}), 500


@main_bp.route('/compile_report', methods=['POST'])
@login_required
def compile_report():
    """Compile report from selected records."""
    try:
        # Use async task for long-running report generation
        task = compile_report_task.delay(current_user.id)
        
        return jsonify({
            'message': 'Report compilation started',
            'task_id': task.id
        }), 202
        
    except Exception as e:
        logger.error(f"Error starting report compilation: {e}")
        return jsonify({'error': 'Failed to start report compilation'}), 500


@main_bp.route('/report/<int:report_id>')
@login_required
def display_report(report_id):
    """Display report."""
    report = Report.query.options(
        joinedload(Report.associated_records)
    ).get_or_404(report_id)
    
    # Check permissions
    if not current_user.can_view_report(report):
        abort(403)
    
    # Get graph files
    graph_images = {}
    graph_base_path = Path('static/graphs')

    for folder in ['phylo_tree', 'dot_plot', 'heat_map', 'stacked_bar', 'gc_line', 'gc_skew', 'nuc_pie']:
        folder_path = graph_base_path / folder
        if folder_path.exists():
            images = [
                str(img.relative_to('static'))
                for img in folder_path.glob(f'*{report_id}*')
            ]
            if images:
                graph_images[folder] = images

    graph_folders = list(graph_images.keys())

    return render_template(
        'report/display.html',
        report=report,
        graph_images=graph_images,
        graph_folders=graph_folders
    )


@main_bp.route('/report/<int:report_id>/download')
@login_required
def download_report(report_id):
    """Download report as PDF."""
    report = Report.query.get_or_404(report_id)
    
    # Check permissions
    if not current_user.can_print_report(report):
        abort(403)
    
    try:
        # Generate PDF
        pdf_buffer = PDFReportGenerator.generate_pdf(report)
        
        return send_file(
            pdf_buffer,
            mimetype='application/pdf',
            as_attachment=True,
            download_name=f'report_{report_id}_{datetime.now().strftime("%Y%m%d")}.pdf'
        )
        
    except Exception as e:
        logger.error(f"Error generating PDF for report {report_id}: {e}")
        flash('Failed to generate PDF report', 'error')
        return redirect(url_for('main.display_report', report_id=report_id))


# Blueprint for admin routes
admin_bp = Blueprint('admin', __name__, url_prefix='/admin')

@admin_bp.route('/')
@login_required
@role_required('admin')
def dashboard():
    """Admin dashboard."""
    users = User.query.all()
    reports_count = Report.query.count()
    records_count = Record.query.count()
    
    stats = {
        'total_users': len(users),
        'total_reports': reports_count,
        'total_records': records_count,
        'users_by_role': {}
    }
    
    for user in users:
        role = user.role
        stats['users_by_role'][role] = stats['users_by_role'].get(role, 0) + 1
    
    return render_template('admin/dashboard.html', users=users, stats=stats)


@admin_bp.route('/user/<int:user_id>', methods=['GET', 'POST'])
@login_required
@role_required('admin')
def edit_user(user_id):
    """Edit user details and permissions."""
    user = User.query.get_or_404(user_id)
    form = UserEditForm(obj=user)
    
    if form.validate_on_submit():
        try:
            user.username = form.username.data
            user.email = form.email.data
            user.role = form.role.data
            
            # Update permissions
            user.view_reports = form.view_reports.data
            user.delete_reports = form.delete_reports.data
            user.print_reports = form.print_reports.data
            user.change_reports = form.change_reports.data
            
            # Update password if provided
            if form.password.data:
                user.password = generate_password_hash(form.password.data)
            
            db.session.commit()
            flash(f'User {user.username} updated successfully', 'success')
            return redirect(url_for('admin.dashboard'))
            
        except SQLAlchemyError as e:
            logger.error(f"Error updating user {user_id}: {e}")
            db.session.rollback()
            flash('Failed to update user', 'error')
    
    return render_template('admin/edit_user.html', form=form, user=user)


# API Blueprint for RESTful endpoints
api_bp = Blueprint('api', __name__, url_prefix='/api/v1')

@api_bp.route('/reports', methods=['GET'])
@login_required
def get_reports():
    """Get user's reports with pagination."""
    page = request.args.get('page', 1, type=int)
    per_page = request.args.get('per_page', 10, type=int)
    
    query = Report.query.filter_by(employee_id=current_user.id)
    pagination = query.paginate(page=page, per_page=per_page, error_out=False)
    
    reports = [{
        'id': report.id,
        'created_at': report.created_at.isoformat(),
        'organisms': report.organisms,
        'nucleotide_ids': report.nucleotide_ids
    } for report in pagination.items]
    
    return jsonify({
        'reports': reports,
        'total': pagination.total,
        'pages': pagination.pages,
        'current_page': page
    })


@api_bp.route('/reports/<int:report_id>', methods=['DELETE'])
@login_required
def delete_report(report_id):
    """Delete a report."""
    report = Report.query.get_or_404(report_id)
    
    # Check permissions
    if not current_user.can_delete_report(report):
        return jsonify({'error': 'Permission denied'}), 403
    
    try:
        # Delete associated files
        for path in [report.phylo_tree, report.bar_chart] + report.dot_line_graph + report.heat_map:
            if path and os.path.exists(path):
                os.remove(path)
        
        db.session.delete(report)
        db.session.commit()
        
        return jsonify({'message': 'Report deleted successfully'})
        
    except Exception as e:
        logger.error(f"Error deleting report {report_id}: {e}")
        db.session.rollback()
        return jsonify({'error': 'Failed to delete report'}), 500


# Celery tasks
@celery.task(bind=True)
def compile_report_task(self, user_id):
    """Async task to compile report."""
    try:
        # Get records
        records = Record.query.filter_by(employee_id=user_id).all()
        
        if len(records) < 2:
            raise ValueError("At least 2 records required for report")
        
        # Extract data
        nucleotide_ids = [r.nucleotide_id for r in records]
        organisms = [r.organism for r in records]
        nucleotides = [r.nucleotides for r in records]
        
        # Create report
        report = Report(
            nucleotide_ids=nucleotide_ids,
            organisms=organisms,
            nucleotides=nucleotides,
            employee_id=user_id
        )
        db.session.add(report)
        db.session.commit()
        
        # Generate graphs
        graph_paths = ReportGenerator.generate_graphs(
            report.id, nucleotides, organisms
        )
        
        # Update report with graph paths
        report.phylo_tree = graph_paths['phylo_tree']
        report.dot_line_graph = graph_paths['dot_plots']
        report.heat_map = graph_paths['heat_maps']
        report.bar_chart = graph_paths['bar_chart']
        report.gc_line_graphs = graph_paths['gc_lines']
        report.gc_skew_graphs = graph_paths['gc_skews']
        report.nuc_pie_charts = graph_paths['nuc_pies']
        
        # Associate records
        report.associated_records.extend(records)
        db.session.commit()
        
        return {
            'status': 'success',
            'report_id': report.id
        }
        
    except Exception as e:
        logger.error(f"Report compilation failed: {e}")
        self.update_state(
            state='FAILURE',
            meta={'error': str(e)}
        )
        raise


# Application factory
def create_app(config_class=Config):
    """Create Flask application."""
    app = Flask(__name__)
    app.config.from_object(config_class)
    
    # Initialize extensions
    db.init_app(app)
    login_manager.init_app(app)
    from models import User
    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))
    csrf.init_app(app)
    cache.init_app(app)
    limiter.init_app(app)
    
    # Register blueprints
    app.register_blueprint(auth_bp)
    app.register_blueprint(main_bp)
    app.register_blueprint(admin_bp)
    app.register_blueprint(api_bp)
    
    # Configure login
    login_manager.login_view = 'auth.login'
    login_manager.login_message = 'Please log in to access this page.'
    
    # Create Celery
    global celery
    celery = make_celery(app)
    app.celery = celery
    
    # Error handlers
    @app.errorhandler(404)
    def not_found_error(error):
        return render_template('errors/404.html'), 404
    
    @app.errorhandler(403)
    def forbidden_error(error):
        return render_template('errors/403.html'), 403
    
    @app.errorhandler(500)
    def internal_error(error):
        db.session.rollback()
        return render_template('errors/500.html'), 500
    
    # Security headers
    @app.after_request
    def set_security_headers(response):
        response.headers['X-Content-Type-Options'] = 'nosniff'
        response.headers['X-Frame-Options'] = 'DENY'
        response.headers['X-XSS-Protection'] = '1; mode=block'
        response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'
        return response
    
    return app


# Initialize extensions (imported by models.py)
# (instances created above)

if __name__ == "__main__":
    app = create_app()
    app.run(debug=False)
