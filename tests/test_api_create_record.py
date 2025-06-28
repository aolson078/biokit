import sys
import types
sys.path.append('.')
import bio_algos.utilities as utilities

# Provide minimal helpers before importing routes
if not hasattr(utilities, 'GC_content'):
    utilities.GC_content = lambda seq: (seq.count('G') + seq.count('C')) / len(seq)
if not hasattr(utilities, 'fetch_records'):
    utilities.fetch_records = lambda query: []

import models
if not hasattr(models, 'is_manager'):
    def is_manager(f):
        return f
    models.is_manager = is_manager
if not hasattr(models, 'is_admin'):
    def is_admin(f):
        return f
    models.is_admin = is_admin

# Create lightweight flask_bio_app module with required members
flask_module = types.ModuleType('flask_bio_app')
from marshmallow import Schema, fields, validate, ValidationError
from flask import request, jsonify
def validate_json(schema_class):
    def decorator(f):
        def wrapper(*args, **kwargs):
            schema = schema_class()
            try:
                data = schema.load(request.get_json())
                request.validated_data = data
            except ValidationError as err:
                return jsonify({'errors': err.messages}), 400
            return f(*args, **kwargs)
        return wrapper
    return decorator

class CreateRecordSchema(Schema):
    nucleotide_id = fields.Str(required=True, validate=validate.Length(min=1, max=20))
    organism = fields.Str(required=True, validate=validate.Length(min=1, max=80))
    gene_info = fields.Str(required=True, validate=validate.Length(min=1, max=100))
    nucleotides = fields.Str(required=True, validate=validate.Length(min=1))

flask_module.validate_json = validate_json
flask_module.CreateRecordSchema = CreateRecordSchema
sys.modules['flask_bio_app'] = flask_module

import pytest
import routes

class FakeUser:
    def __init__(self, user_id=1):
        self.is_authenticated = True
        self.id = user_id

    def get_id(self):
        return str(self.id)

@pytest.fixture
def client(monkeypatch):
    app = routes.app
    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///:memory:'
    with app.app_context():
        routes.db.create_all()
    # Provide minimal implementations for missing utility functions
    if not hasattr(routes, 'dna_to_rna'):
        routes.dna_to_rna = lambda seq: seq.replace('T', 'U')
    if not hasattr(routes, 'transcribe_dna'):
        routes.transcribe_dna = lambda seq: seq
    # Stub complex analysis functions to avoid heavy dependencies
    routes.select_target_sequence = lambda rna_seq: ('A'*21, 0.5)
    routes.create_rna_strands = lambda target: ('A'*21, 'U'*21)
    routes.calculate_similarity = lambda a,b: 1.0
    routes.calculate_molecular_weight = lambda s: 1.0
    routes.calculate_melting_temp = lambda s: 1.0
    routes.design_siRNA = lambda s: ['A'*21]
    routes.predict_efficiency = lambda s: 1.0
    monkeypatch.setattr('flask_login.utils._get_user', lambda: FakeUser())
    with app.test_client() as client:
        yield client

def test_create_record_valid(client):
    payload = {
        'nucleotide_id': 'ABC123',
        'organism': 'Test org',
        'gene_info': 'Some gene',
        'nucleotides': 'ATGCGTAC'
    }
    resp = client.post('/api/create_record', json=payload)
    assert resp.status_code == 201
    data = resp.get_json()
    assert 'record_id' in data

def test_create_record_invalid(client):
    payload = {
        'organism': 'Test org',
        'gene_info': 'Some gene',
        'nucleotides': 'ATGCGTAC'
    }
    resp = client.post('/api/create_record', json=payload)
    assert resp.status_code == 400
    data = resp.get_json()
    assert 'errors' in data
