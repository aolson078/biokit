# BioKit – Bioinformatic Toolkit

BioKit is a Flask web application offering researchers a simple interface for processing nucleotide sequences and generating reports. It was developed by Alex Olson, Sola Yun and Will Hutcheon as part of a W.A.S. software project.

## Features
- Guided NCBI/FASTA nucleotide discovery and analysis workspace
- Upload sequences and store results
- Automatic phylogenetic tree generation
- Role-based access (employee, manager, admin)
- Asynchronous job queue powered by Celery
- JSON API endpoints for automation
- Docker configuration for easy deployment
- GC content line graphs for individual sequences
- GC content heatmaps comparing multiple sequences
- GC skew plots highlighting replication biases
- Cumulative GC skew plots for detecting replication origins
- Nucleotide composition pie charts for each sequence
- Sequence entropy line graphs showing local complexity
- K-mer frequency histograms for exploring motif patterns
- Automatic inclusion of entropy and k-mer visualizations in compiled reports

## Quick start
1. Create a virtual environment and install the dependencies
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```
2. Run the application
   ```bash
   python manage.py
   ```
3. Visit `http://127.0.0.1:5000` and log in using the default credentials:
   - **email**: `1@gmail.com`
   - **password**: `11111`

### Background tasks
Make sure a Redis server is running (``redis-server``) and start a Celery worker to handle long running jobs:
```bash
celery -A tasks.celery worker --loglevel=info
```
Use `/compile_report_async` to queue a report and check `/task_status/<task_id>` for progress.

### API endpoints
- `POST /api/create_record` – add a sequence record
- `POST /compile_report_async` – generate a report asynchronously
- `GET /api/report/<id>` – download a report

### Running tests
Execute the unit tests with `pytest`:
```bash
pytest
```

### Docker
Run the complete stack with Docker:
```bash
docker-compose up --build
```
This launches the web server, Celery worker and Redis broker.

The sequence workspace is shipped behind two flags. Keep both off until the
database fingerprint, backup/restore evidence, and migrations are approved.
After migration, enable routes first with
`SEQUENCE_WORKSPACE_ENABLED=true`; enable the home-page cutover separately
with `SEQUENCE_WORKSPACE_CTA_ENABLED=true`.

For local Vite development:

```bash
docker compose -f docker-compose.yml -f docker-compose.dev.yml up --build
```

Flask remains the page and API origin at `http://localhost:5000`; Vite serves
development assets at port 5173. Configure `NCBI_EMAIL` to enable public NCBI
search. Manual DNA/RNA and FASTA import remains available without it.

The authenticated routes are `/sequence-wizard` for the novice flow and
`/sequence-workbench` for the advanced workspace.

![Use case](https://github.com/aolson078/biokit/assets/69769089/97b19bdb-c369-4d3a-a883-24ca40f4b959)
