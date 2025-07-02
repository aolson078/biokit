# BioKit – Bioinformatic Toolkit

BioKit is a Flask web application offering researchers a simple interface for processing nucleotide sequences and generating reports. It was developed by Alex Olson, Sola Yun and Will Hutcheon as part of a W.A.S. software project.

## Features
- Upload sequences and store results
- Automatic phylogenetic tree generation
- Role-based access (employee, manager, admin)
- Asynchronous job queue powered by Celery
- JSON API endpoints for automation
- Docker configuration for easy deployment
- GC content line graphs for individual sequences
- GC content heatmaps comparing multiple sequences
- GC skew plots highlighting replication biases
- Nucleotide composition pie charts for each sequence

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

![Use case](https://github.com/aolson078/biokit/assets/69769089/97b19bdb-c369-4d3a-a883-24ca40f4b959)
