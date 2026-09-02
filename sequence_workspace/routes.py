"""Authenticated HTTP and job boundaries for the sequence workspace."""

from __future__ import annotations

from datetime import datetime, timezone
from functools import wraps
import hashlib
import json
from pathlib import Path
from typing import Any, Callable
from uuid import uuid4

from flask import Blueprint, abort, current_app, jsonify, render_template, request, send_file
from flask_login import current_user, login_required
from flask_wtf.csrf import generate_csrf
from sqlalchemy.exc import IntegrityError

from .catalog import ANALYSIS_CATALOG, public_catalog
from .domain import (
    MAX_IMPORT_BYTES,
    MAX_SEQUENCE_LENGTH,
    SequenceDomainError,
    normalize_record,
    parse_import,
    run_analysis,
    sequence_summary,
)
from .ncbi import NcbiClient, NcbiError

JOB_TTL_SECONDS = 86_400
RESULT_TTL_SECONDS = 86_400


def api_error(code: str, message: str, status: int, *, retryable: bool = False):
    return jsonify({"error": {"code": code, "message": message, "retryable": retryable}}), status


def workspace_required(function: Callable[..., Any]):
    @wraps(function)
    def decorated(*args: Any, **kwargs: Any):
        if not current_app.config.get("SEQUENCE_WORKSPACE_ENABLED", False):
            return api_error("workspace_disabled", "The sequence workspace is not enabled.", 404)
        return function(*args, **kwargs)
    return decorated


def create_sequence_blueprint(db: Any, cache: Any) -> Blueprint:
    blueprint = Blueprint("sequence_workspace", __name__)

    @blueprint.errorhandler(SequenceDomainError)
    def handle_sequence_error(error: SequenceDomainError):
        return api_error(
            error.code,
            error.message,
            503 if error.retryable else 400,
            retryable=error.retryable,
        )

    @blueprint.route("/sequence-wizard")
    @login_required
    def wizard_page():
        return _workspace_page("wizard")

    @blueprint.route("/sequence-workbench")
    @login_required
    def workbench_page():
        return _workspace_page("workbench")

    @blueprint.route("/api/v1/nucleotides/search", methods=["POST"])
    @login_required
    @workspace_required
    def search_nucleotides():
        payload = request.get_json(silent=True) or {}
        try:
            client = _ncbi_client(cache)
            return jsonify(client.search(
                str(payload.get("query", "")),
                str(payload.get("organism", "")) or None,
                int(payload.get("page", 1)),
                int(payload.get("page_size", 10)),
            ))
        except (ValueError, TypeError):
            return api_error("invalid_page", "Page values must be integers.", 400)
        except NcbiError as error:
            return api_error(error.code, error.message, 503 if error.retryable or error.code == "ncbi_not_configured" else 400, retryable=error.retryable)

    @blueprint.route("/api/v1/nucleotides/ncbi/<accession_version>")
    @login_required
    @workspace_required
    def fetch_nucleotide(accession_version: str):
        try:
            return jsonify(_ncbi_client(cache).fetch_record(accession_version))
        except NcbiError as error:
            status = 503 if error.retryable or error.code == "ncbi_not_configured" else 404 if error.code == "record_not_found" else 400
            return api_error(error.code, error.message, status, retryable=error.retryable)

    @blueprint.route("/api/v1/nucleotides/import", methods=["POST"])
    @login_required
    @workspace_required
    def import_nucleotides():
        payload = request.get_json(silent=True) or {}
        try:
            records = parse_import(
                payload.get("text", ""),
                molecule_hint=payload.get("molecule_type"),
                source_name=payload.get("source_name"),
            )
            return jsonify({"records": records, "limits": {"records": 25, "bytes": MAX_IMPORT_BYTES, "sequence_length": MAX_SEQUENCE_LENGTH}})
        except SequenceDomainError as error:
            return api_error(error.code, error.message, 400, retryable=error.retryable)

    @blueprint.route("/api/v1/analyses")
    @login_required
    @workspace_required
    def analyses_catalog():
        return jsonify({"analyses": public_catalog()})

    @blueprint.route("/api/v1/analyses/summary", methods=["POST"])
    @login_required
    @workspace_required
    def summarize_sequence():
        payload = request.get_json(silent=True) or {}
        try:
            return jsonify({"analysis_id": "summary", "result": sequence_summary(payload.get("record", {}))})
        except SequenceDomainError as error:
            return api_error(error.code, error.message, 503 if error.retryable else 400, retryable=error.retryable)

    @blueprint.route("/api/v1/analyses/<analysis_id>", methods=["POST"])
    @login_required
    @workspace_required
    def execute_analysis(analysis_id: str):
        payload = request.get_json(silent=True) or {}
        records = payload.get("records") or ([payload["record"]] if payload.get("record") else [])
        parameters = payload.get("parameters") or {}
        try:
            if not isinstance(parameters, dict):
                raise SequenceDomainError("invalid_parameters", "Analysis parameters must be an object.")
            normalized_records = [normalize_record(record) for record in records]
            entry = ANALYSIS_CATALOG.get(analysis_id)
            if entry is None:
                raise SequenceDomainError("unknown_analysis", "Unknown analysis.")
            cache_key = _result_cache_key(current_user.id, analysis_id, normalized_records, parameters)
            cached = _cache_get(cache, cache_key)
            if cached is not None:
                return jsonify({"analysis_id": analysis_id, "result": cached, "cached": True})
            execution = entry["execution"]
            should_queue = execution == "async" or (
                execution == "adaptive" and any(len(record["sequence"]) > 20_000 for record in normalized_records)
            )
            if should_queue:
                job_id = _queue_job(cache, analysis_id, normalized_records, parameters)
                return jsonify({"job_id": job_id, "status": "queued"}), 202
            result = run_analysis(analysis_id, normalized_records, parameters)
            _cache_set(cache, cache_key, result, timeout=RESULT_TTL_SECONDS)
            return jsonify({"analysis_id": analysis_id, "result": result, "cached": False})
        except SequenceDomainError as error:
            return api_error(error.code, error.message, 400)

    @blueprint.route("/api/v1/records", methods=["GET", "POST"])
    @login_required
    @workspace_required
    def records_collection():
        from models import Record

        if request.method == "GET":
            records = Record.query.filter_by(employee_id=current_user.id).order_by(Record.id.desc()).all()
            return jsonify({"records": [_record_json(record) for record in records]})

        payload = request.get_json(silent=True) or {}
        source = payload.get("source")
        user_label = _normalize_label(payload.get("user_label"))
        try:
            if source == "ncbi":
                accession = str(payload.get("source_accession", "")).strip().upper()
                existing = Record.query.filter_by(employee_id=current_user.id, source_accession=accession).first()
                if existing:
                    return jsonify({"created": False, "record": _record_json(existing)}), 200
                candidate = _ncbi_client(cache).fetch_record(accession)
                source_accession = candidate["source_accession"]
            elif source == "manual":
                candidate = normalize_record({
                    **payload,
                    "source": "manual",
                    "client_id": payload.get("client_id") or f"local:{uuid4()}",
                })
                source_accession = None
            else:
                raise SequenceDomainError("invalid_source", "Source must be ncbi or manual.")

            record = Record(
                nucleotide_id=source_accession or f"manual:{uuid4()}",
                source=source,
                source_accession=source_accession,
                source_title=candidate.get("source_title"),
                user_label=user_label,
                source_retrieved_at=datetime.now(timezone.utc),
                source_updated_at=_parse_optional_datetime(candidate.get("source_updated_at")),
                molecule_type=candidate["molecule_type"],
                sequence_alphabet=candidate["sequence_alphabet"],
                sequence_length=len(candidate["sequence"]),
                organism=candidate.get("organism"),
                gene_info=candidate.get("source_title"),
                nucleotides=candidate["sequence"],
                employee_id=current_user.id,
            )
            db.session.add(record)
            db.session.commit()
            return jsonify({"created": True, "record": _record_json(record)}), 201
        except IntegrityError:
            db.session.rollback()
            if source == "ncbi" and payload.get("source_accession"):
                existing = Record.query.filter_by(
                    employee_id=current_user.id,
                    source_accession=str(payload["source_accession"]).strip().upper(),
                ).first()
                if existing:
                    return jsonify({"created": False, "record": _record_json(existing)}), 200
            return api_error("record_conflict", "The record could not be saved because it conflicts with existing data.", 409)
        except (SequenceDomainError, NcbiError) as error:
            db.session.rollback()
            return api_error(error.code, error.message, 400 if not error.retryable else 503, retryable=error.retryable)

    @blueprint.route("/api/v1/records/<int:record_id>", methods=["GET", "PATCH"])
    @login_required
    @workspace_required
    def record_item(record_id: int):
        from models import Record

        record = Record.query.filter_by(id=record_id, employee_id=current_user.id).first()
        if record is None:
            return api_error("record_not_found", "Record was not found.", 404)
        if request.method == "GET":
            return jsonify({"record": _record_json(record)})
        payload = request.get_json(silent=True) or {}
        if set(payload) != {"user_label"}:
            return api_error("invalid_patch", "Only user_label can be changed.", 400)
        try:
            record.user_label = _normalize_label(payload.get("user_label"))
        except SequenceDomainError as error:
            return api_error(error.code, error.message, 400)
        db.session.commit()
        return jsonify({"record": _record_json(record)})

    @blueprint.route("/api/v1/comparisons", methods=["POST"])
    @login_required
    @workspace_required
    def create_comparison():
        from models import Record

        payload = request.get_json(silent=True) or {}
        analysis_id = payload.get("analysis_id", "stacked_composition")
        record_ids = payload.get("record_ids") or []
        if (
            not isinstance(record_ids, list)
            or len(record_ids) < 2
            or len(record_ids) > 25
            or not all(type(value) is int for value in record_ids)
            or len(record_ids) != len(set(record_ids))
        ):
            return api_error("invalid_record_ids", "Select between 2 and 25 saved records.", 400)
        entry = ANALYSIS_CATALOG.get(analysis_id)
        if entry is None or entry["max_records"] < 2:
            return api_error("invalid_comparison", "Select a multi-record analysis.", 400)
        owned_count = Record.query.filter(Record.employee_id == current_user.id, Record.id.in_(record_ids)).count()
        if owned_count != len(set(record_ids)):
            return api_error("record_not_found", "One or more selected records are unavailable.", 404)
        job_id = _queue_job(cache, analysis_id, [], payload.get("parameters") or {}, record_ids=record_ids)
        return jsonify({"job_id": job_id, "status": "queued"}), 202

    @blueprint.route("/api/v1/jobs/<job_id>")
    @login_required
    @workspace_required
    def job_status(job_id: str):
        envelope = _cache_get(cache, _job_key(job_id))
        if not envelope or envelope.get("user_id") != current_user.id:
            return api_error("job_not_found", "Job was not found or has expired.", 404)
        return jsonify({key: value for key, value in envelope.items() if key != "user_id"})

    @blueprint.route("/api/v1/jobs/<job_id>/artifacts/<artifact_id>")
    @login_required
    @workspace_required
    def job_artifact(job_id: str, artifact_id: str):
        envelope = _cache_get(cache, _job_key(job_id))
        if not envelope or envelope.get("user_id") != current_user.id:
            return api_error("job_not_found", "Job was not found or has expired.", 404)
        artifact = (envelope.get("artifacts") or {}).get(artifact_id)
        if not artifact:
            return api_error("artifact_not_found", "Artifact was not found.", 404)
        root = (Path(current_app.instance_path) / "analysis_artifacts" / job_id).resolve()
        path = (root / artifact["name"]).resolve()
        if root not in path.parents or not path.is_file():
            abort(404)
        if artifact.get("content_type") not in {"image/png", "application/json", "text/plain"}:
            abort(404)
        return send_file(path, mimetype=artifact["content_type"], download_name=artifact["name"])

    return blueprint


def register_sequence_tasks(celery: Any, cache: Any) -> dict[str, Any]:
    task_name = "biokit.sequence_analysis"
    existing = celery.tasks.get(task_name)
    if existing:
        return {"analysis": existing}

    @celery.task(bind=True, name=task_name)
    def sequence_analysis_task(self: Any, job_id: str, user_id: int, analysis_id: str, records: list[dict[str, Any]], parameters: dict[str, Any], record_ids: list[int] | None = None):
        envelope = _cache_get(cache, _job_key(job_id)) or {}
        envelope.update({"user_id": user_id, "job_id": job_id, "analysis_id": analysis_id, "status": "running"})
        _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
        try:
            if record_ids:
                from models import Record

                owned = Record.query.filter(Record.employee_id == user_id, Record.id.in_(record_ids)).all()
                if len(owned) != len(set(record_ids)):
                    raise SequenceDomainError("record_not_found", "One or more selected records are unavailable.")
                by_id = {record.id: record for record in owned}
                records = [_record_candidate(by_id[record_id]) for record_id in record_ids]
            cache_key = _result_cache_key(user_id, analysis_id, records, parameters)
            cached = _cache_get(cache, cache_key)
            if cached is not None:
                envelope.update({"status": "completed", "result": cached, "cached": True, "completed_at": _now_iso()})
                _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
                return cached
            result = run_analysis(analysis_id, records, parameters)
            _cache_set(cache, cache_key, result, timeout=RESULT_TTL_SECONDS)
            envelope.update({"status": "completed", "result": result, "cached": False, "completed_at": _now_iso()})
            _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
            return result
        except SequenceDomainError as error:
            envelope.update({"status": "failed", "error": {"code": error.code, "message": error.message, "retryable": False}, "completed_at": _now_iso()})
            _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
            raise

    return {"analysis": sequence_analysis_task}


def _workspace_page(mode: str):
    if not current_app.config.get("SEQUENCE_WORKSPACE_ENABLED", False):
        abort(404)
    dev_server = current_app.config.get("VITE_DEV_SERVER_URL")
    scripts: list[str] = []
    styles: list[str] = []
    if not dev_server:
        manifest_path = Path(current_app.static_folder) / "sequence-workspace" / "manifest.json"
        if not manifest_path.is_file():
            return render_template("sequence_workspace_unavailable.html"), 503
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        entry = manifest.get("src/main.tsx") or manifest.get("index.html")
        if not entry:
            return render_template("sequence_workspace_unavailable.html"), 503
        scripts.append(f"/static/sequence-workspace/{entry['file']}")
        styles.extend(f"/static/sequence-workspace/{style}" for style in entry.get("css", []))
    return render_template(
        "sequence_workspace.html",
        mode=mode,
        csrf_token=generate_csrf(),
        workspace_user_key=str(current_user.get_id()),
        dev_server=dev_server.rstrip("/") if dev_server else None,
        scripts=scripts,
        styles=styles,
    )


def _ncbi_client(cache: Any) -> NcbiClient:
    return NcbiClient(
        cache,
        email=current_app.config.get("NCBI_EMAIL", ""),
        api_key=current_app.config.get("NCBI_API_KEY"),
    )


def _queue_job(cache: Any, analysis_id: str, records: list[dict[str, Any]], parameters: dict[str, Any], *, record_ids: list[int] | None = None) -> str:
    job_id = str(uuid4())
    envelope = {
        "user_id": current_user.id,
        "job_id": job_id,
        "analysis_id": analysis_id,
        "status": "queued",
        "created_at": _now_iso(),
        "artifacts": {},
    }
    _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
    task = current_app.extensions["sequence_tasks"]["analysis"]
    try:
        task.apply_async(args=[job_id, current_user.id, analysis_id, records, parameters, record_ids], task_id=job_id)
    except Exception as error:
        envelope.update({
            "status": "failed",
            "error": {"code": "job_service_unavailable", "message": "Background analysis is temporarily unavailable.", "retryable": True},
            "completed_at": _now_iso(),
        })
        _cache_set(cache, _job_key(job_id), envelope, timeout=JOB_TTL_SECONDS)
        raise SequenceDomainError(
            "job_service_unavailable",
            "Background analysis is temporarily unavailable.",
            retryable=True,
        ) from error
    return job_id


def _record_json(record: Any, *, include_sequence: bool = True) -> dict[str, Any]:
    value = {
        "client_id": f"saved:{record.id}",
        "id": record.id,
        "nucleotide_id": record.nucleotide_id,
        "source": record.source,
        "source_accession": record.source_accession,
        "source_title": record.source_title,
        "user_label": record.user_label,
        "organism": record.organism,
        "molecule_type": record.molecule_type,
        "sequence_alphabet": record.sequence_alphabet,
        "length": record.sequence_length or (len(record.nucleotides) if record.nucleotides else 0),
        "source_retrieved_at": record.source_retrieved_at.isoformat() if record.source_retrieved_at else None,
        "source_updated_at": record.source_updated_at.isoformat() if record.source_updated_at else None,
        "created_at": record.created_at.isoformat() if record.created_at else None,
    }
    if include_sequence:
        value["sequence"] = record.nucleotides
    return value


def _cache_get(cache: Any, key: str) -> Any:
    try:
        return cache.get(key)
    except Exception as error:
        raise SequenceDomainError(
            "cache_unavailable",
            "The cache or background-job service is temporarily unavailable.",
            retryable=True,
        ) from error


def _cache_set(cache: Any, key: str, value: Any, *, timeout: int) -> None:
    try:
        cache.set(key, value, timeout=timeout)
    except Exception as error:
        raise SequenceDomainError(
            "cache_unavailable",
            "The cache or background-job service is temporarily unavailable.",
            retryable=True,
        ) from error


def _record_candidate(record: Any) -> dict[str, Any]:
    return {
        "client_id": f"saved:{record.id}",
        "source": record.source,
        "source_accession": record.source_accession,
        "source_title": record.source_title,
        "organism": record.organism,
        "molecule_type": record.molecule_type,
        "sequence_alphabet": record.sequence_alphabet,
        "sequence": record.nucleotides,
    }


def _normalize_label(value: Any) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise SequenceDomainError("invalid_user_label", "User label must be text.")
    label = value.strip()
    if not label:
        return None
    if len(label) > 120:
        raise SequenceDomainError("invalid_user_label", "User label cannot exceed 120 characters.")
    return label


def _parse_optional_datetime(value: Any):
    if not value or not isinstance(value, str):
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def _result_cache_key(user_id: int, analysis_id: str, records: list[dict[str, Any]], parameters: dict[str, Any]) -> str:
    entry = ANALYSIS_CATALOG[analysis_id]
    sequence_hashes = [hashlib.sha256(record["sequence"].encode("ascii")).hexdigest() for record in records]
    material = json.dumps({
        "schema": 1,
        "user_id": user_id,
        "analysis_id": analysis_id,
        "algorithm_version": entry["algorithm_version"],
        "molecules": [(record["molecule_type"], record["sequence_alphabet"]) for record in records],
        "sequence_hashes": sequence_hashes,
        "parameters": parameters,
    }, sort_keys=True, separators=(",", ":"))
    return f"sequence:result:{hashlib.sha256(material.encode('utf-8')).hexdigest()}"


def _job_key(job_id: str) -> str:
    return f"sequence:job:{job_id}"


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
