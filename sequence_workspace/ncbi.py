"""Structured NCBI E-utilities client with bounded retries and stale caching."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json
import random
import time
from typing import Any, Callable
from urllib.error import HTTPError, URLError

from Bio import Entrez

from .domain import MAX_SEQUENCE_LENGTH, SequenceDomainError, normalize_sequence

Entrez.max_tries = 1
SEARCH_FRESH_SECONDS = 3600
SEARCH_RETAIN_SECONDS = 86_400
RECORD_FRESH_SECONDS = 86_400
RECORD_RETAIN_SECONDS = 604_800


class NcbiError(SequenceDomainError):
    pass


class NcbiClient:
    def __init__(self, cache: Any, *, email: str, api_key: str | None = None, tool: str = "biokit"):
        if not email or email == "your-email@example.com":
            raise NcbiError("ncbi_not_configured", "NCBI search requires a configured contact email.")
        self.cache = cache
        self.email = email
        self.api_key = api_key
        self.tool = tool

    def search(self, query: str, organism: str | None = None, page: int = 1, page_size: int = 10) -> dict[str, Any]:
        normalized_query = " ".join(query.split())
        normalized_organism = " ".join((organism or "").split())
        if not normalized_query or len(normalized_query) > 200:
            raise NcbiError("invalid_query", "Search query must contain 1 to 200 characters.")
        if len(normalized_organism) > 100:
            raise NcbiError("invalid_organism", "Organism filter cannot exceed 100 characters.")
        if page < 1 or page > 5 or page_size != 10:
            raise NcbiError("invalid_page", "Search supports pages 1 through 5 with 10 results per page.")

        key = f"sequence:ncbi:search:v1:{normalized_query.lower()}:{normalized_organism.lower()}:{page}:{page_size}"
        return self._cached_request(
            key,
            SEARCH_FRESH_SECONDS,
            SEARCH_RETAIN_SECONDS,
            lambda: self._search_live(normalized_query, normalized_organism, page, page_size),
        )

    def fetch_record(self, accession_version: str) -> dict[str, Any]:
        accession = accession_version.strip().upper()
        if not accession or len(accession) > 64 or not all(character.isalnum() or character in "._" for character in accession):
            raise NcbiError("invalid_accession", "Invalid accession.version value.")
        key = f"sequence:ncbi:record:v1:{accession}"
        return self._cached_request(
            key,
            RECORD_FRESH_SECONDS,
            RECORD_RETAIN_SECONDS,
            lambda: self._fetch_live(accession),
        )

    def _search_live(self, query: str, organism: str, page: int, page_size: int) -> dict[str, Any]:
        term = f"({query})"
        if organism:
            term += f" AND ({organism}[Organism])"
        esearch = self._json_request("esearch.fcgi", {
            "db": "nucleotide",
            "term": term,
            "retstart": (page - 1) * page_size,
            "retmax": page_size,
            "retmode": "json",
        })
        search_result = esearch.get("esearchresult", {})
        uids = search_result.get("idlist", [])
        summaries = self._summary(uids)
        return {
            "query": query,
            "organism": organism or None,
            "page": page,
            "page_size": page_size,
            "total": min(int(search_result.get("count", 0)), 50),
            "results": summaries,
        }

    def _fetch_live(self, accession: str) -> dict[str, Any]:
        esearch = self._json_request("esearch.fcgi", {
            "db": "nucleotide",
            "term": f"{accession}[Accession]",
            "retmax": 1,
            "retmode": "json",
        })
        uids = esearch.get("esearchresult", {}).get("idlist", [])
        if not uids:
            raise NcbiError("record_not_found", "No NCBI nucleotide record matched that accession.")
        summaries = self._summary(uids)
        summary = summaries[0] if summaries else {}
        if summary.get("length", 0) > MAX_SEQUENCE_LENGTH:
            raise NcbiError("sequence_too_long", f"NCBI record exceeds the {MAX_SEQUENCE_LENGTH:,} nucleotide limit.")
        fasta = self._text_request("efetch.fcgi", {
            "db": "nucleotide",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
        })
        lines = [line.strip() for line in fasta.splitlines() if line.strip()]
        sequence = "".join(line for line in lines if not line.startswith(">"))
        molecule = summary.get("molecule_type", "unknown")
        normalized = normalize_sequence(
            sequence,
            molecule_hint=molecule,
            authoritative_molecule=molecule in {"dna", "rna"},
        )
        return {
            "status": "ready",
            "client_id": summary.get("accession_version") or accession,
            "source": "ncbi",
            "source_accession": summary.get("accession_version") or accession,
            "source_title": summary.get("title") or (lines[0][1:] if lines and lines[0].startswith(">") else accession),
            "organism": summary.get("organism"),
            "source_updated_at": summary.get("updated_at"),
            "molecule_type": normalized.molecule_type,
            "sequence_alphabet": normalized.sequence_alphabet,
            "sequence": normalized.sequence,
            "length": len(normalized.sequence),
        }

    def _summary(self, uids: list[str]) -> list[dict[str, Any]]:
        if not uids:
            return []
        payload = self._json_request("esummary.fcgi", {
            "db": "nucleotide",
            "id": ",".join(uids),
            "retmode": "json",
        })
        result = payload.get("result", {})
        summaries: list[dict[str, Any]] = []
        for uid in result.get("uids", uids):
            item = result.get(str(uid), {})
            title = item.get("title") or item.get("caption") or str(uid)
            summaries.append({
                "uid": str(uid),
                "accession_version": item.get("accessionversion") or item.get("caption") or str(uid),
                "title": title,
                "organism": item.get("organism"),
                "length": int(item.get("slen") or 0),
                "molecule_type": _molecule_type(item),
                "updated_at": item.get("updatedate"),
            })
        return summaries

    def _json_request(self, endpoint: str, parameters: dict[str, Any]) -> dict[str, Any]:
        return json.loads(self._request(endpoint, parameters, lambda response: response.read().decode("utf-8")))

    def _text_request(self, endpoint: str, parameters: dict[str, Any]) -> str:
        return self._request(endpoint, parameters, lambda response: response.read().decode("utf-8"))

    def _request(self, endpoint: str, parameters: dict[str, Any], reader: Callable[[Any], Any]) -> Any:
        query = dict(parameters)
        query.update({"email": self.email, "tool": self.tool})
        if self.api_key:
            query["api_key"] = self.api_key
        operation = getattr(Entrez, endpoint.removesuffix(".fcgi"))
        last_error: Exception | None = None
        for attempt in range(2):
            try:
                with operation(**query) as response:
                    return reader(response)
            except HTTPError as error:
                if error.code not in {429, 500, 502, 503, 504}:
                    if error.code == 404:
                        raise NcbiError("record_not_found", "NCBI record was not found.") from error
                    raise NcbiError("ncbi_request_rejected", "NCBI rejected the request.") from error
                last_error = error
            except (TimeoutError, URLError, ConnectionError) as error:
                last_error = error
            if attempt == 0:
                time.sleep(random.uniform(0.25, 0.75))
        raise NcbiError("ncbi_unavailable", "NCBI is temporarily unavailable.", retryable=True) from last_error

    def _cached_request(
        self,
        key: str,
        fresh_seconds: int,
        retain_seconds: int,
        loader: Callable[[], dict[str, Any]],
    ) -> dict[str, Any]:
        now = datetime.now(timezone.utc)
        try:
            cached = self.cache.get(key)
        except Exception as error:
            raise NcbiError(
                "cache_unavailable",
                "NCBI search caching is temporarily unavailable.",
                retryable=True,
            ) from error
        if cached and _parse_time(cached["fresh_until"]) > now:
            return _with_cache_metadata(cached["value"], cached, stale=False)
        try:
            value = loader()
        except NcbiError as error:
            if cached and error.retryable:
                return _with_cache_metadata(cached["value"], cached, stale=True, upstream_error_code=error.code)
            raise
        entry = {
            "value": value,
            "cached_at": now.isoformat(),
            "fresh_until": (now + timedelta(seconds=fresh_seconds)).isoformat(),
            "retain_until": (now + timedelta(seconds=retain_seconds)).isoformat(),
        }
        try:
            self.cache.set(key, entry, timeout=retain_seconds)
        except Exception as error:
            raise NcbiError(
                "cache_unavailable",
                "NCBI search caching is temporarily unavailable.",
                retryable=True,
            ) from error
        return _with_cache_metadata(value, entry, stale=False)


def _molecule_type(item: dict[str, Any]) -> str:
    value = " ".join(str(item.get(key, "")) for key in ("moltype", "biomol", "subtype")).lower()
    if "rna" in value:
        return "rna"
    if "dna" in value or "genomic" in value:
        return "dna"
    return "unknown"


def _parse_time(value: str) -> datetime:
    return datetime.fromisoformat(value)


def _with_cache_metadata(
    value: dict[str, Any],
    entry: dict[str, Any],
    *,
    stale: bool,
    upstream_error_code: str | None = None,
) -> dict[str, Any]:
    result = dict(value)
    result["cache"] = {
        "stale": stale,
        "cached_at": entry["cached_at"],
        "expires_at": entry["fresh_until"],
        "upstream_error_code": upstream_error_code,
    }
    return result
