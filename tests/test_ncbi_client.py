from datetime import datetime, timedelta, timezone
from urllib.error import HTTPError

import pytest

from sequence_workspace.ncbi import NcbiClient, NcbiError
import sequence_workspace.ncbi as ncbi_module


class MemoryCache:
    def __init__(self):
        self.values = {}

    def get(self, key):
        return self.values.get(key)

    def set(self, key, value, timeout=None):
        self.values[key] = value


def test_search_maps_esearch_and_esummary(monkeypatch):
    client = NcbiClient(MemoryCache(), email="test@example.com")
    responses = iter([
        {"esearchresult": {"idlist": ["1"], "count": "1"}},
        {"result": {"uids": ["1"], "1": {"accessionversion": "NM_1.2", "title": "Example", "organism": "Testus", "slen": 42, "moltype": "mRNA"}}},
    ])
    monkeypatch.setattr(client, "_json_request", lambda *_args, **_kwargs: next(responses))
    result = client.search("example")
    assert result["results"] == [{
        "uid": "1",
        "accession_version": "NM_1.2",
        "title": "Example",
        "organism": "Testus",
        "length": 42,
        "molecule_type": "rna",
        "updated_at": None,
    }]
    assert result["cache"]["stale"] is False


def test_retryable_failure_serves_labeled_stale_cache(monkeypatch):
    cache = MemoryCache()
    client = NcbiClient(cache, email="test@example.com")
    key = "sequence:ncbi:search:v1:example::1:10"
    now = datetime.now(timezone.utc)
    cache.values[key] = {
        "value": {"results": [{"accession_version": "OLD.1"}]},
        "cached_at": (now - timedelta(hours=2)).isoformat(),
        "fresh_until": (now - timedelta(hours=1)).isoformat(),
        "retain_until": (now + timedelta(hours=22)).isoformat(),
    }
    monkeypatch.setattr(client, "_search_live", lambda *_args: (_ for _ in ()).throw(NcbiError("ncbi_unavailable", "Unavailable", retryable=True)))
    result = client.search("example")
    assert result["cache"]["stale"] is True
    assert result["cache"]["upstream_error_code"] == "ncbi_unavailable"


def test_missing_configuration_disables_ncbi_only():
    with pytest.raises(NcbiError) as error:
        NcbiClient(MemoryCache(), email="")
    assert error.value.code == "ncbi_not_configured"


def test_entrez_request_retries_once_for_429(monkeypatch):
    client = NcbiClient(MemoryCache(), email="test@example.com")
    attempts = []

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self):
            return b'{"ok": true}'

    def fake_esearch(**_query):
        attempts.append(1)
        if len(attempts) == 1:
            raise HTTPError("https://example.invalid", 429, "rate limited", {}, None)
        return Response()

    monkeypatch.setattr(ncbi_module.Entrez, "esearch", fake_esearch)
    monkeypatch.setattr(ncbi_module.time, "sleep", lambda _seconds: None)
    assert client._json_request("esearch.fcgi", {"db": "nucleotide"}) == {"ok": True}
    assert len(attempts) == 2


def test_cache_failure_is_typed_and_retryable():
    class BrokenCache:
        def get(self, _key):
            raise ConnectionError("unavailable")

    client = NcbiClient(BrokenCache(), email="test@example.com")
    with pytest.raises(NcbiError) as error:
        client.search("example")
    assert error.value.code == "cache_unavailable"
    assert error.value.retryable is True
