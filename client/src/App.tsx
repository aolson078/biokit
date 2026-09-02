import { useEffect, useMemo, useReducer, useRef, useState } from "react";
import type { FormEvent, ReactNode } from "react";
import { apiError, apiRequest } from "./api";
import { initialState, reducer, type WorkspaceState } from "./state";
import type { AnalysisCatalogEntry, CandidateRecord, SearchResult, WorkspaceMode } from "./types";

type AppProps = { initialMode: WorkspaceMode; userKey: string; csrfToken: string };
type CatalogResponse = { analyses: AnalysisCatalogEntry[] };
type SavedResponse = { records: CandidateRecord[] };

const STORAGE_SCHEMA = 1;
const MAX_IMPORT_BYTES = 2 * 1024 * 1024;

export function App({ initialMode, userKey, csrfToken }: AppProps) {
  const storageKey = `biokit.sequence-workspace.v${STORAGE_SCHEMA}.${userKey}`;
  const [state, dispatch] = useReducer(reducer, initialMode, (mode) => restoreState(storageKey, mode));
  const [searchQuery, setSearchQuery] = useState("");
  const [organism, setOrganism] = useState("");
  const [importText, setImportText] = useState("");
  const [moleculeHint, setMoleculeHint] = useState<"" | "dna" | "rna">("");
  const [analysisParameters, setAnalysisParameters] = useState({ min_length: 90, frame: 1, genetic_code: 1 });
  const [comparisonIds, setComparisonIds] = useState<number[]>([]);
  const [comparisonAnalysisId, setComparisonAnalysisId] = useState("stacked_composition");
  const activeRequests = useRef(new Set<AbortController>());
  const latestSearch = useRef<AbortController | null>(null);
  const requestedAnalysisId = useRef(new URLSearchParams(window.location.search).get("analysis")).current;

  const activeRecord = state.activeRecordId ? state.recordsById[state.activeRecordId] : null;
  const enabledAnalyses = useMemo(
    () => state.catalog.filter((entry) => entry.id !== "summary" && entry.min_records === 1),
    [state.catalog],
  );
  const hasUnsavedRecords = state.recordOrder.some((id) => !state.recordsById[id].id);

  useEffect(() => {
    const controllers = activeRequests.current;
    const controller = new AbortController();
    controllers.add(controller);
    const catalogRequest = apiRequest<CatalogResponse>("/api/v1/analyses", csrfToken, { signal: controller.signal })
      .then((catalog) => {
        if (!controller.signal.aborted) dispatch({ type: "catalogReceived", catalog: catalog.analyses });
      })
      .catch((error) => {
        if (!controller.signal.aborted) dispatch({ type: "noticeShown", notice: apiError(error).message });
      });
    const recordsRequest = apiRequest<SavedResponse>("/api/v1/records", csrfToken, { signal: controller.signal })
      .then((saved) => {
        if (!controller.signal.aborted) dispatch({ type: "savedRecordsReceived", records: saved.records });
      })
      .catch((error) => {
        if (!controller.signal.aborted) dispatch({ type: "noticeShown", notice: apiError(error).message });
      });
    void Promise.allSettled([catalogRequest, recordsRequest]).finally(() => controllers.delete(controller));

    return () => {
      for (const active of controllers) active.abort();
      controllers.clear();
    };
  }, [csrfToken]);

  useEffect(() => {
    const recoverable = {
      schema: STORAGE_SCHEMA,
      mode: state.mode,
      wizardStep: state.wizardStep,
      recordsById: state.recordsById,
      recordOrder: state.recordOrder,
      activeRecordId: state.activeRecordId,
      selectedAnalysisId: state.selectedAnalysisId,
      analysis: state.analysis,
    };
    sessionStorage.setItem(storageKey, JSON.stringify(recoverable));
  }, [state.mode, state.wizardStep, state.recordsById, state.recordOrder, state.activeRecordId, state.selectedAnalysisId, state.analysis, storageKey]);

  useEffect(() => {
    function handlePopState() {
      const mode = window.location.pathname.endsWith("workbench") ? "workbench" : "wizard";
      dispatch({ type: "modeChanged", mode });
    }
    window.addEventListener("popstate", handlePopState);
    return () => window.removeEventListener("popstate", handlePopState);
  }, []);

  useEffect(() => {
    if (!hasUnsavedRecords) return;
    function warnBeforeUnload(event: BeforeUnloadEvent) {
      event.preventDefault();
      event.returnValue = "";
    }
    window.addEventListener("beforeunload", warnBeforeUnload);
    return () => window.removeEventListener("beforeunload", warnBeforeUnload);
  }, [hasUnsavedRecords]);

  useEffect(() => {
    const frame = window.requestAnimationFrame(() => {
      const target = document.querySelector<HTMLElement>("main .panel");
      if (target) {
        target.tabIndex = -1;
        target.focus();
      }
    });
    return () => window.cancelAnimationFrame(frame);
  }, [state.mode, state.wizardStep]);

  async function handleSearch(event: FormEvent) {
    event.preventDefault();
    latestSearch.current?.abort();
    const controller = new AbortController();
    latestSearch.current = controller;
    activeRequests.current.add(controller);
    dispatch({ type: "searchRequested" });
    try {
      const response = await apiRequest<{ results: SearchResult[] }>("/api/v1/nucleotides/search", csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify({ query: searchQuery, organism: organism || undefined, page: 1, page_size: 10 }),
      });
      if (latestSearch.current === controller) dispatch({ type: "searchSucceeded", results: response.results });
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "searchFailed", error: apiError(error) });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function handleOpenSearchResult(result: SearchResult) {
    const controller = requestController();
    try {
      const record = await apiRequest<CandidateRecord>(`/api/v1/nucleotides/ncbi/${encodeURIComponent(result.accession_version)}`, csrfToken, { signal: controller.signal });
      dispatch({ type: "recordsAcquired", records: [record] });
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "noticeShown", notice: apiError(error).message });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function handleImport(event: FormEvent) {
    event.preventDefault();
    if (new TextEncoder().encode(importText).byteLength > MAX_IMPORT_BYTES) {
      dispatch({ type: "noticeShown", notice: "Import content cannot exceed 2 MiB." });
      return;
    }
    const controller = requestController();
    try {
      const response = await apiRequest<{ records: Array<CandidateRecord | { status: "error"; error: { message: string } }> }>("/api/v1/nucleotides/import", csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify({ text: importText, molecule_type: moleculeHint || undefined, source_name: "Manual import" }),
      });
      const ready = response.records.filter((record): record is CandidateRecord => record.status !== "error");
      const failed = response.records.length - ready.length;
      if (ready.length) dispatch({ type: "recordsAcquired", records: ready });
      dispatch({ type: "noticeShown", notice: failed ? `${ready.length} imported; ${failed} record(s) need attention.` : `${ready.length} record(s) imported.` });
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "noticeShown", notice: apiError(error).message });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function handleImportFile(file: File) {
    if (file.size > MAX_IMPORT_BYTES) {
      dispatch({ type: "noticeShown", notice: "Import files cannot exceed 2 MiB." });
      return;
    }
    setImportText(await file.text());
  }

  async function handleSelectRecord(recordId: string) {
    const record = state.recordsById[recordId];
    if (!record) return;
    dispatch({ type: "recordSelected", recordId });
    if (requestedAnalysisId && state.catalog.some(
      (entry) => entry.id === requestedAnalysisId
        && entry.status === "verified"
        && entry.molecule_types.includes(record.molecule_type),
    )) {
      dispatch({ type: "analysisSelected", analysisId: requestedAnalysisId });
    }
    dispatch({ type: "summaryRequested" });
    const controller = requestController();
    try {
      const response = await apiRequest<{ result: unknown }>("/api/v1/analyses/summary", csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify({ record }),
      });
      dispatch({ type: "summarySucceeded", result: response.result });
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "summaryFailed", error: apiError(error) });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function handleRunAnalysis() {
    if (!activeRecord || !state.selectedAnalysisId) return;
    const analysisId = state.selectedAnalysisId;
    dispatch({ type: "analysisRequested", analysisId });
    const controller = requestController();
    try {
      const response = await apiRequest<{ result?: unknown; job_id?: string }>(`/api/v1/analyses/${analysisId}`, csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify({ record: activeRecord, parameters: analysisParameters }),
      });
      if (response.job_id) {
        dispatch({ type: "analysisQueued", analysisId, jobId: response.job_id });
        await pollJob(response.job_id, analysisId, controller.signal);
      } else {
        dispatch({ type: "analysisSucceeded", analysisId, result: response.result });
      }
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "analysisFailed", analysisId, error: apiError(error) });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function pollJob(jobId: string, analysisId: string, signal: AbortSignal) {
    while (!signal.aborted) {
      const job = await apiRequest<{ status: string; result?: unknown; error?: { code: string; message: string; retryable: boolean } }>(`/api/v1/jobs/${jobId}`, csrfToken, { signal });
      if (job.status === "completed") {
        dispatch({ type: "analysisSucceeded", analysisId, result: job.result });
        return;
      }
      if (job.status === "failed") {
        dispatch({ type: "analysisFailed", analysisId, error: job.error ?? { code: "job_failed", message: "Analysis failed.", retryable: false } });
        return;
      }
      await abortableDelay(1000, signal);
    }
  }

  async function handleSave() {
    if (!activeRecord) return;
    const controller = requestController();
    try {
      const payload = activeRecord.source === "ncbi"
        ? { source: "ncbi", source_accession: activeRecord.source_accession, user_label: activeRecord.user_label }
        : { ...activeRecord, source: "manual" };
      const response = await apiRequest<{ record: CandidateRecord }>("/api/v1/records", csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify(payload),
      });
      dispatch({ type: "recordSaved", clientId: activeRecord.client_id, record: response.record });
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "noticeShown", notice: apiError(error).message });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  async function handleCompare() {
    if (comparisonIds.length < 2) return;
    const controller = requestController();
    try {
      const response = await apiRequest<{ job_id: string }>("/api/v1/comparisons", csrfToken, {
        method: "POST",
        signal: controller.signal,
        body: JSON.stringify({ record_ids: comparisonIds, analysis_id: comparisonAnalysisId }),
      });
      dispatch({ type: "analysisQueued", analysisId: comparisonAnalysisId, jobId: response.job_id });
      await pollJob(response.job_id, comparisonAnalysisId, controller.signal);
    } catch (error) {
      if (!controller.signal.aborted) dispatch({ type: "analysisFailed", analysisId: comparisonAnalysisId, error: apiError(error) });
    } finally {
      activeRequests.current.delete(controller);
    }
  }

  function requestController() {
    const controller = new AbortController();
    activeRequests.current.add(controller);
    return controller;
  }

  function changeMode(mode: WorkspaceMode) {
    const path = mode === "wizard" ? "/sequence-wizard" : "/sequence-workbench";
    window.history.pushState({}, "", path);
    dispatch({ type: "modeChanged", mode });
  }

  function changeAnalysisParameters(value: { min_length: number; frame: number; genetic_code: number }) {
    setAnalysisParameters(value);
    if (state.selectedAnalysisId) {
      dispatch({ type: "analysisSelected", analysisId: state.selectedAnalysisId });
    }
  }

  function resetWorkspace() {
    sessionStorage.removeItem(storageKey);
    dispatch({ type: "reset", mode: state.mode });
  }

  return (
    <div className="workspace-shell">
      <header className="workspace-header">
        <a className="brand" href="/">BioKit</a>
        <nav aria-label="Sequence workspace views">
          <button className={state.mode === "wizard" ? "active" : ""} onClick={() => changeMode("wizard")}>Guided wizard</button>
          <button className={state.mode === "workbench" ? "active" : ""} onClick={() => changeMode("workbench")}>Advanced workbench</button>
        </nav>
        <div className="header-actions">
          <button onClick={resetWorkspace}>Reset</button>
          <a href="/logout" onClick={() => sessionStorage.removeItem(storageKey)}>Logout</a>
        </div>
      </header>

      <main>
        {state.notice && <div className="notice" role="status">{state.notice}<button aria-label="Dismiss message" onClick={() => dispatch({ type: "noticeShown", notice: null })}>×</button></div>}
        {state.mode === "wizard" ? (
          <Wizard
            state={state}
            activeRecord={activeRecord}
            enabledAnalyses={enabledAnalyses}
            acquisition={<Acquisition state={state} searchQuery={searchQuery} organism={organism} importText={importText} moleculeHint={moleculeHint} onSearchQuery={setSearchQuery} onOrganism={setOrganism} onImportText={setImportText} onImportFile={handleImportFile} onMoleculeHint={setMoleculeHint} onSearch={handleSearch} onImport={handleImport} onOpenResult={handleOpenSearchResult} />}
            onSelect={handleSelectRecord}
            onAnalysisSelected={(analysisId) => dispatch({ type: "analysisSelected", analysisId })}
            onRun={handleRunAnalysis}
            onSave={handleSave}
            onBack={() => dispatch({ type: "wizardBacked" })}
            onAdvance={() => dispatch({ type: "wizardAdvanced" })}
            onWorkbench={() => changeMode("workbench")}
            parameters={analysisParameters}
            setParameters={changeAnalysisParameters}
          />
        ) : (
          <Workbench
            state={state}
            activeRecord={activeRecord}
            enabledAnalyses={enabledAnalyses}
            acquisition={<Acquisition state={state} searchQuery={searchQuery} organism={organism} importText={importText} moleculeHint={moleculeHint} onSearchQuery={setSearchQuery} onOrganism={setOrganism} onImportText={setImportText} onImportFile={handleImportFile} onMoleculeHint={setMoleculeHint} onSearch={handleSearch} onImport={handleImport} onOpenResult={handleOpenSearchResult} />}
            onSelect={handleSelectRecord}
            onAnalysisSelected={(analysisId) => dispatch({ type: "analysisSelected", analysisId })}
            onRun={handleRunAnalysis}
            onSave={handleSave}
            comparisonIds={comparisonIds}
            setComparisonIds={setComparisonIds}
            comparisonAnalysisId={comparisonAnalysisId}
            setComparisonAnalysisId={setComparisonAnalysisId}
            onCompare={handleCompare}
            parameters={analysisParameters}
            setParameters={changeAnalysisParameters}
          />
        )}
      </main>
    </div>
  );
}

type AcquisitionProps = {
  state: WorkspaceState;
  searchQuery: string;
  organism: string;
  importText: string;
  moleculeHint: "" | "dna" | "rna";
  onSearchQuery(value: string): void;
  onOrganism(value: string): void;
  onImportText(value: string): void;
  onImportFile(file: File): void;
  onMoleculeHint(value: "" | "dna" | "rna"): void;
  onSearch(event: FormEvent): void;
  onImport(event: FormEvent): void;
  onOpenResult(result: SearchResult): void;
};

function Acquisition(props: AcquisitionProps) {
  return (
    <section className="panel acquisition" aria-labelledby="acquire-title">
      <h2 id="acquire-title">Find or import a nucleotide sequence</h2>
      <p className="muted">Search public NCBI records by gene, organism, or accession. You can also paste DNA/RNA or load a FASTA file.</p>
      <div className="acquire-grid">
        <form onSubmit={props.onSearch}>
          <h3>Search NCBI</h3>
          <label>Gene, description, or accession<input required maxLength={200} value={props.searchQuery} onChange={(event) => props.onSearchQuery(event.target.value)} placeholder="BRCA1, human insulin, NM_007294.4" /></label>
          <label>Organism <span>(optional)</span><input maxLength={100} value={props.organism} onChange={(event) => props.onOrganism(event.target.value)} placeholder="Homo sapiens" /></label>
          <button className="primary" disabled={props.state.search.status === "pending"}>{props.state.search.status === "pending" ? "Searching…" : "Search NCBI"}</button>
          {props.state.search.status === "failed" && <p className="error" role="alert">{props.state.search.error.message}</p>}
        </form>
        <form onSubmit={props.onImport}>
          <h3>Paste or upload FASTA</h3>
          <label>Sequence data<textarea required rows={7} value={props.importText} onChange={(event) => props.onImportText(event.target.value)} placeholder={">example\nATGGCC…"} /></label>
          <label>Ambiguous sequence type<select value={props.moleculeHint} onChange={(event) => props.onMoleculeHint(event.target.value as "" | "dna" | "rna")}><option value="">Detect from T or U</option><option value="dna">DNA</option><option value="rna">RNA</option></select></label>
          <label className="file-control">Load FASTA file<input type="file" accept=".fa,.fasta,.txt,text/plain" onChange={(event) => { const file = event.target.files?.[0]; if (file) void props.onImportFile(file); }} /></label>
          <button className="primary">Import records</button>
        </form>
      </div>
      {props.state.search.status === "succeeded" && (
        <div className="results" aria-label="NCBI search results">
          {props.state.search.data.length === 0 ? <p>No matching records.</p> : props.state.search.data.map((result) => (
            <article key={result.accession_version}>
              <div><strong>{result.accession_version}</strong><p>{result.title}</p><small>{result.organism || "Organism unavailable"} · {result.length.toLocaleString()} nt</small></div>
              <button onClick={() => void props.onOpenResult(result)}>Open record</button>
            </article>
          ))}
        </div>
      )}
    </section>
  );
}

type SharedFlowProps = {
  state: WorkspaceState;
  activeRecord: CandidateRecord | null;
  enabledAnalyses: AnalysisCatalogEntry[];
  acquisition: ReactNode;
  onSelect(id: string): void;
  onAnalysisSelected(id: string): void;
  onRun(): void;
  onSave(): void;
  parameters: { min_length: number; frame: number; genetic_code: number };
  setParameters(value: { min_length: number; frame: number; genetic_code: number }): void;
};

function Wizard(props: SharedFlowProps & { onBack(): void; onAdvance(): void; onWorkbench(): void }) {
  const stepNames = ["Find", "Choose", "Review", "Analyze", "Results"];
  return (
    <div className="wizard-layout">
      <section className="wizard-intro"><p className="eyebrow">Guided sequence analysis</p><h1>Find a sequence and understand it step by step</h1><ol className="steps">{stepNames.map((name, index) => <li key={name} className={props.state.wizardStep === index + 1 ? "current" : props.state.wizardStep > index + 1 ? "complete" : ""}><span>{index + 1}</span>{name}</li>)}</ol></section>
      {props.state.wizardStep === 1 && props.acquisition}
      {props.state.wizardStep === 2 && <RecordChooser state={props.state} activeId={props.activeRecord?.client_id ?? null} onSelect={props.onSelect} />}
      {props.state.wizardStep === 3 && <section className="panel"><h2>Review the sequence</h2><RecordIdentity record={props.activeRecord} /><ResultBlock state={props.state.summary} /><div className="flow-actions"><button onClick={props.onBack}>Back</button><button className="primary" disabled={props.state.summary.status !== "succeeded"} onClick={props.onAdvance}>Choose an analysis</button></div></section>}
      {props.state.wizardStep === 4 && <AnalysisChooser {...props} />}
      {props.state.wizardStep === 5 && <section className="panel"><h2>Analysis result</h2><AnalysisOutput state={props.state} /><div className="flow-actions"><button onClick={props.onBack}>Back</button><button onClick={props.onSave}>Save sequence</button><button className="primary" onClick={props.onWorkbench}>Continue in workbench</button></div></section>}
    </div>
  );
}

function Workbench(props: SharedFlowProps & { comparisonIds: number[]; setComparisonIds(ids: number[]): void; comparisonAnalysisId: string; setComparisonAnalysisId(id: string): void; onCompare(): void }) {
  return (
    <div className="workbench-layout">
      <section className="workbench-title"><p className="eyebrow">Advanced workspace</p><h1>Explore several sequences without losing context</h1></section>
      <div className="workbench-grid">
        <div>{props.acquisition}<RecordChooser state={props.state} activeId={props.activeRecord?.client_id ?? null} onSelect={props.onSelect} /></div>
        <div><section className="panel"><h2>Active sequence</h2><RecordIdentity record={props.activeRecord} />{props.activeRecord && <button onClick={props.onSave}>Save sequence</button>}<ResultBlock state={props.state.summary} /></section><AnalysisChooser {...props} /><section className="panel"><h2>Current result</h2><AnalysisOutput state={props.state} /></section></div>
        <section className="panel saved-tray"><h2>Saved comparison tray</h2>{props.state.savedRecords.length === 0 ? <p className="muted">Save at least two records to compare them.</p> : props.state.savedRecords.map((record) => <label className="saved-row" key={record.id}><input type="checkbox" checked={record.id ? props.comparisonIds.includes(record.id) : false} onChange={(event) => { if (!record.id) return; props.setComparisonIds(event.target.checked ? [...props.comparisonIds, record.id] : props.comparisonIds.filter((id) => id !== record.id)); }} /><span><strong>{record.user_label || record.source_accession || record.nucleotide_id || "Saved sequence"}</strong><small>{record.length.toLocaleString()} nt</small></span></label>)}<label>Comparison tool<select value={props.comparisonAnalysisId} onChange={(event) => props.setComparisonAnalysisId(event.target.value)}><option value="stacked_composition">Comparative composition</option><option value="pairwise_alignment">Pairwise alignment</option><option disabled>Dot plot — awaiting verification</option><option disabled>Heat map — awaiting verification</option><option disabled>Phylogenetic tree — unavailable</option></select></label><button className="primary" disabled={props.comparisonIds.length < 2 || (props.comparisonAnalysisId === "pairwise_alignment" && props.comparisonIds.length !== 2)} onClick={props.onCompare}>Run comparison</button></section>
      </div>
    </div>
  );
}

function RecordChooser({ state, activeId, onSelect }: { state: WorkspaceState; activeId: string | null; onSelect(id: string): void }) {
  return <section className="panel"><h2>Choose a record</h2>{state.recordOrder.length === 0 ? <p className="muted">Search or import records first.</p> : <div className="record-list">{state.recordOrder.map((id) => { const record = state.recordsById[id]; return <button key={id} className={activeId === id ? "selected" : ""} onClick={() => void onSelect(id)}><strong>{record.source_accession || record.source_title || "Manual sequence"}</strong><span>{record.organism || record.molecule_type.toUpperCase()} · {record.length.toLocaleString()} nt</span></button>; })}</div>}</section>;
}

function RecordIdentity({ record }: { record: CandidateRecord | null }) {
  if (!record) return <p className="muted">Choose a record to continue.</p>;
  return <dl className="identity"><div><dt>Record</dt><dd>{record.source_accession || record.source_title || "Manual sequence"}</dd></div><div><dt>Molecule</dt><dd>{record.molecule_type.toUpperCase()}</dd></div><div><dt>Length</dt><dd>{record.length.toLocaleString()} nt</dd></div><div><dt>Source</dt><dd>{record.source.toUpperCase()}</dd></div></dl>;
}

function AnalysisChooser(props: SharedFlowProps) {
  const selected = props.enabledAnalyses.find((entry) => entry.id === props.state.selectedAnalysisId);
  const selectedCompatible = Boolean(selected && props.activeRecord && selected.molecule_types.includes(props.activeRecord.molecule_type));
  return <section className="panel"><h2>Choose an analysis</h2><div className="analysis-list">{props.enabledAnalyses.map((entry) => { const compatible = props.activeRecord ? entry.molecule_types.includes(props.activeRecord.molecule_type) : false; const enabled = entry.status === "verified" && compatible; return <button key={entry.id} disabled={!enabled} className={props.state.selectedAnalysisId === entry.id ? "selected" : ""} onClick={() => props.onAnalysisSelected(entry.id)}><strong>{entry.name}</strong><span>{entry.description}</span>{entry.status !== "verified" && <small>{entry.unavailable_reason}</small>}</button>; })}</div>{selected?.id === "orf" && <label>Minimum ORF length<input type="number" min={30} max={5000} value={props.parameters.min_length} onChange={(event) => props.setParameters({ ...props.parameters, min_length: Number(event.target.value) })} /></label>}{selected?.id === "translation" && <div className="parameter-grid"><label>Frame<select value={props.parameters.frame} onChange={(event) => props.setParameters({ ...props.parameters, frame: Number(event.target.value) })}><option value={1}>1</option><option value={2}>2</option><option value={3}>3</option></select></label><label>Genetic code<select value={props.parameters.genetic_code} onChange={(event) => props.setParameters({ ...props.parameters, genetic_code: Number(event.target.value) })}><option value={1}>Standard</option><option value={2}>Vertebrate mitochondrial</option><option value={3}>Yeast mitochondrial</option><option value={11}>Bacterial</option></select></label></div>}<div className="flow-actions"><button className="primary" disabled={!selectedCompatible || selected?.status !== "verified" || props.state.analysis.status === "pending" || props.state.analysis.status === "queued"} onClick={props.onRun}>Run analysis</button></div></section>;
}

function ResultBlock({ state }: { state: WorkspaceState["summary"] }) {
  if (state.status === "idle") return null;
  if (state.status === "pending") return <p role="status">Calculating summary…</p>;
  if (state.status === "failed") return <p className="error" role="alert">{state.error.message}</p>;
  return <pre className="result-json">{JSON.stringify(state.data, null, 2)}</pre>;
}

function AnalysisOutput({ state }: { state: WorkspaceState }) {
  if (state.analysis.status === "idle") return <p className="muted">Run an analysis to see its result.</p>;
  if (state.analysis.status === "pending") return <p role="status">Starting analysis…</p>;
  if (state.analysis.status === "queued") return <p role="status">Analysis queued. This page will update when it completes.</p>;
  if (state.analysis.status === "failed") return <p className="error" role="alert">{state.analysis.error?.message}</p>;
  return <pre className="result-json">{JSON.stringify(state.analysis.result, null, 2)}</pre>;
}

function restoreState(storageKey: string, mode: WorkspaceMode): WorkspaceState {
  const fallback = initialState(mode);
  try {
    const raw = sessionStorage.getItem(storageKey);
    if (!raw) return fallback;
    const stored = JSON.parse(raw) as Partial<WorkspaceState> & { schema?: number };
    if (stored.schema !== STORAGE_SCHEMA) {
      sessionStorage.removeItem(storageKey);
      return { ...fallback, notice: "A previous workspace version could not be restored and was reset." };
    }
    return {
      ...fallback,
      mode,
      wizardStep: stored.wizardStep ?? 1,
      recordsById: stored.recordsById ?? {},
      recordOrder: stored.recordOrder ?? [],
      activeRecordId: stored.activeRecordId ?? null,
      selectedAnalysisId: stored.selectedAnalysisId ?? null,
      analysis: stored.analysis ?? fallback.analysis,
    };
  } catch {
    sessionStorage.removeItem(storageKey);
    return { ...fallback, notice: "The temporary workspace could not be restored and was reset." };
  }
}

function abortableDelay(milliseconds: number, signal: AbortSignal): Promise<void> {
  return new Promise((resolve, reject) => {
    const timeout = window.setTimeout(resolve, milliseconds);
    signal.addEventListener("abort", () => { window.clearTimeout(timeout); reject(new DOMException("Aborted", "AbortError")); }, { once: true });
  });
}
