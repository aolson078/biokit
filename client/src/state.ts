import type { AnalysisCatalogEntry, AnalysisResult, ApiError, CandidateRecord, RequestState, SearchResult, WorkspaceMode } from "./types";

export type WorkspaceState = {
  mode: WorkspaceMode;
  wizardStep: 1 | 2 | 3 | 4 | 5;
  recordsById: Record<string, CandidateRecord>;
  recordOrder: string[];
  activeRecordId: string | null;
  catalog: AnalysisCatalogEntry[];
  selectedAnalysisId: string | null;
  search: RequestState<SearchResult[]>;
  summary: RequestState<unknown>;
  analysis: AnalysisResult;
  savedRecords: CandidateRecord[];
  notice: string | null;
};

export type WorkspaceAction =
  | { type: "modeChanged"; mode: WorkspaceMode }
  | { type: "searchRequested" }
  | { type: "searchSucceeded"; results: SearchResult[] }
  | { type: "searchFailed"; error: ApiError }
  | { type: "recordsAcquired"; records: CandidateRecord[] }
  | { type: "recordSelected"; recordId: string }
  | { type: "summaryRequested" }
  | { type: "summarySucceeded"; result: unknown }
  | { type: "summaryFailed"; error: ApiError }
  | { type: "catalogReceived"; catalog: AnalysisCatalogEntry[] }
  | { type: "savedRecordsReceived"; records: CandidateRecord[] }
  | { type: "analysisSelected"; analysisId: string }
  | { type: "analysisRequested"; analysisId: string }
  | { type: "analysisQueued"; analysisId: string; jobId: string }
  | { type: "analysisSucceeded"; analysisId: string; result: unknown }
  | { type: "analysisFailed"; analysisId: string; error: ApiError }
  | { type: "recordSaved"; clientId: string; record: CandidateRecord }
  | { type: "wizardAdvanced" }
  | { type: "wizardBacked" }
  | { type: "noticeShown"; notice: string | null }
  | { type: "reset"; mode: WorkspaceMode };

export function initialState(mode: WorkspaceMode): WorkspaceState {
  return {
    mode,
    wizardStep: 1,
    recordsById: {},
    recordOrder: [],
    activeRecordId: null,
    catalog: [],
    selectedAnalysisId: null,
    search: { status: "idle" },
    summary: { status: "idle" },
    analysis: { analysisId: "", status: "idle" },
    savedRecords: [],
    notice: null,
  };
}

export function reducer(state: WorkspaceState, action: WorkspaceAction): WorkspaceState {
  switch (action.type) {
    case "modeChanged":
      return { ...state, mode: action.mode };
    case "searchRequested":
      return { ...state, search: { status: "pending" } };
    case "searchSucceeded":
      return { ...state, search: { status: "succeeded", data: action.results } };
    case "searchFailed":
      return { ...state, search: { status: "failed", error: action.error } };
    case "recordsAcquired": {
      const recordsById = { ...state.recordsById };
      const recordOrder = [...state.recordOrder];
      for (const record of action.records) {
        recordsById[record.client_id] = record;
        if (!recordOrder.includes(record.client_id)) recordOrder.push(record.client_id);
      }
      return { ...state, recordsById, recordOrder, wizardStep: state.mode === "wizard" ? 2 : state.wizardStep };
    }
    case "recordSelected":
      if (!state.recordsById[action.recordId]) return state;
      return {
        ...state,
        activeRecordId: action.recordId,
        summary: { status: "idle" },
        analysis: { analysisId: "", status: "idle" },
        selectedAnalysisId: null,
        wizardStep: state.mode === "wizard" ? 3 : state.wizardStep,
      };
    case "summaryRequested":
      return { ...state, summary: { status: "pending" } };
    case "summarySucceeded":
      return { ...state, summary: { status: "succeeded", data: action.result } };
    case "summaryFailed":
      return { ...state, summary: { status: "failed", error: action.error } };
    case "catalogReceived":
      return { ...state, catalog: action.catalog };
    case "savedRecordsReceived":
      return { ...state, savedRecords: action.records };
    case "analysisSelected":
      return { ...state, selectedAnalysisId: action.analysisId, analysis: { analysisId: action.analysisId, status: "idle" } };
    case "analysisRequested":
      return { ...state, analysis: { analysisId: action.analysisId, status: "pending" } };
    case "analysisQueued":
      return { ...state, analysis: { analysisId: action.analysisId, status: "queued", jobId: action.jobId } };
    case "analysisSucceeded":
      return { ...state, analysis: { analysisId: action.analysisId, status: "completed", result: action.result }, wizardStep: state.mode === "wizard" ? 5 : state.wizardStep };
    case "analysisFailed":
      return { ...state, analysis: { analysisId: action.analysisId, status: "failed", error: action.error } };
    case "recordSaved": {
      const recordsById = { ...state.recordsById, [action.clientId]: { ...state.recordsById[action.clientId], id: action.record.id } };
      const savedRecords = [action.record, ...state.savedRecords.filter((record) => record.id !== action.record.id)];
      return { ...state, recordsById, savedRecords, notice: "Sequence saved." };
    }
    case "wizardAdvanced":
      return { ...state, wizardStep: Math.min(5, state.wizardStep + 1) as WorkspaceState["wizardStep"] };
    case "wizardBacked":
      return { ...state, wizardStep: Math.max(1, state.wizardStep - 1) as WorkspaceState["wizardStep"] };
    case "noticeShown":
      return { ...state, notice: action.notice };
    case "reset":
      return initialState(action.mode);
    default:
      return assertNever(action);
  }
}

function assertNever(value: never): never {
  throw new Error(`Unhandled workspace action: ${JSON.stringify(value)}`);
}
