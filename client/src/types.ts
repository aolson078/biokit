export type WorkspaceMode = "wizard" | "workbench";
export type MoleculeType = "dna" | "rna" | "unknown";
export type SequenceAlphabet = "dna" | "rna" | "neutral" | "unknown";

export type CandidateRecord = {
  client_id: string;
  id?: number;
  nucleotide_id?: string;
  status?: "ready" | "error";
  source: "ncbi" | "manual" | "legacy";
  source_accession?: string | null;
  source_title?: string | null;
  user_label?: string | null;
  organism?: string | null;
  molecule_type: MoleculeType;
  sequence_alphabet: SequenceAlphabet;
  sequence: string;
  length: number;
  source_retrieved_at?: string | null;
  source_updated_at?: string | null;
  cache?: {
    stale: boolean;
    cached_at: string;
    expires_at: string;
    upstream_error_code?: string | null;
  };
};

export type SearchResult = {
  accession_version: string;
  title: string;
  organism?: string | null;
  length: number;
  molecule_type: MoleculeType;
  updated_at?: string | null;
};

export type AnalysisCatalogEntry = {
  id: string;
  name: string;
  description: string;
  status: "verified" | "exploratory" | "unavailable";
  execution: "sync" | "async" | "adaptive";
  algorithm_version: string;
  min_records: number;
  max_records: number;
  molecule_types: MoleculeType[];
  max_sequence_length: number;
  unavailable_reason: string;
  compatible_inputs: {
    molecule_types: MoleculeType[];
    min_records: number;
    max_records: number;
  };
  parameters: Array<Record<string, unknown>>;
  limits: Record<string, number>;
  units: string;
  interpretation: string;
  limitations: string;
};

export type ApiError = { code: string; message: string; retryable: boolean };
export type JobState = "queued" | "running" | "completed" | "failed" | "expired";
export type RequestState<T> =
  | { status: "idle" }
  | { status: "pending" }
  | { status: "succeeded"; data: T }
  | { status: "failed"; error: ApiError };

export type AnalysisResult = {
  analysisId: string;
  status: "idle" | "pending" | "queued" | "completed" | "failed";
  jobId?: string;
  result?: unknown;
  error?: ApiError;
};
