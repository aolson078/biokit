import { describe, expect, it } from "vitest";
import { initialState, reducer } from "./state";
import type { CandidateRecord } from "./types";

const record: CandidateRecord = {
  client_id: "local:one",
  source: "manual",
  source_title: "Example",
  molecule_type: "dna",
  sequence_alphabet: "dna",
  sequence: "ATGAAATAA",
  length: 9,
};

describe("sequence workspace reducer", () => {
  it("advances acquisition and invalidates only downstream state after reselection", () => {
    let state = reducer(initialState("wizard"), { type: "recordsAcquired", records: [record] });
    expect(state.wizardStep).toBe(2);

    state = reducer(state, { type: "recordSelected", recordId: record.client_id });
    state = reducer(state, { type: "summarySucceeded", result: { length: 9 } });
    state = reducer(state, { type: "analysisSelected", analysisId: "orf" });
    state = reducer(state, { type: "analysisSucceeded", analysisId: "orf", result: { count: 1 } });
    expect(state.wizardStep).toBe(5);

    state = reducer(state, { type: "recordSelected", recordId: record.client_id });
    expect(state.recordsById[record.client_id]).toEqual(record);
    expect(state.summary.status).toBe("idle");
    expect(state.analysis.status).toBe("idle");
    expect(state.selectedAnalysisId).toBeNull();
    expect(state.wizardStep).toBe(3);
  });

  it("supports bounded wizard backtracking", () => {
    let state = initialState("wizard");
    state = reducer(state, { type: "wizardBacked" });
    expect(state.wizardStep).toBe(1);
    state = reducer(state, { type: "wizardAdvanced" });
    state = reducer(state, { type: "wizardAdvanced" });
    state = reducer(state, { type: "wizardBacked" });
    expect(state.wizardStep).toBe(2);
  });
});
