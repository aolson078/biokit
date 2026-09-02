import { cleanup, render, screen } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { App } from "./App";

function jsonResponse(body: unknown): Response {
  return new Response(JSON.stringify(body), { status: 200, headers: { "Content-Type": "application/json" } });
}

describe("sequence workspace shell", () => {
  beforeEach(() => {
    sessionStorage.clear();
    vi.stubGlobal("fetch", vi.fn((input: RequestInfo | URL) => {
      const url = String(input);
      return Promise.resolve(jsonResponse(url.includes("/records") ? { records: [] } : { analyses: [] }));
    }));
  });

  afterEach(() => {
    cleanup();
    vi.unstubAllGlobals();
  });

  it("renders the novice flow with accessible acquisition choices", () => {
    render(<App initialMode="wizard" userKey="user-1" csrfToken="token" />);
    expect(screen.getByRole("heading", { name: /find a sequence and understand it/i })).toBeTruthy();
    expect(screen.getByRole("heading", { name: "Search NCBI" })).toBeTruthy();
    expect(screen.getByRole("heading", { name: "Paste or upload FASTA" })).toBeTruthy();
  });

  it("discards another schema version and explains the reset", () => {
    sessionStorage.setItem("biokit.sequence-workspace.v1.user-1", JSON.stringify({ schema: 999 }));
    render(<App initialMode="wizard" userKey="user-1" csrfToken="token" />);
    expect(screen.getByRole("status").textContent).toContain("could not be restored");
    expect(sessionStorage.getItem("biokit.sequence-workspace.v1.user-1")).not.toContain("999");
  });
});
