import type { ApiError } from "./types";

export class WorkspaceApiError extends Error {
  readonly detail: ApiError;

  constructor(detail: ApiError) {
    super(detail.message);
    this.detail = detail;
  }
}

export async function apiRequest<T>(path: string, csrfToken: string, init: RequestInit = {}): Promise<T> {
  const headers = new Headers(init.headers);
  headers.set("Accept", "application/json");
  if (init.body) headers.set("Content-Type", "application/json");
  if ((init.method ?? "GET").toUpperCase() !== "GET") headers.set("X-CSRFToken", csrfToken);
  const response = await fetch(path, { ...init, headers, credentials: "same-origin" });
  const payload = await response.json().catch(() => ({}));
  if (!response.ok) {
    const detail: ApiError = payload.error ?? { code: "request_failed", message: "The request failed.", retryable: false };
    if (response.status === 401) window.location.assign(`/login?next=${encodeURIComponent(window.location.pathname)}`);
    throw new WorkspaceApiError(detail);
  }
  return payload as T;
}

export function apiError(error: unknown): ApiError {
  if (error instanceof WorkspaceApiError) return error.detail;
  if (error instanceof DOMException && error.name === "AbortError") {
    return { code: "request_cancelled", message: "Request cancelled.", retryable: true };
  }
  return { code: "unexpected_error", message: "An unexpected error occurred.", retryable: false };
}
