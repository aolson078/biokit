import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { App } from "./App";
import "./styles.css";

const rootElement = document.getElementById("sequence-workspace-root");
if (!rootElement) throw new Error("Sequence workspace root is missing.");

const mode = rootElement.dataset.mode === "workbench" ? "workbench" : "wizard";
const userKey = rootElement.dataset.userKey ?? "unknown";
const csrfToken = document.querySelector<HTMLMetaElement>('meta[name="csrf-token"]')?.content ?? "";

createRoot(rootElement).render(
  <StrictMode>
    <App initialMode={mode} userKey={userKey} csrfToken={csrfToken} />
  </StrictMode>,
);
