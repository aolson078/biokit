import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";
import { resolve } from "node:path";

export default defineConfig({
  plugins: [react()],
  test: {
    environment: "jsdom",
  },
  server: {
    port: 5173,
    strictPort: true,
    cors: { origin: ["http://localhost:5000", "http://127.0.0.1:5000"] },
  },
  build: {
    outDir: resolve(import.meta.dirname, "../static/sequence-workspace"),
    emptyOutDir: true,
    manifest: "manifest.json",
    rollupOptions: { input: resolve(import.meta.dirname, "src/main.tsx") },
  },
});
